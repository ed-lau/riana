# -*- coding: utf-8 -*-

"""Shared project orchestration — the SDRF/manifest spine (M6a).

The per-fraction integrate loop was inlined in :func:`riana.cli.integrate` and
*mirrored* in :mod:`riana.gui.tasks`; the GUI "mirror, don't refactor" note
flagged that duplication. This module is the single extraction both the CLI and
GUI call (PROJECT_REVIEW.md §3, Track A), plus the **fit-time recombination**
that replaces parsing the timepoint out of the ``sample`` string with the
header-authoritative :class:`~riana.records.RunIdentity` carried on the
``riana_manifest.tsv``.

Two entry points:

- :func:`integrate_project` — drive the SDRF → mzTab → per-run integration:
  one ``<mzml_stem>_riana.txt`` per run with the full identity frozen into its
  provenance header, and one ``integrate`` row per run appended to the manifest.
- :func:`fit_project` — read the manifest's ``integrate`` rows, group runs into
  kinetic curves (one per ``(experiment, condition)``), merge fractions of the
  same ``(condition, biological_replicate, labeling_time)`` at peptide level,
  and fit each curve with the x-axis taken from the identity. Different
  biological replicates at the same timepoint stay as **independent replicate
  points** on the curve (locked decision; user-confirmed 2026-06-07).
"""

from __future__ import annotations

import dataclasses
import logging
import os
from concurrent import futures
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from riana.config import FitConfig, IntegrationConfig
from riana.exceptions import DataError
from riana.io.manifest import (
    _STAGES,
    MANIFEST_FILENAME,
    ManifestRow,
    append_manifest,
    read_manifest,
)
from riana.io.mzml import IndexedMzML, list_mzml_files, mzml_stem
from riana.io.mztab import read_mztab
from riana.io.sdrf import SdrfTable
from riana.io.writers import hash_config, make_provenance, write_dataframe_tsv
from riana.records import (
    CURVE_KEY_COLUMNS,
    GROUP_KEY_COLUMNS,
    PSMRecord,
    RunIdentity,
)

_LOGGER = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# integrate
# --------------------------------------------------------------------------- #
@dataclasses.dataclass(frozen=True)
class RunTask:
    """One run to integrate: its identity, mzML path, and per-run PSMs.

    The unit of the integrate *plan*. Picklable (``PSMRecord`` / ``RunIdentity``
    are frozen dataclasses), so the GUI can carry tasks back from a planning
    worker and submit each to its *own* shared pool — the executor-agnostic half
    of the plan/dispatch split (Track E), which is how the GUI gets cross-file
    parallelism without nesting a process pool inside ``integrate_project``.
    """

    file_idx: int
    stem: str
    identity: RunIdentity
    mzml_path: str
    psms: list[PSMRecord]


def plan_integration(
    config: IntegrationConfig,
    sdrf: SdrfTable,
    mzml_dir: str | os.PathLike[str],
    mztab_path: str | os.PathLike[str],
) -> list[RunTask]:
    """Resolve the SDRF + mzTab into per-run :class:`RunTask`s (no integration).

    The cheap planning half of :func:`integrate_project`: parse the ID file
    (DDA mzTab or — when the SDRF declares ``acquisition == "DIA"`` — a DIA-NN
    ``report.parquet``), group PSMs by run, attach each run's SDRF identity, and
    resolve its mzML — failing fast on a missing identity / mzML before any
    (long) integration starts. The CLI (:func:`integrate_project`) and the GUI
    build the *same* tasks here, so the identity/manifest logic can't drift
    between the surfaces.
    """
    if sdrf.acquisition == "DIA":
        if sdrf.is_multiplexed:
            raise DataError(
                "sample-axis multiplexing (dimethyl/SILAC) is not supported for DIA "
                "acquisition yet — channel routing needs the DDA mzTab per-PSM mods."
            )
        from riana.io.diann import read_diann

        all_psms, file_index_map = read_diann(
            mztab_path, sdrf.sample_map
        )
    elif sdrf.is_multiplexed:
        # Sample-axis multiplex: route each PSM to its channel identity by its own
        # label mod (distinct sample/condition/enrichment), not one identity per file.
        all_psms, file_index_map = read_mztab(
            mztab_path, channel_map=sdrf.multiplex_channel_map
        )
    else:
        all_psms, file_index_map = read_mztab(
            mztab_path, sdrf.sample_map
        )
    if not all_psms:
        raise DataError(f"no PSMs parsed from {mztab_path}")

    # Match-between-runs (mzTab/DDA only): fill curve holes by transferring
    # confident precursors into runs that missed them. The mzTab is whole-
    # experiment, so the cross-run donor assembly is free here; transferred rows
    # (scan=-1, evidence="mbr") then flow through the same RT-anchored extraction
    # DIA uses. DIA-NN already propagates across runs, so MBR no-ops there.
    if config.mbr:
        if sdrf.acquisition == "DIA":
            _LOGGER.info(
                "MBR requested but acquisition is DIA (DIA-NN already "
                "propagates across runs) — skipping MBR."
            )
        else:
            from riana.core.mbr import augment_with_mbr

            all_psms = augment_with_mbr(all_psms, config)

    mzml_index = _index_mzml_dir(mzml_dir)
    by_file_idx = _group_by_file_idx(all_psms)

    tasks: list[RunTask] = []
    for file_idx in sorted(by_file_idx):
        stem = file_index_map.get(file_idx, "")
        if stem not in mzml_index:
            raise DataError(
                f"no mzML for run {stem!r} in {mzml_dir} "
                f"(have {sorted(mzml_index)[:5]}...)."
            )
        file_psms = by_file_idx[file_idx]
        if sdrf.is_multiplexed:
            # One task per CHANNEL: each PSM already carries its channel identity
            # (distinct sample / condition / enrichment) from the reader; the physical
            # mzML is shared, but each channel is a separate run at its own m/z. The
            # output stem is the channel's SDRF `source name` (unique per channel).
            by_identity: dict[RunIdentity, list[PSMRecord]] = {}
            for p in file_psms:
                by_identity.setdefault(p.identity, []).append(p)
            for identity, ch_psms in by_identity.items():
                fraction = _assign_pep_ids(ch_psms)
                tasks.append(
                    RunTask(file_idx, identity.sample, identity,
                            mzml_index[stem], fraction)
                )
        else:
            identity = sdrf.sample_map.get(stem)
            if identity is None:  # read_mztab already guards this; defensive.
                raise DataError(f"ms_run {file_idx} ({stem!r}) has no SDRF identity.")
            fraction = _assign_pep_ids(file_psms)
            tasks.append(RunTask(file_idx, stem, identity, mzml_index[stem], fraction))
    return tasks


def finalize_run(
    config: IntegrationConfig,
    task: RunTask,
    df: pd.DataFrame,
    out_dir: str | os.PathLike[str],
    mztab_path: str | os.PathLike[str],
) -> ManifestRow:
    """Write one run's ``<stem>_riana.txt`` (identity-stamped) + its manifest row.

    The finalize half of the plan/dispatch split — pure, fast file I/O, so the
    GUI can call it on the main side after gathering each run's integrated frame
    from its pool.
    """
    out_file = Path(out_dir) / f"{task.stem}_riana.txt"
    provenance = make_provenance(
        dataclasses.asdict(config),
        id_source=str(mztab_path),
        extra=identity_to_extra(task.identity),
    )
    write_dataframe_tsv(out_file, df, provenance, include_index=True)
    return ManifestRow(
        stage="integrate",
        output_path=str(out_file),
        identity=task.identity,
        config_hash=provenance.config_hash,
        git_sha=provenance.git_sha,
    )


def integrate_project(
    config: IntegrationConfig,
    sdrf: SdrfTable,
    mzml_dir: str | os.PathLike[str],
    mztab_path: str | os.PathLike[str],
    out_dir: str | os.PathLike[str],
    *,
    manifest_path: str | os.PathLike[str] | None = None,
    max_workers: int = 1,
    resume: bool = False,
    logger: logging.Logger | None = None,
) -> list[ManifestRow]:
    """Integrate every run in *mztab_path*, keyed by the *sdrf* identity.

    Writes one ``<mzml_stem>_riana.txt`` per run (full identity in its provenance
    header) and appends an ``integrate`` row per run to the manifest. Returns the
    rows for every run (kept + freshly written), in ``file_idx`` order.

    **Crash-resilient:** workers only *integrate* (return the frame); the **main
    process writes each ``_riana.txt`` + appends its manifest row as that run
    completes** — so an interruption keeps the runs already finished (the manifest
    is always sorted on write, so the final result is independent of completion
    order). All file/manifest I/O is in this process, so there are no concurrent
    manifest writers.

    Args:
        max_workers: number of runs to integrate concurrently. Each run holds one
            mzML in memory (one-mzML-per-worker ceiling; pick 2–4 for the big
            animal time series). ``1`` (default) is the serial path the tests pin.
        resume: when true, skip any run whose ``<stem>_riana.txt`` already exists
            **and** whose manifest ``integrate`` row matches the current
            ``config_hash`` (settings) — re-running an interrupted run continues
            where it stopped. Assumes the *same* SDRF / mzTab inputs (a settings
            change re-runs everything; the user owns input identity). Default
            false keeps the always-fresh behavior.
    """
    log = logger or _LOGGER
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if manifest_path is None:
        manifest_path = out_dir / MANIFEST_FILENAME

    tasks = plan_integration(config, sdrf, mzml_dir, mztab_path)
    current_hash = hash_config(dataclasses.asdict(config))
    to_run, kept = _resume_partition(
        tasks, manifest_path, out_dir, current_hash, resume, log)

    # Integrate the remaining runs, writing each run's output + recording its
    # manifest row AS IT COMPLETES so an interruption keeps finished runs.
    new_rows: dict[int, ManifestRow] = {}
    for task, df in _integrate_results(config, to_run, max_workers, log):
        # A failure here escapes the generator and re-triggers the pool's
        # ``shutdown(wait=True)`` freeze (see :func:`_integrate_results`); keep
        # writing the rest of the batch if one run's output can't be finalized.
        try:
            row = finalize_run(config, task, df, out_dir, mztab_path)
            append_manifest(manifest_path, [row])
        except Exception as exc:  # noqa: BLE001 — isolate one run
            log.error("integrate: writing run %s FAILED, skipped (%s: %s)",
                      task.stem, type(exc).__name__, exc)
            continue
        new_rows[task.file_idx] = row
        n_mbr = int((df["evidence"] == "mbr").sum()) if "evidence" in df.columns else 0
        dropped = int(df.attrs.get("n_mbr_dropped", 0))
        suffix = (f"  ({n_mbr:,} MBR kept, {dropped:,} gated/no-apex)"
                  if config.mbr and (n_mbr or dropped) else "")
        log.info("wrote %s%s", row.output_path, suffix)

    by_idx = {**kept, **new_rows}
    return [by_idx[t.file_idx] for t in tasks if t.file_idx in by_idx]


def _resume_partition(
    tasks: list[RunTask],
    manifest_path: str | os.PathLike[str],
    out_dir: Path,
    current_hash: str,
    resume: bool,
    log: logging.Logger,
) -> tuple[list[RunTask], dict[int, ManifestRow]]:
    """Split *tasks* into ``(to_run, {file_idx: kept_row})``.

    With ``resume``, a run whose ``<stem>_riana.txt`` exists and whose manifest
    ``integrate`` row matches *current_hash* is kept (skipped); everything else
    runs. Without it, every task runs.
    """
    kept: dict[int, ManifestRow] = {}
    if not resume or not Path(manifest_path).exists():
        return list(tasks), kept
    prev = {r.output_path: r for r in read_manifest(manifest_path)
            if r.stage == "integrate"}
    to_run: list[RunTask] = []
    for task in tasks:
        out_file = out_dir / f"{task.stem}_riana.txt"
        row = prev.get(str(out_file))
        if row is not None and row.config_hash == current_hash and out_file.exists():
            kept[task.file_idx] = row
            log.info("resume: keeping %s (already integrated)", task.stem)
        else:
            to_run.append(task)
    log.info("resume: %d kept, %d to integrate", len(kept), len(to_run))
    return to_run, kept


def _integrate_results(
    config: IntegrationConfig,
    tasks: list[RunTask],
    max_workers: int,
    log: logging.Logger,
):
    """Yield ``(task, df)`` as each run finishes — serial or over a ProcessPool.

    Workers only integrate (return the frame); the caller does all file/manifest
    writes in the main process, so writes are incremental and uncontended, and
    the GUI (which dispatches over its *own* pool) never nests process pools.
    """
    if not tasks:
        return
    n_parallel = max(1, min(max_workers, len(tasks)))
    if n_parallel > 1:
        log.info("integrating %d runs, %d at a time", len(tasks), n_parallel)
        with futures.ProcessPoolExecutor(max_workers=n_parallel) as pool:
            future_to_task = {
                pool.submit(
                    _integrate_one_run, config, t.mzml_path, t.psms, t.stem): t
                for t in tasks
            }
            for fut in futures.as_completed(future_to_task):
                task = future_to_task[fut]
                # Isolate per-run failures. Letting ``fut.result()`` raise here
                # escapes the ``with pool:`` block, whose ``__exit__`` calls
                # ``shutdown(wait=True)`` — which then blocks on *every* other
                # submitted task before the exception can surface. On a 384-run
                # job that is an indefinite freeze (no output) from a single bad
                # file. Log + skip instead so one run can't sink the batch.
                try:
                    df = fut.result()
                except Exception as exc:  # noqa: BLE001 — isolate one run
                    log.error("integrate: run %s FAILED, skipped (%s: %s)",
                              task.stem, type(exc).__name__, exc)
                    continue
                yield task, df
    else:
        for task in tasks:
            log.info("integrating run %s (%s)",
                     task.stem, _identity_brief(task.identity))
            try:
                df = _integrate_one_run(
                    config, task.mzml_path, task.psms, task.stem)
            except Exception as exc:  # noqa: BLE001 — isolate one run
                log.error("integrate: run %s FAILED, skipped (%s: %s)",
                          task.stem, type(exc).__name__, exc)
                continue
            yield task, df


def _integrate_one_run(
    config: IntegrationConfig,
    mzml_path: str,
    psms: list[PSMRecord],
    stem: str,
) -> pd.DataFrame:
    """Integrate one run's PSMs against its mzML — the parallelizable unit.

    Module-level (not a closure) so it pickles cleanly for the
    :class:`ProcessPoolExecutor` path. Opens the mzML inside the worker so only
    the path crosses the process boundary, keeping one mzML in memory per worker.
    """
    from riana.core.integration import integrate_run

    with IndexedMzML(mzml_path) as mzml:
        return integrate_run(config, psms, mzml, file_label=stem)


def identity_to_extra(identity: RunIdentity) -> dict[str, str]:
    """Flatten a :class:`RunIdentity` into provenance-header ``extra`` lines.

    Freezing identity into the header makes a ``_riana.txt`` self-identifying and
    SDRF-independent at fit time (locked decision #2). ``None`` values are
    omitted so a header line is never a fabricated ``None``.
    """
    fields = {
        "experiment": identity.experiment,
        "sample": identity.sample,
        "data_file": identity.data_file,
        "biological_replicate": str(identity.biological_replicate),
        "technical_replicate": str(identity.technical_replicate),
        "fraction": str(identity.fraction),
        "labeling_time": identity.labeling_time,
        "labeling_time_unit": identity.labeling_time_unit,
        "mixing_proportion": identity.mixing_proportion,
        "condition": identity.condition,
        "acquisition": identity.acquisition,
        "precursor_enrichment": identity.precursor_enrichment,
        "experiment_type": identity.experiment_type,
    }
    return {k: str(v) for k, v in fields.items() if v is not None and v != ""}


# --------------------------------------------------------------------------- #
# fit recombination
# --------------------------------------------------------------------------- #
def recombine_for_fit(
    integrate_rows: list[ManifestRow],
    *,
    fraction_collapse: str = "sum",
) -> dict[tuple[str, str], pd.DataFrame]:
    """Assemble per-curve fit frames from ``integrate`` manifest rows.

    Groups runs into one curve per ``(experiment, condition)``; within a curve,
    collapses fractions of the same ``(biological_replicate, labeling_time)`` at
    peptide level (``fraction_collapse`` policy — see :func:`_merge_fractions`),
    keeping different biological replicates as independent rows. Each returned
    frame carries a numeric ``labeling_time`` column (the fit x-axis) and is ready
    for :func:`riana.core.fitting.fit_run` with ``time_column="labeling_time"``.

    Returns ``{(experiment, condition): frame}``.
    """
    curves: dict[tuple[str, str], list[pd.DataFrame]] = {}
    for row in integrate_rows:
        if row.stage != "integrate":
            continue
        ident = row.identity
        x = ident.independent_value
        if x is None:
            raise DataError(
                f"run {ident.data_file!r} has no labeling time / mixing "
                "proportion; cannot place it on a curve."
            )
        df = pd.read_table(row.output_path, comment="#")
        df = df.assign(
            labeling_time=float(x),
            biological_replicate=int(ident.biological_replicate),
            fraction=int(ident.fraction),
        )
        curves.setdefault(ident.group_key, []).append(df)

    return {key: _merge_fractions(frames, policy=fraction_collapse)
            for key, frames in curves.items()}


def fit_project(
    config: FitConfig,
    manifest_path: str | os.PathLike[str],
    aa_coefficients: Mapping[str, float],
    *,
    logger: logging.Logger | None = None,
    ria_override: float | None = None,
    **fit_kwargs,
) -> pd.DataFrame:
    """Fit every kinetic curve indexed by the manifest, tagged by condition.

    Reads ``stage == "integrate"`` rows, recombines them into curves, and runs
    :func:`riana.core.fitting.fit_run` per curve with the identity-supplied
    x-axis. Results across curves are concatenated with ``experiment`` /
    ``condition`` columns so side-by-side groups are distinguishable.
    """
    from riana.core.fitting import build_fractions_long, fit_run

    log = logger or _LOGGER
    integrate_rows = read_manifest(manifest_path, stage="integrate")
    if not integrate_rows:
        raise DataError(f"no integrate rows in manifest {manifest_path}")
    curves = recombine_for_fit(
        integrate_rows, fraction_collapse=config.fraction_collapse)

    # Experiment-type → model dispatch (the decided behavior): a mixing-proportion
    # run is fit with the calibration recovery line, a labeling-time run with the
    # kinetic model. recombine_for_fit already routes the right x-axis
    # (RunIdentity.independent_value → mixing_proportion for calibration); here we
    # pick the matching model per curve, overriding a kinetic default. An explicit
    # ``--model calibration`` is respected on turnover data too (it just stays).
    exp_type = {r.identity.group_key: r.identity.experiment_type
                for r in integrate_rows if r.stage == "integrate"}
    # Per-experiment precursor enrichment (RIA) from the manifest — physically one
    # value per experiment (the SDRF characteristics[precursor enrichment]). Used
    # to shape the labeled envelope per curve unless an explicit ria_override is
    # given (CLI --ria / GUI spin). Load-bearing for ¹⁸O (AC16=0.0583 vs
    # iPSC=0.0897 differ sharply; the old config.ria_max default of 0.06 would
    # mis-shape the envelope and bias every FS/k).
    exp_ria = {r.identity.group_key: r.identity.precursor_enrichment
               for r in integrate_rows if r.stage == "integrate"}

    results: list[pd.DataFrame] = []
    long_frames: list[pd.DataFrame] = []
    for group_key, frame in sorted(curves.items()):
        experiment, condition = group_key
        curve_config = config
        if exp_type.get(group_key) == "calibration" and config.model != "calibration":
            curve_config = dataclasses.replace(curve_config, model="calibration")
            log.info(
                "calibration run → fitting the FS-vs-mixing-proportion recovery "
                "line (experiment=%s condition=%s)", experiment, condition or "-",
            )
        # Resolve the curve's RIA: explicit override > manifest enrichment > default.
        manifest_ria = exp_ria.get(group_key)
        if ria_override is not None:
            curve_ria, ria_src = ria_override, "--ria"
        elif manifest_ria is not None:
            curve_ria, ria_src = manifest_ria, "manifest precursor_enrichment"
        else:
            curve_ria, ria_src = config.ria_max, "default"
        if curve_ria != curve_config.ria_max:
            curve_config = dataclasses.replace(curve_config, ria_max=curve_ria)
        log.info(
            "fitting curve experiment=%s condition=%s (%d rows; RIA=%.4f from %s)",
            experiment, condition or "-", len(frame), curve_ria, ria_src,
        )
        try:
            result = fit_run(
                curve_config, [frame], aa_coefficients,
                time_column="labeling_time", **fit_kwargs,
            )
        except ValueError as exc:
            # One sparse curve (nothing surviving --q-value / --depth) should not
            # abort a multi-condition run — skip it with a clear warning.
            log.warning(
                "skipping curve experiment=%s condition=%s: %s",
                experiment, condition or "-", exc,
            )
            continue
        # M5: tag the per-timepoint long table with the curve identity so
        # side-by-side conditions stay distinguishable in riana_fit_fractions.txt.
        long = result.attrs.get("fractions_long")
        if long is not None and not long.empty:
            long = long.copy()
            for col, val in zip(GROUP_KEY_COLUMNS, group_key):
                long[col] = val
            long_frames.append(long)
        result = result.copy()
        # Drop the per-curve fractions_long DataFrame from .attrs before concat:
        # pandas' concat → __finalize__ reconciles .attrs by equality-comparing
        # values across frames, which raises on differently-shaped DataFrames once
        # there is >1 curve (the two-condition Δk path). It is already captured in
        # long_frames; out.attrs["fractions_long"] is rebuilt after the concat.
        result.attrs.pop("fractions_long", None)
        for col, val in zip(GROUP_KEY_COLUMNS, group_key):
            result[col] = val
        results.append(result)

    if not results:
        raise DataError(
            f"no fittable curves in manifest {manifest_path} "
            "(check --q-value / --depth and that runs have ≥ depth timepoints)."
        )
    out = pd.concat(results)
    # Each row is one fitted peptide curve, uniquely identified by
    # CURVE_KEY_COLUMNS (experiment, condition, concat) — the per-curve key the
    # rollup keys on. Guard against an accidental collision (e.g. two curves with
    # the same concat that failed to be tagged with distinct conditions).
    curve_id = out.reset_index()[list(CURVE_KEY_COLUMNS)]
    if curve_id.duplicated().any():
        n = int(curve_id.duplicated().sum())
        raise DataError(
            f"{n} duplicate {tuple(CURVE_KEY_COLUMNS)} rows in the fit output — "
            "a curve-key collision (two conditions tagged the same?)."
        )
    out.attrs["fractions_long"] = (
        pd.concat(long_frames, ignore_index=True)
        if long_frames
        else build_fractions_long([])
    )
    return out


# --------------------------------------------------------------------------- #
# manifest stage rows (fit / protein) — the project-directory chain
# --------------------------------------------------------------------------- #
def aggregate_identity(result_df: pd.DataFrame) -> RunIdentity:
    """A coarse, experiment-level identity for an aggregate fit / protein output.

    Fit and rollup each emit one file spanning many curves, so there is no
    per-run identity to record on its manifest row; carry the common
    ``experiment`` (and ``condition`` when single), leaving the run-specific
    fields default. This is for manifest *indexing/provenance* only — downstream
    stages group off the file's own ``experiment`` / ``condition`` columns, not
    this identity.
    """
    def _single(col: str) -> str:
        if col in result_df.columns:
            vals = sorted({str(v) for v in result_df[col].dropna().unique()})
            return vals[0] if len(vals) == 1 else ""
        return ""

    return RunIdentity(
        experiment=_single("experiment"), sample="", data_file="",
        condition=_single("condition"),
    )


def record_stage_rows(
    manifest_path: str | os.PathLike[str],
    stage: str,
    output_paths: list[str | os.PathLike[str]],
    result_df: pd.DataFrame,
    provenance,
) -> list[ManifestRow]:
    """Append one *stage* row per output file to the manifest (idempotent upsert).

    Used by ``fit`` (``stage="fit"``) and ``rollup`` (``stage="rollup"``) to
    register their outputs in the same project manifest ``integrate`` wrote — so
    a single ``--manifest`` drives the whole ``integrate → fit → rollup`` chain.
    """
    identity = aggregate_identity(result_df)
    rows = [
        ManifestRow(
            stage=stage, output_path=str(p), identity=identity,
            config_hash=provenance.config_hash, git_sha=provenance.git_sha,
        )
        for p in output_paths
    ]
    append_manifest(manifest_path, rows)
    return rows


def resolve_manifest_write(
    manifest_path: str | os.PathLike[str],
    requested_out_dir: str | os.PathLike[str],
    stage: str,
) -> tuple[Path, Path]:
    """Decide where a ``--manifest`` *stage* writes, and seed a forked manifest.

    On the manifest path the project's outputs live next to the manifest. This
    resolves the requested output dir into ``(out_dir, target_manifest)``:

    - **Default / same folder** — ``requested_out_dir`` empty, ``"."``, or the
      manifest's own folder → *update in place*: write next to the manifest and
      record the stage rows in it (the existing behaviour).
    - **A different folder** → *fork*: write into that folder, into a **new**
      manifest there seeded with the *upstream* stages' rows (``integrate`` for a
      ``fit`` fork; ``integrate`` + ``fit`` for a ``rollup`` fork) with their
      paths resolved to **absolute** so they still point at the already-computed
      originals — and leave the **source manifest untouched**. The forked folder
      is thus a self-contained project (a later ``rollup`` / GUI "Load results"
      works on it), and a pristine input run is never mutated just by pointing
      ``-o`` elsewhere.

    ``append_manifest`` upserts, so re-forking into an existing folder is
    idempotent.
    """
    src = Path(manifest_path).resolve()
    proj_dir = src.parent
    req = str(requested_out_dir)
    if req in ("", ".") or Path(req).resolve() == proj_dir:
        return proj_dir, src  # update the project in place

    out_dir = Path(req)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / MANIFEST_FILENAME
    if target.resolve() == src:  # a different spelling of the same folder
        return proj_dir, src

    # Seed the fork with the upstream stages' rows, paths -> absolute so they keep
    # pointing at the already-computed originals. Manifest paths are stored
    # **relative to the working directory** (that is how the reader at ``run_fit`` /
    # ``pd.read_table(row.output_path)`` resolves them), so anchor a relative path at
    # the cwd — NOT at ``proj_dir`` (the manifest's folder), which double-counts the
    # ``runs/<project>/`` prefix already inside the stored path.
    upstream = _STAGES[:_STAGES.index(stage)]
    base = Path.cwd()
    seed = [
        dataclasses.replace(row, output_path=str(_abs_under(row.output_path, base)))
        for st in upstream
        for row in read_manifest(src, stage=st)
    ]
    # A fork writes its own manifest; merging into a *different* project's manifest
    # already in the folder would make a hybrid that can't be loaded. Refuse that —
    # an empty folder or a prior fork of this same source (its upstream rows match
    # the seed) is fine.
    if target.exists() and _fork_target_conflicts(target, seed):
        raise DataError(
            f"output folder {out_dir} already contains a {MANIFEST_FILENAME} from a "
            f"different project; choose an empty output folder to fork into (a fork "
            f"writes its own manifest and must not merge with another project's)."
        )
    if seed:
        append_manifest(target, seed)
    return out_dir, target


def _fork_target_conflicts(target: Path, seed: list[ManifestRow]) -> bool:
    """True if *target*'s manifest belongs to a different project than *seed*.

    Compares the target's upstream-stage rows (the stages we're seeding) to the
    seed by resolved absolute output path: any the seed doesn't have means the
    folder holds another project's outputs (or the manifest is unreadable).
    """
    seed_keys = {(r.stage, r.output_path) for r in seed}  # seed paths are absolute
    seed_stages = {r.stage for r in seed}
    base = Path.cwd()  # same anchor the seed uses, so relative rows compare correctly
    try:
        existing = read_manifest(target)
    except DataError:
        return True  # a file named like a manifest but unreadable — don't touch it
    for row in existing:
        if row.stage not in seed_stages:
            continue  # ignore the target's own downstream rows
        if (row.stage, str(_abs_under(row.output_path, base))) not in seed_keys:
            return True
    return False


def _abs_under(path_str: str, base: Path) -> Path:
    """Resolve *path_str* to absolute; a relative path is taken under *base*."""
    p = Path(path_str)
    return p if p.is_absolute() else (base / p).resolve()


def fit_outputs_from_manifest(
    manifest_path: str | os.PathLike[str],
) -> tuple[str, str]:
    """Locate the ``riana_fit_peptides.txt`` / ``riana_fit_fractions.txt`` paths
    from a manifest's ``stage == "fit"`` rows (for ``rollup --manifest``)."""
    rows = read_manifest(manifest_path, stage="fit")
    if not rows:
        raise DataError(
            f"no fit rows in manifest {manifest_path} — run `riana fit "
            "--manifest` first."
        )
    peptides = next(
        (r.output_path for r in rows
         if r.output_path.endswith("fit_peptides.txt")), None)
    fractions = next(
        (r.output_path for r in rows
         if r.output_path.endswith("fit_fractions.txt")), None)
    if peptides is None or fractions is None:
        raise DataError(
            f"manifest {manifest_path} fit rows are missing the peptides and/or "
            "fractions output (expected riana_fit_peptides.txt + "
            "riana_fit_fractions.txt)."
        )
    return peptides, fractions


# --------------------------------------------------------------------------- #
# internals
# --------------------------------------------------------------------------- #
#: The point that fractions collapse onto. ``concat`` carries the charge, so this
#: is one ``(peptidoform, charge, biological replicate, labeling time)``.
_FRACTION_POINT_KEYS = ["concat", "biological_replicate", "labeling_time"]


def _merge_fractions(
    frames: list[pd.DataFrame], policy: str = "sum"
) -> pd.DataFrame:
    """Collapse LC fractions / technical replicates of the same point into one row.

    A "point" is one ``(concat, biological_replicate, labeling_time)`` (``concat``
    carries the charge); fractions of the same sample at the same timepoint are one
    kinetic observation, while different biological replicates stay separate. The
    intensities are combined *before* a single FS is solved downstream — fractions
    are never fit as independent points, and their FS values are never averaged.

    ``policy``:

    - ``"sum"`` (default): sum each ``isoN`` channel across fractions (the
      high-abundance fraction dominates, so a minor fraction's noise is naturally
      down-weighted), and **intensity-weight the per-channel mass/QC columns**
      (``iso{N}_obs_mz`` / ``iso{N}_ppm_error`` by that channel's intensity;
      ``apex_snr`` by the row's total intensity; ``n_scans`` → max). Without the
      weighting the summed envelope would carry one arbitrary fraction's masses,
      which would corrupt the mass-defect (``fs_ds``) estimate.
    - ``"anchor"``: keep only the single highest-total-intensity fraction's row,
      discarding the rest (legacy parity).

    When every run is single-fraction (the common case) this is a structural no-op.
    """
    rdf = pd.concat(frames, ignore_index=True)
    iso_cols = sorted(
        (c for c in rdf.columns if _is_iso_col(c)), key=lambda c: int(c[3:])
    )
    point_keys = _FRACTION_POINT_KEYS
    if not iso_cols or rdf.empty:
        return rdf
    # Nothing to collapse if every point already appears once.
    if not rdf.duplicated(subset=point_keys).any():
        return rdf

    if policy == "anchor":
        return _merge_fractions_anchor(rdf, iso_cols, point_keys)
    return _merge_fractions_sum(rdf, iso_cols, point_keys)


def _merge_fractions_anchor(
    rdf: pd.DataFrame, iso_cols: list[str], point_keys: list[str]
) -> pd.DataFrame:
    """Keep only the single highest-total-intensity fraction per point."""
    work = rdf.copy()
    work["__rowtot"] = work[iso_cols].to_numpy(dtype=float).sum(axis=1)
    keep = work.groupby(point_keys, sort=False)["__rowtot"].idxmax()
    return work.loc[keep, rdf.columns].reset_index(drop=True)


def _sum_min1(s: pd.Series) -> float:
    """A sum that keeps an all-NaN group NaN (``min_count=1``) instead of folding it
    to ``0.0``. Under ``--iso auto`` a channel a short peptidoform never extracted is
    NaN-padded ("not a channel", distinct from an integrated ``0.0``); plain ``sum``
    (``min_count=0``) would fabricate an observed-zero channel that the FS solver then
    scores, biasing FS/k low. Fixed-iso mode is unaffected — its empty channels are a
    real ``0.0`` (present), and a partially-padded channel still sums to its seen value.
    """
    return s.sum(min_count=1)


def _merge_fractions_sum(
    rdf: pd.DataFrame, iso_cols: list[str], point_keys: list[str]
) -> pd.DataFrame:
    """Sum ``isoN`` across fractions; intensity-weight the mass/QC columns."""
    work = rdf.copy()
    # Per-channel mass/error columns are intensity-weighted by their own channel;
    # build weighted numerators to sum alongside the channel intensities.
    wnum: dict[str, str] = {}  # source column -> weighted-numerator helper column
    weight_of: dict[str, str] = {}  # source column -> weight (channel) column
    for ic in iso_cols:
        n = ic[3:]
        for col in (f"iso{n}_obs_mz", f"iso{n}_ppm_error"):
            if col in work.columns:
                helper = f"__wnum_{col}"
                work[helper] = (work[col].to_numpy(dtype=float)
                                * work[ic].to_numpy(dtype=float))
                wnum[col] = helper
                weight_of[col] = ic
    row_tot = work[iso_cols].to_numpy(dtype=float).sum(axis=1)
    has_snr = "apex_snr" in work.columns
    if has_snr:
        work["__rowtot"] = row_tot
        work["__wnum_apex_snr"] = work["apex_snr"].to_numpy(dtype=float) * row_tot

    # NaN-preserving sum for the iso channels (see _sum_min1): an all-NaN --iso auto
    # pad must stay NaN, not become a fabricated observed-0.0. The weighted-numerator
    # helpers below stay plain "sum" — their result is divided by the (now NaN)
    # channel sum and guarded by `denom > 0`, so an all-NaN channel yields NaN anyway.
    agg: dict[str, object] = {c: _sum_min1 for c in iso_cols}
    agg.update({h: "sum" for h in wnum.values()})
    if has_snr:
        agg["__rowtot"] = "sum"
        agg["__wnum_apex_snr"] = "sum"
    if "n_scans" in work.columns:
        agg["n_scans"] = "max"
    for c in work.columns:
        if c in point_keys or c in agg or c.startswith("__"):
            continue
        agg[c] = "min" if c == "percolator q-value" else "first"

    merged = work.groupby(point_keys, as_index=False, sort=False).agg(agg)

    # Resolve the intensity-weighted means: numerator / channel-intensity sum.
    for col, helper in wnum.items():
        denom = merged[weight_of[col]].to_numpy(dtype=float)
        merged[col] = np.where(
            denom > 0, merged[helper].to_numpy(dtype=float) / denom, np.nan)
    if has_snr:
        denom = merged["__rowtot"].to_numpy(dtype=float)
        merged["apex_snr"] = np.where(
            denom > 0,
            merged["__wnum_apex_snr"].to_numpy(dtype=float) / denom, np.nan)
    return merged[rdf.columns]


def _is_iso_col(name: str) -> bool:
    return name.startswith("iso") and name[3:].isdigit()


def _index_mzml_dir(mzml_dir: str | os.PathLike[str]) -> dict[str, str]:
    """``{mzml_stem: full path}`` for the mzMLs in *mzml_dir*."""
    out: dict[str, str] = {}
    for name in list_mzml_files(mzml_dir):
        out[mzml_stem(name)] = os.path.join(str(mzml_dir), name)
    return out


def _group_by_file_idx(psms: list[PSMRecord]) -> dict[int, list[PSMRecord]]:
    by: dict[int, list[PSMRecord]] = {}
    for p in psms:
        by.setdefault(p.file_idx, []).append(p)
    return by


def _assign_pep_ids(psms: list[PSMRecord]) -> list[PSMRecord]:
    """Sort one run's PSMs by scan and assign the per-run sequential pep_id."""
    rows = sorted(psms, key=lambda r: r.scan)
    return [dataclasses.replace(r, pep_id=i) for i, r in enumerate(rows)]


def _identity_brief(identity: RunIdentity) -> str:
    if identity.experiment_type == "calibration":
        return f"mix={identity.mixing_proportion} cond={identity.condition or '-'}"
    return (
        f"t={identity.labeling_time}{identity.labeling_time_unit} "
        f"rep={identity.biological_replicate} cond={identity.condition or '-'}"
    )
