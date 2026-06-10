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

import pandas as pd

from riana.config import FitConfig, IntegrationConfig
from riana.exceptions import DataError
from riana.io.manifest import MANIFEST_FILENAME, ManifestRow, append_manifest, read_manifest
from riana.io.mzml import IndexedMzML, list_mzml_files, mzml_stem
from riana.io.mztab import read_mztab
from riana.io.sdrf import SdrfTable
from riana.io.writers import make_provenance, write_dataframe_tsv
from riana.records import PSMRecord, RunIdentity

_LOGGER = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# integrate
# --------------------------------------------------------------------------- #
def integrate_project(
    config: IntegrationConfig,
    sdrf: SdrfTable,
    mzml_dir: str | os.PathLike[str],
    mztab_path: str | os.PathLike[str],
    out_dir: str | os.PathLike[str],
    *,
    manifest_path: str | os.PathLike[str] | None = None,
    max_workers: int = 1,
    logger: logging.Logger | None = None,
) -> list[ManifestRow]:
    """Integrate every run in *mztab_path*, keyed by the *sdrf* identity.

    Writes one ``<mzml_stem>_riana.txt`` per run (full identity in its provenance
    header) and appends an ``integrate`` row per run to the manifest. Returns the
    rows it wrote.

    Args:
        max_workers: number of runs to integrate concurrently. Each run holds one
            mzML in memory, so this is the file-parallelism knob bounded by the
            one-mzML-per-worker memory ceiling (Track A; pick 2–4 for the big
            animal time series, where a dozen ~350 MB mzMLs are embarrassingly
            parallel). ``1`` (default) keeps the serial, deterministic path the
            tests pin. Output files and manifest rows are always written in
            ``file_idx`` order regardless of completion order, so the result is
            independent of *max_workers*.
    """
    log = logger or _LOGGER
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if manifest_path is None:
        manifest_path = out_dir / MANIFEST_FILENAME

    all_psms, file_index_map = read_mztab(
        mztab_path, sdrf.sample_map, ignored_mods=config.ignored_mods
    )
    if not all_psms:
        raise DataError(f"no PSMs parsed from {mztab_path}")

    mzml_index = _index_mzml_dir(mzml_dir)
    by_file_idx = _group_by_file_idx(all_psms)

    # Resolve identity + mzML path for every run up front so a missing SDRF
    # entry / mzML fails before any (potentially long) integration starts.
    plan: list[tuple[int, str, RunIdentity, str, list[PSMRecord]]] = []
    for file_idx in sorted(by_file_idx):
        stem = file_index_map.get(file_idx, "")
        identity = sdrf.sample_map.get(stem)
        if identity is None:  # read_mztab already guards this; defensive.
            raise DataError(f"ms_run {file_idx} ({stem!r}) has no SDRF identity.")
        if stem not in mzml_index:
            raise DataError(
                f"no mzML for run {stem!r} in {mzml_dir} "
                f"(have {sorted(mzml_index)[:5]}...)."
            )
        fraction = _assign_pep_ids(by_file_idx[file_idx])
        plan.append((file_idx, stem, identity, mzml_index[stem], fraction))

    n_parallel = max(1, min(max_workers, len(plan)))
    if n_parallel > 1:
        log.info("integrating %d runs, %d at a time", len(plan), n_parallel)
        dfs: dict[int, pd.DataFrame] = {}
        with futures.ProcessPoolExecutor(max_workers=n_parallel) as pool:
            future_to_idx = {
                pool.submit(_integrate_one_run, config, mzml_path, fraction, stem):
                (file_idx, stem)
                for file_idx, stem, _identity, mzml_path, fraction in plan
            }
            for fut in futures.as_completed(future_to_idx):
                file_idx, stem = future_to_idx[fut]
                dfs[file_idx] = fut.result()
                log.info("integrated run %s", stem)
    else:
        dfs = {}
        for file_idx, stem, identity, mzml_path, fraction in plan:
            log.info("integrating run %s (%s)", stem, _identity_brief(identity))
            dfs[file_idx] = _integrate_one_run(config, mzml_path, fraction, stem)

    # Write outputs + manifest in file_idx order — deterministic regardless of
    # which worker finished first.
    rows: list[ManifestRow] = []
    for file_idx, stem, identity, _mzml_path, _fraction in plan:
        out_file = out_dir / f"{stem}_riana.txt"
        provenance = make_provenance(
            dataclasses.asdict(config),
            id_source=str(mztab_path),
            extra=identity_to_extra(identity),
        )
        write_dataframe_tsv(out_file, dfs[file_idx], provenance, include_index=True)
        rows.append(
            ManifestRow(
                stage="integrate",
                output_path=str(out_file),
                identity=identity,
                config_hash=provenance.config_hash,
                git_sha=provenance.git_sha,
            )
        )
        log.info("wrote %s", out_file)

    append_manifest(manifest_path, rows)
    log.info("appended %d integrate rows to %s", len(rows), manifest_path)
    return rows


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
) -> dict[tuple[str, str], pd.DataFrame]:
    """Assemble per-curve fit frames from ``integrate`` manifest rows.

    Groups runs into one curve per ``(experiment, condition)``; within a curve,
    merges fractions of the same ``(biological_replicate, labeling_time)`` at
    peptide level (summing the ``isoN`` channels), keeping different biological
    replicates as independent rows. Each returned frame carries a numeric
    ``labeling_time`` column (the fit x-axis) and is ready for
    :func:`riana.core.fitting.fit_run` with ``time_column="labeling_time"``.

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
        curves.setdefault((ident.experiment, ident.condition), []).append(df)

    return {key: _merge_fractions(frames) for key, frames in curves.items()}


def fit_project(
    config: FitConfig,
    manifest_path: str | os.PathLike[str],
    aa_coefficients: Mapping[str, float],
    *,
    logger: logging.Logger | None = None,
    **fit_kwargs,
) -> pd.DataFrame:
    """Fit every kinetic curve indexed by the manifest, tagged by condition.

    Reads ``stage == "integrate"`` rows, recombines them into curves, and runs
    :func:`riana.core.fitting.fit_run` per curve with the identity-supplied
    x-axis. Results across curves are concatenated with ``experiment`` /
    ``condition`` columns so side-by-side groups are distinguishable.
    """
    from riana.core.fitting import fit_run

    log = logger or _LOGGER
    integrate_rows = read_manifest(manifest_path, stage="integrate")
    if not integrate_rows:
        raise DataError(f"no integrate rows in manifest {manifest_path}")
    curves = recombine_for_fit(integrate_rows)

    results: list[pd.DataFrame] = []
    for (experiment, condition), frame in sorted(curves.items()):
        log.info(
            "fitting curve experiment=%s condition=%s (%d rows)",
            experiment, condition or "-", len(frame),
        )
        try:
            result = fit_run(
                config, [frame], aa_coefficients,
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
        result = result.copy()
        result["experiment"] = experiment
        result["condition"] = condition
        results.append(result)

    if not results:
        raise DataError(
            f"no fittable curves in manifest {manifest_path} "
            "(check --q-value / --depth and that runs have ≥ depth timepoints)."
        )
    return pd.concat(results)


# --------------------------------------------------------------------------- #
# internals
# --------------------------------------------------------------------------- #
def _merge_fractions(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Sum ``isoN`` channels across fractions of the same point.

    A "point" is one ``(concat, biological_replicate, labeling_time)``: fractions
    of the same sample at the same timepoint are summed (more signal, one
    envelope to solve), while different biological replicates stay separate. When
    every run is single-fraction (the common case) this is a structural no-op.
    """
    rdf = pd.concat(frames, ignore_index=True)
    iso_cols = sorted(
        (c for c in rdf.columns if _is_iso_col(c)), key=lambda c: int(c[3:])
    )
    point_keys = ["concat", "biological_replicate", "labeling_time"]
    if not iso_cols or rdf.empty:
        return rdf
    # Nothing to merge if every point already appears once.
    if not rdf.duplicated(subset=point_keys).any():
        return rdf

    agg: dict[str, object] = {c: "sum" for c in iso_cols}
    for c in rdf.columns:
        if c in point_keys or c in iso_cols:
            continue
        agg[c] = "min" if c == "percolator q-value" else "first"
    merged = rdf.groupby(point_keys, as_index=False, sort=False).agg(agg)
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
