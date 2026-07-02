# -*- coding: utf-8 -*-

"""Qt-free worker functions for the GUI's ``ProcessPoolExecutor`` (M4 Phase 2).

These run in worker *processes*, so they must be importable at module level and
take only picklable arguments — which is why they live apart from the widgets
(no PySide6 import here) and take an mzML *path* rather than an open
:class:`~riana.io.mzml.IndexedMzML` (the reader holds a ``threading.local`` and
locks and is built inside the worker).

Each function mirrors the exact sequence :func:`riana.cli.integrate` runs and
calls the *same* core entry points (:func:`riana.core.integration.integrate_run`
and friends), so the GUI and CLI cannot diverge numerically. Because they are
pure and Qt-free, they are unit-tested headlessly with no display.
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd

from riana.algorithms.calibration import DriftSummary
from riana.config import FitConfig, IntegrationConfig
from riana.core.integration import (
    PeptideTrace,
    extract_peptide_trace,
    integrate_run,
)
from riana.exceptions import IntegrationError
from riana.io.mzml import IndexedMzML
from riana.io.percolator import read_percolator
from riana.records import PSMRecord


def read_psms(
    id_path: str,
    sample: str,
) -> list[PSMRecord]:
    """Parse the Percolator id file into typed records (worker side).

    Thin wrapper over :func:`riana.io.percolator.read_percolator` so the GUI can
    do the (potentially slow) parse off the Qt event loop. Same call the CLI
    makes in :func:`riana.cli.integrate`.
    """
    return read_percolator(id_path, sample=sample)


def integrate_fraction(
    config: IntegrationConfig,
    psms: Sequence[PSMRecord],
    mzml_path: str,
    file_label: str,
) -> tuple[pd.DataFrame, DriftSummary | None]:
    """Integrate one fraction's PSMs against ``mzml_path``.

    Opens the indexed reader in-process, calls the shared
    :func:`riana.core.integration.integrate_run`, and returns the per-fraction
    frame plus its :class:`~riana.algorithms.calibration.DriftSummary`. The
    drift is pulled out of ``df.attrs`` here (in-process) and returned
    explicitly, since ``DataFrame.attrs`` is not guaranteed to survive the
    pickle back to the GUI process.
    """
    with IndexedMzML(mzml_path) as mzml:
        df = integrate_run(config, list(psms), mzml, file_label=file_label)
    drift = df.attrs.get("drift_summary")
    return df, drift


def extract_trace(
    config: IntegrationConfig,
    psm: PSMRecord,
    mzml_path: str,
    scan_span: tuple[int, int] | None = None,
) -> PeptideTrace | None:
    """Extract one peptide-charge's XICs + integration window for the plot.

    Returns ``None`` when the peptide has no extractable intensity profile (an
    honest "nothing to show" rather than propagating
    :class:`~riana.exceptions.IntegrationError` across the process boundary).
    """
    try:
        with IndexedMzML(mzml_path) as mzml:
            return extract_peptide_trace(config, psm, mzml, scan_span=scan_span)
    except IntegrationError:
        return None


def run_fit(
    config: FitConfig,
    riana_paths: Sequence[str],
    coefficients: str | None,
) -> pd.DataFrame:
    """Read the integrate-output time series and run the kinetic fit (worker side).

    Mirrors the sequence :func:`riana.cli.fit` runs — load the per-AA coefficient
    table, read one DataFrame per timepoint file, then call the shared
    :func:`riana.core.fitting.fit_run`. Returns the per-peptide result frame
    (indexed by ``concat``, carrying the per-peptide ``t`` / ``fs`` lists used by
    the GUI's fitted-curve plot). Errors (bad coefficients, nothing surviving
    ``--depth``, o18 guard) propagate for the tab to surface.
    """
    from riana.core.fitting import fit_run, load_aa_coefficients

    coeffs = load_aa_coefficients(coefficients) if coefficients else {}
    dfs = [pd.read_table(p, comment="#") for p in riana_paths]
    return fit_run(config, dfs, coeffs)


def plan_sdrf_integration(
    config: IntegrationConfig,
    sdrf_path: str,
    mzml_dir: str,
    mztab_path: str,
) -> list:
    """Read the SDRF and plan the per-run integrate tasks (worker side).

    Mirrors :func:`riana.cli.integrate`'s SDRF planning (``read_sdrf`` →
    :func:`riana.core.pipeline.plan_integration`) off the Qt event loop, since
    the mzTab parse is the slow part. Returns the list of
    :class:`riana.core.pipeline.RunTask` (picklable) the Integrate tab then
    dispatches over its *own* shared pool — so the GUI gets cross-file
    parallelism without nesting a process pool inside ``integrate_project``.
    """
    from riana.core.pipeline import plan_integration
    from riana.io.sdrf import read_sdrf

    sdrf = read_sdrf(sdrf_path)
    return plan_integration(config, sdrf, mzml_dir, mztab_path)


def run_fit_manifest(
    config: FitConfig,
    manifest_path: str,
    coefficients: str | None,
) -> pd.DataFrame:
    """Fit every curve indexed by a manifest (worker side; the SDRF fit path).

    Mirrors :func:`riana.cli.fit`'s ``--manifest`` branch — load the coefficient
    table, then :func:`riana.core.pipeline.fit_project`, which groups the
    manifest's ``integrate`` rows into curves by ``(experiment, condition)`` with
    the timepoint from the SDRF identity. ``fit_project`` fits serially per curve
    (each curve's ``fit_run`` uses an in-process thread pool), so there is no
    nested process pool to worry about. The per-timepoint long table rides on
    ``df.attrs["fractions_long"]``.
    """
    from riana.core.fitting import load_aa_coefficients
    from riana.core.pipeline import fit_project

    coeffs = load_aa_coefficients(coefficients) if coefficients else {}
    # The Model tab's RIA spin is an explicit user value → pass it as the override
    # so it wins over the manifest's per-experiment enrichment (preserves the
    # GUI's current behavior; CLI without --ria defers to the manifest instead).
    return fit_project(config, manifest_path, coeffs, ria_override=config.ria_max)


def run_rollup(
    fit_dir: str,
    model: str,
    kp: float,
    kr: float,
    rp: float,
    parsimony: str,
    min_peptides: int,
    min_points: int,
    min_r2: float | None = None,
    k_cv_max: float = 0.2,
    rescue_r2: float = 0.6,
    method: str = "weighted",
    workers: int = 1,
    phi_limit: float = -4.0,
    reference_condition: str | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Read the ``riana fit`` outputs in *fit_dir* and roll peptides up to proteins.

    Mirrors :func:`riana.cli.rollup`: read ``riana_fit_peptides.txt`` +
    ``riana_fit_fractions.txt`` and call the shared
    :func:`riana.core.protein.rollup_proteins`, so the GUI and CLI cannot
    diverge. Errors (bad model/parsimony, missing files) propagate for the tab
    to surface.

    ``workers`` (process-level parallelism) is forwarded to ``rollup_proteins``.
    **The tab must dispatch this on a main-process thread when ``workers > 1``**
    (``run_in_executor(None, …)``), not the shared ``ProcessPoolExecutor`` — a
    pool worker spawning its own pool is a nested pool, which breaks
    (``BrokenProcessPool``). On a thread the pool is created from the main
    process. ``phi_limit`` / ``reference_condition`` apply only to
    ``model="linear simple"``.

    Returns ``(protein_table, points)`` — ``points`` is the
    ``{(experiment, condition, protein): (t_list, fs_list)}`` collapsed-refit
    map, pulled out of ``DataFrame.attrs`` here (in-process) and returned
    explicitly since ``attrs`` is not guaranteed to survive the pickle back to
    the GUI process (same posture as :func:`integrate_fraction`'s drift).
    """
    from pathlib import Path

    from riana.core.protein import rollup_proteins

    fd = Path(fit_dir)
    peptides = pd.read_table(fd / "riana_fit_peptides.txt", comment="#")
    fractions = pd.read_table(fd / "riana_fit_fractions.txt", comment="#")
    result = rollup_proteins(
        peptides, fractions, model=model, method=method,
        kinetic_kwargs=dict(k_p=kp, k_r=kr, r_p=rp),
        parsimony=parsimony, min_peptides=int(min_peptides),
        min_points=int(min_points), min_r2=min_r2,
        k_cv_max=float(k_cv_max), rescue_r2=float(rescue_r2),
        workers=int(workers), phi_limit=float(phi_limit),
        reference_condition=reference_condition,
    )
    return result, result.attrs.get("protein_points", {})


# --------------------------------------------------------------------------- #
# Load prior results for display (the GUI "Load results" button) — no recompute.
# These reconstruct the exact in-memory shapes `run_fit_manifest` / `run_rollup`
# return, from the on-disk outputs, so the tabs' existing table + curve code
# renders a loaded result identically to a freshly-run one.
# --------------------------------------------------------------------------- #
def load_fit_results(manifest_path: str) -> tuple[pd.DataFrame, dict]:
    """Load a prior ``riana fit`` result for display (no recompute).

    Reconstructs the frame ``fit_project`` returns — concat-indexed, carrying the
    per-timepoint ``t`` / ``fs`` / ``evidence`` / ``metox`` / ``fs_ds`` / ``dmass`` /
    ``dspacing`` list-cells the Model-tab curve needs — from the on-disk
    ``riana_fit_peptides.txt`` (the scalar summary) + ``riana_fit_fractions.txt``
    (the per-timepoint substrate), located from the manifest's ``stage="fit"``
    rows. Also returns the peptides file's provenance header, so the caller can
    plot the fitted curve with the **model that produced it** (only ``k_deg`` is in
    the table; the curve shape needs the model).
    """
    from riana.core.pipeline import fit_outputs_from_manifest
    from riana.io.writers import read_provenance_header

    pep_path, frac_path = fit_outputs_from_manifest(manifest_path)
    peptides = pd.read_table(pep_path, comment="#")
    fractions = pd.read_table(frac_path, comment="#")
    return _reconstruct_fit_wide(peptides, fractions), read_provenance_header(pep_path)


def _reconstruct_fit_wide(
    peptides: pd.DataFrame, fractions: pd.DataFrame
) -> pd.DataFrame:
    """Attach the per-timepoint list-cells from *fractions* onto *peptides*.

    The written peptides file is a scalar summary; the per-timepoint detail lives
    in the fractions file. Grouping the fractions by the curve key and collecting
    each channel back into a list rebuilds the ``t`` / ``fs`` / … list-cells the
    curve plot reads, indexed by ``concat`` exactly like the fitted frame.
    """
    keys = [k for k in ("concat", "experiment", "condition")
            if k in peptides.columns and k in fractions.columns]

    def _iso_cols(prefix: str) -> list[str]:
        return sorted((c for c in fractions.columns if c.startswith(prefix)),
                      key=lambda c: int(c.rsplit("iso", 1)[1]))

    dmass_cols, dspacing_cols = _iso_cols("dmass_iso"), _iso_cols("dspacing_iso")
    # metox round-trips as a real bool via pandas; coerce defensively if it came
    # back as text ("False" is truthy, which would mis-colour every point).
    if "metox" in fractions.columns and fractions["metox"].dtype == object:
        fractions = fractions.copy()
        fractions["metox"] = (fractions["metox"].astype(str).str.strip()
                              .str.lower().isin(("true", "1")))

    ordered = fractions.sort_values("labeling_time", kind="mergesort")
    per: dict[tuple, dict] = {}
    for key_vals, grp in ordered.groupby(keys, sort=False):
        rec = {"t": grp["labeling_time"].tolist(), "fs": grp["fs"].tolist()}
        for col in ("evidence", "metox", "fs_ds"):
            if col in grp.columns:
                rec[col] = grp[col].tolist()
        if dmass_cols:
            rec["dmass"] = grp[dmass_cols].to_numpy().tolist()
        if dspacing_cols:
            rec["dspacing"] = grp[dspacing_cols].to_numpy().tolist()
        per[key_vals if isinstance(key_vals, tuple) else (key_vals,)] = rec

    out = peptides.copy()
    pep_keys = list(zip(*[peptides[k] for k in keys]))  # key tuple per peptides row
    for col in ("t", "fs", "evidence", "metox", "fs_ds", "dmass", "dspacing"):
        out[col] = [per.get(kt, {}).get(col, []) for kt in pep_keys]
    return out.set_index("concat")


def load_rollup_results(manifest_path: str) -> tuple[pd.DataFrame, dict, dict]:
    """Load a prior ``riana rollup`` result for display (no recompute).

    Returns ``(protein_table, points, header)`` — the ``riana_rollup_proteins.txt``
    table, the ``{(experiment, condition, protein): (t_list, fs_list)}`` collapsed
    refit points reconstructed from ``riana_rollup_fractions.txt`` (the same
    ``protein_points`` the worker returns, for the refit / φ-space curve), and the
    proteins file's provenance header (for the model). Located from the manifest's
    ``stage="rollup"`` rows.
    """
    from pathlib import Path

    from riana.exceptions import DataError
    from riana.io.manifest import read_manifest
    from riana.io.writers import read_provenance_header

    rows = read_manifest(manifest_path, stage="rollup")
    prot = next((r.output_path for r in rows
                 if r.output_path.endswith("rollup_proteins.txt")), None)
    frac = next((r.output_path for r in rows
                 if r.output_path.endswith("rollup_fractions.txt")), None)
    if prot is None:
        raise DataError(
            f"no rollup output in manifest {manifest_path} — run rollup first.")
    proteins = pd.read_table(prot, comment="#")

    points: dict[tuple, tuple[list, list]] = {}
    if frac is not None and Path(frac).exists():
        fr = pd.read_table(frac, comment="#")
        keys = [k for k in ("experiment", "condition", "protein")
                if k in fr.columns]
        for key_vals, grp in fr.sort_values(
                "labeling_time", kind="mergesort").groupby(keys, sort=False):
            k = key_vals if isinstance(key_vals, tuple) else (key_vals,)
            points[k] = (grp["labeling_time"].tolist(), grp["fs"].tolist())
    return proteins, points, read_provenance_header(prot)


def load_integrate_results(out_dir: str) -> pd.DataFrame:
    """Load a prior ``riana integrate`` result for display (no recompute).

    Reads the per-run ``<stem>_riana.txt`` files recorded as the ``stage="integrate"``
    rows of the manifest in *out_dir* and concatenates them into the one frame the
    Integrate tab shows (same as a fresh run's ``pd.concat`` of the per-run frames).
    The manifest is the project locator here — the Integrate tab has no manifest
    field, so its Output dir doubles as "where this project's results live". Rows
    whose output file is missing are skipped.
    """
    from pathlib import Path

    from riana.io.manifest import MANIFEST_FILENAME, read_manifest

    mf = Path(out_dir) / MANIFEST_FILENAME
    if not mf.is_file():
        return pd.DataFrame()
    rows = read_manifest(mf, stage="integrate")
    frames = [pd.read_table(r.output_path, comment="#")
              for r in rows if Path(r.output_path).exists()]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
