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
    ignored_mods: Sequence[float] = (),
) -> list[PSMRecord]:
    """Parse the Percolator id file into typed records (worker side).

    Thin wrapper over :func:`riana.io.percolator.read_percolator` so the GUI can
    do the (potentially slow) parse off the Qt event loop. Same call the CLI
    makes in :func:`riana.cli.integrate`.
    """
    return read_percolator(id_path, sample=sample, ignored_mods=ignored_mods)


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
    return fit_project(config, manifest_path, coeffs)


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
    alt_k: float = 0.025,
    alt_se: float = 0.05,
    threads: int = 1,
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
        alt_k=float(alt_k), alt_se=float(alt_se), threads=int(threads),
        workers=int(workers), phi_limit=float(phi_limit),
        reference_condition=reference_condition,
    )
    return result, result.attrs.get("protein_points", {})
