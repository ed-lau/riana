"""Core integration entry point — Phase A: faithful 0.9.0 fixed-window port.

This is the *port* of :mod:`riana.riana_integrate`. The contract for Phase A
is strict: same numerical output as 0.9.0 within 1e-3 relative tolerance
(PROJECT_REVIEW.md "Verification" §), so the Week 3 benchmarks have a clean
diff to attribute Phase C peak-detection / Phase D mass-accuracy changes
against. Phase A delivers:

- :func:`integrate_run` — single-fraction entry point. Takes typed inputs
  (``IntegrationConfig`` + ``list[PSMRecord]`` + ``IndexedMzML``) and
  returns the ``pandas.DataFrame`` that the legacy pipeline writes as
  ``<sample>_riana.txt``.
- A consciously thin scope: the multi-fraction loop lives in the caller
  (Week 4 CLI / ``bench_id_path``) so this function stays focused.
- Streaming mzML access via :meth:`IndexedMzML.peaks` — the memory peak on
  a 2 GB mzML drops from "all spectra" to one MS1 at a time
  (PROJECT_REVIEW.md M3 verification §).

Phase C will swap the inline `_extract_per_psm` body for the
``algorithms/peaks.py`` + ``algorithms/baseline.py`` pipeline; Phase D
populates the mass-accuracy fields. Phase A keeps the legacy semantics
bit-near-identical so each later change is its own attributable diff.
"""

from __future__ import annotations

import re
from concurrent import futures
from dataclasses import replace
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.signal

from riana import constants
from riana.config import IntegrationConfig
from riana.exceptions import IntegrationError
from riana.io.mzml import IndexedMzML
from riana.io.percolator import filter_by_q_value, fraction_psms
from riana.records import PSMRecord


# --- Column schema -----------------------------------------------------------

# The PSM metadata columns the legacy ``_riana.txt`` writes, in order.
# Names match the Percolator id_df columns so downstream readers
# (bench_m0_ma_recovery, bench_aa_coefficients) keep parsing.
_PSM_COLUMNS = [
    "file_idx", "scan", "charge",
    "spectrum precursor m/z", "spectrum neutral mass", "peptide mass",
    "percolator score", "percolator q-value", "percolator PEP",
    "distinct matches/spectrum",
    "sequence", "protein id", "flanking aa",
    "concat", "sample", "pep_id", "evidence",
]


def integrate_run(
    config: IntegrationConfig,
    psms: Sequence[PSMRecord],
    mzml: IndexedMzML,
    *,
    file_label: str | None = None,
) -> pd.DataFrame:
    """Integrate one fraction's PSMs against an indexed mzML.

    Args:
        config: pinned configuration (mass tol, isotopomers, threads, ...).
        psms: PSMs already filtered to a single fraction. The caller is
            responsible for that — emitted by :func:`io.percolator.fraction_psms`
            or :func:`io.mztab.read_mztab` after a per-file filter.
        mzml: streaming mzML reader for the fraction's run.
        file_label: written into the ``file`` column. Defaults to the mzML
            basename (without ``.mzML.gz``).

    Returns:
        DataFrame whose columns match the legacy ``<sample>_riana.txt``
        schema: PSM metadata + ``iso{N}`` (or ``mod{m}_iso{N}`` for
        forced mods ≠ 0) area columns + ``file``.
    """
    # Filter by q-value, assign per-fraction pep_id 0..N-1 sorted by scan.
    kept = filter_by_q_value(psms, config.q_value)
    if config.unique_only:
        kept = [p for p in kept if "," not in p.protein_id]
    if not kept:
        raise IntegrationError(
            "No PSMs survive q-value/unique filtering. "
            "Relax --q_value or --unique."
        )
    # The legacy pipeline assigns pep_id per fraction; the caller is already
    # at one fraction so it's safe to reset here.
    kept = sorted(kept, key=lambda p: p.scan)
    kept = [replace(p, pep_id=i) for i, p in enumerate(kept)]

    # Per-concat scan span — for the use_range branch we need min/max scan
    # across every PSM with the same (sequence, charge).
    concat_scans: dict[str, tuple[int, int]] = {}
    for p in kept:
        lo, hi = concat_scans.get(p.concat, (p.scan, p.scan))
        concat_scans[p.concat] = (min(lo, p.scan), max(hi, p.scan))

    forced = tuple(config.forced_mods or (0.0,))
    isos = tuple(config.isotopomers)

    # Threaded per-PSM extraction → list of per-PSM DataFrames (one row per
    # MS1 scan, columns are mod{m}_iso{n}). Same worker pattern as legacy.
    def _do(idx: int) -> pd.DataFrame:
        return _extract_per_psm(
            kept[idx], concat_scans, forced, isos, config, mzml
        )

    with futures.ThreadPoolExecutor(max_workers=config.threads) as ex:
        intensity_dfs = list(ex.map(_do, range(len(kept))))

    # Trapezoidal integration per (PSM × mod × iso) on the per-PSM rt trace.
    iso_cols = [f"mod{_g(m)}_iso{n}" for m in forced for n in isos]
    integrated_rows: list[list] = []
    for idf in intensity_dfs:
        row: list = [idf["pep_id"].iloc[0]]
        for col in iso_cols:
            row.append(float(np.trapezoid(idf[col].to_numpy(), x=idf["rt"].to_numpy())))
        integrated_rows.append(row)

    integrated_df = pd.DataFrame(integrated_rows, columns=["pep_id"] + iso_cols)
    # Legacy strips ``mod0_`` from output column names (the default
    # no-forced-mod case); other forced mods stay namespaced. Reproduces the
    # ac16 baseline header exactly.
    integrated_df.columns = [re.sub(r"^mod0_", "", c) for c in integrated_df.columns]

    # Build the PSM-metadata frame the legacy pipeline emits, then merge.
    psm_df = _psm_metadata_df(kept)
    out = pd.merge(psm_df, integrated_df, on="pep_id", how="left")
    out["file"] = file_label if file_label is not None else _mzml_basename(mzml)
    return out


# --- internals ---------------------------------------------------------------


def _extract_per_psm(
    psm: PSMRecord,
    concat_scans: dict[str, tuple[int, int]],
    forced: tuple[float, ...],
    isos: tuple[int, ...],
    config: IntegrationConfig,
    mzml: IndexedMzML,
) -> pd.DataFrame:
    """Per-MS1-scan isotopomer intensities for one PSM (Phase A: fixed window).

    Matches the legacy ``get_isotopomer_intensity`` numerics:
    ``prec_iso_am = (peptide_mass + z*proton)/z + forced/z + iso*Δ/z``;
    ``±N ppm`` half-width around it; centroids summed within window.
    """
    proton = constants.PROTON_MASS
    iso_added_mass = config.mass_difference

    peptide_mass = psm.peptide_mass
    charge = float(psm.charge)
    peptide_prec = (peptide_mass + charge * proton) / charge

    # Choose the MS1 scans to integrate over. use_range=True is the 0.9.0
    # default: span all PSM scans of the same (sequence, charge), then
    # widen by ±r_time. The searchsorted-then-minus-1 trick lands on the
    # MS1 *preceding* the (possibly MS2) PSM scan — the precursor cycle.
    if config.use_range:
        min_scan, max_scan = concat_scans[psm.concat]
        rt_lo = mzml.rt_idx[np.searchsorted(mzml.scan_idx, min_scan, side="left") - 1]
        rt_hi = mzml.rt_idx[np.searchsorted(mzml.scan_idx, max_scan, side="left") - 1]
        nearby = mzml.scan_idx[
            (mzml.rt_idx - rt_lo > -config.r_time)
            & (mzml.rt_idx - rt_hi < config.r_time)
        ]
    else:
        rt_center = mzml.rt_idx[np.searchsorted(mzml.scan_idx, psm.scan, side="left") - 1]
        nearby = mzml.scan_idx[np.abs(mzml.rt_idx - rt_center) <= config.r_time]

    # Lazy peak fetch + centroid sum per (mod, iso) per MS1.
    rows: list[tuple[str, float, float]] = []
    ppm_tol = float(config.mass_tol_ppm) * 1e-6
    for scan in nearby:
        scan_int = int(scan)
        mz_arr, intens_arr = mzml.peaks(scan_int)
        rt = float(mzml.rt_idx[mzml.scan_idx == scan_int].item())
        for mod in forced:
            prec_shifted = peptide_prec + (mod / charge)
            for iso in isos:
                target = prec_shifted + (iso * iso_added_mass / charge)
                delta = target * ppm_tol
                mask = np.abs(mz_arr - target) <= delta
                summed = float(np.sum(intens_arr[mask])) if mask.any() else 0.0
                rows.append((f"mod{_g(mod)}_iso{iso}", rt, summed))

    if not rows:
        if config.use_range:
            min_scan, max_scan = concat_scans[psm.concat]
            raise IntegrationError(
                f"No intensity profile for peptide {peptide_prec} "
                f"(scans {min_scan}-{max_scan}, iso {list(isos)}). "
                "Try widening --r_time or --mass_tol."
            )
        raise IntegrationError(
            f"No intensity profile for peptide {peptide_prec} "
            f"(iso {list(isos)}). Try widening --r_time or --mass_tol."
        )

    # Pivot (mod_iso × rt) → wide table. Matching legacy: rt rounded to 6 dp
    # so the pivot dedupes by per-MS1 timestamp.
    long = pd.DataFrame(rows, columns=["mod_iso", "rt", "int"])
    long["rt"] = long["rt"].round(6)
    wide = long.pivot(index="rt", columns="mod_iso", values="int")
    iso_cols = [f"mod{_g(m)}_iso{n}" for m in forced for n in isos]
    wide = wide[iso_cols]

    wide["rt"] = wide.index
    wide = wide.reset_index(drop=True)
    wide["pep_id"] = psm.pep_id
    wide["concat"] = psm.concat

    # Smoothing — Phase A preserves the legacy polyorder=1 path so the port
    # is bit-near-identical. Phase B upgrades the default polyorder to 2.
    if config.smoothing is not None and config.smoothing > 0:
        for col in iso_cols:
            if wide[col].sum() > 0:
                wide[col] = scipy.signal.savgol_filter(
                    wide[col].to_numpy(),
                    window_length=config.smoothing,
                    polyorder=1,
                    mode="nearest",
                )
    return wide


def _psm_metadata_df(psms: Sequence[PSMRecord]) -> pd.DataFrame:
    """Build the legacy id_df schema from typed records."""
    rows = []
    for p in psms:
        rows.append({
            "file_idx": p.file_idx,
            "scan": p.scan,
            "charge": p.charge,
            "spectrum precursor m/z": p.precursor_mz,
            "spectrum neutral mass": p.neutral_mass,
            "peptide mass": p.peptide_mass,
            "percolator score": p.percolator_score,
            "percolator q-value": p.percolator_q_value,
            "percolator PEP": p.percolator_pep,
            "distinct matches/spectrum": p.distinct_matches,
            "sequence": p.sequence,
            "protein id": p.protein_id,
            "flanking aa": p.flanking_aa,
            "concat": p.concat,
            "sample": p.sample,
            "pep_id": p.pep_id,
            "evidence": p.evidence,
        })
    return pd.DataFrame(rows, columns=_PSM_COLUMNS)


def _mzml_basename(mzml: IndexedMzML) -> str:
    """``data/.../foo.mzML.gz`` -> ``foo``."""
    name = mzml.path.name
    if name.endswith(".mzML.gz"):
        return name[: -len(".mzML.gz")]
    if name.endswith(".mzML"):
        return name[: -len(".mzML")]
    return mzml.path.stem


def _g(value: float) -> str:
    """``{value:g}`` — strips trailing zeros so ``0.0`` ↔ ``0`` for column names.

    Reproduces the legacy quirk that ``forced_mods=[0]`` (argparse default, an
    int) emits ``mod0_iso0`` while a user-supplied ``-F 0`` (parsed as 0.0)
    would have emitted ``mod0.0_iso0``. We pick the int-shaped form for both
    because the committed v0.9.0 baselines were run with default ``-F``.
    """
    return f"{value:g}"
