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

import logging
import re
from concurrent import futures
from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.signal

from riana import constants
from riana.algorithms import baseline as ba
from riana.algorithms import calibration as cal
from riana.algorithms import peaks as pk
from riana.algorithms import smoothing as sm
from riana.config import IntegrationConfig
from riana.exceptions import IntegrationError
from riana.io.mzml import IndexedMzML
from riana.io.percolator import filter_by_q_value, fraction_psms
from riana.records import Chromatogram, PSMRecord

_LOGGER = logging.getLogger(__name__)


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

    # Threaded per-PSM extraction → list of per-PSM (DataFrame,
    # mass_accuracy_dict). Same worker pattern as legacy.
    def _do(idx: int):
        return _extract_per_psm(
            kept[idx], concat_scans, forced, isos, config, mzml
        )

    with futures.ThreadPoolExecutor(max_workers=config.threads) as ex:
        extract_results = list(ex.map(_do, range(len(kept))))
    intensity_dfs = [r[0] for r in extract_results]
    mass_accuracies = [r[1] for r in extract_results]

    # Trapezoidal integration per (PSM × mod × iso) on the per-PSM rt trace.
    # Phase C: with peak_rt="apex" / integration_half_width="auto" the window
    # comes from the apex finder on the iso0 XIC (+ a co-elution check vs iso1);
    # baseline-subtraction uses baseline_method (default "none"; the narrow
    # window excludes background). Per-peptide fallback to fixed on failure.
    iso_cols = [f"mod{_g(m)}_iso{n}" for m in forced for n in isos]
    mz_cols = [f"{c}_obs_mz" for c in iso_cols]
    ppm_cols = [f"{c}_ppm_error" for c in iso_cols]
    integrated_rows: list[list] = []
    n_detected = n_fallback = 0
    iso0_col = f"mod{_g(forced[0])}_iso0"
    iso1_col = f"mod{_g(forced[0])}_iso1" if 1 in isos else None
    iso0_ppm_errors: list[float] = []
    for idf, ma, psm in zip(intensity_dfs, mass_accuracies, kept):
        row: list = [idf["pep_id"].iloc[0]]
        rt_arr = idf["rt"].to_numpy(dtype=np.float64)

        boundary = _peak_boundary(idf, psm, mzml, config, rt_arr, iso0_col, iso1_col)
        if boundary is not None:
            n_detected += 1
        else:
            n_fallback += 1

        for col in iso_cols:
            trace = idf[col].to_numpy(dtype=np.float64)
            if boundary is not None:
                lo, hi = boundary.lo, boundary.hi
                # Baseline subtraction applies to ALL isotopomers, because
                # the noise floor is approximately constant across iso0..iso5
                # (same instrument, similar m/z window). Subtraction is tiny
                # for iso0 (noise << apex) and meaningful for iso{N>0} where
                # noise is comparable to apex — without it iso5 would be
                # over-estimated by integrating the noise floor across the
                # window. Per-trace noise-floor (default) is robust to
                # endpoint spikes (the failure mode of the first Phase C cut).
                corrected = _baseline_corrected_slice(
                    rt_arr, trace, lo, hi, config.baseline_method
                )
                # Clamp the integrated area (not per-scan): a single noisy
                # boundary scan can produce small negatives in `corrected`,
                # but the *integrated* area still reflects the underlying
                # peak — clipping per-scan amputates whole isotopomers.
                area = max(0.0, float(np.trapezoid(corrected, x=rt_arr[lo : hi + 1])))
            else:
                area = float(np.trapezoid(trace, x=rt_arr))
            row.append(area)
        # Append mass-accuracy columns in the same order as iso_cols.
        for col in iso_cols:
            obs_mz, _ppm = ma.get(col, (None, None))
            row.append(np.nan if obs_mz is None else obs_mz)
        for col in iso_cols:
            _obs_mz, ppm_error = ma.get(col, (None, None))
            row.append(np.nan if ppm_error is None else ppm_error)
        # iso0's ppm_error is the strongest per-PSM drift indicator (M+1 etc.
        # can be contaminated by other peptides at high D2O); use it for the
        # per-fraction drift footer.
        iso0_ppm = ma.get(iso0_col, (None, None))[1]
        if iso0_ppm is not None and not np.isnan(iso0_ppm):
            iso0_ppm_errors.append(float(iso0_ppm))
        integrated_rows.append(row)

    if config.peak_rt in ("apex", "consensus") or config.integration_half_width == "auto":
        _LOGGER.info(
            "integrate_run: %d PSMs detected, %d fell back to fixed-window "
            "(detection or co-elution failure).",
            n_detected, n_fallback,
        )

    integrated_df = pd.DataFrame(
        integrated_rows,
        columns=["pep_id"] + iso_cols + mz_cols + ppm_cols,
    )
    # Legacy strips ``mod0_`` from output column names (the default
    # no-forced-mod case); other forced mods stay namespaced. Reproduces the
    # ac16 baseline header exactly. The same stripping applies to the new
    # mass-accuracy columns (``mod0_isoN_obs_mz`` -> ``isoN_obs_mz``).
    integrated_df.columns = [re.sub(r"^mod0_", "", c) for c in integrated_df.columns]

    # Per-fraction drift summary — log a warning if the median ppm error
    # exceeds --ppm-alert.
    drift = cal.drift_summary(np.asarray(iso0_ppm_errors, dtype=np.float64))
    integrate_run.last_drift = drift  # type: ignore[attr-defined]
    if drift.n > 0 and abs(drift.median_ppm) > config.ppm_alert:
        _LOGGER.warning(
            "integrate_run: per-fraction iso0 ppm-error median %+0.2f ppm "
            "(MAD %.2f, n=%d) exceeds --ppm-alert %+0.1f ppm. Suggested "
            "calibration shift: %+0.2f ppm.",
            drift.median_ppm, drift.mad_ppm, drift.n,
            config.ppm_alert, drift.suggested_shift_ppm,
        )

    # Build the PSM-metadata frame the legacy pipeline emits, then merge.
    psm_df = _psm_metadata_df(kept)
    out = pd.merge(psm_df, integrated_df, on="pep_id", how="left")
    out["file"] = file_label if file_label is not None else _mzml_basename(mzml)
    # Stash the drift summary on the DataFrame as attrs so callers/writers can
    # emit the footer without re-computing.
    out.attrs["drift_summary"] = drift
    return out


@dataclass(frozen=True, slots=True)
class PeptideTrace:
    """Per-isotopomer XICs for one peptide-charge plus the integrated window.

    The data behind the M4 Phase 2 GUI chromatogram view. Built by
    :func:`extract_peptide_trace`; Qt-free and picklable (only stdlib + the
    frozen :class:`riana.records.Chromatogram`) so it can cross a
    ``ProcessPoolExecutor`` boundary back to the GUI process.
    """

    #: ``sequence_charge`` identifier of the peptide-charge.
    concat: str
    #: iso index → the extracted-ion chromatogram for that isotopomer (the
    #: first forced-mod cluster — ``mod0`` in the no-SILAC default).
    chromatograms: dict[int, Chromatogram]
    #: ``(lo_rt, hi_rt)`` of the integrated window in RT minutes, or ``None``
    #: when the whole extraction was integrated (the ms2 / fixed-window path).
    window: tuple[float, float] | None


def extract_peptide_trace(
    config: IntegrationConfig,
    psm: PSMRecord,
    mzml: IndexedMzML,
    *,
    scan_span: tuple[int, int] | None = None,
) -> PeptideTrace:
    """Extract one peptide-charge's per-isotopomer XICs + integration window.

    Reuses the same extraction (:func:`_extract_per_psm`) and boundary detection
    (:func:`_peak_boundary`) the integrator runs, so the trace and shaded window
    the GUI shows match what was integrated for that PSM. Gives the otherwise
    orphaned :class:`riana.records.Chromatogram` record a producer.

    Args:
        config: the same :class:`~riana.config.IntegrationConfig` used to run.
        psm: the selected peptide-spectrum match.
        mzml: streaming reader for the fraction's run.
        scan_span: the ``(min_scan, max_scan)`` the integrator spans for this
            peptide-charge across the fraction (``config.use_range``). Pass the
            fraction's real span so the trace matches the integration exactly;
            defaults to the PSM's own scan.

    Returns:
        :class:`PeptideTrace` — chromatograms keyed by iso index + the window.
    """
    forced = tuple(config.forced_mods or (0.0,))
    isos = tuple(config.isotopomers)
    concat_scans = {psm.concat: scan_span or (psm.scan, psm.scan)}

    idf, _ma = _extract_per_psm(psm, concat_scans, forced, isos, config, mzml)
    rt_arr = idf["rt"].to_numpy(dtype=np.float64)

    iso0_col = f"mod{_g(forced[0])}_iso0"
    iso1_col = f"mod{_g(forced[0])}_iso1" if 1 in isos else None
    boundary = _peak_boundary(idf, psm, mzml, config, rt_arr, iso0_col, iso1_col)
    window = (
        (float(rt_arr[boundary.lo]), float(rt_arr[boundary.hi]))
        if boundary is not None
        else None
    )

    # Recover the scan per (rounded) rt — _extract_per_psm pivots on rt rounded
    # to 6 dp, and those rt values originate from mzml.rt_idx, so the map is exact.
    rt_to_scan = {
        round(float(r), 6): int(s) for s, r in zip(mzml.scan_idx, mzml.rt_idx)
    }
    scans_for_rt = tuple(rt_to_scan.get(round(float(r), 6), -1) for r in rt_arr)

    proton = constants.PROTON_MASS
    peptide_prec = (psm.peptide_mass + psm.charge * proton) / psm.charge
    mod0 = forced[0]
    chromatograms: dict[int, Chromatogram] = {}
    for iso in isos:
        col = f"mod{_g(mod0)}_iso{iso}"
        intensity = idf[col].to_numpy(dtype=np.float64)
        target = (
            peptide_prec
            + (mod0 / psm.charge)
            + (iso * config.mass_difference / psm.charge)
        )
        chromatograms[iso] = Chromatogram(
            isotopomer=iso,
            target_mz=float(target),
            mass_tol_ppm=float(config.mass_tol_ppm),
            scans=scans_for_rt,
            rt=tuple(float(r) for r in rt_arr),
            intensity=tuple(float(v) for v in intensity),
        )
    return PeptideTrace(
        concat=psm.concat, chromatograms=chromatograms, window=window
    )


# --- internals ---------------------------------------------------------------


def _extract_per_psm(
    psm: PSMRecord,
    concat_scans: dict[str, tuple[int, int]],
    forced: tuple[float, ...],
    isos: tuple[int, ...],
    config: IntegrationConfig,
    mzml: IndexedMzML,
) -> tuple[pd.DataFrame, dict[str, tuple[float | None, float | None]]]:
    """Per-MS1-scan isotopomer intensities + mass-accuracy for one PSM.

    Returns ``(wide_df, mass_accuracy)`` where ``mass_accuracy`` maps each
    ``mod{m}_iso{n}`` column to ``(obs_mz, ppm_error)``. Both are ``None``
    when no centroid matched in the window (an honest "not observed"
    rather than a fabricated zero — same posture
    :class:`riana.records.IsotopomerPeak` takes).

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
    # widen by ±extraction_half_width. The searchsorted-then-minus-1 lands on
    # MS1 *preceding* the (possibly MS2) PSM scan — the precursor cycle.
    if config.use_range:
        min_scan, max_scan = concat_scans[psm.concat]
        rt_lo = mzml.rt_idx[np.searchsorted(mzml.scan_idx, min_scan, side="left") - 1]
        rt_hi = mzml.rt_idx[np.searchsorted(mzml.scan_idx, max_scan, side="left") - 1]
        nearby = mzml.scan_idx[
            (mzml.rt_idx - rt_lo > -config.extraction_half_width)
            & (mzml.rt_idx - rt_hi < config.extraction_half_width)
        ]
    else:
        rt_center = mzml.rt_idx[np.searchsorted(mzml.scan_idx, psm.scan, side="left") - 1]
        nearby = mzml.scan_idx[np.abs(mzml.rt_idx - rt_center) <= config.extraction_half_width]

    # Lazy peak fetch + centroid sum per (mod, iso) per MS1. As a Phase D
    # addition we also accumulate intensity-weighted observed m/z per
    # (mod, iso) across the integration window so the writer can emit the
    # iso{N}_obs_mz / iso{N}_ppm_error mass-accuracy columns.
    rows: list[tuple[str, float, float]] = []
    iso_cols_local = [f"mod{_g(m)}_iso{n}" for m in forced for n in isos]
    targets: dict[str, float] = {}
    mz_weighted: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    intens_total: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    ppm_tol = float(config.mass_tol_ppm) * 1e-6
    for scan in nearby:
        scan_int = int(scan)
        mz_arr, intens_arr = mzml.peaks(scan_int)
        rt = float(mzml.rt_idx[mzml.scan_idx == scan_int].item())
        for mod in forced:
            prec_shifted = peptide_prec + (mod / charge)
            for iso in isos:
                col = f"mod{_g(mod)}_iso{iso}"
                target = prec_shifted + (iso * iso_added_mass / charge)
                targets[col] = target
                delta = target * ppm_tol
                mask = np.abs(mz_arr - target) <= delta
                if mask.any():
                    matched_mz = mz_arr[mask]
                    matched_i = intens_arr[mask]
                    summed = float(matched_i.sum())
                    mz_weighted[col] += float((matched_mz * matched_i).sum())
                    intens_total[col] += summed
                else:
                    summed = 0.0
                rows.append((col, rt, summed))

    if not rows:
        if config.use_range:
            min_scan, max_scan = concat_scans[psm.concat]
            raise IntegrationError(
                f"No intensity profile for peptide {peptide_prec} "
                f"(scans {min_scan}-{max_scan}, iso {list(isos)}). "
                "Try widening --extraction_half_width or --mass_tol."
            )
        raise IntegrationError(
            f"No intensity profile for peptide {peptide_prec} "
            f"(iso {list(isos)}). Try widening --extraction_half_width or --mass_tol."
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

    # Smoothing. When the user asks for fixed-window + smoothing, polyorder=1
    # is preserved for bit-near-identity with the 0.9.0 baseline (the only
    # gate that can catch a regression in the legacy semantics). For the
    # detected/Phase-C science path we use the §2c fix polyorder=2 via the
    # algorithms.smoothing wrapper.
    if config.smoothing is not None and config.smoothing > 0:
        if config.peak_rt == "ms2" and config.integration_half_width != "auto":
            for col in iso_cols:
                if wide[col].sum() > 0:
                    wide[col] = scipy.signal.savgol_filter(
                        wide[col].to_numpy(),
                        window_length=config.smoothing,
                        polyorder=1,
                        mode="nearest",
                    )
        else:
            for col in iso_cols:
                wide[col] = sm.savgol(
                    wide[col].to_numpy(),
                    window=config.smoothing,
                    polyorder=config.smoothing_polyorder,
                )

    # Mass-accuracy aggregation: intensity-weighted observed m/z over the
    # whole integration window, per (mod, iso). None when no centroid landed
    # in window — honest "not observed" over a fabricated zero.
    mass_accuracy: dict[str, tuple[float | None, float | None]] = {}
    for col in iso_cols_local:
        weight = intens_total[col]
        target = targets[col]
        if weight > 0:
            obs_mz = mz_weighted[col] / weight
            ppm_error = (obs_mz - target) / target * 1e6
            mass_accuracy[col] = (obs_mz, ppm_error)
        else:
            mass_accuracy[col] = (None, None)
    return wide, mass_accuracy


def _peak_boundary(
    idf: pd.DataFrame,
    psm: PSMRecord,
    mzml: IndexedMzML,
    config: IntegrationConfig,
    rt_arr: np.ndarray,
    iso0_col: str,
    iso1_col: str | None,
) -> "pk.PeakBoundary | None":
    """Locate the integration [lo, hi] per ``peak_rt`` / ``integration_half_width``.

    ``None`` ⇒ integrate the whole extraction (ms2 + fixed width = 0.9.0).
    Otherwise: apex ± fixed width (``peak_rt="apex"``), or detected boundaries
    with an iso1 co-elution gate (``integration_half_width="auto"``).
    """
    auto = config.integration_half_width == "auto"
    if config.peak_rt == "ms2" and not auto:
        return None  # 0.9.0 fixed: integrate the whole ±extraction_half_width
    iso0_trace = idf[iso0_col].to_numpy(dtype=np.float64)
    if iso0_trace.sum() <= 0:
        return None
    psm_rt = float(
        mzml.rt_idx[np.searchsorted(mzml.scan_idx, psm.scan, side="left") - 1]
    )
    if auto:
        iso0_b = pk.detect_peak(
            rt_arr, iso0_trace, scan_prior_rt=psm_rt,
            rel_height=config.width_rel_height, prominence_k=config.prominence_k,
        )
        if iso0_b is None:
            return None
        if iso1_col is not None:
            iso1_trace = idf[iso1_col].to_numpy(dtype=np.float64)
            iso1_b = (
                pk.detect_peak(
                    rt_arr, iso1_trace, scan_prior_rt=psm_rt,
                    rel_height=config.width_rel_height,
                    prominence_k=config.prominence_k,
                )
                if iso1_trace.sum() > 0
                else None
            )
            if not pk.coelution_ok(iso0_b, iso1_b, rt_arr):
                return None
        return iso0_b
    # peak_rt in {"apex","consensus"}: fixed-width window around a detected apex.
    if config.peak_rt == "consensus":
        # Median apex over m0..m{n-1} (co-elution consensus): labelling-
        # independent and rejects a contaminated channel regardless of intensity.
        base = iso0_col[:-1]  # "mod0_iso0" -> "mod0_iso"
        cons_cols = [f"{base}{k}" for k in range(config.apex_n_consensus)
                     if f"{base}{k}" in idf.columns]
        traces = [idf[c].to_numpy(dtype=np.float64) for c in cons_cols]
        res = pk.consensus_apex(
            rt_arr, traces, scan_prior_rt=psm_rt,
            prominence_k=config.prominence_k,
            apex_search_half_width=config.apex_search_half_width,
            selection=config.apex_selection,
        )
        apex = res[0] if res is not None else None
    else:  # "apex"
        apex = pk.find_apex(
            rt_arr, iso0_trace, scan_prior_rt=psm_rt,
            prominence_k=config.prominence_k,
            apex_search_half_width=config.apex_search_half_width,
            selection=config.apex_selection,
        )
    if apex is None:
        return None
    half = float(config.integration_half_width)
    apex_rt = float(rt_arr[apex])
    win = np.where(np.abs(rt_arr - apex_rt) <= half)[0]
    return pk.PeakBoundary(
        apex_idx=apex, lo=int(win[0]), hi=int(win[-1]),
        apex_intensity=float(iso0_trace[apex]), prominence=float("nan"),
    )


def _baseline_corrected_slice(
    rt_arr: np.ndarray,
    trace: np.ndarray,
    lo: int,
    hi: int,
    method: str,
) -> np.ndarray:
    """Return ``trace[lo:hi+1] - baseline`` according to ``method``.

    For ``"none"`` we skip subtraction (just slice). ``"noise_floor"`` is a
    flat p10-of-(full-trace) constant; ``"snip"`` / ``"asls"`` route through
    :mod:`algorithms.baseline`. (Skyline-style ``"linear"`` was discarded —
    see PROJECT_REVIEW; it over-subtracts on narrow on-peak boundaries.)

    The corrected slice can contain small negatives when a single
    boundary scan happens to sit above the trace just inside. Those
    integrate to small negative or near-zero areas — the caller clips the
    *final* integrated area at zero. Clipping per-scan (the original Phase
    C behavior) zeroed entire isotopomers from one noisy boundary scan.
    """
    sliced = trace[lo : hi + 1]
    if method == "none":
        return sliced
    if method == "noise_floor":
        bl = ba.noise_floor(trace)
    elif method == "snip":
        bl = ba.snip(trace)
    elif method == "asls":
        bl = ba.asls(trace)
    else:  # defensive; __post_init__ already gated
        return sliced
    return sliced - bl[lo : hi + 1]


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
