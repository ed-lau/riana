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
from riana.exceptions import DataError, IntegrationError
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
    "sequence", "protein id", "mod sites", "flanking aa",
    "concat", "sample", "pep_id", "evidence",
]


@dataclass(frozen=True, slots=True)
class ScanRtCheck:
    """Result of reconciling a run's PSM ``spectra_ref`` scans against mzML RT.

    Produced by :func:`check_scan_rt_consistency`. ``n_checked`` counts the PSMs
    that carried a usable mzTab ``retention_time`` (> 0); the Percolator path
    supplies none, so the guard no-ops there (``n_checked == 0`` ⇒ :attr:`ok`).
    """

    #: PSMs with a usable reported RT that were reconciled.
    n_checked: int
    #: Median |mzML scan-RT − reported RT| over the checked PSMs, in RT minutes.
    median_offset_min: float
    #: Fraction of checked PSMs within :attr:`tol_min` of their reported RT.
    frac_within_tol: float
    #: The tolerance the median was gated against (RT minutes).
    tol_min: float

    @property
    def ok(self) -> bool:
        """True when there was nothing to check or the median is within tol."""
        return self.n_checked == 0 or self.median_offset_min <= self.tol_min


def check_scan_rt_consistency(
    psms: Sequence[PSMRecord],
    mzml: IndexedMzML,
    tol_min: float,
) -> ScanRtCheck:
    """Reconcile each PSM's ``spectra_ref`` scan → mzML MS1 RT vs its mzTab RT.

    On the quantms mzTab path every PSM carries the ``retention_time`` the
    search/quant pipeline reported (seconds). Looking the PSM's ``scan`` up in
    *this* mzML's MS1 index must land near that RT; a large per-run **median**
    offset means the ``spectra_ref`` scans do not belong to this mzML — the
    quantms filename-prefix scan-scramble (mzML basenames that are prefixes of
    one another) or a wrong mzML↔mzTab pairing (PROJECT_REVIEW Track A intake
    guard). The failure was previously silent.

    The **median** is the gate (not per-PSM): it is robust to the ~10% of PSMs
    that legitimately mismatch and to the run-dependent ProteomicsLFQ alignment
    offset (≤~0.9 min measured on real output), while a scrambled run sits tens
    of minutes off — a ~25× separation, so the exact tolerance barely matters.

    PSMs without a usable ``retention_time`` (≤ 0 — the Percolator path) are
    skipped; with none, the guard no-ops (``n_checked == 0``, :attr:`ScanRtCheck.ok`).
    The caller decides policy (error vs. log).
    """
    if not psms:
        return ScanRtCheck(0, float("nan"), float("nan"), tol_min)
    scans = np.array([p.scan for p in psms], dtype=np.int64)
    rep_rt_min = np.array(
        [p.retention_time / 60.0 for p in psms], dtype=np.float64
    )  # mzTab RT is seconds; mzML rt_idx is minutes.
    valid = rep_rt_min > 0
    n = int(valid.sum())
    if n == 0:
        return ScanRtCheck(0, float("nan"), float("nan"), tol_min)
    offset = np.abs(mzml.rt_for_scans(scans[valid]) - rep_rt_min[valid])
    return ScanRtCheck(
        n_checked=n,
        median_offset_min=float(np.median(offset)),
        frac_within_tol=float(np.mean(offset <= tol_min)),
        tol_min=tol_min,
    )


def resolve_rt_anchored_scans(
    psms: Sequence[PSMRecord],
    mzml: IndexedMzML,
    *,
    bounds_tol_min: float = 1.0,
    file_label: str | None = None,
) -> list[PSMRecord]:
    """Map RT-anchored (DIA) PSMs to the nearest MS1 scan in *mzml*.

    The DIA intake (:mod:`riana.io.diann`) has no MS2 scan to anchor on, so it
    emits ``scan = -1`` and carries DIA-NN's apex into ``retention_time``
    (seconds). Here each such PSM's RT is resolved to the nearest MS1 scan in
    *this* mzML, so the downstream scan-based extraction runs unchanged. PSMs
    that already carry a real scan (the DDA path) are returned untouched — this
    is a no-op there.

    The resolved scan is at most one MS1 cycle off the true apex (the
    ``searchsorted - 1`` precursor-cycle lookup), which the apex finder
    (``apex_search_half_width`` ≫ a cycle) absorbs by re-centring on the MS1
    apex. As a sanity check we count PSMs whose reported RT falls outside the
    mzML's MS1 RT span (± ``bounds_tol_min``): a large fraction means the report
    was paired with the wrong mzML, and we raise.

    This replaces the scan↔RT scramble guard for DIA, which is circular there
    (the scan is *derived* from the RT, so it always reconciles).
    """
    label = file_label or _mzml_basename(mzml)
    rt_idx = mzml.rt_idx
    scan_idx = mzml.scan_idx
    if len(rt_idx) == 0:
        raise DataError(f"{label}: mzML has no MS1 scans to anchor DIA RTs to.")

    out = list(psms)
    todo = [i for i, p in enumerate(out) if p.scan < 0 and p.retention_time > 0]
    if not todo:
        return out

    rt_min = np.array([out[i].retention_time / 60.0 for i in todo], dtype=np.float64)
    # rt_idx is ascending (MS1 RTs increase) → nearest via searchsorted.
    pos = np.clip(np.searchsorted(rt_idx, rt_min), 1, len(rt_idx) - 1)
    take_left = (rt_min - rt_idx[pos - 1]) <= (rt_idx[pos] - rt_min)
    nearest = np.where(take_left, pos - 1, pos)
    resolved_scans = scan_idx[nearest]

    rt_lo, rt_hi = float(rt_idx.min()), float(rt_idx.max())
    oob = int(np.sum((rt_min < rt_lo - bounds_tol_min) | (rt_min > rt_hi + bounds_tol_min)))
    frac_oob = oob / len(todo)
    if frac_oob > 0.5:
        raise DataError(
            f"{label}: {frac_oob:.0%} of DIA PSM apex RTs fall outside this "
            f"mzML's MS1 RT span [{rt_lo:.1f}, {rt_hi:.1f}] min (± "
            f"{bounds_tol_min:.1f}). The DIA-NN report is most likely paired "
            "with the wrong mzML — check the SDRF `comment[data file]` ↔ mzML "
            "names."
        )
    if oob:
        _LOGGER.warning(
            "%s: %d/%d DIA PSM apex RTs outside the mzML MS1 span "
            "[%.1f, %.1f] min — integrated at the nearest edge scan.",
            label, oob, len(todo), rt_lo, rt_hi,
        )

    for k, i in enumerate(todo):
        out[i] = replace(out[i], scan=int(resolved_scans[k]))
    _LOGGER.info(
        "%s: resolved %d DIA PSMs to MS1 scans by apex RT.", label, len(todo)
    )
    return out


def integrate_run(
    config: IntegrationConfig,
    psms: Sequence[PSMRecord],
    mzml: IndexedMzML,
    *,
    file_label: str | None = None,
) -> pd.DataFrame:
    """Integrate one fraction's PSMs against an indexed mzML.

    Args:
        config: pinned configuration (mass tol, isotopomers, ...).
        psms: PSMs already filtered to a single fraction. The caller is
            responsible for that — emitted by :func:`io.percolator.fraction_psms`
            or :func:`io.mztab.read_mztab` after a per-file filter.
        mzml: streaming mzML reader for the fraction's run.
        file_label: written into the ``file`` column. Defaults to the mzML
            basename (without ``.mzML.gz``).

    Returns:
        DataFrame whose columns match the legacy ``<sample>_riana.txt``
        schema: PSM metadata + ``iso{N}`` area columns + ``file``.
    """
    # Filter by q-value, assign per-fraction pep_id 0..N-1 sorted by scan.
    if mzml.ms1_centroid is False:
        _LOGGER.warning(
            "%s: MS1 is PROFILE, not centroid. The ±%d ppm window assumes "
            "centroid data (one line per isotopomer ~ mass accuracy); on profile "
            "data it under-captures the peak (FWHM ~30 ppm at 700 m/z). Centroid "
            "the mzML, or widen --mass_tol to the profile peak width.",
            file_label or _mzml_basename(mzml), config.mass_tol_ppm,
        )

    # DIA (RT-anchored) intake: io.diann emits scan=-1 because DIA has no MS2
    # scan, carrying DIA-NN's apex in retention_time. Resolve each to the
    # nearest MS1 scan in THIS mzML (+ an RT-in-bounds check) so the scan-based
    # extraction below runs unchanged; the apex finder then re-centres on the
    # true MS1 apex. No-op on the DDA path (every PSM already has a real scan).
    rt_anchored = any(p.scan < 0 for p in psms)
    if rt_anchored:
        psms = resolve_rt_anchored_scans(psms, mzml, file_label=file_label)

    # Intake scan↔RT guard (Track A): before integrating, verify the mzTab
    # spectra_ref scans actually index THIS mzML. Runs on the full PSM set
    # (best statistics) and no-ops on the Percolator path (no retention_time).
    # Skipped for DIA: the scan was just *derived* from the RT, so the check is
    # circular (resolve_rt_anchored_scans does the RT-in-bounds check instead).
    if config.check_scan_rt and not rt_anchored:
        label = file_label or _mzml_basename(mzml)
        check = check_scan_rt_consistency(psms, mzml, config.scan_rt_tol_min)
        if check.n_checked > 0 and not check.ok:
            raise DataError(
                f"{label}: scan↔RT reconciliation FAILED — median offset "
                f"{check.median_offset_min:.2f} min over {check.n_checked} PSMs "
                f"exceeds tol {check.tol_min:.1f} min "
                f"({check.frac_within_tol:.0%} within tol). The mzTab spectra_ref "
                "scans do not line up with this mzML's retention times: most "
                "likely the quantms filename-prefix scan-scramble (zero-pad / "
                "de-prefix the mzML basenames before the quantms run) or a wrong "
                "mzML↔mzTab pairing. Override with --no-rt-check only if this run "
                "is knowingly correct."
            )
        if check.n_checked > 0:
            _LOGGER.info(
                "%s: scan↔RT reconciled — median %.2f min, %.0f%% within "
                "%.1f min (n=%d).",
                label, check.median_offset_min, 100 * check.frac_within_tol,
                check.tol_min, check.n_checked,
            )

    # Integrate ALL peptides, shared included — protein attribution (unique /
    # isoform parsimony) is a summarize-time decision in `riana rollup`, not an
    # integrate-time filter (a shared peptide is still a valid per-peptide
    # measurement; only its protein roll-up is ambiguous).
    kept = filter_by_q_value(psms, config.q_value)
    if not kept:
        raise IntegrationError(
            "No PSMs survive q-value filtering. Relax --q_value."
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

    # Per-concat apex anchor — the scan of the **best-q (most confident) PSM**
    # for each peptide-charge. When a peptide is identified at several scans, the
    # apex search (and the use_range=False extraction centre) keys on this single
    # anchor rather than each PSM's own scan, so every row of the concat locates
    # the same peak near the confident ID instead of drifting per-PSM. Paired with
    # ``apex_search_half_width`` it stops the apex roaming to a co-eluting isobar.
    concat_anchor: dict[str, int] = {}
    _best_q: dict[str, float] = {}
    for p in kept:
        if p.concat not in _best_q or p.percolator_q_value < _best_q[p.concat]:
            _best_q[p.concat] = p.percolator_q_value
            concat_anchor[p.concat] = p.scan

    isos = tuple(config.isotopomers)

    # Per-PSM extraction → list of per-PSM (DataFrame, mass_accuracy_dict).
    # Serial: the hot path (pyteomics mzML parse + IsoSpec) is GIL-bound, so
    # threading it gave no speedup (in fact slower); cross-run parallelism is the
    # ProcessPool lever (`integrate --workers`), see core/pipeline.
    def _do(idx: int):
        return _extract_per_psm(
            kept[idx], concat_scans, isos, config, mzml,
            anchor_scan=concat_anchor[kept[idx].concat],
        )

    extract_results = [_do(i) for i in range(len(kept))]
    intensity_dfs = [r[0] for r in extract_results]
    mass_accuracies = [r[1] for r in extract_results]

    # Trapezoidal integration per (PSM × mod × iso) on the per-PSM rt trace.
    # Phase C: with peak_rt="apex" / integration_half_width="auto" the window
    # comes from the apex finder on the iso0 XIC (+ a co-elution check vs iso1);
    # baseline-subtraction uses baseline_method (default "none"; the narrow
    # window excludes background). Per-peptide fallback to fixed on failure.
    iso_cols = [f"iso{n}" for n in isos]
    mz_cols = [f"{c}_obs_mz" for c in iso_cols]
    ppm_cols = [f"{c}_ppm_error" for c in iso_cols]
    integrated_rows: list[list] = []
    n_detected = n_fallback = 0
    iso0_col = "iso0"
    iso1_col = "iso1" if 1 in isos else None
    iso0_ppm_errors: list[float] = []
    for idf, ma, psm in zip(intensity_dfs, mass_accuracies, kept):
        row: list = [idf["pep_id"].iloc[0]]
        rt_arr = idf["rt"].to_numpy(dtype=np.float64)

        boundary = _peak_boundary(
            idf, psm, mzml, config, rt_arr, iso0_col, iso1_col,
            anchor_scan=concat_anchor[psm.concat],
        )
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
    #: iso index → the extracted-ion chromatogram for that isotopomer.
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
    isos = tuple(config.isotopomers)
    concat_scans = {psm.concat: scan_span or (psm.scan, psm.scan)}

    idf, _ma = _extract_per_psm(psm, concat_scans, isos, config, mzml)
    rt_arr = idf["rt"].to_numpy(dtype=np.float64)

    iso0_col = "iso0"
    iso1_col = "iso1" if 1 in isos else None
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
    chromatograms: dict[int, Chromatogram] = {}
    for iso in isos:
        col = f"iso{iso}"
        intensity = idf[col].to_numpy(dtype=np.float64)
        target = (
            peptide_prec
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
    isos: tuple[int, ...],
    config: IntegrationConfig,
    mzml: IndexedMzML,
    *,
    anchor_scan: int | None = None,
) -> tuple[pd.DataFrame, dict[str, tuple[float | None, float | None]]]:
    """Per-MS1-scan isotopomer intensities + mass-accuracy for one PSM.

    Returns ``(wide_df, mass_accuracy)`` where ``mass_accuracy`` maps each
    ``iso{n}`` column to ``(obs_mz, ppm_error)``. Both are ``None``
    when no centroid matched in the window (an honest "not observed"
    rather than a fabricated zero — same posture
    :class:`riana.records.IsotopomerPeak` takes).

    Matches the legacy ``get_isotopomer_intensity`` numerics:
    ``prec_iso_am = (peptide_mass + z*proton)/z + iso*Δ/z``;
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
        center_scan = anchor_scan if anchor_scan is not None else psm.scan
        rt_center = mzml.rt_idx[np.searchsorted(mzml.scan_idx, center_scan, side="left") - 1]
        nearby = mzml.scan_idx[np.abs(mzml.rt_idx - rt_center) <= config.extraction_half_width]

    # Lazy peak fetch + centroid sum per (mod, iso) per MS1. As a Phase D
    # addition we also accumulate intensity-weighted observed m/z per
    # (mod, iso) across the integration window so the writer can emit the
    # iso{N}_obs_mz / iso{N}_ppm_error mass-accuracy columns.
    rows: list[tuple[str, float, float]] = []
    iso_cols_local = [f"iso{n}" for n in isos]
    targets: dict[str, float] = {}
    mz_weighted: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    intens_total: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    ppm_tol = float(config.mass_tol_ppm) * 1e-6
    for scan in nearby:
        scan_int = int(scan)
        mz_arr, intens_arr = mzml.peaks(scan_int)
        rt = float(mzml.rt_idx[mzml.scan_idx == scan_int].item())
        for iso in isos:
            col = f"iso{iso}"
            target = peptide_prec + (iso * iso_added_mass / charge)
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
    iso_cols = [f"iso{n}" for n in isos]
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
    *,
    anchor_scan: int | None = None,
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
    prior_scan = anchor_scan if anchor_scan is not None else psm.scan
    psm_rt = float(
        mzml.rt_idx[np.searchsorted(mzml.scan_idx, prior_scan, side="left") - 1]
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
        base = iso0_col[:-1]  # "iso0" -> "iso"
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
            "mod sites": p.mod_sites,
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
