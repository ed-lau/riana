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
class ScanPrecursorCheck:
    """Result of reconciling a run's PSM scans against THIS mzML's precursors.

    Produced by :func:`check_scan_precursor_consistency`. A real file/search
    mismatch shows up as :attr:`frac_matched` near 0 (the scans are absent or
    point at different precursors); RT alignment leaves it ~1.0.
    """

    #: Sampled PSMs with a usable scan + reported precursor m/z.
    n_checked: int
    #: Median |reported − mzML selected-ion m/z| over matched scans, in ppm.
    median_ppm: float
    #: Fraction of sampled scans whose mzML precursor m/z matches the report.
    frac_matched: float
    #: Minimum acceptable :attr:`frac_matched` (below ⇒ mismatch).
    min_frac: float

    @property
    def ok(self) -> bool:
        """True when there was nothing to check or enough scans matched."""
        return self.n_checked == 0 or self.frac_matched >= self.min_frac


# Scan↔precursor guard knobs. A legitimate file matches ~100% within a few ppm
# (the mzTab ``exp_mass_to_charge`` IS the mzML selected-ion m/z); a wrong file /
# scan-scramble matches ~0%, so the 0.5 floor separates them with a wide margin.
# 20 ppm absorbs any precursor-refinement rounding. Sampling ~300 PSMs keeps the
# guard at ~0.1 s/run (vs. iterating every MS2).
_PRECURSOR_SAMPLE_N = 300
_PRECURSOR_TOL_PPM = 10.0  # default; overridable via config.scan_precursor_tol_ppm
_PRECURSOR_MIN_FRAC = 0.5


def check_scan_precursor_consistency(
    psms: Sequence[PSMRecord],
    mzml: IndexedMzML,
    *,
    sample_n: int = _PRECURSOR_SAMPLE_N,
    tol_ppm: float = _PRECURSOR_TOL_PPM,
    min_frac: float = _PRECURSOR_MIN_FRAC,
) -> ScanPrecursorCheck:
    """Verify the mzTab ``spectra_ref`` scans index THIS mzML, by precursor m/z.

    For a sample of PSMs, each one's MS2 scan in this mzML must carry a precursor
    (selected-ion) m/z matching the PSM's reported ``precursor_mz``. Matching is
    by **mass** (ppm), so — unlike the scan↔RT guard it replaces — it is immune
    to OpenMS RT alignment, while still catching a real file/search mismatch
    (wrong file → scans absent or pointing at different precursors → low
    :attr:`~ScanPrecursorCheck.frac_matched`). PSMs without a usable scan or
    reported precursor m/z (the Percolator path) are skipped; with none, the
    guard no-ops (``n_checked == 0``).
    """
    cand = [p for p in psms if p.scan > 0 and p.precursor_mz > 0]
    if not cand:
        return ScanPrecursorCheck(0, float("nan"), float("nan"), min_frac)
    if len(cand) > sample_n:
        sample = [cand[i]
                  for i in np.linspace(0, len(cand) - 1, sample_n).astype(int)]
    else:
        sample = cand
    dppm: list[float] = []
    matched = 0
    for p in sample:
        obs = mzml.precursor_mz(p.scan)
        if obs is not None and obs > 0:
            d = abs(p.precursor_mz - obs) / obs * 1e6
            dppm.append(d)
            if d <= tol_ppm:
                matched += 1
        # obs is None ⇒ scan absent from this mzML ⇒ checked-but-unmatched
    return ScanPrecursorCheck(
        n_checked=len(sample),
        median_ppm=float(np.median(dppm)) if dppm else float("nan"),
        frac_matched=matched / len(sample),
        min_frac=min_frac,
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

    # RT-anchored intake (DIA, and MBR transfers): io.diann and the MBR pass emit
    # scan=-1 because there is no MS2 scan to anchor on, carrying the apex in
    # retention_time. Resolve each to the nearest MS1 scan in THIS mzML (+ an
    # RT-in-bounds check) so the scan-based extraction below runs unchanged; the
    # apex finder then re-centres on the true MS1 apex. `directly_scanned` is
    # captured BEFORE resolution (order-preserving) so the guard below can tell
    # real MS2 IDs from RT-anchored rows; no-op when every PSM has a real scan.
    directly_scanned = [p.scan >= 0 for p in psms]
    if not all(directly_scanned):
        psms = resolve_rt_anchored_scans(psms, mzml, file_label=file_label)

    # Intake scan↔precursor guard (Track A): verify the mzTab spectra_ref scans
    # actually index THIS mzML, by precursor m/z — mass-based, so immune to the
    # OpenMS RT alignment that made the old scan↔RT guard misfire on legitimate
    # aligned data (median offsets a few min, while the scans were correct). Runs
    # ONLY on the originally directly-scanned PSMs (real MS2 IDs); RT-anchored
    # rows (DIA, MBR) are excluded — their scan was *derived* from the RT
    # (resolve_rt_anchored_scans does the RT-in-bounds check for them instead).
    # Pure-DIA → empty subset → no-op; also no-ops on the Percolator path (no
    # reported precursor m/z).
    if config.check_scan_id:
        direct_psms = [p for p, d in zip(psms, directly_scanned) if d]
        if direct_psms:
            label = file_label or _mzml_basename(mzml)
            check = check_scan_precursor_consistency(
                direct_psms, mzml, tol_ppm=config.scan_precursor_tol_ppm)
            if check.n_checked > 0 and not check.ok:
                raise DataError(
                    f"{label}: scan↔precursor reconciliation FAILED — only "
                    f"{check.frac_matched:.0%} of {check.n_checked} sampled PSMs match a "
                    f"precursor m/z in this mzML (median |Δ| {check.median_ppm:.1f} ppm). "
                    "The mzTab spectra_ref scans do not point at the matching precursors "
                    "in this mzML — a real file/search mismatch: (1) the search ran on "
                    "different source files than these mzML; (2) the quantms filename-"
                    "prefix scan-scramble (zero-pad / de-prefix the basenames before the "
                    "quantms run); (3) a wrong mzML↔mzTab pairing. (Mass-based — immune "
                    "to RT alignment.) Override with --no-id-check only if this run is "
                    "knowingly correct."
                )
            if check.n_checked > 0:
                _LOGGER.info(
                    "%s: scan↔precursor reconciled — %.0f%% of %d sampled PSMs match "
                    "(median |Δ| %.1f ppm).",
                    label, 100 * check.frac_matched, check.n_checked, check.median_ppm,
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

    # Isotopomer channel set. Fixed path: the shared ``config.isotopomers`` tuple.
    # Adaptive path (``--iso auto``): each peptidoform's channels come from its
    # IsoSpec init∪final envelope (``_adaptive_concat_channels``, cached per
    # concat). The run's output columns span the widest peptide (run-max); short
    # peptides extract only their own channels and are NaN-padded in the row loop.
    if config.adaptive_iso:
        concat_channels = _adaptive_concat_channels(kept, config)
        run_max = max((len(v) for v in concat_channels.values()), default=1)
        isos = tuple(range(run_max))
    else:
        concat_channels = None
        isos = tuple(config.isotopomers)

    # Decode every MS1 once up front. Overlapping per-PSM RT windows otherwise
    # re-fetch+re-decode the same spectra hundreds of times in `peaks()` — the
    # dominant cost on dense runs (~90% of wall time; the redundant-decode
    # blow-up scaled the integrate to ~2 h/run on 7.8k-PSM fractions). The cache
    # is numerically transparent (identical arrays) and small for centroided MS1
    # (~0.4 GB for a 15k-scan run). After this the extraction's `peaks()` calls
    # are RAM lookups; the bottleneck reverts to the inherent per-centroid math.
    mzml.preload_peaks()

    # Per-PSM extraction → list of per-PSM (DataFrame, mass_accuracy_dict).
    # Serial: the hot path (pyteomics mzML parse + IsoSpec) is GIL-bound, so
    # threading it gave no speedup (in fact slower); cross-run parallelism is the
    # ProcessPool lever (`integrate --workers`), see core/pipeline.
    def _do(idx: int):
        psm = kept[idx]
        cm = concat_channels[psm.concat] if concat_channels is not None else None
        psm_isos = tuple(range(len(cm))) if cm is not None else isos
        return _extract_per_psm(
            psm, concat_scans, psm_isos, config, mzml,
            anchor_scan=concat_anchor[psm.concat],
            channel_masses=cm,
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
    dropped_pep_ids: set[int] = set()
    n_mbr_dropped = 0
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
        apex_snr = boundary.snr if boundary is not None else float("nan")
        # Nonzero scans in the integrated window — how many real MS1 points back
        # the peak. A sparse XIC (1–2 nonzero scans) can't define a reliable peak;
        # it is also why apex_snr goes inf (MAD=0). Emitted as a diagnostic + the
        # interpretable half of the MBR quality gate.
        iso0_window = (
            idf[iso0_col].to_numpy(dtype=np.float64)[boundary.lo : boundary.hi + 1]
            if boundary is not None
            else idf[iso0_col].to_numpy(dtype=np.float64)
        )
        n_scans = int(np.count_nonzero(iso0_window))
        # MBR two-part graceful failure: drop a transferred precursor with no
        # detectable apex, too few nonzero scans (--mbr-min-scans — sparse, can't
        # trust the peak), or an apex that fails the SNR floor (--mbr-min-snr; an
        # inf SNR = MAD=0 = no noise floor also FAILS — it is not a defined SNR).
        # The real-data apex_snr/n_scans stratification sets both
        # (reports/2026-06-17_mbr_v1_design.md). Real q-value PSMs keep the
        # fixed-window fallback (their ID is independent evidence the peptide is
        # here); R²>0.95 is the downstream backstop.
        if psm.evidence == "mbr" and (
            boundary is None
            or (config.mbr_min_scans > 0 and n_scans < config.mbr_min_scans)
            or (config.mbr_min_snr > 0
                and not (np.isfinite(apex_snr) and apex_snr >= config.mbr_min_snr))
        ):
            dropped_pep_ids.add(int(psm.pep_id))
            n_mbr_dropped += 1
            continue
        if boundary is not None:
            n_detected += 1
        else:
            n_fallback += 1

        for col in iso_cols:
            # Adaptive N_ISO: ``iso_cols`` spans the run-wide max, but a short
            # peptidoform only extracted its own (narrower) channel set — its
            # higher channels are absent from ``idf``. NaN-pad them: an honest
            # "not a channel for this peptide" (the solver trims trailing NaN),
            # distinct from an integrated 0.0 ("looked, found nothing").
            if col not in idf.columns:
                row.append(np.nan)
                continue
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
        # apex_snr (prominence / local-noise; NaN on the fixed-window / consensus
        # path, inf on a MAD=0 sparse trace) + n_scans (nonzero scans in the
        # window): diagnostics (apex_snr also informs the calibration noise-floor
        # question) and the two-part MBR quality gate.
        row.append(apex_snr)
        row.append(n_scans)
        integrated_rows.append(row)

    if config.peak_rt in ("apex", "consensus") or config.integration_half_width == "auto":
        _LOGGER.info(
            "integrate_run: %d PSMs detected, %d fell back to fixed-window "
            "(detection or co-elution failure).",
            n_detected, n_fallback,
        )

    integrated_df = pd.DataFrame(
        integrated_rows,
        columns=["pep_id"] + iso_cols + mz_cols + ppm_cols + ["apex_snr", "n_scans"],
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
    # Drop MBR transfers with no detectable apex (see the extraction loop): they
    # carry no integrated signal, so they must not reach the output as NaN rows.
    if dropped_pep_ids:
        out = out[~out["pep_id"].isin(dropped_pep_ids)].reset_index(drop=True)
        _LOGGER.info(
            "integrate_run: dropped %d MBR transfer(s) with no detectable apex.",
            n_mbr_dropped,
        )
    out["file"] = file_label if file_label is not None else _mzml_basename(mzml)
    # Stash the drift summary on the DataFrame as attrs so callers/writers can
    # emit the footer without re-computing.
    out.attrs["drift_summary"] = drift
    # MBR gate drops happen here (in a worker on the --workers path), so the per-run
    # INFO log never reaches the main logfile; stash the count for the main process
    # (finalize_run) to surface.
    out.attrs["n_mbr_dropped"] = n_mbr_dropped
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


def _adaptive_concat_channels(
    psms: Sequence[PSMRecord], config: IntegrationConfig,
) -> dict[str, tuple[float, ...]]:
    """Per-concat adaptive isotopomer channel masses for ``--iso auto``.

    One IsoSpec init+final envelope per *distinct* peptidoform (memoized in
    :mod:`riana.algorithms.isotope_dist`), so the cost is bounded by the
    peptidoform count, not the PSM count. Returns ``{concat: (neutral channel
    masses)}`` with iso0 == the precursor m0.

    Derives ``(bare_seq, mods, pep_mass)`` exactly as the fit path does
    (:mod:`riana.core.fitting`) — bracket-stripped sequence + parsed ``[UNIMOD:N]``
    ids + the record's recomputed ``peptide_mass`` — so the integrate envelope and
    the fit envelope are built from the same atoms. A peptidoform whose envelope
    can't be built (non-canonical residue B/X/U) falls back to a fixed m0..m5
    analytic set so it still integrates.
    """
    from riana.algorithms import isotope_dist as idist
    from riana.algorithms.mass_calc import parse_unimod_ids
    from riana.utils import strip_concat

    out: dict[str, tuple[float, ...]] = {}
    for p in psms:
        if p.concat in out:
            continue
        try:
            bare = strip_concat(p.sequence)
            mods = tuple(parse_unimod_ids(p.sequence))
            out[p.concat] = idist.adaptive_channel_masses(
                bare, p.peptide_mass,
                ria_max=config.ria_max, mods=mods,
                abundance_floor=config.iso_abundance_floor,
                iso_max=config.iso_max,
            )
        except (KeyError, ValueError):
            md = config.mass_difference
            out[p.concat] = tuple(p.peptide_mass + i * md for i in range(6))
    return out


#: MS1 centroid count at/above which the binary-searched window beats the O(n)
#: boolean mask in :func:`_window_sum`. The searchsorted path is ~flat in n
#: (O(log n)) but carries a fixed two-dispatch overhead, so on sparse spectra the
#: mask is faster; the measured crossover is ~7k centroids (single-protein BSA is
#: ~2k → mask; a complex-proteome Orbitrap MS1 is 10-40k → searchsorted, up to
#: ~2.7× at 30k). Both branches are byte-identical, so this only picks the faster
#: one — never changes a number.
_SEARCHSORTED_MIN_PEAKS = 8000


def _window_sum(
    mz_arr: np.ndarray,
    intens_arr: np.ndarray,
    target: float,
    delta: float,
    n_peaks: int,
) -> tuple[float, float]:
    """``(Σ intensity, Σ intensity·m/z)`` of the centroids within
    ``|m/z − target| ≤ delta`` in one MS1 scan.

    Two byte-identical implementations, selected by centroid count
    (:data:`_SEARCHSORTED_MIN_PEAKS`):

    * **Sparse** — the plain O(n) ``np.abs(mz_arr − target) ≤ delta`` mask.
    * **Dense** — ``mz_arr`` is ascending (centroid mzML convention, verified once
      per scan in :meth:`io.mzml.IndexedMzML.preload_peaks`), so the ±delta window
      is a contiguous run: two binary searches (``np.searchsorted``, O(log n))
      bracket it instead of scanning every centroid — the residual per-PSM cost
      after the MS1-decode precache, dominant on dense runs where ``use_range``
      spans the whole concat scan span.

    The dense branch is byte-identical to the mask by construction: the search only
    *narrows* the candidates (padded one index each side so ULP rounding of
    ``target ± delta`` can never drop a boundary peak the exact test keeps), then
    the exact predicate makes the final selection — over a contiguous ascending
    slice, so the matched sequence and its float reductions match bit-for-bit.
    Returns ``(0.0, 0.0)`` for an empty window (incl. an empty scan).
    """
    if n_peaks < _SEARCHSORTED_MIN_PEAKS:
        mask = np.abs(mz_arr - target) <= delta
        if not mask.any():
            return 0.0, 0.0
        matched_mz = mz_arr[mask]
        matched_i = intens_arr[mask]
        return float(matched_i.sum()), float((matched_mz * matched_i).sum())

    lo = max(0, int(np.searchsorted(mz_arr, target - delta, side="left")) - 1)
    hi = min(n_peaks, int(np.searchsorted(mz_arr, target + delta, side="right")) + 1)
    if hi <= lo:
        return 0.0, 0.0
    cand_mz = mz_arr[lo:hi]
    mask = np.abs(cand_mz - target) <= delta
    if not mask.any():
        return 0.0, 0.0
    matched_mz = cand_mz[mask]
    matched_i = intens_arr[lo:hi][mask]
    return float(matched_i.sum()), float((matched_mz * matched_i).sum())


def _extract_per_psm(
    psm: PSMRecord,
    concat_scans: dict[str, tuple[int, int]],
    isos: tuple[int, ...],
    config: IntegrationConfig,
    mzml: IndexedMzML,
    *,
    anchor_scan: int | None = None,
    channel_masses: tuple[float, ...] | None = None,
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
        _win = (
            (mzml.rt_idx - rt_lo > -config.extraction_half_width)
            & (mzml.rt_idx - rt_hi < config.extraction_half_width)
        )
    else:
        center_scan = anchor_scan if anchor_scan is not None else psm.scan
        rt_center = mzml.rt_idx[np.searchsorted(mzml.scan_idx, center_scan, side="left") - 1]
        _win = np.abs(mzml.rt_idx - rt_center) <= config.extraction_half_width
    # Carry each window scan's RT straight from the index (``rt_idx[_win]``)
    # instead of re-finding it per scan with ``scan_idx == scan`` (an O(n_MS1)
    # linear scan that dominated once the decode was cached). Same values.
    nearby = mzml.scan_idx[_win]
    nearby_rt = mzml.rt_idx[_win]

    # Lazy peak fetch + centroid sum per (mod, iso) per MS1. As a Phase D
    # addition we also accumulate intensity-weighted observed m/z per
    # (mod, iso) across the integration window so the writer can emit the
    # iso{N}_obs_mz / iso{N}_ppm_error mass-accuracy columns.
    rows: list[tuple[str, float, float]] = []
    iso_cols_local = [f"iso{n}" for n in isos]
    # Per-isotopomer extraction target m/z, computed once. Fixed path: analytic
    # ``m0 + iso·Δ/z``. Adaptive path (``channel_masses`` given, aligned to
    # ``isos`` by position): the IsoSpec **init (unlabeled)** averaged-isotopolog
    # accurate mass per channel — iso0 is the precursor m0. The ppm_error below is
    # taken vs this target, so on the adaptive path it reads the mass-defect drift
    # from the θ=0 position (the orthogonal mass-defect-θ substrate, Track B).
    iso_target_mz: dict[int, float] = {}
    for k, iso in enumerate(isos):
        if channel_masses is not None:
            iso_target_mz[iso] = (channel_masses[k] + charge * proton) / charge
        else:
            iso_target_mz[iso] = peptide_prec + (iso * iso_added_mass / charge)
    targets: dict[str, float] = {}
    mz_weighted: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    intens_total: dict[str, float] = {c: 0.0 for c in iso_cols_local}
    ppm_tol = float(config.mass_tol_ppm) * 1e-6
    for scan, rt in zip(nearby.tolist(), nearby_rt.tolist()):
        scan_int = int(scan)
        mz_arr, intens_arr = mzml.peaks(scan_int)
        n_peaks = mz_arr.shape[0]
        for iso in isos:
            col = f"iso{iso}"
            target = iso_target_mz[iso]
            targets[col] = target
            delta = target * ppm_tol
            # Binary-searched centroid window (byte-identical to the old O(n)
            # boolean mask; mz_arr is ascending — see _window_sum). Adding the
            # 0.0 returned for an empty window is a no-op, so the accumulation
            # matches the old "only on match" form bit-for-bit.
            summed, mz_w = _window_sum(mz_arr, intens_arr, target, delta, n_peaks)
            mz_weighted[col] += mz_w
            intens_total[col] += summed
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
    apex_snr = float("nan")  # set by the "apex" path; nan for consensus
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
        res = pk.find_apex(
            rt_arr, iso0_trace, scan_prior_rt=psm_rt,
            prominence_k=config.prominence_k,
            apex_search_half_width=config.apex_search_half_width,
            selection=config.apex_selection,
        )
        apex = res[0] if res is not None else None
        apex_snr = res[1] if res is not None else float("nan")
    if apex is None:
        return None
    half = float(config.integration_half_width)
    apex_rt = float(rt_arr[apex])
    win = np.where(np.abs(rt_arr - apex_rt) <= half)[0]
    return pk.PeakBoundary(
        apex_idx=apex, lo=int(win[0]), hi=int(win[-1]),
        apex_intensity=float(iso0_trace[apex]), prominence=float("nan"),
        snr=apex_snr,
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
