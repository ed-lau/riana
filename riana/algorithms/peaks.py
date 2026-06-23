"""Peak detection on chromatographic isotopomer XICs (Skyline-shaped).

Per PROJECT_REVIEW.md §2c.

The contract:
1. Detect candidate peaks in the iso0 XIC inside the PSM RT window via
   :func:`scipy.signal.find_peaks`, with a prominence floor from a local
   noise estimate (median-absolute-deviation of the trace).
2. Pick the apex nearest the PSM-scan RT prior — that's our strongest
   indicator of where the peptide elutes.
3. Determine the integration boundaries via :func:`scipy.signal.peak_widths`
   at ``rel_height=0.05`` (Skyline's "95% of apex below peak").
4. Co-elution check: iso0 and iso1 apices must agree within 2× the median
   MS1 cycle time. If they don't, the caller falls back to fixed-window
   integration for that peptide (peak detection's failure mode is
   integrating a co-eluting decoy).
5. Quality scoring (S/N + symmetry) feeds the IsotopomerPeak record so
   downstream curation can flag suspicious peaks.

All functions are pure (numpy in, dataclass/tuple out) so Phase B unit
tests can verify each piece independently, and Phase C composes them
inside the integrator's per-PSM loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.signal


@dataclass(frozen=True, slots=True)
class PeakBoundary:
    """Outcome of peak detection for one isotopomer trace."""

    #: Index of the chosen apex within the input trace.
    apex_idx: int
    #: Inclusive left boundary index.
    lo: int
    #: Inclusive right boundary index.
    hi: int
    #: Intensity at the apex (post-smoothing if applied).
    apex_intensity: float
    #: Prominence of the apex (above the local baseline at its base).
    prominence: float
    #: Apex signal-to-noise = prominence / (1.4826·MAD of the trace) — the units
    #: of the admission gate (a peak passes when ``snr >= prominence_k``). NaN when
    #: not computed (the consensus apex / fixed-window fallback paths). Emitted as
    #: the ``apex_snr`` output column and used by the MBR SNR floor.
    snr: float = float("nan")


def detect_peak(
    rt: np.ndarray,
    intensity: np.ndarray,
    *,
    scan_prior_rt: float,
    rel_height: float = 0.05,
    prominence_k: float = 3.0,
) -> PeakBoundary | None:
    """Detect the dominant peak nearest the PSM-scan RT prior.

    Args:
        rt: per-scan retention times (monotonic).
        intensity: per-scan summed intensities.
        scan_prior_rt: RT of the PSM scan — the strongest prior on where
            this peptide elutes.
        rel_height: width measurement threshold; ``0.05`` matches Skyline's
            "95% of apex below peak" default.
        prominence_k: multiplier on local-MAD noise for the prominence
            floor.  Higher ⇒ stricter peak admission.

    Returns:
        :class:`PeakBoundary` if a peak passes the prominence floor, else
        ``None`` (caller falls back to fixed-window).
    """
    intensity = np.asarray(intensity, dtype=np.float64)
    rt = np.asarray(rt, dtype=np.float64)
    if intensity.size < 3:
        return None
    if not np.any(intensity > 0):
        return None

    # Robust noise estimate: 1.4826 * MAD is the gaussian-σ-equivalent
    # of the median absolute deviation, less brittle than std().
    med = float(np.median(intensity))
    mad = float(np.median(np.abs(intensity - med)))
    noise = 1.4826 * mad
    prominence_floor = max(prominence_k * noise, 1.0)

    apex_idxs, props = scipy.signal.find_peaks(
        intensity, prominence=prominence_floor
    )
    if apex_idxs.size == 0:
        return None

    # Choose the apex nearest the PSM-scan RT prior.
    best = int(apex_idxs[np.argmin(np.abs(rt[apex_idxs] - scan_prior_rt))])
    best_local = int(np.where(apex_idxs == best)[0][0])

    widths, _, lo_f, hi_f = scipy.signal.peak_widths(
        intensity, [best], rel_height=1.0 - rel_height
    )
    lo = int(np.clip(np.floor(lo_f[0]), 0, intensity.size - 1))
    hi = int(np.clip(np.ceil(hi_f[0]), 0, intensity.size - 1))
    if lo > hi:
        lo, hi = hi, lo

    prom = float(props["prominences"][best_local])
    return PeakBoundary(
        apex_idx=best,
        lo=lo,
        hi=hi,
        apex_intensity=float(intensity[best]),
        prominence=prom,
        snr=(prom / noise) if noise > 0 else float("inf"),
    )


def find_apex(
    rt: np.ndarray,
    intensity: np.ndarray,
    *,
    scan_prior_rt: float,
    prominence_k: float = 3.0,
    apex_search_half_width: float = 0.0,
    selection: str = "nearest",
) -> "tuple[int, float] | None":
    """Apex index **and SNR** of a prominent peak on one channel's XIC.

    Gate candidate maxima by a local-MAD prominence floor, restrict them to a
    search window around the PSM-scan RT prior, then choose one:

    - ``selection="nearest"`` — the candidate closest to the MS2 RT prior (DDA
      fired the scan on *this* precursor, so proximity disambiguates a
      co-eluting neighbour).
    - ``selection="tallest"`` — the most intense candidate (the shipped default;
      iso0 is often suppressed at high D₂O, so nearest-by-RT can latch onto a
      noise bump — the cross-proportion A/B preferred tallest).

    ``apex_search_half_width`` bounds how far the apex may sit from the PSM RT
    (``<= 0`` ⇒ whole trace, the wide default). A tight bound keeps the apex near
    the labelling-independent PSM RT, which matters for cross-proportion
    stability (and for MBR, where the transferred RT is a trustworthy anchor).

    Returns ``(apex_idx, snr)`` where ``snr = prominence / (1.4826·MAD)`` — the
    admission-gate units (a candidate clears the gate at ``snr >= prominence_k``)
    — or ``None`` when nothing clears the floor inside the window.
    """
    intensity = np.asarray(intensity, dtype=np.float64)
    rt = np.asarray(rt, dtype=np.float64)
    if intensity.size < 3 or not np.any(intensity > 0):
        return None
    med = float(np.median(intensity))
    mad = float(np.median(np.abs(intensity - med)))
    noise = 1.4826 * mad
    prominence_floor = max(prominence_k * noise, 1.0)
    apex_idxs, props = scipy.signal.find_peaks(intensity, prominence=prominence_floor)
    if apex_idxs.size == 0:
        return None
    proms = props["prominences"]
    if apex_search_half_width > 0:
        keep = np.abs(rt[apex_idxs] - scan_prior_rt) <= apex_search_half_width
        apex_idxs, proms = apex_idxs[keep], proms[keep]
        if apex_idxs.size == 0:
            return None
    if selection == "tallest":
        sel = int(np.argmax(intensity[apex_idxs]))
    else:
        sel = int(np.argmin(np.abs(rt[apex_idxs] - scan_prior_rt)))
    snr = float(proms[sel] / noise) if noise > 0 else float("inf")
    return int(apex_idxs[sel]), snr


def consensus_apex(
    rt: np.ndarray,
    channel_traces,
    *,
    scan_prior_rt: float,
    prominence_k: float = 3.0,
    apex_search_half_width: float = 0.0,
    selection: str = "nearest",
) -> "tuple[int, float] | None":
    """Co-elution consensus apex: the **median** of the per-channel apexes.

    Every isotopomer of a peptide shares one true apex (co-elution), so the
    robust estimate is the median of the per-channel apex RTs. It ignores a
    channel whose apex drifts (contamination) *regardless of that channel's
    intensity* — unlike summing the channels, which a strong contaminant would
    pull — and a suppressed channel simply contributes no apex. So it is
    labelling-independent (there is always a strong channel somewhere in
    m0..m3) and cross-proportion stable (the median RT does not move as
    labelling shifts which channel dominates). :func:`find_apex` is the
    single-channel special case (one trace).

    Args:
        channel_traces: iterable of per-channel intensity arrays (e.g. iso0..3),
            all on the shared ``rt`` axis.

    Returns ``(apex_idx, spread)`` where ``spread`` is the MAD of the surviving
    per-channel apex RTs — a co-elution quality score (small = clean
    agreement) — or ``None`` if no channel yields an apex.
    """
    rt = np.asarray(rt, dtype=np.float64)
    apex_rts = []
    for trace in channel_traces:
        a = find_apex(
            rt, trace, scan_prior_rt=scan_prior_rt, prominence_k=prominence_k,
            apex_search_half_width=apex_search_half_width, selection=selection,
        )
        if a is not None:
            apex_rts.append(rt[a[0]])
    if not apex_rts:
        return None
    arr = np.asarray(apex_rts, dtype=np.float64)
    consensus_rt = float(np.median(arr))
    spread = float(np.median(np.abs(arr - consensus_rt)))
    apex_idx = int(np.argmin(np.abs(rt - consensus_rt)))
    return apex_idx, spread


def coelution_ok(
    iso0: PeakBoundary,
    iso1: PeakBoundary | None,
    rt: np.ndarray,
    *,
    cycle_tolerance: int = 2,
) -> bool:
    """True iff iso0/iso1 apices agree within ``cycle_tolerance`` MS1 cycles.

    Isotopomers of the same peptide must co-elute. A failed check is the
    signal that the iso0 "peak" is actually a co-eluting decoy whose iso1
    lives elsewhere — better to fall back to fixed-window than integrate
    contaminants.

    Args:
        iso0: detected iso0 boundary (the anchor).
        iso1: detected iso1 boundary, or ``None`` (no detectable iso1).
        rt: per-scan RT array (shared by both isotopomers, length-matched).
        cycle_tolerance: max apex separation in MS1 cycles.

    Returns:
        ``True`` if the apices agree, or if iso1 has no detectable peak
        (genuine low-iso1 envelope at high D₂O ⇒ don't reject on absence).
    """
    if iso1 is None:
        return True
    cycle = float(np.median(np.diff(rt)))
    return bool(abs(rt[iso1.apex_idx] - rt[iso0.apex_idx]) <= cycle_tolerance * cycle)


def snr(intensity: np.ndarray, baseline: np.ndarray) -> float:
    """Signal-to-noise: apex above baseline / MAD of (intensity − baseline)."""
    intensity = np.asarray(intensity, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    residual = intensity - baseline
    noise = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
    if noise <= 0:
        return float("inf") if residual.max() > 0 else 0.0
    return float((intensity.max() - baseline[intensity.argmax()]) / noise)


def symmetry(rt: np.ndarray, intensity: np.ndarray, apex_idx: int) -> float:
    """Peak symmetry ratio — area before apex / area after apex.

    1.0 is perfectly symmetric. Distance from 1.0 in either direction
    indicates fronting (<1) or tailing (>1). Reported as ``min(r, 1/r)``
    so the score lives in (0, 1].
    """
    intensity = np.asarray(intensity, dtype=np.float64)
    rt = np.asarray(rt, dtype=np.float64)
    if not 0 < apex_idx < intensity.size - 1:
        return 0.0
    left = float(np.trapezoid(intensity[: apex_idx + 1], x=rt[: apex_idx + 1]))
    right = float(np.trapezoid(intensity[apex_idx:], x=rt[apex_idx:]))
    if left <= 0 or right <= 0:
        return 0.0
    r = left / right
    return min(r, 1.0 / r)


def quality(snr_value: float, symmetry_value: float) -> float:
    """Composite quality — geometric mean of saturating S/N and symmetry.

    Both inputs are bounded into (0, 1] before mixing so a single bad
    factor cannot drag the composite implausibly low. S/N saturates at 10
    (anything above is "clearly real"), symmetry is already in (0, 1].
    """
    snr_capped = max(0.0, min(snr_value, 10.0)) / 10.0
    sym = max(0.0, min(symmetry_value, 1.0))
    if snr_capped <= 0 or sym <= 0:
        return 0.0
    return float(np.sqrt(snr_capped * sym))
