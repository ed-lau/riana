"""Mass-domain refinement and per-fraction calibration drift summary.

Two responsibilities:

1. Per-isotopomer mass-accuracy outputs (PROJECT_REVIEW.md "Mass-accuracy
   output" §): intensity-weighted observed m/z and the ppm error vs. the
   theoretical target. Populates ``IsotopomerPeak.obs_mz`` /
   ``ppm_error`` so a Phase D ``_riana.txt`` column comes from a
   well-defined number, not a heuristic.

2. Per-fraction drift summary — median ppm error, MAD ppm, and a
   suggested calibration shift — written to the per-fraction footer of
   ``_riana.txt`` and surfaced as a logger warning when the median
   exceeds the ``--ppm-alert`` threshold (Phase D wire-up).

Pure numpy functions; Phase B unit-tests them on synthetic peak data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class DriftSummary:
    """Per-fraction mass-accuracy drift statistics for the footer."""

    n: int
    median_ppm: float
    mad_ppm: float
    suggested_shift_ppm: float


def weighted_obs_mz(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    target_mz: float,
    half_width_ppm: float,
) -> tuple[float | None, float | None]:
    """Intensity-weighted observed m/z and ppm error vs. target.

    Uses the same ``±N ppm`` half-window the integrator extracts on. If no
    centroid lands in window, both outputs are ``None`` — an honest
    "not observed" rather than a fabricated zero (the same posture
    :class:`IsotopomerPeak` takes for unset mass-accuracy fields).

    Args:
        mz_array: centroid m/z values from one MS1 spectrum.
        intensity_array: parallel centroid intensities.
        target_mz: theoretical m/z for the isotopomer.
        half_width_ppm: ± half-width (mirrors :attr:`IntegrationConfig.mass_tol_ppm`).

    Returns:
        ``(obs_mz, ppm_error)``. ``ppm_error`` is signed: positive ⇒
        observed mass higher than theoretical.
    """
    if mz_array.size == 0:
        return None, None
    delta = target_mz * half_width_ppm * 1e-6
    mask = np.abs(mz_array - target_mz) <= delta
    if not mask.any():
        return None, None
    mz_in = mz_array[mask]
    i_in = intensity_array[mask]
    weight = float(i_in.sum())
    if weight <= 0:
        return None, None
    obs_mz = float(np.sum(mz_in * i_in) / weight)
    ppm_error = float((obs_mz - target_mz) / target_mz * 1e6)
    return obs_mz, ppm_error


def drift_summary(ppm_errors: np.ndarray) -> DriftSummary:
    """Aggregate per-PSM ppm errors into the fraction-level footer.

    NaN entries are treated as "no observation" and excluded.

    The suggested shift is simply the median: applying it as a fixed
    correction zeroes the median. The user's actual mass calibration
    decision is theirs — Riana just reports.
    """
    arr = np.asarray(ppm_errors, dtype=np.float64)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return DriftSummary(n=0, median_ppm=float("nan"),
                            mad_ppm=float("nan"),
                            suggested_shift_ppm=0.0)
    med = float(np.median(arr))
    mad_ppm = float(1.4826 * np.median(np.abs(arr - med)))
    return DriftSummary(
        n=int(arr.size),
        median_ppm=med,
        mad_ppm=mad_ppm,
        suggested_shift_ppm=-med,
    )
