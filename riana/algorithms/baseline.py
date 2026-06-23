"""Chromatographic baseline subtraction.

Pure ``numpy`` functions, one per method. Phase C composes them inside
``core/integration.py`` after peak boundary detection.

Methods (matching ``bench_baseline.py``'s expected ``--method`` set):

- ``noise_floor`` — flat baseline at a low quantile (default p10) of the
  full extracted trace. Robust to single-scan noise spikes and to where
  exactly the detected boundaries land. **This is the recommended default
  for our calibration data** — see PROJECT_REVIEW.md Week 3 notes.
- ``snip`` — Statistics-sensitive Non-linear Iterative Peak clipping
  (``pybaselines``). Robust to broad humps.
- ``asls`` — Asymmetric Least Squares (``pybaselines``). Smooth,
  parameterised baseline; good for noisy traces with structured drift.

``pybaselines`` is imported lazily inside the snip/asls calls so the rest
of the pipeline (and the unit tests on ``noise_floor``) doesn't pay the
import cost.

NB: ``local_linear`` (Skyline-style baseline between the integration
boundaries) was implemented and **discarded** 2026-06-05 — it is the correct
Skyline algorithm but over-subtracts here because our narrow / fixed-apex
boundaries sit *on* the peak, not at the chromatographic floor Skyline's
boundaries land on. It failed on all three calibration lines (worst with
narrow windows). See PROJECT_REVIEW for the full rationale.
"""

from __future__ import annotations

import numpy as np


def noise_floor(
    intensity: np.ndarray,
    *,
    quantile: float = 0.10,
) -> np.ndarray:
    """Flat noise-floor baseline from a low quantile of the trace.

    The instrument noise floor at a given m/z window is approximately
    constant across the RT range we extract (±extraction_half_width around
    the PSM scan). A low quantile of the **full extracted trace** captures
    that floor: the peak occupies only part of the window, so the lower tail
    of the intensity distribution is the off-peak noise.

    Returns a per-scan baseline array (same length as ``intensity``) so
    callers can use the same subtraction code path as the other methods.
    Robust to single noise spikes (a flat quantile, not anchored on any scan).

    TODO (off-peak estimate): because the quantile is over the *full
    extraction*, its quality depends on the (hidden) ``extraction_half_width``
    — too narrow an extraction lets the peak dominate and p10 creeps onto the
    flank, over-subtracting. The principled fix is to estimate the floor from
    the *off-peak* scans only (outside ``[lo, hi]``), decoupling it from peak
    occupancy and extraction width. Deferred — ``none`` is the default and the
    narrow apex window already removes background by exclusion (see
    PROJECT_REVIEW).

    Args:
        intensity: per-scan summed intensities for the full extracted RT
            window.
        quantile: which quantile of the trace to treat as noise floor.
            Default 0.10 — robust to peak occupying up to ~90% of window
            (rare; usually the peak is < 50% of the window).
    """
    arr = np.asarray(intensity, dtype=np.float64)
    if arr.size == 0:
        return arr.copy()
    floor = float(np.quantile(arr, quantile))
    return np.full(arr.size, floor, dtype=np.float64)


def snip(intensity: np.ndarray, n_iter: int = 40) -> np.ndarray:
    """SNIP baseline via :mod:`pybaselines`.

    Args:
        intensity: per-scan summed intensities.
        n_iter: number of iterations (the SNIP window grows up to this).
    """
    from pybaselines import Baseline  # lazy import

    arr = np.asarray(intensity, dtype=np.float64)
    if arr.size == 0:
        return arr.copy()
    fitter = Baseline(x_data=None)
    base, _ = fitter.snip(arr, max_half_window=n_iter)
    return np.asarray(base, dtype=np.float64)


def asls(
    intensity: np.ndarray,
    lam: float = 1e6,
    p: float = 0.01,
) -> np.ndarray:
    """Asymmetric Least Squares baseline via :mod:`pybaselines`.

    Args:
        intensity: per-scan summed intensities.
        lam: smoothness penalty (Eilers & Boelens defaults are 1e6-1e9).
        p: asymmetry weight in (0, 1); lower ⇒ more aggressive on the
            below-baseline side.
    """
    from pybaselines import Baseline  # lazy import

    arr = np.asarray(intensity, dtype=np.float64)
    if arr.size == 0:
        return arr.copy()
    fitter = Baseline(x_data=None)
    base, _ = fitter.asls(arr, lam=lam, p=p)
    return np.asarray(base, dtype=np.float64)
