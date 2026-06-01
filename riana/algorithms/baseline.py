"""Chromatographic baseline subtraction.

Pure ``numpy`` functions, one per method. Phase C composes them inside
``core/integration.py`` after peak boundary detection.

Methods (matching ``bench_baseline.py``'s expected ``--method`` set):

- ``noise_floor`` — flat baseline at a low quantile (default p10) of the
  full extracted trace. Robust to single-scan noise spikes and to where
  exactly the detected boundaries land. **This is the recommended default
  for our calibration data** — see PROJECT_REVIEW.md Week 3 notes.
- ``local_linear`` — Skyline classic. Linear interpolation between the
  intensities at the detected peak boundaries (lo, hi). Works on Skyline-
  shaped data where boundary detection lands on true noise floor; brittle
  on our peak_widths(rel_height=0.05) boundaries because a single noise
  spike at one endpoint inflates the whole baseline.
- ``snip`` — Statistics-sensitive Non-linear Iterative Peak clipping
  (``pybaselines``). Robust to broad humps.
- ``asls`` — Asymmetric Least Squares (``pybaselines``). Smooth,
  parameterised baseline; good for noisy traces with structured drift.

``pybaselines`` is imported lazily inside the snip/asls calls so the rest
of the pipeline (and the unit tests on ``local_linear`` / ``noise_floor``)
doesn't pay the import cost.
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
    constant across the RT range we extract (±r_time around the PSM
    scan). A low quantile of the **full extracted trace** captures that
    floor: the peak occupies only part of the window, so the lower tail
    of the intensity distribution is the off-peak noise.

    Returns a per-scan baseline array (same length as ``intensity``) so
    callers can use the same subtraction code path as the other methods.
    Robust to single noise spikes — unlike :func:`local_linear` which
    anchors on two specific endpoint scans.

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


def local_linear(
    rt: np.ndarray,
    intensity: np.ndarray,
    lo: int,
    hi: int,
) -> np.ndarray:
    """Linear baseline between the intensities at ``rt[lo]`` and ``rt[hi]``.

    Returns a per-scan baseline array the same length as ``intensity``,
    suitable for direct subtraction. Inside ``[lo, hi]`` it interpolates;
    outside, it extends the endpoint values flat (no extrapolation noise).

    Args:
        rt: per-scan retention times.
        intensity: per-scan summed intensities.
        lo: index of the left peak boundary (inclusive).
        hi: index of the right peak boundary (inclusive).
    """
    rt = np.asarray(rt, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    n = intensity.size
    if not (0 <= lo <= hi < n):
        raise ValueError(
            f"boundaries (lo={lo}, hi={hi}) out of range for n={n}"
        )
    out = np.empty(n, dtype=np.float64)
    if lo == hi:
        out[:] = intensity[lo]
        return out

    slope = (intensity[hi] - intensity[lo]) / (rt[hi] - rt[lo])
    # Inside the peak: linear interpolation. Outside: flat-extend the
    # endpoint values so subtraction outside the peak doesn't introduce
    # baseline-extrapolation noise.
    out[:lo] = intensity[lo]
    out[lo : hi + 1] = intensity[lo] + slope * (rt[lo : hi + 1] - rt[lo])
    out[hi + 1 :] = intensity[hi]
    return out


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
