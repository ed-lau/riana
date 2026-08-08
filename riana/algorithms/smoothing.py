"""Smoothing for chromatographic intensity traces.

A pure-function module on ``numpy`` arrays — Phase B lays it down without
wiring; Phase C composes it inside ``core/integration.py``.

The 0.9.0 path used Savitzky–Golay with ``polyorder=1`` (PROJECT_REVIEW.md
§2c point 3). That is mathematically a moving average, which distorts peak
heights and shifts the area integral with window size (M2's
``bench_smoothing.py`` reproduced this: per-AA coefficients drift up to
0.22 at S=9 on iPSC). The roadmap's fix is to require ``polyorder ≥ 2``,
which preserves second-derivative features and is the standard Skyline
choice.
"""

from __future__ import annotations

from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.signal


def savgol(
    intensity: npt.NDArray[np.float64],
    window: int,
    polyorder: int = 2,
) -> npt.NDArray[np.float64]:
    """Savitzky–Golay smoothing with ``polyorder ≥ 2`` enforced.

    Mirrors the legacy ``mode='nearest'`` boundary handling so a like-for-like
    swap-in produces the same shape but corrected curvature. An all-zero or
    sub-window trace is returned unchanged — the 0.9.0 caller already gated
    the call on ``np.sum(intensity) > 0`` for the same reason.

    Args:
        intensity: 1-D intensity trace.
        window: odd window length, ≥ ``polyorder + 1``.
        polyorder: polynomial order; defaults to 2 (the §2c fix).

    Raises:
        ValueError: window even, < 3, ≤ polyorder, or polyorder < 2.
    """
    if polyorder < 2:
        raise ValueError(
            f"polyorder must be ≥ 2 (the §2c fix); got {polyorder}. "
            "polyorder=1 is mathematically a moving average and distorts areas."
        )
    if window < 3 or window % 2 == 0:
        raise ValueError(f"window must be an odd integer ≥ 3, got {window}")
    if window <= polyorder:
        raise ValueError(
            f"window ({window}) must exceed polyorder ({polyorder})"
        )

    arr = np.asarray(intensity, dtype=np.float64)
    if arr.size < window or not np.any(arr > 0):
        return arr.copy()
    return cast(
        "npt.NDArray[np.float64]",
        scipy.signal.savgol_filter(
            arr, window_length=window, polyorder=polyorder, mode="nearest"
        ),
    )
