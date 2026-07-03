"""Byte-identity + guard tests for the searchsorted masking vectorization.

``core.integration._window_sum`` replaced the per-(scan, isotopomer) O(n)
``np.abs(mz_arr - target) <= delta`` boolean mask with a binary-searched window;
it MUST reproduce that mask's float reductions bit-for-bit, and
``IndexedMzML.preload_peaks`` must reject a non-ascending m/z array (the rewrite's
precondition) rather than silently under-summing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import riana.core.integration as integ
from riana.core.integration import _window_sum
from riana.exceptions import DataError
from riana.io import mzml as mzml_mod
from riana.io.mzml import IndexedMzML

SAMPLE1 = Path("tests/data/sample1")
BSA = SAMPLE1 / "20180216_BSA.mzML.gz"


@pytest.fixture(params=[0, 10**9], ids=["searchsorted", "boolmask"])
def both_branches(request, monkeypatch):
    """Force _window_sum through each of its two byte-identical branches, so the
    small test arrays exercise the searchsorted path (threshold 0) as well as the
    mask path (threshold ∞), not just whichever their size selects."""
    monkeypatch.setattr(integ, "_SEARCHSORTED_MIN_PEAKS", request.param)
    return request.param


def _ref(mz, intens, target, delta):
    """The pre-change O(n) boolean-mask reduction — the ground truth."""
    mask = np.abs(mz - target) <= delta
    if not mask.any():
        return 0.0, 0.0
    mm, mi = mz[mask], intens[mask]
    return float(mi.sum()), float((mm * mi).sum())


def _assert_identical(mz, intens, target, delta):
    got = _window_sum(mz, intens, target, delta, mz.shape[0])
    exp = _ref(mz, intens, target, delta)
    assert got == exp, f"target={target!r} delta={delta!r}: {got} != {exp}"


@pytest.mark.slow
def test_window_sum_matches_boolean_mask_on_real_scan(both_branches):
    """Across exact peak hits, inter-peak midpoints, and off-ends on a real
    BSA MS1 scan, the searchsorted window equals the boolean mask bit-for-bit."""
    with IndexedMzML(BSA) as mz:
        scan = int(mz.scan_idx[len(mz.scan_idx) // 2])
        a, i = mz.peaks(scan)
    assert a.shape[0] > 100 and bool(np.all(a[:-1] <= a[1:]))
    ppm = 50e-6
    targets = [float(v) for v in a[::7]]                                   # exact hits
    targets += [float((a[k] + a[k + 1]) / 2)                              # midpoints
                for k in range(0, a.shape[0] - 1, 13)]
    targets += [float(a[0] - 5), float(a[-1] + 5),                        # off both ends
                float(a[0]), float(a[-1])]                                # first / last peak
    for t in targets:
        for delta in (t * ppm, t * ppm * 10, 0.0):                       # incl. delta=0
            _assert_identical(a, i, t, float(delta))


def test_window_sum_boundary_and_empty(both_branches):
    """Peaks exactly at ``target ± delta`` are inclusive; empty / single-peak
    scans (the ``n_peaks`` diff-guard boundary) return zero without crashing."""
    target, delta = 500.0, 0.01
    mz = np.array([target - delta, target - delta / 2, target,
                   target + delta, target + 2 * delta], dtype=np.float64)
    intens = np.array([1.0, 2.0, 4.0, 8.0, 16.0], dtype=np.float64)
    _assert_identical(mz, intens, target, delta)
    # The four in-window peaks (both boundaries inclusive), not the +2δ one.
    assert _window_sum(mz, intens, target, delta, mz.shape[0])[0] == 15.0

    empty = np.array([], dtype=np.float64)
    assert _window_sum(empty, empty, target, delta, 0) == (0.0, 0.0)

    one = np.array([target], dtype=np.float64)
    _assert_identical(one, np.array([3.0]), target, delta)
    _assert_identical(one, np.array([3.0]), target + 1.0, delta)  # single peak, no hit


@pytest.mark.slow
def test_preload_rejects_non_ascending_mz(monkeypatch):
    """The sortedness guard fires loudly rather than silently under-summing."""
    real = mzml_mod._spec_peaks

    def descending(spec):
        m, i = real(spec)
        return m[::-1].copy(), i  # reverse a real (ascending) scan → non-ascending

    monkeypatch.setattr(mzml_mod, "_spec_peaks", descending)
    with IndexedMzML(BSA) as mz:
        with pytest.raises(DataError, match="non-ascending"):
            mz.preload_peaks()
