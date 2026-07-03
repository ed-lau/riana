"""Tests for the linearized turnover model + cross-sample Δk (`linear simple`)."""

import numpy as np
import pandas as pd
import pytest

from riana.core.linear_model import (
    LINEAR_COLUMNS,
    fit_linear_deltak,
    to_phi,
    truncate_plateau,
)


def _curve(exp, prot, cond, k, times, *, noise=0.0, seed=0):
    """A synthetic (condition, t, θ) point set for one protein from θ = 1−e^(−kt)."""
    rng = np.random.default_rng(seed)
    rows = []
    for t in times:
        theta = 1.0 - np.exp(-k * t) + (rng.normal(0, noise) if noise else 0.0)
        rows.append({"experiment": exp, "protein": prot, "condition": cond,
                     "labeling_time": float(t), "theta": float(theta)})
    return rows


def test_to_phi_clamps_and_transforms():
    phi = to_phi(np.array([0.0, 0.5, 0.999, 1.5, -0.2]))
    assert phi[0] == pytest.approx(np.log(0.999), abs=1e-9)  # floor 0.001
    assert phi[1] == pytest.approx(np.log(0.5))
    assert np.isfinite(phi).all()  # θ≥1 / θ<0 clamped, never -inf/NaN
    assert phi[3] == pytest.approx(np.log(1 - 0.999))  # ceiling


def test_truncate_plateau_drops_saturated_tail():
    day = np.array([0, 1, 2, 4, 8, 16], dtype=float)
    phi = np.array([0.0, -1.0, -2.0, -3.5, -5.0, -6.0])  # saturates past day 4
    t, p = truncate_plateau(day, phi, phi_limit=-4.0)
    assert list(t) == [0, 1, 2, 4]  # day 8/16 (φ=-5/-6) dropped
    assert (p > -4.0).all()


def test_recovers_known_slopes_and_delta_k():
    times = [0, 1, 2, 3, 4, 6, 8, 10]
    pts = pd.DataFrame(
        _curve("e", "P1", "control", 0.05, times, noise=0.005, seed=1)
        + _curve("e", "P1", "atrium", 0.10, times, noise=0.005, seed=2)
    )
    out = fit_linear_deltak(pts, reference_condition="control")
    assert list(out.columns) == LINEAR_COLUMNS
    k = out.set_index("condition")["k_deg"]
    assert k["control"] == pytest.approx(0.05, abs=0.01)
    assert k["atrium"] == pytest.approx(0.10, abs=0.02)
    # Δk = k_atrium − k_control ≈ +0.05 (control as reference), significant.
    dk = out["delta_k"].dropna().iloc[0]
    assert dk == pytest.approx(0.05, abs=0.02)
    assert (out["delta_k_p"].dropna() < 0.05).all()


def test_delta_k_sign_and_per_protein_bh():
    times = [0, 1, 2, 3, 4, 6, 8, 10]
    rows = []
    # P1: atrium faster (Δk>0). P2: equal (Δk≈0, not significant).
    rows += _curve("e", "P1", "control", 0.05, times, noise=0.004, seed=1)
    rows += _curve("e", "P1", "atrium", 0.12, times, noise=0.004, seed=2)
    rows += _curve("e", "P2", "control", 0.07, times, noise=0.004, seed=3)
    rows += _curve("e", "P2", "atrium", 0.07, times, noise=0.004, seed=4)
    out = fit_linear_deltak(pd.DataFrame(rows), reference_condition="control")
    by = out.dropna(subset=["delta_k"]).drop_duplicates("protein").set_index("protein")
    assert by.loc["P1", "delta_k"] > 0.02
    assert abs(by.loc["P2", "delta_k"]) < 0.02
    # BH-adjusted p present for both proteins; P1 stays significant.
    assert by["delta_k_p_adj"].notna().all()
    assert by.loc["P1", "delta_k_p_adj"] < 0.05


def test_single_condition_protein_has_no_delta_k():
    times = [0, 1, 2, 3, 4, 6, 8]
    out = fit_linear_deltak(pd.DataFrame(
        _curve("e", "P1", "control", 0.06, times, noise=0.003, seed=1)))
    assert (out["condition"] == "control").all()
    assert out["delta_k"].isna().all()  # no contrast without a second condition
    assert out["k_deg"].iloc[0] == pytest.approx(0.06, abs=0.01)


def test_named_pair_contrast_with_more_than_two_conditions():
    """Option B: a named reference/test pair contrasts exactly those two from the
    joint fit even when a 3rd condition is present (which still gets its k row)."""
    times = [0, 1, 2, 3, 4, 6, 8, 10]
    rows = (
        _curve("e", "P1", "control", 0.05, times, noise=0.004, seed=1)
        + _curve("e", "P1", "drugA", 0.10, times, noise=0.004, seed=2)
        + _curve("e", "P1", "drugB", 0.20, times, noise=0.004, seed=3)
    )
    out = fit_linear_deltak(pd.DataFrame(rows),
                            reference_condition="control", test_condition="drugA")
    # All three conditions still get a per-condition k.
    k = out.set_index("condition")["k_deg"]
    assert set(k.index) == {"control", "drugA", "drugB"}
    assert k["control"] == pytest.approx(0.05, abs=0.01)
    assert k["drugB"] == pytest.approx(0.20, abs=0.03)
    # The single Δk is the NAMED pair (drugA − control ≈ +0.05), not drugB.
    dk = out["delta_k"].dropna().unique()
    assert len(dk) == 1 and dk[0] == pytest.approx(0.05, abs=0.02)
    assert (out["delta_k_p"].dropna() < 0.05).all()


def test_named_pair_reference_sets_sign():
    """Swapping reference/test flips the Δk sign (Δk = k(test) − k(reference))."""
    times = [0, 1, 2, 3, 4, 6, 8, 10]
    df = pd.DataFrame(
        _curve("e", "P1", "control", 0.05, times, noise=0.004, seed=1)
        + _curve("e", "P1", "drug", 0.10, times, noise=0.004, seed=2)
        + _curve("e", "P1", "other", 0.03, times, noise=0.004, seed=3)
    )
    up = fit_linear_deltak(df, reference_condition="control", test_condition="drug")
    down = fit_linear_deltak(df, reference_condition="drug", test_condition="control")
    dk_up = up["delta_k"].dropna().iloc[0]
    dk_down = down["delta_k"].dropna().iloc[0]
    assert dk_up == pytest.approx(-dk_down, abs=1e-9)
    assert dk_up > 0  # k(drug) > k(control)


def test_named_pair_absent_in_a_protein_gives_nan():
    """A protein missing one of the named pair gets its k rows but no Δk."""
    times = [0, 1, 2, 3, 4, 6, 8]
    out = fit_linear_deltak(
        pd.DataFrame(
            _curve("e", "P1", "control", 0.05, times, noise=0.003, seed=1)
            + _curve("e", "P1", "drug", 0.10, times, noise=0.003, seed=2)
            + _curve("e", "P2", "control", 0.06, times, noise=0.003, seed=3)
            + _curve("e", "P2", "other", 0.06, times, noise=0.003, seed=4)
        ),
        reference_condition="control", test_condition="drug",
    )
    assert out[out["protein"] == "P1"]["delta_k"].notna().all()  # has the pair
    assert out[out["protein"] == "P2"]["delta_k"].isna().all()   # lacks 'drug'


def test_unknown_contrast_condition_raises():
    """A typo'd reference/test condition fails loudly rather than all-NaN silently."""
    from riana.exceptions import DataError
    times = [0, 1, 2, 3, 4, 6]
    df = pd.DataFrame(
        _curve("e", "P1", "control", 0.05, times)
        + _curve("e", "P1", "drug", 0.10, times)
    )
    with pytest.raises(DataError, match="not among the conditions"):
        fit_linear_deltak(df, reference_condition="control", test_condition="tpyo")
    with pytest.raises(DataError, match="not among the conditions"):
        fit_linear_deltak(df, reference_condition="nope")


def test_same_reference_and_test_condition_raises():
    """A self-contrast (reference == test) is rejected — it would degenerate into a
    spurious slope-vs-zero test (delta_k = −k, p ≈ 0), not a Δk between conditions."""
    from riana.exceptions import DataError
    times = [0, 1, 2, 3, 4, 6, 8]
    df = pd.DataFrame(
        _curve("e", "P1", "control", 0.05, times, noise=0.003, seed=1)
        + _curve("e", "P1", "drug", 0.10, times, noise=0.003, seed=2)
    )
    with pytest.raises(DataError, match="same"):
        fit_linear_deltak(df, reference_condition="control", test_condition="control")


def test_plateau_truncation_changes_fast_curve_slope():
    """A fast curve with a saturated tail: truncation must keep the slope honest
    rather than letting the floor-noise flatten it."""
    # k=0.5 saturates fast; later points sit at the θ ceiling.
    times = [0, 1, 2, 3, 5, 8, 12, 20, 30]
    pts = pd.DataFrame(_curve("e", "P1", "control", 0.5, times, seed=1))
    k_trunc = fit_linear_deltak(pts, phi_limit=-4.0).set_index("condition")["k_deg"]
    k_notrunc = fit_linear_deltak(
        pts, phi_limit=-50.0  # effectively no truncation
    ).set_index("condition")["k_deg"]
    # With the saturated tail dropped, k is closer to the true 0.5 and larger
    # than the flattened no-truncation estimate.
    assert k_trunc["control"] > k_notrunc["control"]
