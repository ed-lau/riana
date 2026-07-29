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


# --- weighting (WLS by default since 2026-07) --------------------------------
# phi = log(1-theta) is a LOG of FS-scale noise, so Var(phi) = sigma^2/(1-theta)^2:
# the phi-residuals are heteroscedastic and an unweighted fit is anti-conservative
# and biased low in the fast tail. See reports/2026-07-13_linear_model_wls.md.

def test_unknown_weight_scheme_raises():
    from riana.exceptions import DataError
    pts = pd.DataFrame(_curve("e", "P1", "control", 0.05, [1, 2, 4]))
    with pytest.raises(DataError, match="weights must be one of"):
        fit_linear_deltak(pts, weights="bogus")


def _two_cond(seed=3, **kw):
    return pd.DataFrame(
        _curve("e", "P1", "control", 0.08, [0, 1, 2, 4, 8], noise=0.02, seed=seed)
        + _curve("e", "P1", "atrium", 0.12, [0, 1, 2, 4, 8], noise=0.02, seed=seed + 1))


def test_wls_var_falls_back_to_wls_without_theta_var():
    """wls-var on a points table with no per-point variance must reproduce wls exactly."""
    pts = _two_cond()
    a = fit_linear_deltak(pts, weights="wls", reference_condition="control")
    b = fit_linear_deltak(pts, weights="wls-var", reference_condition="control")
    for col in ("k_deg", "delta_k"):
        np.testing.assert_allclose(a[col].to_numpy(), b[col].to_numpy(), equal_nan=True)


def test_wls_var_equals_wls_when_variance_is_constant():
    """A constant Var̂(θ) cancels out of the RELATIVE weights, so wls-var ≡ wls."""
    pts = _two_cond()
    pts["theta_var"] = 4e-4        # constant across all points
    pts["theta_df"] = 8.0
    a = fit_linear_deltak(pts, weights="wls", reference_condition="control")
    b = fit_linear_deltak(pts, weights="wls-var", reference_condition="control")
    np.testing.assert_allclose(a["k_deg"].to_numpy(), b["k_deg"].to_numpy(), rtol=1e-9)


def test_wls_var_reweights_under_heteroscedastic_variance():
    """With genuinely varying per-point variance, wls-var re-weights → k differs from wls
    (but tracks it); a tiny-sample fitFDist falls back to the fixed d0."""
    rng = np.random.default_rng(7)
    rows, var = [], []
    for i in range(40):
        for cond, k in (("control", 0.08), ("atrium", 0.12)):
            for t in [1, 2, 4, 8, 12]:
                sd = 0.06 if rng.random() < 0.5 else 0.012      # half the points 5x noisier
                rows.append({"experiment": "e", "protein": f"P{i}", "condition": cond,
                             "labeling_time": float(t), "theta": float(1 - np.exp(-k * t) + rng.normal(0, sd))})
                var.append(sd ** 2)
    pts = pd.DataFrame(rows); pts["theta_var"] = var; pts["theta_df"] = 6.0
    a = fit_linear_deltak(pts, weights="wls", reference_condition="control").set_index("protein")["k_deg"]
    b = fit_linear_deltak(pts, weights="wls-var", reference_condition="control").set_index("protein")["k_deg"]
    assert (a - b).abs().median() > 0        # the per-point weighting bites
    assert a.corr(b) > 0.9                    # but the two track each other


def test_fit_fdist_fallback_on_tiny_sample():
    from riana.core.linear_model import _fit_fdist
    d0, s0 = _fit_fdist(np.array([1e-3, 2e-3]), np.array([6.0, 6.0]), d0_fallback=2.0)
    assert d0 == 2.0 and s0 > 0               # < 8 variances → the fixed-2 fallback


def test_t0_point_is_excluded_and_does_not_move_k():
    """t=0 has ZERO leverage on a through-origin slope, so dropping it cannot change
    k — but its theta is pinned by the floor clamp, so its residual is artificially ~0
    and would deflate the residual variance. It must be excluded from the fit."""
    with_t0 = pd.DataFrame(
        _curve("e", "P1", "control", 0.08, [0, 1, 2, 4, 8], noise=0.02, seed=3))
    # drop the t=0 ROW so every other point is bit-identical (not a re-draw)
    no_t0 = with_t0[with_t0["labeling_time"] > 0].reset_index(drop=True)
    a = fit_linear_deltak(with_t0)
    b = fit_linear_deltak(no_t0)
    # identical k (t=0 has zero leverage on a through-origin slope) ...
    assert a["k_deg"].iloc[0] == pytest.approx(b["k_deg"].iloc[0], rel=1e-9)
    # ... and t=0 is not counted among the fitted points
    assert int(a["n_points"].iloc[0]) == len(no_t0)


def _fast_curve_panel(k_true, n_rep, times, sigma, seed):
    """n_rep synthetic proteins, each a two-condition curve at the SAME k, with
    homoscedastic noise on the FS scale (the real error structure)."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_rep):
        for cond in ("control", "atrium"):
            for t in times:
                th = 1.0 - np.exp(-k_true * t) + rng.normal(0, sigma)
                rows.append({"experiment": "e", "protein": f"P{i}",
                             "condition": cond, "labeling_time": float(t),
                             "theta": float(th)})
    return pd.DataFrame(rows)


def test_wls_removes_the_fast_tail_bias_that_ols_has():
    """Under FS-scale noise the unweighted fit reads a FAST curve's k low (~-14% on
    real data); the delta-method weighting is essentially unbiased."""
    k_true, times = 0.30, [1, 2, 3, 4, 6, 8, 10, 15, 20, 25, 30]
    pts = _fast_curve_panel(k_true, 120, times, sigma=0.056, seed=17)
    k_ols = fit_linear_deltak(pts, weights="ols")["k_deg"].median()
    k_wls = fit_linear_deltak(pts, weights="wls")["k_deg"].median()
    assert k_ols < 0.9 * k_true          # OLS is biased LOW on a fast curve
    assert k_wls == pytest.approx(k_true, rel=0.10)   # WLS recovers it
    assert abs(k_wls - k_true) < abs(k_ols - k_true)  # and is strictly better


def test_wls_is_the_default():
    k_true, times = 0.30, [1, 2, 3, 4, 6, 8, 10, 15, 20, 25, 30]
    pts = _fast_curve_panel(k_true, 60, times, sigma=0.056, seed=23)
    default = fit_linear_deltak(pts)["k_deg"].median()
    wls = fit_linear_deltak(pts, weights="wls")["k_deg"].median()
    ols = fit_linear_deltak(pts, weights="ols")["k_deg"].median()
    assert default == pytest.approx(wls, rel=1e-9)
    assert default != pytest.approx(ols, rel=1e-6)


def test_ols_keeps_t0_while_wls_drops_it():
    """Regression: ``weights="ols"`` must reproduce the pre-2026-07 fit, which KEPT
    the t=0 point; only ``wls`` drops it (the clamped ≈0 residual at t=0 deflates σ̂²).
    t=0 has zero leverage on a through-origin slope, so k is unchanged either way — the
    difference is purely which points enter the fit. Plateau truncation is independent
    of the scheme (see ``test_truncate_plateau_drops_saturated_tail``)."""
    times = [0, 1, 2, 3, 4, 6, 8, 10]     # includes t=0
    pts = pd.DataFrame(
        _curve("e", "P", "control", 0.05, times, noise=0.01, seed=1)
        + _curve("e", "P", "atrium", 0.10, times, noise=0.01, seed=2))
    ols = fit_linear_deltak(pts, weights="ols", reference_condition="control")
    wls = fit_linear_deltak(pts, weights="wls", reference_condition="control")
    n_ols = int(ols["n_points"].iloc[0])
    n_wls = int(wls["n_points"].iloc[0])
    assert n_ols > n_wls                                     # ols retains t=0, wls drops it
    # k is unchanged (t=0 has no leverage on a through-origin slope)
    assert ols["k_deg"].iloc[0] == pytest.approx(wls["k_deg"].iloc[0], abs=5e-3)
