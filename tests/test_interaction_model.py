"""Tests for the experiment × condition interaction (random-animal mixed model)."""

import numpy as np
import pandas as pd
import pytest

import riana.core.interaction_model as interaction_model
from riana.core.interaction_model import (
    INTERACTION_COLUMNS,
    containment_df,
    fit_interaction,
    satterthwaite_ddf,
)


def split_plot_points(
    k_ref: dict[str, float],
    test_ratio: dict[str, float],
    *,
    proteins=("P1",),
    days=(1.0, 3.0, 7.0, 10.0, 14.0),
    n_rep: int = 3,
    sd_animal: float = 0.006,
    sd_theta: float = 0.02,
    seed: int = 0,
) -> pd.DataFrame:
    """Synthetic split-plot points: each animal (condition × day × replicate) is
    sampled in every experiment, with one random slope shared by its experiments
    and FS-scale noise. ``test_ratio[s]`` = k_test / k_reference in experiment s."""
    rng = np.random.default_rng(seed)
    rows = []
    for prot in proteins:
        for cond in ("control", "KO"):
            for t in days:
                for br in range(1, n_rep + 1):
                    u = rng.normal(0.0, sd_animal)
                    for s, k0 in k_ref.items():
                        k = k0 * (test_ratio[s] if cond == "KO" else 1.0) + u
                        theta = 1.0 - np.exp(-k * t) + rng.normal(0.0, sd_theta)
                        rows.append(dict(protein=prot, experiment=s, condition=cond,
                                         biological_replicate=br, labeling_time=t,
                                         theta=theta))
    return pd.DataFrame(rows)


# The dataset the nlme reference values below were computed on, and how:
#   d = split_plot_points(PARITY_K, PARITY_RATIO, seed=7)   (no point is filtered out)
#   R (nlme 3.1.169): design columns as in riana.core.interaction_model (strata sorted,
#   contr.sum), w = exp(2 * pmax(fitted(lm(phi ~ 0 + X)), -4)),
#   m <- lme(phi ~ 0 + X, random = ~ 0 + day | animal, weights = varFixed(~ 1/w))
PARITY_K = {"A": 0.10, "B": 0.06, "C": 0.08}
PARITY_RATIO = {"A": 0.70, "B": 1.00, "C": 0.90}
NLME = {
    "ref_slope": [-0.09859358796, -0.05712798687, -0.07702201042],
    "shift": 0.01201385383, "shift_se": 0.002360497819, "shift_df": 29,
    "add_F": 88.98039809, "add_df2": 56,
    "sigma_animal": 0.005234727122,
}
# Satterthwaite denominator df on the same dataset, from R (lme4 2.0.1, lmerTest 3.2.1):
#   lmer(phi*sw ~ 0 + X*sw + (0 + day*sw | animal), REML = TRUE)  (= the √w-scaled fit)
#   contestMD(m, L) with L = the ratio test's delta-method contrast (2 × 6), and with
#   L = the shift column alone (1 × 6, the one-df branch).
# Script: runs/interaction_satterthwaite/{export.py, lmertest_ref.R} (gitignored).
LMERTEST = {"ratio_ddf": 66.280716, "shift_ddf": 23.139048}


def _protein(out: pd.DataFrame, prot: str = "P1") -> pd.Series:
    return out[out.protein == prot].iloc[0]


def test_containment_df_matches_nlme_rule():
    """nlme getFixDF: between-animal columns n_animals − p_between; within-animal
    N − n_animals − p_within + 1 (no intercept). ATF6 layout: 30 animals × 4."""
    rng = np.random.default_rng(0)
    animal = np.repeat(np.arange(30), 4)
    day = np.repeat(rng.choice([1.0, 3.0, 7.0, 10.0, 14.0], 30), 4)
    ko = np.repeat(np.arange(30) % 2, 4).astype(float)
    s = np.tile(np.arange(4), 30)
    C = np.vstack([np.eye(3), -np.ones(3)])
    X = np.column_stack([day * (s == j) for j in range(4)] + [day * ko]
                        + [day * ko * C[s, j] for j in range(3)])
    assert containment_df(X, animal).tolist() == [84, 84, 84, 84, 29, 84, 84, 84]


def test_matches_nlme_reference():
    """The fit reproduces nlme lme(random = ~0 + day | animal, weights = varFixed)
    on the stored dataset — estimates, SE, df, the additive interaction F, σ_animal."""
    d = split_plot_points(PARITY_K, PARITY_RATIO, seed=7)
    out = fit_interaction(d, reference_condition="control", test_condition="KO")
    assert list(out.columns) == INTERACTION_COLUMNS
    assert out.experiment.tolist() == ["A", "B", "C"]
    k_ref = out.k_reference.to_numpy()
    assert k_ref == pytest.approx([-v for v in NLME["ref_slope"]], rel=2e-3)
    row = _protein(out)
    assert row.mean_delta_k == pytest.approx(-NLME["shift"], rel=2e-3)
    assert row.mean_delta_k_se == pytest.approx(NLME["shift_se"], rel=5e-2)
    assert row.mean_delta_k_df == NLME["shift_df"]
    assert row.interaction_additive_F == pytest.approx(NLME["add_F"], rel=2e-2)
    assert row.interaction_additive_df2 == NLME["add_df2"]
    assert row.sigma_animal == pytest.approx(NLME["sigma_animal"], rel=1e-2)
    # ratio-scale test: S − 1 numerator df, Satterthwaite denominator df (lmerTest)
    assert row.interaction_df1 == 2
    assert row.interaction_df2 == pytest.approx(LMERTEST["ratio_ddf"], rel=1e-4)
    assert row.n_experiments == 3 and row.n_animals == 30
    # per-experiment Δk and ratio are consistent with k_reference / k_test
    assert out.delta_k.to_numpy() == pytest.approx(
        (out.k_test - out.k_reference).to_numpy(), abs=1e-12)
    assert out.k_ratio.to_numpy() == pytest.approx(
        (out.k_test / out.k_reference).to_numpy(), rel=1e-12)


def test_uniform_proportional_effect_is_not_a_ratio_interaction():
    """Test condition slows every experiment by the same 30 %: the absolute Δk
    differs (fast vs slow experiments) but the ratio does not."""
    d = split_plot_points({"fast": 0.15, "slow": 0.04, "mid": 0.08},
                          {"fast": 0.7, "slow": 0.7, "mid": 0.7},
                          sd_theta=0.01, seed=3)
    row = _protein(fit_interaction(d, reference_condition="control", test_condition="KO"))
    assert row.interaction_additive_p < 1e-3
    assert row.interaction_p > 0.05


def test_experiment_specific_effect_is_a_ratio_interaction():
    d = split_plot_points({"A": 0.08, "B": 0.08, "C": 0.08},
                          {"A": 0.6, "B": 1.0, "C": 1.0}, sd_theta=0.01, seed=4)
    out = fit_interaction(d, reference_condition="control", test_condition="KO")
    row = _protein(out)
    assert row.interaction_p < 1e-3
    ratio = out.set_index("experiment").k_ratio
    assert ratio["A"] == pytest.approx(0.6, abs=0.05)
    assert ratio["B"] == pytest.approx(1.0, abs=0.05)


def test_unpaired_experiments_are_dropped_and_single_experiment_proteins_skipped():
    d = split_plot_points({"A": 0.08, "B": 0.06, "C": 0.07},
                          {"A": 0.8, "B": 0.8, "C": 0.8}, proteins=("P1", "P2", "P3"),
                          seed=5)
    # P2: experiment C has no test-condition points -> fit on A and B only.
    d = d[~((d.protein == "P2") & (d.experiment == "C") & (d.condition == "KO"))]
    # P3: observed in one experiment only -> nothing to compare, no rows.
    d = d[~((d.protein == "P3") & (d.experiment != "A"))]
    out = fit_interaction(d, reference_condition="control", test_condition="KO")
    assert set(out.protein) == {"P1", "P2"}
    assert sorted(out[out.protein == "P2"].experiment) == ["A", "B"]
    assert _protein(out, "P2").n_experiments == 2
    assert _protein(out, "P2").interaction_df1 == 1
    # one BH family per test: every tested protein gets an adjusted p
    per = out.drop_duplicates("protein")
    assert per.interaction_p_adj.notna().all()
    assert per.interaction_additive_p_adj.notna().all()
    assert per.mean_delta_k_p_adj.notna().all()


def test_other_conditions_are_ignored():
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.8}, seed=6)
    extra = d[d.condition == "KO"].assign(condition="drug")
    a = fit_interaction(d, reference_condition="control", test_condition="KO")
    b = fit_interaction(pd.concat([d, extra]), reference_condition="control",
                        test_condition="KO")
    pd.testing.assert_frame_equal(a, b)


def test_workers_give_identical_results():
    d = split_plot_points({"A": 0.08, "B": 0.06, "C": 0.07},
                          {"A": 0.7, "B": 0.9, "C": 0.8},
                          proteins=("P1", "P2", "P3", "P4"), seed=8)
    serial = fit_interaction(d, reference_condition="control", test_condition="KO")
    pooled = fit_interaction(d, reference_condition="control", test_condition="KO",
                             workers=2)
    pd.testing.assert_frame_equal(serial, pooled)


def test_input_errors():
    from riana.exceptions import DataError
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.8}, seed=1)
    with pytest.raises(DataError, match="missing columns"):
        fit_interaction(d.drop(columns="biological_replicate"),
                        reference_condition="control", test_condition="KO")
    with pytest.raises(DataError, match="same"):
        fit_interaction(d, reference_condition="KO", test_condition="KO")
    with pytest.raises(DataError, match="not among the conditions"):
        fit_interaction(d, reference_condition="control", test_condition="drug")


def test_no_qualifying_protein_returns_empty_schema():
    d = split_plot_points({"A": 0.08}, {"A": 0.8}, seed=2)       # one experiment only
    out = fit_interaction(d, reference_condition="control", test_condition="KO")
    assert out.empty and list(out.columns) == INTERACTION_COLUMNS


def test_satterthwaite_one_df_branch_matches_lmertest(monkeypatch):
    """The 1-df branch (a single contrast row), on the inputs the ratio test used."""
    seen = []
    orig = interaction_model.satterthwaite_ddf

    def spy(L, X, y, z, groups, s2a, s2e):
        seen.append((X, y, z, groups, s2a, s2e))
        return orig(L, X, y, z, groups, s2a, s2e)

    monkeypatch.setattr(interaction_model, "satterthwaite_ddf", spy)
    fit_interaction(split_plot_points(PARITY_K, PARITY_RATIO, seed=7),
                    reference_condition="control", test_condition="KO")
    X, y, z, groups, s2a, s2e = seen[0]
    L = np.zeros((1, X.shape[1]))
    L[0, 3] = 1.0                                     # the shift column (S = 3)
    assert orig(L, X, y, z, groups, s2a, s2e) == pytest.approx(
        LMERTEST["shift_ddf"], rel=1e-4)


def test_satterthwaite_undefined_returns_nan():
    """No variance at all: Σ is singular, so there is no df (the fit then falls
    back to the between-animal containment df)."""
    X = np.column_stack([np.ones(6), np.arange(6.0)])
    ddf = satterthwaite_ddf(np.eye(2), X, np.arange(6.0), np.ones(6),
                            np.repeat([0, 1, 2], 2), 0.0, 0.0)
    assert np.isnan(ddf)


def _with_individual(d: pd.DataFrame, fmt: str) -> pd.DataFrame:
    return d.assign(individual=[
        fmt.format(c=c, t=t, r=r) for c, t, r in
        zip(d.condition, d.labeling_time, d.biological_replicate)])


def test_individual_one_to_one_with_composite_gives_identical_fit():
    d = split_plot_points(PARITY_K, PARITY_RATIO, seed=7)
    base = fit_interaction(d, reference_condition="control", test_condition="KO")
    named = fit_interaction(_with_individual(d, "{c}_d{t:g}_r{r}"),
                            reference_condition="control", test_condition="KO")
    pd.testing.assert_frame_equal(base, named)
    assert base.attrs["animal_key"] == "composite"
    assert named.attrs["animal_key"] == "individual"


def test_individual_ids_restarting_per_condition_are_flagged(caplog):
    """IDs 1..n reused in both conditions merge two animals: warn."""
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.8}, seed=6)
    with caplog.at_level("WARNING", logger="riana.core.interaction_model"):
        out = fit_interaction(_with_individual(d, "d{t:g}_r{r}"),
                              reference_condition="control", test_condition="KO")
    assert "appear under both conditions" in caplog.text
    assert _protein(out).n_animals == 15               # 30 animals collapsed to 15


def test_longitudinal_individuals_fit_as_one_animal_over_days(caplog):
    """An animal sampled at every day (serial sampling): one random slope across its
    days, so fewer, larger animal groups."""
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.7, "B": 0.9}, seed=9)
    with caplog.at_level("INFO", logger="riana.core.interaction_model"):
        out = fit_interaction(_with_individual(d, "{c}_m{r}"),
                              reference_condition="control", test_condition="KO")
    assert "longitudinal" in caplog.text
    row = _protein(out)
    assert row.n_animals == 6
    assert np.isfinite(row.interaction_p)


def test_partial_individuals_fall_back_to_composite(caplog):
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.8}, seed=6)
    named = _with_individual(d, "{c}_d{t:g}_r{r}")
    named.loc[named.labeling_time == 1.0, "individual"] = np.nan
    with caplog.at_level("WARNING", logger="riana.core.interaction_model"):
        out = fit_interaction(named, reference_condition="control", test_condition="KO")
    assert "have no individual" in caplog.text
    assert out.attrs["animal_key"] == "individual+composite"
    assert _protein(out).n_animals == 30


def test_singular_fit_is_kept_as_the_weighted_least_squares_solution():
    """REML optimum at σ²_animal = 0: statsmodels stops there unconverged with wrong
    fixed effects. The protein is kept with the singular fit lme4 returns — the
    per-experiment WLS slopes (= the shipped linear model's k) — and the ratio test's
    Satterthwaite df is N − p (lmerTest at a singular fit)."""
    from riana.core.linear_model import fit_linear_deltak
    d = split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.9},
                          sd_animal=0.0, seed=0)
    out = fit_interaction(d, reference_condition="control", test_condition="KO")
    row = _protein(out)
    assert row.sigma_animal == 0.0
    assert row.interaction_df2 == pytest.approx(out.n_points.sum() - 4)   # N − 2S
    lin = fit_linear_deltak(d, reference_condition="control").set_index(
        ["experiment", "condition"]).k_deg
    for _, r in out.iterrows():
        assert r.k_reference == pytest.approx(lin[(r.experiment, "control")], rel=1e-9)
        assert r.k_test == pytest.approx(lin[(r.experiment, "KO")], rel=1e-9)


# --- weights: the interaction follows the per-experiment fit's --linear-weights ---------


def _with_var(d: pd.DataFrame, seed: int = 3) -> pd.DataFrame:
    """Per-point variances spread over an order of magnitude, as the collapse writes."""
    rng = np.random.default_rng(seed)
    return d.assign(theta_var=10 ** rng.uniform(-4.5, -3.5, len(d)),
                    theta_df=rng.uniform(4.0, 20.0, len(d)))


def test_wls_var_weights_match_the_per_experiment_fit():
    """At a singular fit the joint model reduces to the per-experiment weighted fits, so
    its slopes must equal fit_linear_deltak's under the same scheme — same moderation
    prior, same per-point weights (1 − θ̂)²/Ṽar."""
    from riana.core.linear_model import fit_linear_deltak
    d = _with_var(split_plot_points({"A": 0.08, "B": 0.06}, {"A": 0.8, "B": 0.9},
                                    sd_animal=0.0, seed=0))
    out = fit_interaction(d, reference_condition="control", test_condition="KO",
                          weights="wls-var")
    assert (out.sigma_animal == 0.0).all()
    lin = fit_linear_deltak(d, reference_condition="control", weights="wls-var").set_index(
        ["experiment", "condition"]).k_deg
    wls = fit_linear_deltak(d, reference_condition="control").set_index(
        ["experiment", "condition"]).k_deg
    assert not np.allclose(lin.to_numpy(), wls.to_numpy(), rtol=1e-6)   # weights matter here
    for _, r in out.iterrows():
        assert r.k_reference == pytest.approx(lin[(r.experiment, "control")], rel=1e-9)
        assert r.k_test == pytest.approx(lin[(r.experiment, "KO")], rel=1e-9)


def test_wls_var_discounts_imprecise_points():
    """Biased points that carry a large variance (the early points of a precursor-corrected
    series) pull the wls fit but not the wls-var fit."""
    d = split_plot_points(PARITY_K, PARITY_RATIO, seed=7)
    rng = np.random.default_rng(5)
    d = d.assign(theta_var=1e-4 * rng.uniform(0.8, 1.2, len(d)), theta_df=10.0)
    clean = _protein(fit_interaction(d, reference_condition="control",
                                     test_condition="KO", weights="wls-var"))
    bad = (d.condition == "KO") & (d.labeling_time <= 3.0)
    dc = d.copy()
    dc.loc[bad, "theta"] += 0.3
    dc.loc[bad, "theta_var"] = 0.1
    wls, var = (_protein(fit_interaction(dc, reference_condition="control",
                                         test_condition="KO", weights=w))
                for w in ("wls", "wls-var"))
    assert abs(wls.mean_delta_k - clean.mean_delta_k) > 0.005
    assert abs(var.mean_delta_k - clean.mean_delta_k) < 0.1 * abs(wls.mean_delta_k
                                                                  - clean.mean_delta_k)


def test_wls_var_without_variances_and_ols_fall_back_to_wls(caplog):
    d = split_plot_points(PARITY_K, PARITY_RATIO, seed=7)
    base = fit_interaction(d, reference_condition="control", test_condition="KO")
    with caplog.at_level("INFO", logger="riana.core.interaction_model"):
        no_var = fit_interaction(d, reference_condition="control", test_condition="KO",
                                 weights="wls-var")
        ols = fit_interaction(d, reference_condition="control", test_condition="KO",
                              weights="ols")
    pd.testing.assert_frame_equal(base, no_var)
    pd.testing.assert_frame_equal(base, ols)
    assert "no theta_var/theta_df" in caplog.text and "audit path" in caplog.text


def test_wls_var_workers_and_bad_scheme():
    from riana.exceptions import DataError
    d = _with_var(split_plot_points({"A": 0.08, "B": 0.06, "C": 0.07},
                                    {"A": 0.7, "B": 0.9, "C": 0.8},
                                    proteins=("P1", "P2", "P3"), seed=8))
    serial = fit_interaction(d, reference_condition="control", test_condition="KO",
                             weights="wls-var")
    pooled = fit_interaction(d, reference_condition="control", test_condition="KO",
                             weights="wls-var", workers=2)
    pd.testing.assert_frame_equal(serial, pooled)
    with pytest.raises(DataError, match="weights"):
        fit_interaction(d, reference_condition="control", test_condition="KO",
                        weights="gls")
