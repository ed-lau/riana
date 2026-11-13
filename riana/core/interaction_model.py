"""Two-factor (experiment × condition) interaction for ``linear simple``: a random-animal
mixed model.

When the experiments are strata of the **same animals** (e.g. four heart chambers sampled
from every mouse, split with ``integrate --experiment-column``), the per-experiment Δk in
``riana_rollup_proteins.txt`` answers "does the condition change k in this tissue?" but not
"does the condition affect the tissues differently?". That needs one joint fit across the
experiments. Each animal then contributes one point per experiment, so the fit is a
split-plot: condition and labeling time vary between animals, experiment within them.
Without an animal term the condition effect is pseudo-replicated and the interaction
under-powered (on ATF6, 186 → 1 and 3 → 132 BH hits;
``reports/2026-11-10_random_animal_lme_feasibility.md``).

**Model**, per protein, on the same collapsed points, φ = log(1 − θ) transform and point
filters (plateau, t > 0) as :func:`riana.core.linear_model.fit_linear_deltak`::

    φ = Σ_s β_s·day·[s] + β_shift·day·[test] + Σ_j γ_j·day·[test]·C[s, j] + u_animal·day + ε

- ``β_s`` — the reference-condition slope in experiment *s* (k = −slope);
- ``β_shift`` — the test-condition slope shift averaged over experiments (C = sum contrasts);
- ``γ`` — how that shift differs between experiments (the interaction);
- ``u_animal ~ N(0, σ_animal²)`` — a per-animal slope (nlme ``random = ~0 + day | animal``);
- ``Var(ε) ∝ 1/w``, the weights of the per-experiment linear fit (``--linear-weights``), from a
  first OLS pass as there:

  - ``wls`` (default): ``w = (1 − θ̂)²``, the delta-method weight (nlme
    ``varExp(form = ~fitted(.))`` with δ fixed at −1; ATF6 estimates δ ≈ −1.1);
  - ``wls-var``: ``w = (1 − θ̂)²/Ṽar(θ)``, adding each collapsed point's precision — the
    collapse's ``theta_var`` moderated toward a pooled prior exactly as
    :func:`riana.core.linear_model.fit_linear_deltak` does (robust-fitFDist ``d0`` and
    ``s0²``, estimated once over every point). Opt-in, as for the per-experiment fit. Not
    for the rising-precursor fraction new f: there the per-point precision spans ~10³ over a
    series and these weights make the tests anti-conservative (Monte Carlo,
    ``reports/2026-11-12_deltak_inference_manual.md`` §6.2); ``wls`` is calibrated there.

**Animal** = the ``individual`` column (SDRF ``characteristics[individual]``, threaded
through fit and rollup) where present. Otherwise the (condition, labeling_time,
biological_replicate) composite, which assumes ``biological_replicate`` names the same animal
in every experiment of one condition × timepoint (as an SDRF stratified by tissue does when
each animal keeps one replicate index across its tissues) and that each animal is sampled at
one timepoint. Real IDs also allow longitudinal designs (one animal at several timepoints:
the random slope is then a genuine within-animal slope). Experiments that are *different*
animals (e.g. separate cohorts) do not fit this model.

**Fit.** REML with statsmodels ``MixedLM`` on the √w-transformed data, which is nlme
``lme(..., weights = varFixed(~1/w))``. BFGS, not statsmodels' default L-BFGS: on these small
variances (~1e-5) L-BFGS stops at a wrong point while reporting convergence. Checked against
nlme on ATF6 (estimates within 1.5e-3, F within 1.3 %, df identical). When the REML optimum
is σ²_animal = 0 (a *singular* fit; ~1 % of ATF6 proteins) statsmodels stops there
unconverged with wrong fixed effects, so that fit is computed directly (:func:`_singular_fit`:
the weighted least-squares solution, as lme4 returns); ``sigma_animal`` is then 0.

**Tests reported per protein** (BH across proteins, one family per test):

- ``interaction_*`` (ratio scale, primary) — is ``k_test / k_reference`` the same in every
  experiment? A uniform proportional effect is not an interaction on this scale, though it
  gives a larger absolute Δk in faster-turnover experiments. Wald F on the log ratios (delta
  method). The contrast mixes between-animal (``β_shift``) and within-animal coefficients, so
  no containment df fits it; its denominator df is Satterthwaite's
  (:func:`satterthwaite_ddf`, lmerTest's ``contestMD`` method; checked against lmerTest
  3.2.1). Where that is undefined the smaller, between-animal containment df is used.
- ``interaction_additive_*`` — are the Δk equal across experiments (``γ = 0``)? nlme's
  containment df (:func:`containment_df`, exact within-animal df); the test nlme's ``anova``
  gives for the interaction columns.
- ``mean_delta_k_*`` — Δk averaged over experiments (= −β_shift), the condition main effect,
  at the between-animal containment df. Per-experiment Δk *tests* stay in
  ``riana_rollup_proteins.txt``.
"""

from __future__ import annotations

import logging
import warnings
from concurrent import futures
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats

from riana.core.linear_model import LINEAR_WEIGHT_SCHEMES, _fit_fdist, to_phi
from riana.progress import iter_progress

_LOGGER = logging.getLogger(__name__)

#: Output columns of :func:`fit_interaction`: one row per ``(protein, experiment)``. The
#: per-experiment estimates come from the joint model; the remaining columns are
#: protein-level and repeat on each of the protein's rows.
INTERACTION_COLUMNS = [
    "protein", "experiment", "n_points",
    "k_reference", "k_test", "delta_k", "delta_k_se", "k_ratio",
    "n_experiments", "n_animals", "sigma_animal",
    "mean_delta_k", "mean_delta_k_se", "mean_delta_k_df", "mean_delta_k_p",
    "mean_delta_k_p_adj",
    "interaction_F", "interaction_df1", "interaction_df2", "interaction_p",
    "interaction_p_adj",
    "interaction_additive_F", "interaction_additive_df2", "interaction_additive_p",
    "interaction_additive_p_adj",
]

_ADJUSTED = {
    "mean_delta_k_p_adj": "mean_delta_k_p",
    "interaction_p_adj": "interaction_p",
    "interaction_additive_p_adj": "interaction_additive_p",
}


def containment_df(X: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Denominator df per fixed-effect column, by nlme's containment rule (``getFixDF``),
    for one grouping level and a model without an intercept.

    A column constant within every group is a between-group term and gets
    ``n_groups − p_between``; any other column is within-group and gets
    ``N − n_groups − p_within + 1`` (nlme adds the 1 when there is no intercept).
    """
    X = np.asarray(X, dtype=float)
    codes, uniq = pd.factorize(pd.Series(groups))
    grouped = pd.DataFrame(X).groupby(codes)
    spread = (grouped.max() - grouped.min()).to_numpy().max(axis=0)
    between = np.isclose(spread, 0.0, atol=1e-12)
    n_groups = len(uniq)
    p_between, p_within = int(between.sum()), int((~between).sum())
    return np.where(between, n_groups - p_between,
                    len(X) - n_groups - p_within + 1).astype(float)


def satterthwaite_ddf(
    L: np.ndarray,
    X: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    groups: np.ndarray,
    sigma2_animal: float,
    sigma2_resid: float,
) -> float:
    """Satterthwaite denominator df for the F test of ``L β`` in a REML fit of
    ``y = Xβ + z·u_group + ε`` (``u ~ N(0, σ²_animal)``, ``ε ~ N(0, σ²_resid·I)``).

    lmerTest's ``contestMD`` construction: with ``L V Lᵀ = P diag(d) Pᵀ`` and
    ``V = (Xᵀ Σ⁻¹ X)⁻¹``, each eigen-contrast ``l_m = P_mᵀ L`` gets
    ``ν_m = 2 d_m² / (g_mᵀ A g_m)``, where ``g_m`` is the gradient of ``l_m V l_mᵀ``
    in the variance parameters and ``A`` the inverse observed REML information; the
    ν_m combine as ``2E / (E − q)`` with ``E = Σ ν_m / (ν_m − 2)`` (2 if any ν_m ≤ 2).
    At ``σ²_animal = 0`` (a singular fit) only ``σ²_resid`` counts, giving ``N − p``
    as lmerTest does. Returns ``nan`` where the information is not positive definite.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if not (sigma2_resid > 0 and sigma2_animal >= 0):
        return float("nan")
    # Σ is block-diagonal by animal, so everything is accumulated per block (each a
    # few points) instead of with dense N × N matrices. Per block, with D_a = z zᵀ and
    # D_e = I (∂Σ/∂σ²_animal, ∂Σ/∂σ²_resid), Σ⁻¹ by Sherman–Morrison.
    codes = pd.factorize(pd.Series(groups))[0]
    order = np.argsort(codes, kind="stable")
    bounds = np.flatnonzero(np.diff(codes[order])) + 1
    blocks = []
    p = X.shape[1]
    XtSX = np.zeros((p, p))
    XtSy = np.zeros(p)
    for idx in np.split(order, bounds):
        Xg, yg, zg = X[idx], y[idx], z[idx]
        c = sigma2_animal / (sigma2_resid + sigma2_animal * (zg @ zg))
        Sg = (np.eye(len(idx)) - c * np.outer(zg, zg)) / sigma2_resid
        SXg = Sg @ Xg
        XtSX += Xg.T @ SXg
        XtSy += SXg.T @ yg
        blocks.append((Xg, yg, Sg, SXg, (np.outer(zg, zg), np.eye(len(idx)))))
    try:
        V = np.linalg.inv(XtSX)
    except np.linalg.LinAlgError:
        return float("nan")
    b = V @ XtSy
    # P = Σ⁻¹ − Σ⁻¹X V XᵀΣ⁻¹.  tr(PΣ_kPΣ_l) = Σ_g tr(A_k A_l) − 2 tr(V M_kl)
    # + tr(V B_k V B_l), with A_k = Σ_g⁻¹D_k, M_kl = Σ XᵀA_kA_lΣ⁻¹X, B_k = XᵀΣ⁻¹D_kΣ⁻¹X.
    T1 = np.zeros((2, 2))
    M = np.zeros((2, 2, p, p))
    B = np.zeros((2, p, p))
    w = [[], []]                                   # D_k · Py, per block
    XtSw = np.zeros((2, p))
    for Xg, yg, Sg, SXg, D in blocks:
        rg = Sg @ (yg - Xg @ b)                    # Py on this block
        A = [Sg @ Dk for Dk in D]
        for k in range(2):
            B[k] += SXg.T @ D[k] @ SXg
            wk = D[k] @ rg
            w[k].append(wk)
            XtSw[k] += SXg.T @ wk
            for l in range(2):
                T1[k, l] += np.trace(A[k] @ A[l])
                M[k, l] += Xg.T @ A[k] @ A[l] @ SXg
    # At a singular fit (σ²_animal = 0) only σ²_resid is estimated (lmerTest), which
    # makes every ν_m = N − p.
    par = [0, 1] if sigma2_animal > 0 else [1]
    info = np.empty((len(par), len(par)))
    for i, k in enumerate(par):
        for j, l in enumerate(par):
            trPP = T1[k, l] - 2 * np.trace(V @ M[k, l]) + np.trace(V @ B[k] @ V @ B[l])
            # yᵀPΣ_kPΣ_lPy = (D_k Py)ᵀ P (D_l Py)
            corr = V @ XtSw[l]
            quad = sum(wk @ (Sg @ wl - SXg @ corr)
                       for wk, wl, (_, _, Sg, SXg, _) in zip(w[k], w[l], blocks))
            info[i, j] = -0.5 * trPP + quad
    if not np.all(np.isfinite(info)) or np.linalg.eigvalsh(info).min() <= 0:
        return float("nan")
    A_cov = np.linalg.inv(info)
    dV = [V @ B[k] @ V for k in par]                    # ∂V/∂θ_k
    d, Pe = np.linalg.eigh(L @ V @ L.T)
    keep = d > max(np.finfo(float).eps * d.max(), 0.0)
    nu = []
    for m in np.flatnonzero(keep):
        lm = Pe[:, m] @ L
        g = np.array([lm @ dVk @ lm for dVk in dV])
        nu.append(2.0 * d[m] ** 2 / (g @ A_cov @ g))
    nu = np.asarray(nu)
    if nu.size == 0 or not np.all(np.isfinite(nu)):
        return float("nan")
    if nu.size == 1:
        return float(nu[0])
    if np.all(np.abs(np.diff(nu)) < 1e-8):
        return float(nu.mean())
    if np.any(nu <= 2):
        return 2.0
    E = float(np.sum(nu / (nu - 2)))
    return 2 * E / (E - nu.size)


#: σ²_animal at or below this fraction of σ²_resid is treated as the boundary (0).
_SINGULAR_TOL = 1e-8


def _singular_fit(
    X: np.ndarray, y: np.ndarray, z: np.ndarray, groups: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """The REML fit at σ²_animal = 0 (a "singular" fit): OLS on the scaled data,
    ``σ²_resid = RSS / (N − p)``. Returns ``(β, Cov β, σ²_resid)`` when zero is the
    REML optimum — the score in σ²_animal there, ``½(yᵀPΣ_aPy − tr PΣ_a)``, is not
    positive — else ``None`` (an interior optimum the optimizer failed to reach).
    """
    N, p = X.shape
    if N <= p:
        return None
    try:
        XtX_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    b = XtX_inv @ (X.T @ y)
    r = y - X @ b
    s2e = float(r @ r) / (N - p)
    if not s2e > 0:
        return None
    codes = pd.factorize(pd.Series(groups))[0]
    zr = np.bincount(codes, weights=z * r)                       # z_gᵀ r_g
    Xz = np.column_stack([np.bincount(codes, weights=z * X[:, j]) for j in range(p)])
    tr_PSa = (float(z @ z) - float(np.einsum("gi,ij,gj->", Xz, XtX_inv, Xz))) / s2e
    score = 0.5 * (float(zr @ zr) / s2e ** 2 - tr_PSa)
    if score > 0:
        return None
    return b, s2e * XtX_inv, s2e


def _prepare(
    points: pd.DataFrame,
    reference: str,
    test: str,
    *,
    phi_limit: float,
    theta_floor: float,
    theta_ceiling: float,
    min_points_per_condition: int,
) -> pd.DataFrame:
    """The linear model's point filters, then keep, per protein, only the experiments
    where both conditions still have ``min_points_per_condition`` points."""
    p = points[points["condition"].astype(str).isin([reference, test])].copy()
    p["condition"] = p["condition"].astype(str)
    p["experiment"] = p["experiment"].astype(str)
    phi = to_phi(p["theta"].to_numpy(), theta_floor=theta_floor,
                 theta_ceiling=theta_ceiling)
    t = p["labeling_time"].to_numpy(dtype=float)
    keep = np.isfinite(phi) & np.isfinite(t) & (phi > phi_limit) & (t > 0)
    p = p.loc[keep].assign(phi=phi[keep])
    n = p.groupby(["protein", "experiment", "condition"])["phi"].transform("size")
    p = p[n >= min_points_per_condition]
    paired = p.groupby(["protein", "experiment"])["condition"].transform("nunique")
    p = p[paired == 2].copy()
    composite = ("composite:" + p["condition"] + "|"
                 + p["labeling_time"].astype(float).astype(str)
                 + "|" + p["biological_replicate"].astype(str))
    if "individual" in p.columns:
        from riana.core.protein import individual_label
        named = p["individual"].map(individual_label)
        p["animal"] = np.where(named != "", "individual:" + named, composite)
    else:
        p["animal"] = composite
    return p


def _check_animals(p: pd.DataFrame) -> str:
    """Log how animals were identified; return the key used (for provenance):
    ``"individual"``, ``"composite"`` or ``"individual+composite"``."""
    named = p["animal"].str.startswith("individual:")
    if not named.any():
        return "composite"
    a = p.loc[named]
    per_animal = a.groupby("animal").agg(
        conditions=("condition", "nunique"), times=("labeling_time", "nunique"))
    n_multi_cond = int((per_animal["conditions"] > 1).sum())
    if n_multi_cond:
        _LOGGER.warning(
            "--interaction: %d individual(s) appear under both conditions and are "
            "treated as one animal each. If animal IDs restart within each group "
            "(1, 2, 3 per condition), make them unique in characteristics[individual].",
            n_multi_cond)
    n_long = int((per_animal["times"] > 1).sum())
    if n_long:
        _LOGGER.info("--interaction: %d individual(s) sampled at several labeling "
                     "times (longitudinal).", n_long)
    if named.all():
        return "individual"
    _LOGGER.warning(
        "--interaction: %d of %d points have no individual; those use the "
        "(condition, labeling time, biological replicate) composite.",
        int((~named).sum()), len(p))
    return "individual+composite"


def _fit_protein(
    prot: str,
    g: pd.DataFrame,
    reference: str,
    test: str,
    phi_limit: float,
    min_shared_animals: int,
    weights: str = "wls",
    var_prior: tuple[float, float] | None = None,
) -> list[dict]:
    """Fit one protein; ``[]`` when it cannot be fitted (fewer than two paired
    experiments, too few shared animals or points, or no convergence).

    ``weights`` is ``"wls"`` or ``"wls-var"``; the latter needs ``g``'s
    ``theta_var``/``theta_df`` and ``var_prior = (d0, s0²)``."""
    strata = sorted(g["experiment"].unique())
    S = len(strata)
    if S < 2:
        return []
    sidx = g["experiment"].map({s: i for i, s in enumerate(strata)}).to_numpy()
    day = g["labeling_time"].to_numpy(dtype=float)
    is_test = (g["condition"] == test).to_numpy(dtype=float)
    C = np.vstack([np.eye(S - 1), -np.ones(S - 1)])          # sum contrasts, S × (S−1)
    X = np.column_stack(
        [day * (sidx == s) for s in range(S)]
        + [day * is_test]
        + [day * is_test * C[sidx, j] for j in range(S - 1)])
    y = g["phi"].to_numpy(dtype=float)
    # Integer codes in order of appearance, so the fit does not depend on how the
    # animals are named (statsmodels orders groups by label).
    animal = pd.factorize(g["animal"])[0]
    p = X.shape[1]                                             # 2S
    shared = int((pd.Series(animal).value_counts() >= 2).sum())
    if shared < min_shared_animals or len(y) <= p + 1:
        return []

    from statsmodels.regression.mixed_linear_model import MixedLM

    # The per-experiment fits' weights, from a first OLS pass (the slopes are saturated
    # per cell, so its fitted values are those fits' too): w = (1 − θ̂)² = exp(2·φ̂),
    # and under wls-var also ÷ the moderated per-point variance (Ṽar = (d0·s0² +
    # df·var)/(d0 + df); a point without one gets the prior s0², as in
    # fit_linear_deltak). √w transforms the model to unit-variance errors.
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    sw = np.exp(np.maximum(X @ b0, phi_limit))
    if weights == "wls-var" and var_prior is not None:
        d0, s0 = var_prior
        var = g["theta_var"].to_numpy(dtype=float)
        vdf = g["theta_df"].to_numpy(dtype=float)
        vmod = (d0 * s0 + vdf * var) / (d0 + vdf)
        sw = sw / np.sqrt(np.where(np.isfinite(vmod) & (vmod > 0), vmod, s0))
    Xt, yt, zt = X * sw[:, None], y * sw, day * sw
    res = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = MixedLM(yt, Xt, groups=animal, exog_re=zt[:, None]).fit(
                reml=True, method=["bfgs", "nm"])
    except Exception as exc:  # noqa: BLE001 - statsmodels raises various
        _LOGGER.debug("interaction fit failed for %s: %s", prot, exc)
    s2a = float(np.asarray(res.cov_re)[0, 0]) if res is not None else float("nan")
    if res is not None and res.converged and s2a > _SINGULAR_TOL * float(res.scale):
        b = np.asarray(res.fe_params, dtype=float)
        V = np.asarray(res.cov_params(), dtype=float)[:p, :p]
        s2e = float(res.scale)
    else:
        # statsmodels does not settle on σ²_animal = 0 (it stops there unconverged,
        # with wrong fixed effects): take the singular fit directly when zero is the
        # REML optimum, as lme4 does.
        boundary = _singular_fit(Xt, yt, zt, animal)
        if boundary is None:
            _LOGGER.debug("interaction fit did not converge for %s", prot)
            return []
        b, V, s2e = boundary
        s2a = 0.0
    dfs = containment_df(X, animal)
    ref_slope = b[:S]
    test_slope = ref_slope + b[S] + C @ b[S + 1:]
    sigma_animal = float(np.sqrt(max(s2a, 0.0)))

    def wald(est: np.ndarray, cov: np.ndarray) -> float:
        try:
            return float(est @ np.linalg.solve(cov, est) / len(est))
        except np.linalg.LinAlgError:
            return float("nan")

    # Mean Δk (condition main effect): −β_shift, a between-animal column.
    shift_se = float(np.sqrt(V[S, S]))
    mean_df = float(dfs[S])
    mean_t = b[S] / shift_se if shift_se > 0 else np.nan
    mean_p = float(2 * stats.t.sf(abs(mean_t), mean_df)) if np.isfinite(mean_t) else np.nan

    # Additive interaction: γ = 0 (within-animal columns, one df).
    gi = slice(S + 1, p)
    add_F = wald(b[gi], V[gi, gi])
    add_df2 = float(dfs[S + 1])
    add_p = float(stats.f.sf(add_F, S - 1, add_df2)) if np.isfinite(add_F) else np.nan

    # Ratio-scale interaction: equal log(k_test / k_ref) across experiments.
    rat_F = rat_df2 = rat_p = np.nan
    ddf_fallback = False
    if np.all(ref_slope < 0) and np.all(test_slope < 0):
        J = np.zeros((S, p))                       # ∂ log(test_s / ref_s) / ∂b
        for s in range(S):
            J[s, s] = 1.0 / test_slope[s] - 1.0 / ref_slope[s]
            J[s, S] = 1.0 / test_slope[s]
            J[s, S + 1:] = C[s] / test_slope[s]
        D = np.hstack([np.eye(S - 1), -np.ones((S - 1, 1))])   # r_s − r_last
        DJ = D @ J
        r = np.log(test_slope / ref_slope)
        rat_F = wald(D @ r, DJ @ V @ DJ.T)
        rat_df2 = satterthwaite_ddf(DJ, Xt, yt, zt, animal, s2a, s2e)
        if not np.isfinite(rat_df2) or rat_df2 <= 0:
            # Undefined (information not positive definite): the conservative
            # between-animal containment df.
            used = np.abs(DJ).max(axis=0) > 1e-12 * np.abs(DJ).max()
            rat_df2 = float(dfs[used].min())
            ddf_fallback = True
        rat_p = float(stats.f.sf(rat_F, S - 1, rat_df2)) if np.isfinite(rat_F) else np.nan

    protein_level = {
        "n_experiments": S, "n_animals": int(pd.Series(animal).nunique()),
        "sigma_animal": sigma_animal,
        "mean_delta_k": -float(b[S]), "mean_delta_k_se": shift_se,
        "mean_delta_k_df": mean_df, "mean_delta_k_p": mean_p,
        "interaction_F": rat_F, "interaction_df1": float(S - 1),
        "interaction_df2": rat_df2, "interaction_p": rat_p,
        "interaction_additive_F": add_F, "interaction_additive_df2": add_df2,
        "interaction_additive_p": add_p,
        "_ddf_fallback": ddf_fallback,
    }
    rows = []
    for s, name in enumerate(strata):
        l_d = np.zeros(p)
        l_d[S] = 1.0
        l_d[S + 1:] = C[s]
        rows.append({
            "protein": prot, "experiment": name,
            "n_points": int((sidx == s).sum()),
            "k_reference": -float(ref_slope[s]), "k_test": -float(test_slope[s]),
            "delta_k": -float(l_d @ b), "delta_k_se": float(np.sqrt(l_d @ V @ l_d)),
            "k_ratio": float(test_slope[s] / ref_slope[s]),
            **protein_level,
        })
    return rows


# Process-pool plumbing: the protein frames are pickled to each worker once via the
# initializer; only an integer index crosses per task (mirrors core/protein.py).
_WORKER_STATE: dict = {}


def _init_worker(groups, args) -> None:
    _WORKER_STATE["groups"] = groups
    _WORKER_STATE["args"] = args


def _fit_worker(index: int) -> list[dict]:
    prot, g = _WORKER_STATE["groups"][index]
    return _fit_protein(prot, g, *_WORKER_STATE["args"])


def fit_interaction(
    points: pd.DataFrame,
    *,
    reference_condition: str,
    test_condition: str,
    phi_limit: float = -4.0,
    theta_floor: float = 0.001,
    theta_ceiling: float = 0.999,
    min_points_per_condition: int = 2,
    min_shared_animals: int = 2,
    weights: str = "wls",
    workers: int = 1,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """Per-protein experiment × condition interaction (random-animal mixed model).

    Args:
        points: long table with ``experiment``, ``condition``, ``protein``,
            ``biological_replicate``, ``labeling_time``, ``theta`` and optionally
            ``individual`` — the collapsed rollup points
            (``riana_rollup_fractions.txt`` with ``fs`` as ``theta``).
        reference_condition / test_condition: the two conditions compared; Δk and
            the ratio are ``test`` relative to ``reference``. Other conditions are
            ignored.
        phi_limit, theta_floor, theta_ceiling, min_points_per_condition: as in
            :func:`riana.core.linear_model.fit_linear_deltak`. Per protein, an
            experiment enters the fit only if both conditions keep
            ``min_points_per_condition`` points; a protein needs two such experiments.
        min_shared_animals: a protein needs at least this many animals contributing
            to two or more experiments (otherwise σ_animal is not identifiable).
        weights: the per-experiment fit's scheme (``--linear-weights``). ``"wls"``
            (default) or ``"wls-var"`` (needs ``theta_var``/``theta_df``; without them
            it falls back to ``"wls"``, as the per-experiment fit does). ``"ols"``, the
            per-experiment fit's audit-only unweighted path, has no interaction
            counterpart: it is logged and ``"wls"`` is used.
        workers: processes (statsmodels is GIL-bound; threads would not help).
        progress_callback: ``callback(done, total)`` per protein.

    Returns:
        One row per ``(protein, experiment)`` in :data:`INTERACTION_COLUMNS`; empty
        when no protein qualifies. ``attrs["animal_key"]`` records how animals were
        identified (``"individual"``, ``"composite"`` or ``"individual+composite"``).

    Raises:
        DataError: missing columns, an unknown ``weights``, a condition not present in
            the data, or the two conditions are the same.
    """
    from riana.exceptions import DataError

    if weights not in LINEAR_WEIGHT_SCHEMES:
        raise DataError(
            f"interaction weights must be one of {LINEAR_WEIGHT_SCHEMES}, got {weights!r}.")
    if weights == "ols":
        _LOGGER.info("--interaction: 'ols' is the per-experiment fit's audit path; the "
                     "interaction uses 'wls'.")
        weights = "wls"

    need = {"experiment", "condition", "protein", "biological_replicate",
            "labeling_time", "theta"}
    missing = need - set(points.columns)
    if missing:
        raise DataError(f"interaction points missing columns {sorted(missing)}.")
    reference, test = str(reference_condition), str(test_condition)
    if reference == test:
        raise DataError(
            f"reference and test condition are the same ({reference!r}).")
    present = set(points["condition"].astype(str).unique())
    for role, cond in (("reference", reference), ("test", test)):
        if cond not in present:
            raise DataError(
                f"{role} condition {cond!r} is not among the conditions in the data "
                f"({sorted(present)}).")

    # wls-var: the moderation prior from every point, as the per-experiment fit
    # estimates it from every collapsed cell (same points → same d0, s0²).
    var_prior = None
    if weights == "wls-var":
        if {"theta_var", "theta_df"} <= set(points.columns):
            v = points["theta_var"].to_numpy(dtype=float)
            vdf = points["theta_df"].to_numpy(dtype=float)
            ok = np.isfinite(v) & (v > 0) & np.isfinite(vdf) & (vdf > 0)
            var_prior = _fit_fdist(v[ok], vdf[ok])
            _LOGGER.info("--interaction wls-var: robust fitFDist d0=%.2f, s0^2=%.3g "
                         "(%d points with a variance)", var_prior[0], var_prior[1],
                         int(ok.sum()))
        else:
            _LOGGER.info("--interaction: wls-var requested but the points carry no "
                         "theta_var/theta_df — using wls.")
            weights = "wls"

    p = _prepare(points, reference, test, phi_limit=phi_limit,
                 theta_floor=theta_floor, theta_ceiling=theta_ceiling,
                 min_points_per_condition=min_points_per_condition)
    animal_key = _check_animals(p) if not p.empty else "composite"
    groups = [(str(k), g) for k, g in p.groupby("protein", sort=True)
              if g["experiment"].nunique() >= 2]
    args = (reference, test, phi_limit, min_shared_animals, weights, var_prior)
    if workers > 1 and len(groups) > 1:
        chunk = max(1, len(groups) // (workers * 8))
        with futures.ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker, initargs=(groups, args),
        ) as ex:
            computed = list(iter_progress(
                ex.map(_fit_worker, range(len(groups)), chunksize=chunk),
                len(groups), progress_callback))
    else:
        computed = list(iter_progress(
            (_fit_protein(prot, g, *args) for prot, g in groups),
            len(groups), progress_callback))

    n_fallback = sum(bool(rows[0]["_ddf_fallback"]) for rows in computed if rows)
    if n_fallback:
        _LOGGER.info("--interaction: Satterthwaite df undefined for %d protein(s); "
                     "their ratio-scale test uses the between-animal df.", n_fallback)
    out = pd.DataFrame([r for rows in computed for r in rows],
                       columns=[c for c in INTERACTION_COLUMNS if c not in _ADJUSTED])
    for adj, raw in _ADJUSTED.items():
        out[adj] = np.nan
    out.attrs["animal_key"] = animal_key
    if out.empty:
        return out[INTERACTION_COLUMNS]

    # BH across proteins — one p per protein for each test (one family each).
    from statsmodels.stats.multitest import multipletests
    per_protein = out.drop_duplicates("protein").set_index("protein")
    for adj, raw in _ADJUSTED.items():
        pv = per_protein[raw].dropna()
        if not pv.empty:
            padj = dict(zip(pv.index, multipletests(pv.to_numpy(), method="fdr_bh")[1]))
            out[adj] = out["protein"].map(padj)
    return out[INTERACTION_COLUMNS]
