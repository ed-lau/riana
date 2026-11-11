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
- ``Var(ε) ∝ 1/(1 − θ̂)²`` — the shipped delta-method weight from a first OLS pass (nlme
  ``varExp(form = ~fitted(.))`` with δ fixed at −1; ATF6 estimates δ ≈ −1.1).

**Animal** = (condition, labeling_time, biological_replicate), so ``biological_replicate``
must identify the same animal in every experiment of one condition × timepoint (as an SDRF
stratified by tissue does when each animal keeps one replicate index across its tissues).
Experiments that are *different* animals (e.g. separate cohorts) do not fit this model.

**Fit.** REML with statsmodels ``MixedLM`` on the √w-transformed data, which is nlme
``lme(..., weights = varFixed(~1/w))``. BFGS, not statsmodels' default L-BFGS: on these small
variances (~1e-5) L-BFGS stops at a wrong point while reporting convergence. Denominator df
follow nlme's containment rule (:func:`containment_df`). Checked against nlme on ATF6
(estimates within 1.5e-3, F within 1.3 %, df identical).

**Tests reported per protein** (BH across proteins, one family per test):

- ``interaction_*`` (ratio scale, primary) — is ``k_test / k_reference`` the same in every
  experiment? A uniform proportional effect is not an interaction on this scale, though it
  gives a larger absolute Δk in faster-turnover experiments. Wald F on the log ratios (delta
  method). The contrast draws on between-animal (``β_shift``) and within-animal coefficients,
  which nlme would not give one df; the smaller (between-animal) df is used, which is
  conservative.
- ``interaction_additive_*`` — are the Δk equal across experiments (``γ = 0``)? Exact
  within-animal df; the test nlme's ``anova`` gives for the interaction columns.
- ``mean_delta_k_*`` — Δk averaged over experiments (= −β_shift), the condition main effect,
  at between-animal df. Per-experiment Δk *tests* stay in ``riana_rollup_proteins.txt``.
"""

from __future__ import annotations

import logging
import warnings
from concurrent import futures
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats

from riana.core.linear_model import to_phi
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
    p["animal"] = (p["condition"] + "|" + p["labeling_time"].astype(float).astype(str)
                   + "|" + p["biological_replicate"].astype(str))
    return p


def _fit_protein(
    prot: str,
    g: pd.DataFrame,
    reference: str,
    test: str,
    phi_limit: float,
    min_shared_animals: int,
) -> list[dict]:
    """Fit one protein; ``[]`` when it cannot be fitted (fewer than two paired
    experiments, too few shared animals or points, or no convergence)."""
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
    animal = g["animal"].to_numpy()
    p = X.shape[1]                                             # 2S
    shared = int((pd.Series(animal).value_counts() >= 2).sum())
    if shared < min_shared_animals or len(y) <= p + 1:
        return []

    from statsmodels.regression.mixed_linear_model import MixedLM

    # The shipped weight w = (1 − θ̂)² = exp(2·φ̂) from a first OLS pass (the slopes are
    # saturated per cell, so these are the per-experiment fits' weights); √w transforms
    # the model to unit-variance errors.
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    sw = np.exp(np.maximum(X @ b0, phi_limit))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = MixedLM(y * sw, X * sw[:, None], groups=animal,
                          exog_re=(day * sw)[:, None]).fit(reml=True, method=["bfgs", "nm"])
    except Exception as exc:  # noqa: BLE001 - statsmodels raises various
        _LOGGER.debug("interaction fit failed for %s: %s", prot, exc)
        return []
    if not res.converged:
        _LOGGER.debug("interaction fit did not converge for %s", prot)
        return []

    b = np.asarray(res.fe_params, dtype=float)
    V = np.asarray(res.cov_params(), dtype=float)[:p, :p]
    dfs = containment_df(X, animal)
    ref_slope = b[:S]
    test_slope = ref_slope + b[S] + C @ b[S + 1:]
    sigma_animal = float(np.sqrt(max(float(np.asarray(res.cov_re)[0, 0]), 0.0)))

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
        used = np.abs(DJ).max(axis=0) > 1e-12 * np.abs(DJ).max()
        rat_df2 = float(dfs[used].min())
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
    workers: int = 1,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """Per-protein experiment × condition interaction (random-animal mixed model).

    Args:
        points: long table with ``experiment``, ``condition``, ``protein``,
            ``biological_replicate``, ``labeling_time``, ``theta`` — the collapsed
            rollup points (``riana_rollup_fractions.txt`` with ``fs`` as ``theta``).
        reference_condition / test_condition: the two conditions compared; Δk and
            the ratio are ``test`` relative to ``reference``. Other conditions are
            ignored.
        phi_limit, theta_floor, theta_ceiling, min_points_per_condition: as in
            :func:`riana.core.linear_model.fit_linear_deltak`. Per protein, an
            experiment enters the fit only if both conditions keep
            ``min_points_per_condition`` points; a protein needs two such experiments.
        min_shared_animals: a protein needs at least this many animals contributing
            to two or more experiments (otherwise σ_animal is not identifiable).
        workers: processes (statsmodels is GIL-bound; threads would not help).
        progress_callback: ``callback(done, total)`` per protein.

    Returns:
        One row per ``(protein, experiment)`` in :data:`INTERACTION_COLUMNS`; empty
        when no protein qualifies.

    Raises:
        DataError: missing columns, or a condition not present in the data, or the
            two conditions are the same.
    """
    from riana.exceptions import DataError

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

    p = _prepare(points, reference, test, phi_limit=phi_limit,
                 theta_floor=theta_floor, theta_ceiling=theta_ceiling,
                 min_points_per_condition=min_points_per_condition)
    groups = [(str(k), g) for k, g in p.groupby("protein", sort=True)
              if g["experiment"].nunique() >= 2]
    args = (reference, test, phi_limit, min_shared_animals)
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

    out = pd.DataFrame([r for rows in computed for r in rows],
                       columns=[c for c in INTERACTION_COLUMNS if c not in _ADJUSTED])
    for adj, raw in _ADJUSTED.items():
        out[adj] = np.nan
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
