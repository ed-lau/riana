"""Linearized turnover model + cross-sample Δk (the ``linear simple`` rollup model).

The nonlinear rollup models (``simple``/``guan``/``fornasiero``) fit one ``k_deg``
per ``(experiment, condition, protein)`` independently via ``curve_fit``. This
module is the **linearized** alternative: it transforms fraction-new θ into
clearance φ = ``log(1 − θ)`` — which is **linear in time through the origin**
(φ = −k·t) for the simple one-exponent model — and fits the two conditions of a
protein **jointly** in one OLS with a ``day:condition`` interaction. That joint
fit is what yields a **cross-sample Δk** (the difference in slope) with a proper
shared-variance p-value, then Benjamini-Hochberg across proteins.

Adapted from `data/notebook/03_R_linearmodel_reference.Rmd` (the R `lm` +
`emmeans::emtrends` + `p.adjust(BH)` workup); the study-specific compartment /
strain / JCAST machinery there does not apply. The pairwise slope contrast is
emmeans/emtrends for a fixed-effects OLS — a linear combination of the fitted
coefficients — so we get it from statsmodels (`OLS.t_test`) rather than
reimplementing the contrast algebra (PROJECT_REVIEW.md §3, Track C).

**Weighting (WLS by default, since 2026-07) — load-bearing.** θ carries roughly
homoscedastic *measurement* noise on the FS scale, but φ is a **log** of it, so by the
delta method ``Var(φ) = σ_θ²/(1 − θ)²`` — the φ-residuals are strongly
**heteroscedastic**, their SD blowing up as θ → 1 (measured on real data: residual SD
0.060 → 0.603 across θ bins; regressing ``log|resid|`` on ``−log(1−θ)`` gives slope
**0.80** [0.77, 0.83], vs 1.0 predicted for pure FS-scale noise and 0 for φ-scale — i.e.
predominantly FS-scale, likely with a smaller φ-scale floor). An **unweighted** OLS
treats that 10× SD range as equal, which makes it badly anti-conservative: on RIANA's own
design a true null is rejected **~28 %** of the time at α = 0.05, 95 % CIs cover ~54 %, and
k is biased **−14 %** in the fast tail (real data, vs the nonlinear MLE).

So the fit is **weighted** by the delta-method inverse variance ``(1 − θ)²``, taken from
the **fitted** value (one IRLS step; see :func:`fit_linear_deltak`) rather than the observed
θ — weighting by the *observed* θ makes each weight a function of that point's own error and
biases k low. The weighted estimator is the delta-method linearization of the exact MLE
(which is plain nonlinear LS on the FS scale — the ``simple`` model) and recovers **~97 % of
its efficiency**, while keeping the closed-form joint covariance the Δk contrast needs.
``weights="ols"`` restores the old unweighted fit for audit. Full workup:
``reports/2026-07-13_linear_model_wls.md``.

**Plateau truncation (φ-limit) — required, and linear-only.** φ = `log(1 − θ)`
descends without bound as θ → 1, but θ saturates at a measurement ceiling: once a
fast curve hits θ ≈ 0.95–0.99 the later timepoints are noise around the floor, not
slope, and including them **drags the through-origin slope flat and breaks
linearity**. So per curve we drop points past a configurable ``phi_limit``
(default −4 ≈ θ 0.98; −3 ≈ θ 0.95). The nonlinear models do **not** truncate —
their asymptote parameter fits the plateau — so this is specific to the linear
model.
"""

from __future__ import annotations

import logging
from typing import Callable

import numpy as np
import pandas as pd

from riana.progress import iter_progress

_LOGGER = logging.getLogger(__name__)

#: Output columns of :func:`fit_linear_deltak`, one row per
#: ``(experiment, protein, condition)``.
#:
#: **CI provenance — analytic, NOT bootstrap.** ``ci_lo`` / ``ci_hi`` (per-condition
#: k) and ``delta_k_se`` / ``delta_k_p`` (the Δk contrast) come from the OLS
#: coefficient covariance via statsmodels (``conf_int`` / ``t_test``), i.e. a
#: closed-form t-distribution interval. This is deliberately different from the
#: nonlinear models (``simple`` / ``guan`` / ``fornasiero``), whose CIs are a
#: residual bootstrap in :func:`riana.core.protein._fit_kdeg`. The linearized
#: model has a closed-form covariance, so bootstrapping it would add nothing —
#: ``n_boot`` / ``random_state`` are ignored on this path.
LINEAR_COLUMNS = [
    "experiment", "condition", "protein",
    "n_points", "k_deg", "ci_lo", "ci_hi", "R_squared",
    "delta_k", "delta_k_se", "delta_k_p", "delta_k_p_adj",
]


def to_phi(
    theta: np.ndarray,
    *,
    theta_floor: float = 0.001,
    theta_ceiling: float = 0.999,
) -> np.ndarray:
    """θ (fraction new) → clearance φ = ``log(1 − θ)``.

    θ is clamped to ``[theta_floor, theta_ceiling]`` first so φ is finite (θ = 1
    would be −∞, θ > 1 from measurement noise would be NaN). φ = 0 at θ = 0, which
    is why the downstream fit is through the origin.
    """
    th = np.clip(np.asarray(theta, dtype=float), theta_floor, theta_ceiling)
    return np.log(1.0 - th)


def truncate_plateau(
    day: np.ndarray, phi: np.ndarray, *, phi_limit: float = -4.0
) -> tuple[np.ndarray, np.ndarray]:
    """Drop saturated points: keep only those with φ **above** ``phi_limit``.

    Past the limit, θ is at its measurement ceiling and φ is floor-noise, not
    slope (φ ≈ −4 is θ ≈ 0.98; −3 is θ ≈ 0.95). Filtering by value (rather than
    "first crossing") is robust to non-monotone noise — φ is monotone in t in
    expectation, so this is the saturated tail.
    """
    day = np.asarray(day, dtype=float)
    phi = np.asarray(phi, dtype=float)
    keep = np.isfinite(phi) & np.isfinite(day) & (phi > phi_limit)
    return day[keep], phi[keep]


def _fit_fdist(s2, df, *, trim=0.10, d0_fallback=2.0):
    """Robust Smyth (2004) fitFDist → ``(d0, s0²)`` from per-point variances ``s2`` and
    their df. Winsorizes the log-variance deviations (limma ``robust=TRUE`` / Phipson
    2016) — proteomics variances are heavy-tailed and the outliers bias d0 LOW
    (under-shrinkage). Returns ``d0_fallback`` if the moment-match degenerates."""
    from scipy.special import digamma, polygamma
    from scipy.optimize import brentq
    s2 = np.asarray(s2, float); df = np.asarray(df, float)
    if s2.size < 8:
        return float(d0_fallback), (float(np.mean(s2)) if s2.size else 1.0)
    e = np.log(s2) - digamma(df / 2.0) + np.log(df / 2.0)
    if trim > 0:
        lo, hi = np.quantile(e, [trim, 1.0 - trim]); e = np.clip(e, lo, hi)
    evar = float(np.var(e, ddof=1) - np.mean(polygamma(1, df / 2.0)))
    d0 = float(d0_fallback)
    if evar > 0:
        try:
            d0 = 2.0 * brentq(lambda x: polygamma(1, x) - evar, 1e-4, 1e5)
        except Exception:
            d0 = float(d0_fallback)
    if not (np.isfinite(d0) and d0 > 0):
        d0 = float(d0_fallback)
    s0 = float(np.exp(np.mean(e) + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    if not (np.isfinite(s0) and s0 > 0):
        s0 = float(np.mean(s2))
    return d0, s0


#: Weighting schemes for the linearized fit (:func:`fit_linear_deltak`).
#:
#: ``"wls"`` (default) — **weighted** LS with the delta-method inverse variance ``(1−θ̂)²``,
#: from the FITTED value (one IRLS step). ``"wls-var"`` — opt-in; multiplies that by the
#: **per-point** precision ``1/Var̂(θ)`` (the inverse-variance collapse's ``theta_var``),
#: eBayes-moderated with a robust-fitFDist ``d0`` (fixed-2 fallback). Improves
#: biological-replicate consistency ~10% on rich multi-peptide time series (see
#: ``reports/2026-07-24_linear_wls_per_point_var.md``); falls back to ``"wls"`` if the points
#: carry no per-point variance. ``"ols"`` — the pre-2026-07 unweighted fit, for audit only
#: (anti-conservative; see the module note).
LINEAR_WEIGHT_SCHEMES = ("wls", "ols", "wls-var")


def fit_linear_deltak(
    points: pd.DataFrame,
    *,
    phi_limit: float = -4.0,
    theta_floor: float = 0.001,
    theta_ceiling: float = 0.999,
    min_points: int = 3,
    min_points_per_condition: int = 2,
    reference_condition: str | None = None,
    test_condition: str | None = None,
    weights: str = "wls",
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """Per-protein linearized k per condition + a cross-condition Δk test.

    Args:
        points: long table with columns ``experiment``, ``condition``,
            ``protein``, ``labeling_time``, ``theta`` — one row per collapsed
            ``(condition, timepoint)`` point (the inverse-variance peptide θ the
            weighted rollup already builds, ``riana_rollup_fractions.txt``).
        phi_limit: plateau-truncation threshold; points with φ ≤ this are dropped
            per curve (default −4 ≈ θ 0.98).
        theta_floor / theta_ceiling: θ clamp before the log (keeps φ finite).
        min_points: a protein needs at least this many surviving points total to
            be fit.
        min_points_per_condition: a condition needs at least this many surviving
            points to get a slope (its k). In auto mode (no ``test_condition``) the
            Δk contrast needs **exactly two** qualifying conditions; with a named
            ``reference`` / ``test`` pair it needs both of *those* to qualify.
        reference_condition: the baseline of the Δk contrast — ``delta_k`` is
            ``k(test) − k(reference)``. Defaults to the alphabetically-first
            condition; pass e.g. ``"control"`` to make a treatment read positive
            when faster.
        test_condition: the comparison condition of the Δk contrast. When BOTH
            ``reference_condition`` and ``test_condition`` are given, the contrast is
            that **named pair**, computed from the joint (all-condition) fit — so it
            works even when a protein has **more than two** conditions (an interim
            for multi-group projects, ahead of full all-pairwise/Tukey). The extra
            conditions still contribute their own k rows **and** the shared residual
            variance the contrast's SE pools over, so scope the conditions in the
            SDRF/project deliberately if that pooling is unwanted. When
            ``test_condition`` is ``None`` the legacy auto mode applies: a contrast
            is emitted only for a protein with exactly two qualifying conditions.
        weights: ``"wls"`` (default) or ``"ols"`` — see :data:`LINEAR_WEIGHT_SCHEMES`
            and the "Weighting" note in the module docstring. ``"ols"`` reproduces the
            pre-2026-07 unweighted estimator (anti-conservative; for audit only).

    Returns:
        One row per ``(experiment, protein, condition)`` with ``k_deg`` (= −slope)
        and its CI (from the joint model covariance — emtrends), the joint
        uncentered ``R_squared``, ``n_points``, plus protein-level ``delta_k`` /
        ``delta_k_se`` / ``delta_k_p`` (the pairwise slope contrast,
        ``k(test) − k(reference)``, present for a protein that carries the contrast
        pair — the named pair, or the two conditions in auto mode) and
        ``delta_k_p_adj`` (Benjamini-Hochberg across all proteins that have a Δk).
        Columns are :data:`LINEAR_COLUMNS`.
    """
    import statsmodels.formula.api as smf

    if weights not in LINEAR_WEIGHT_SCHEMES:
        from riana.exceptions import DataError
        raise DataError(
            f"linear-model weights must be one of {LINEAR_WEIGHT_SCHEMES}, got {weights!r}."
        )

    need = {"experiment", "condition", "protein", "labeling_time", "theta"}
    missing = need - set(points.columns)
    if missing:
        from riana.exceptions import DataError
        raise DataError(f"linear-model points missing columns {sorted(missing)}.")

    # wls-var (opt-in): the per-point precision needs the collapse's theta_var/theta_df.
    # Absent them (a pooled/legacy points table) fall back to wls; otherwise estimate the
    # eBayes hyperparameters ONCE across every cell — the robust-fitFDist d0 and pooled s0²
    # applied per protein below.
    _d0 = _s0 = None
    if weights == "wls-var":
        if {"theta_var", "theta_df"} <= set(points.columns):
            _v = points["theta_var"].to_numpy(float); _dfv = points["theta_df"].to_numpy(float)
            _ok = np.isfinite(_v) & (_v > 0) & np.isfinite(_dfv) & (_dfv > 0)
            _d0, _s0 = _fit_fdist(_v[_ok], _dfv[_ok])
            _LOGGER.info("linear wls-var: robust fitFDist d0=%.2f, s0^2=%.3g "
                         "(%d cells with a per-point variance)", _d0, _s0, int(_ok.sum()))
        else:
            _LOGGER.info("linear wls-var requested but points carry no "
                         "theta_var/theta_df — using wls.")
            weights = "wls"

    # A named contrast condition must actually be in the data — otherwise every
    # protein silently misses it and delta_k is all-NaN. Fail loudly with the
    # available choices (the CLI surfaces this as a bad-parameter error).
    present = set(points["condition"].astype(str).unique())
    for role, cond in (("reference", reference_condition), ("test", test_condition)):
        if cond is not None and str(cond) not in present:
            from riana.exceptions import DataError
            raise DataError(
                f"{role} condition {cond!r} is not among the conditions in the data "
                f"({sorted(present)}). Check the spelling / the SDRF condition values."
            )
    if (reference_condition is not None and test_condition is not None
            and str(reference_condition) == str(test_condition)):
        from riana.exceptions import DataError
        raise DataError(
            f"reference and test condition are the same ({reference_condition!r}) — "
            "pick two different conditions. A self-contrast has no Δk: the ±1 "
            "coefficients land on one slope and degenerate into a spurious "
            "slope-vs-zero test (delta_k = −k, p ≈ 0), not a difference."
        )

    rows: list[dict] = []
    grouped = points.groupby(["experiment", "protein"], sort=False)
    for (exp, prot), grp in iter_progress(grouped, grouped.ngroups, progress_callback):
        # φ-transform and per-condition plateau truncation.
        per_cond: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        per_cond_var: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for cond, cg in grp.groupby("condition", sort=True):
            day0 = cg["labeling_time"].to_numpy()
            phi = to_phi(cg["theta"].to_numpy(),
                         theta_floor=theta_floor, theta_ceiling=theta_ceiling)
            t, p = truncate_plateau(day0, phi, phi_limit=phi_limit)
            # Drop t = 0. In a through-origin fit it has ZERO leverage on the slope
            # (x = 0 contributes nothing to Σxy or Σx²), so k is unchanged — but its θ
            # is pinned by the ``theta_floor`` clamp (true θ(0) = 0, so ~half the
            # measurements go negative and clamp to the floor), which makes its residual
            # artificially ≈ 0 against the model's exact 0. That deflates the residual
            # variance and shrinks EVERY standard error. Excluding it only corrects the
            # inference. See reports/2026-07-13_linear_model_wls.md.
            #
            # EXCEPT under ``weights="ols"``: that scheme exists solely to reproduce the
            # pre-2026-07 fit for audit, and that fit KEPT t = 0. Dropping it there would
            # make ``ols`` a hybrid that never shipped, so the audit path retains t = 0.
            keep = np.ones(len(t), dtype=bool) if weights == "ols" else (t > 0)
            t, p = t[keep], p[keep]
            if len(t) >= min_points_per_condition:
                per_cond[str(cond)] = (t, p)
                if weights == "wls-var":
                    # the SAME plateau∘(t>0) mask on the original arrays → align var/df
                    m = (np.isfinite(phi) & np.isfinite(day0)
                         & (phi > phi_limit) & (day0 > 0))
                    per_cond_var[str(cond)] = (cg["theta_var"].to_numpy()[m],
                                               cg["theta_df"].to_numpy()[m])
        n_total = sum(len(t) for t, _ in per_cond.values())
        if not per_cond or n_total < min_points:
            continue

        # Joint through-origin OLS: φ ~ 0 + day:C(condition). Each condition's
        # slope is a direct coefficient (= −k); the interaction shares one
        # residual variance across conditions, which is what gives the Δk test
        # honest pooled standard errors.
        fit_df = pd.DataFrame({
            "day": np.concatenate([t for t, _ in per_cond.values()]),
            "phi": np.concatenate([p for _, p in per_cond.values()]),
            "condition": np.concatenate(
                [[c] * len(t) for c, (t, _) in per_cond.items()]),
        })
        if weights == "wls-var":
            fit_df["theta_var"] = np.concatenate([per_cond_var[c][0] for c in per_cond])
            fit_df["theta_df"] = np.concatenate([per_cond_var[c][1] for c in per_cond])
        formula = "phi ~ 0 + day:C(condition)"
        try:
            res = smf.ols(formula, data=fit_df).fit()
            if weights == "wls":
                # Delta-method inverse variance: Var(φ) = σ²/(1−θ)², so the optimal
                # weight is (1−θ)². Take it from the FITTED value, not the observed θ:
                # w = (1−θ̂)² = exp(2·φ̂) = exp(−2·k̂·t). Using the *observed* θ would make
                # each weight a function of that point's own error — down-weighting the
                # points noise pushed high (the most-negative φ), which flattens the slope
                # and biases k LOW. The fitted value depends on k̂ (all n points) and t
                # (noise-free), so the weights are exogenous and the bias vanishes: the
                # feasible-GLS rule. ONE IRLS step suffices (iterating further lets the
                # weights chase the noise realization and degrades Type-I).
                #
                # φ̂ is floored at ``phi_limit`` so the weight range matches the plateau
                # truncation and cannot underflow to 0 for a fast curve at a late t.
                w = np.exp(2.0 * np.maximum(res.fittedvalues.to_numpy(), phi_limit))
                res = smf.wls(formula, data=fit_df, weights=w).fit()
            elif weights == "wls-var":
                # wls's exogenous (1−θ̂)² transform weight × the per-point precision
                # 1/Ṽar(θ): the collapse's ``theta_var`` eBayes-moderated toward the
                # pooled ``s0²`` by the robust-fitFDist ``d0`` (per-point df = ``theta_df``),
                # Ṽar = (d0·s0² + df·var)/(d0 + df). Only RELATIVE weights matter (WLS
                # estimates its own scale), so the variance units cancel. A point with no
                # usable per-point variance (``theta_var`` NaN) falls back to the pooled
                # prior ``s0²`` (``_s0``) — the df→0 limit of the moderation formula — so it
                # sits on the SAME scale as the moderated points and keeps the plain
                # (1−θ̂)² weight. (A constant 1.0 divisor there would be ~1/s0² ≈ 100–1000×
                # lighter than its finite-variance siblings, silently dropping a real θ
                # measurement from the joint fit.)
                w = np.exp(2.0 * np.maximum(res.fittedvalues.to_numpy(), phi_limit))
                var = fit_df["theta_var"].to_numpy()
                df = fit_df["theta_df"].to_numpy()
                vmod = (_d0 * _s0 + df * var) / (_d0 + df)
                w = w / np.where(np.isfinite(vmod) & (vmod > 0), vmod, _s0)
                res = smf.wls(formula, data=fit_df, weights=w).fit()
        except Exception as exc:  # noqa: BLE001 - statsmodels raises various
            _LOGGER.debug("linear fit failed for %s/%s: %s", exp, prot, exc)
            continue

        slope = res.params
        ci = res.conf_int()
        r2 = float(res.rsquared)  # uncentered for a no-intercept model
        conds = sorted(per_cond)
        coef = {c: f"day:C(condition)[{c}]" for c in conds}

        # Cross-condition Δk, Δk = k(test) − k(reference). With a named
        # reference/test pair, contrast THAT pair from the joint (all-condition)
        # fit — so it works even when the protein has >2 conditions (option B, the
        # multi-group interim). Without a named test, fall back to the legacy auto
        # mode: contrast the two conditions of an exactly-two-condition protein.
        def _contrast(ref: str, other: str) -> tuple[float, float, float]:
            names = list(slope.index)
            c_vec = np.zeros(len(names))
            c_vec[names.index(coef[other])] = 1.0   # slope_test
            c_vec[names.index(coef[ref])] = -1.0    # − slope_reference
            tt = res.t_test(c_vec)
            # k = −slope, so Δk = k_test − k_ref = −(slope_test − slope_ref).
            return (-float(np.ravel(tt.effect)[0]),
                    float(np.ravel(tt.sd)[0]),
                    float(np.ravel(tt.pvalue)[0]))

        delta_k = delta_se = delta_p = float("nan")
        if reference_condition is not None and test_condition is not None:
            if (reference_condition in coef and test_condition in coef
                    and coef[reference_condition] in slope.index
                    and coef[test_condition] in slope.index):
                delta_k, delta_se, delta_p = _contrast(
                    reference_condition, test_condition)
        elif len(conds) == 2 and all(coef[c] in slope.index for c in conds):
            ref = reference_condition if reference_condition in conds else conds[0]
            other = next(c for c in conds if c != ref)
            delta_k, delta_se, delta_p = _contrast(ref, other)

        for c in conds:
            name = coef[c]
            if name not in slope.index:
                continue
            s = float(slope[name])
            lo_s, hi_s = float(ci.loc[name, 0]), float(ci.loc[name, 1])
            rows.append({
                "experiment": exp, "condition": c, "protein": prot,
                "n_points": int(len(per_cond[c][0])),
                "k_deg": -s,
                # k = −slope, so the CI bounds swap and negate.
                "ci_lo": -hi_s, "ci_hi": -lo_s,
                "R_squared": r2,
                "delta_k": delta_k, "delta_k_se": delta_se,
                "delta_k_p": delta_p,
            })

    out = pd.DataFrame(rows, columns=[c for c in LINEAR_COLUMNS
                                      if c != "delta_k_p_adj"])
    out["delta_k_p_adj"] = np.nan
    if out.empty:
        return out[LINEAR_COLUMNS]

    # Benjamini-Hochberg across proteins (one p per protein, not per row).
    prot_p = (
        out.loc[out["delta_k_p"].notna(), ["experiment", "protein", "delta_k_p"]]
        .drop_duplicates(["experiment", "protein"])
    )
    if not prot_p.empty:
        from statsmodels.stats.multitest import multipletests
        prot_p = prot_p.assign(
            delta_k_p_adj=multipletests(
                prot_p["delta_k_p"].to_numpy(), method="fdr_bh")[1])
        adj = dict(zip(zip(prot_p["experiment"], prot_p["protein"]),
                       prot_p["delta_k_p_adj"]))
        out["delta_k_p_adj"] = [
            adj.get((e, p), np.nan)
            for e, p in zip(out["experiment"], out["protein"])
        ]
    return out[LINEAR_COLUMNS]
