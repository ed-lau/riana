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

    need = {"experiment", "condition", "protein", "labeling_time", "theta"}
    missing = need - set(points.columns)
    if missing:
        from riana.exceptions import DataError
        raise DataError(f"linear-model points missing columns {sorted(missing)}.")

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
        for cond, cg in grp.groupby("condition", sort=True):
            phi = to_phi(cg["theta"].to_numpy(),
                         theta_floor=theta_floor, theta_ceiling=theta_ceiling)
            t, p = truncate_plateau(
                cg["labeling_time"].to_numpy(), phi, phi_limit=phi_limit)
            if len(t) >= min_points_per_condition:
                per_cond[str(cond)] = (t, p)
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
        try:
            res = smf.ols("phi ~ 0 + day:C(condition)", data=fit_df).fit()
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
