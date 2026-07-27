"""`linear simple` Δk — OLS vs WLS: the diagnostic + Monte-Carlo harness.

The workup behind `reports/2026-07-13_linear_model_wls.md` and the switch of the
`linear simple` rollup from unweighted OLS to WLS (`--linear-weights wls`, the default).

THE ARGUMENT, in one line. θ (fraction new) carries roughly homoscedastic *measurement*
noise on the FS scale, but the model fits φ = log(1−θ), a **log** of it. By the delta
method ``Var(φ) = σ_θ²/(1−θ)²`` — so the φ-residuals are strongly **heteroscedastic**, their
SD blowing up as θ → 1. An unweighted OLS treats that ~10× SD range as equal, which makes the
Δk test badly anti-conservative and biases k low in the fast tail.

FOUR ESTIMATORS are compared throughout:
  ols          unweighted (the pre-2026-07 default)
  wls_obs      weights = (1−θ_obs)²   — the naive inverse-variance form. Fixes calibration but
               BIASES k LOW: the weight is a function of that point's OWN error, so points noise
               pushed high get down-weighted, discarding the most-negative φ and flattening the
               slope. (Endogenous weights.)
  wls_fit      weights = (1−θ̂)² = exp(2·φ̂)  — from the FITTED value (one IRLS step). The weight
               then depends only on k̂ (all n points) and t (noise-free), so it is exogenous and
               the bias vanishes: the feasible-GLS rule. THIS IS WHAT PRODUCTION NOW DOES
               (`core.linear_model.fit_linear_deltak(weights="wls")`).
  wls_fit_var  weights = (1−θ̂)²/Var(θ_i)  — the candidate refit. Multiplies wls_fit's transform
               factor by the EXOGENOUS per-point precision (1/Var(θ_i), from the peptide's
               PI-width / inverse-variance collapse). Still needs the one IRLS step for (1−θ̂).
               Only differs from wls_fit when point precision genuinely VARIES — so it is
               exercised in the heteroscedastic regime below, and MUST be fed the *estimated*
               (noisy) variance, as production would. NOT yet wired into production.

The candidate is validated with the `varweight` subcommand, which sweeps a homoscedastic and a
heteroscedastic noise regime × an oracle and a noisy per-point-variance estimate, and ranks
against the WEIGHTED nonlinear MLE (`wmle`, the correct reference under heteroscedastic FS-scale
noise — `mle` is only correct when σ_θ is constant).

All simulation is faithful to RIANA's real chain: clamp θ→[0.001, 0.999] → φ = log(1−θ) →
truncate φ > `--phi-limit` → drop t=0 → joint through-origin fit `phi ~ 0 + day:C(condition)`.

Subcommands
-----------
  diagnose    Is the noise FS-scale or φ-scale?  (real data; the PREMISE the whole thing rests on)
  recover     MC: bias / CI coverage / Type-I across a k grid, all estimators
  power       MC: power curve for the Δk contrast (the δ=0 row is Type-I)
  efficiency  MC: how close is each estimator to the EXACT MLE (nonlinear LS on FS)?
  varweight   MC: (1−θ̂)²/Var(θ) vs (1−θ̂)² across homo/hetero × oracle/noisy-Var — the effect size
  irls        MC: how many IRLS steps, and does dropping t=0 matter?
  real        Real-data impact on a run: k, Δk, significance, agreement with the nonlinear MLE
  all         diagnose + recover + power + efficiency

Usage
-----
    python -m tests.benchmark.bench_linear_weights all
    python -m tests.benchmark.bench_linear_weights all --regime both --plot 36_out  # homo vs hetero
    python -m tests.benchmark.bench_linear_weights recover --regime hetero --spread 0.8 --var-df 16
    python -m tests.benchmark.bench_linear_weights varweight --nsim 500      # the 2×2 effect size
    python -m tests.benchmark.bench_linear_weights real --run runs/lve_atr_clean \
        --reference control --test atrium

``--regime {homo,hetero,both}`` is the homo/hetero switch for recover/power/efficiency (with
``--spread`` the hetero magnitude and ``--var-df`` the Var-estimate fidelity); ``both`` overlays
the two scenarios in one table/plot.

Calibrated defaults come from `runs/lve_atr_clean`: σ_θ ≈ 0.056, days 0..30, and the
empirical protein-k percentiles (5/25/50/75/95 = 0.03/0.05/0.075/0.11/0.216 /day).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import curve_fit
from statsmodels.stats.multitest import multipletests

from riana.core.linear_model import to_phi

# ---- calibrated to lve_atr (see the report) ---------------------------------
DAYS = np.array([0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 25, 30], float)
SIGMA_THETA = 0.056          # per-point noise SD on the θ (FS) scale
PHI_LIMIT = -4.0             # RIANA's plateau truncation (θ ≈ 0.982)
_PI_SPAN = 3.29              # the M5 5–95% PI spans 3.29·σ (= _PI_SPAN_SIGMA in core.protein)
K_GRID = [0.03, 0.05, 0.075, 0.11, 0.216, 0.30]   # lve_atr 5/25/50/75/95 pct + a fast one
METHODS = ("ols", "wls_obs", "wls_fit", "wls_fit_var")


# --------------------------------------------------------------------------- #
# heteroscedastic per-point precision (the regime where per-point Var(θ) can matter)
# --------------------------------------------------------------------------- #
def per_point_sigma(n, sigma, spread, rng):
    """Per-point θ-noise SD. ``spread == 0`` → homoscedastic (every point = ``sigma``).
    ``spread > 0`` → a lognormal spread of precisions (a deep/shallow-peptide mix),
    normalized so E[σ_i²] = sigma² — the AVERAGE noise level is held fixed, so a
    homo-vs-hetero comparison isolates the heteroscedasticity, not just "more noise"."""
    if spread <= 0:
        return np.full(n, sigma)
    s = rng.lognormal(0.0, spread, n) / np.exp(spread ** 2)   # E[s²] = 1
    return sigma * s


def estimate_var(true_sigma, var_df, rng):
    """The per-point variance ESTIMATE the WLS is fed. ``var_df = inf`` → the oracle
    (true σ_i²). Finite ``var_df`` → an unbiased but noisy estimate σ_i²·χ²(df)/df, as
    a PI-width/bootstrap variance from ~df effective observations would be — this is
    what production actually has, and the honest test of whether the gain survives."""
    tv = np.asarray(true_sigma, float) ** 2
    if not np.isfinite(var_df):
        return tv
    return tv * rng.chisquare(var_df, len(tv)) / var_df


# --------------------------------------------------------------------------- #
# regime selection — the homo/hetero split, a first-class knob for every MC command
# --------------------------------------------------------------------------- #
_MCOLOR = {"ols": "C0", "wls_obs": "C1", "wls_fit": "C2", "wls_fit_var": "C3"}
_REGIME_LS = {"homo": "-", "hetero": "--"}


def regimes_from_args(a):
    """``[(label, spread), ...]`` selected by ``--regime``; ``--spread`` sets the hetero
    magnitude. ``both`` runs homoscedastic and heteroscedastic so the scenarios sit side
    by side in one table / plot. recover / power / efficiency all honour this."""
    homo, hetero = ("homo", 0.0), ("hetero", a.spread)
    return {"homo": [homo], "hetero": [hetero], "both": [homo, hetero]}[a.regime]


# --------------------------------------------------------------------------- #
# the estimators
# --------------------------------------------------------------------------- #
def _design(day, cond):
    """Through-origin `day:condition` design (one slope column per condition)."""
    X = np.zeros((len(day), 2))
    X[cond == 0, 0] = day[cond == 0]
    X[cond == 1, 1] = day[cond == 1]
    return X


def fit(day, cond, phi, theta, how, phi_limit=PHI_LIMIT, var=None):
    """Joint through-origin fit; returns per-condition k, a CI for k(A), and the Δk test.

    ``var`` = per-point Var(θ_i), used only by ``wls_fit_var`` (ignored otherwise;
    ``None`` there degrades it to ``wls_fit``)."""
    X = _design(day, cond)
    if (cond == 0).sum() < 2 or (cond == 1).sum() < 2:
        return None
    if how == "ols":
        m = sm.OLS(phi, X).fit()
    elif how == "wls_obs":
        m = sm.WLS(phi, X, weights=(1.0 - theta) ** 2).fit()
    elif how == "wls_fit":
        phat = sm.OLS(phi, X).fit().fittedvalues            # φ̂ = −k̂·t
        w = np.exp(2.0 * np.maximum(phat, phi_limit))       # (1 − θ̂)², floored
        m = sm.WLS(phi, X, weights=w).fit()
    elif how == "wls_fit_var":
        # (1−θ̂)²/Var(θ_i): the delta-method transform factor (from the FITTED θ̂ — same
        # one-step OLS pilot as wls_fit, so the ONLY difference is the /Var factor) times
        # the EXOGENOUS per-point precision. var=None falls back to wls_fit.
        phat = sm.OLS(phi, X).fit().fittedvalues
        w = np.exp(2.0 * np.maximum(phat, phi_limit))
        if var is not None:
            w = w / np.clip(np.asarray(var, float), 1e-12, None)
        m = sm.WLS(phi, X, weights=w).fit()
    else:
        raise ValueError(how)
    k = -m.params                                            # k = −slope
    ci = m.conf_int()
    tt = m.t_test(np.array([-1.0, 1.0]))                     # Δk = k_B − k_A
    return dict(kA=k[0], kB=k[1], kA_lo=-ci[0, 1], kA_hi=-ci[0, 0],
                dk=-float(np.ravel(tt.effect)[0]),
                dk_p=float(np.ravel(tt.pvalue)[0]))


def mle_nls(day, theta):
    """The EXACT MLE under HOMOscedastic FS-scale Gaussian noise: unweighted nonlinear LS
    on θ. (RIANA's nonlinear `simple` model — the unbiased reference when σ_θ is constant;
    suboptimal, though still unbiased, when it is not.)"""
    try:
        return float(curve_fit(lambda x, k: 1 - np.exp(-k * x), day, theta,
                               p0=[0.05], bounds=(0, 5), maxfev=5000)[0][0])
    except Exception:
        return np.nan


def wmle_nls(day, theta, var):
    """The EXACT MLE under HETEROscedastic FS-scale noise: nonlinear LS on θ weighted by
    1/Var(θ_i) (via ``sigma=√var``). The correctly-specified reference the per-point-Var
    WLS is chasing; reduces to :func:`mle_nls` when var is constant."""
    try:
        return float(curve_fit(lambda x, k: 1 - np.exp(-k * x), day, theta,
                               p0=[0.05], bounds=(0, 5), maxfev=5000,
                               sigma=np.sqrt(np.clip(var, 1e-12, None)),
                               absolute_sigma=True)[0][0])
    except Exception:
        return np.nan


def simulate(kA, kB, sigma, rng, days=DAYS, phi_limit=PHI_LIMIT, drop_t0=True,
             spread=0.0, var_df=np.inf):
    """One two-condition dataset through RIANA's exact chain.

    ``spread`` sets the per-point precision spread (0 = homoscedastic); ``var_df`` the
    fidelity of the variance ESTIMATE fed to the WLS (inf = oracle). Returns the
    linear-space arrays + the per-point variance ESTIMATE (for wls_fit_var), then the raw
    (day, cond, θ) and the per-point TRUE variance (for the MLE references).
    """
    day = np.concatenate([days, days])
    cond = np.concatenate([np.zeros(len(days), int), np.ones(len(days), int)])
    ktrue = np.where(cond == 0, kA, kB)
    sig_i = per_point_sigma(len(day), sigma, spread, rng)     # per-point θ-noise SD
    theta_raw = 1 - np.exp(-ktrue * day) + rng.normal(0, sig_i)  # FS-scale noise
    var_est = estimate_var(sig_i, var_df, rng)                # what the WLS is fed
    var_true = sig_i ** 2                                     # for the weighted MLE
    phi = to_phi(theta_raw)                                   # RIANA clamp + log
    keep = phi > phi_limit                                    # RIANA plateau truncation
    if drop_t0:
        keep &= day > 0                                       # zero leverage + clamp artifact
    return (day[keep], cond[keep], phi[keep], np.clip(theta_raw[keep], 0.001, 0.999),
            var_est[keep], day, cond, theta_raw, var_true)


# --------------------------------------------------------------------------- #
# 1. diagnose — is the noise FS-scale (=> weight) or φ-scale (=> OLS is fine)?
# --------------------------------------------------------------------------- #
def cmd_diagnose(a):
    import statsmodels.formula.api as smf
    pts = pd.read_table(Path(a.run) / "riana_rollup_fractions.txt", comment="#")
    pts = pts.rename(columns={"fs": "theta"})
    rec = []
    for (_e, prot), grp in pts.groupby(["experiment", "protein"], sort=False):
        per = {}
        for cond, cg in grp.groupby("condition", sort=True):
            th = cg["theta"].to_numpy()
            phi = to_phi(th)
            day = cg["labeling_time"].to_numpy().astype(float)
            keep = np.isfinite(phi) & (phi > a.phi_limit) & (day > 0)
            if keep.sum() >= 2:
                per[str(cond)] = (day[keep], phi[keep], th[keep])
        if len(per) < 2:
            continue
        d = pd.DataFrame({
            "day": np.concatenate([v[0] for v in per.values()]),
            "phi": np.concatenate([v[1] for v in per.values()]),
            "theta": np.concatenate([v[2] for v in per.values()]),
            "condition": np.concatenate([[c] * len(v[0]) for c, v in per.items()]),
        })
        try:
            r = smf.ols("phi ~ 0 + day:C(condition)", data=d).fit()
        except Exception:
            continue
        d["resid"] = r.resid
        d["protein"] = prot
        rec.append(d)
    d = pd.concat(rec, ignore_index=True)

    low = d[d["theta"] < 0.2]
    sigma_hat = float(low["resid"].std() * (1 - low["theta"]).mean())
    d["bin"] = pd.cut(d["theta"].clip(0, 1), [0, .2, .4, .6, .8, .9, .95, 1.01])
    g = d.groupby("bin", observed=True).agg(
        n=("resid", "size"), resid_sd=("resid", "std"), mean_theta=("theta", "mean"))
    g["pred_FS_scale"] = sigma_hat / (1 - g["mean_theta"])    # rises like 1/(1-θ)
    g["pred_phi_scale"] = float(low["resid"].std())           # flat

    print(f"proteins {d.protein.nunique()}   points {len(d)}   "
          f"calibrated σ_θ ≈ {sigma_hat:.4f}\n")
    print(g.to_string(float_format=lambda x: f"{x:9.4f}"))

    m = d[(d["resid"].abs() > 1e-9) & (d["theta"] < 0.995)].copy()
    m["y"] = np.log(m["resid"].abs())
    m["x"] = -np.log(1 - m["theta"].clip(upper=0.99))
    f = smf.ols("y ~ x", data=m).fit()
    b, (lo, hi) = f.params["x"], f.conf_int().loc["x"]
    print(f"\nlog|resid| ~ −log(1−θ):  slope = {b:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")
    print("  slope ≈ 1 → FS-scale noise (φ heteroscedastic) → WEIGHT the fit")
    print("  slope ≈ 0 → φ-scale noise (homoscedastic)      → OLS is already right")
    if a.plot:
        _plot_diagnose(g, sigma_hat, a.plot)
    return g


# --------------------------------------------------------------------------- #
# 2. recover — bias / coverage / Type-I  (kA == kB, so every rejection is a FALSE positive)
# --------------------------------------------------------------------------- #
def cmd_recover(a):
    rng = np.random.default_rng(a.seed)
    rows = []
    for reg, spread in regimes_from_args(a):
        for k in a.k_grid:
            acc = {m: {"k": [], "cov": [], "rej": []} for m in METHODS}
            for _ in range(a.nsim):
                day, cond, phi, th, var, *_ = simulate(
                    k, k, a.sigma, rng, phi_limit=a.phi_limit,
                    spread=spread, var_df=a.var_df)
                for m in METHODS:
                    r = fit(day, cond, phi, th, m, a.phi_limit, var=var)
                    if r is None:
                        continue
                    acc[m]["k"].append(r["kA"])
                    acc[m]["cov"].append(r["kA_lo"] <= k <= r["kA_hi"])
                    acc[m]["rej"].append(r["dk_p"] < 0.05)
            for m in METHODS:
                kk = np.array(acc[m]["k"])
                rows.append(dict(regime=reg, k_true=k, method=m, n=len(kk),
                                 rel_bias=(kk.mean() - k) / k,
                                 rmse=float(np.sqrt(np.mean((kk - k) ** 2))),
                                 coverage=float(np.mean(acc[m]["cov"])),
                                 type1=float(np.mean(acc[m]["rej"]))))
    rec = pd.DataFrame(rows)
    print("RECOVERY / CI COVERAGE / TYPE-I   (true Δk = 0; nominal coverage .95, Type-I .05)")
    print(f"var-df = {a.var_df:g} (∞ = oracle Var)\n")
    for reg in rec.regime.unique():
        tag = f"  (σ-spread {a.spread})" if reg == "hetero" else ""
        print(f"================ regime: {reg}{tag} ================")
        for m in METHODS:
            print(f"--- {m} ---")
            print(rec[(rec.method == m) & (rec.regime == reg)]
                  [["k_true", "rel_bias", "rmse", "coverage", "type1"]]
                  .to_string(index=False, float_format=lambda x: f"{x:8.3f}"), "\n")
    print("SUMMARY over the k grid")
    print(rec.groupby(["regime", "method"], sort=False).agg(
        mean_rel_bias=("rel_bias", "mean"), mean_coverage=("coverage", "mean"),
        mean_type1=("type1", "mean"), worst_type1=("type1", "max"),
    ).to_string(float_format=lambda x: f"{x:8.3f}"))
    if a.plot:
        _plot_recover(rec, a.plot)
    return rec


# --------------------------------------------------------------------------- #
# 3. power — the Δk contrast. delta = 0 is Type-I; delta > 0 is power.
# --------------------------------------------------------------------------- #
def cmd_power(a):
    rng = np.random.default_rng(a.seed)
    rows = []
    for reg, spread in regimes_from_args(a):
        for d_k in a.deltas:
            acc = {m: [] for m in METHODS}
            for _ in range(a.nsim):
                day, cond, phi, th, var, *_ = simulate(
                    a.k0, a.k0 + d_k, a.sigma, rng, phi_limit=a.phi_limit,
                    spread=spread, var_df=a.var_df)
                for m in METHODS:
                    r = fit(day, cond, phi, th, m, a.phi_limit, var=var)
                    if r:
                        acc[m].append(r["dk_p"] < 0.05)
            rows.append(dict(regime=reg, delta_k=d_k,
                             **{m: float(np.mean(acc[m])) for m in METHODS}))
    pw = pd.DataFrame(rows)
    print(f"POWER of the Δk test at k0={a.k0} (α=0.05). The delta_k=0 row is TYPE-I.\n")
    for reg in pw.regime.unique():
        tag = f"  (σ-spread {a.spread})" if reg == "hetero" else ""
        print(f"--- regime: {reg}{tag} ---")
        print(pw[pw.regime == reg].drop(columns="regime")
              .to_string(index=False, float_format=lambda x: f"{x:8.3f}"), "\n")
    print("Note: OLS's apparent sensitivity at small δ is largely its false-positive rate.")
    if a.plot:
        _plot_power(pw, a.plot)
    return pw


# --------------------------------------------------------------------------- #
# 4. efficiency — vs the EXACT MLE (nonlinear LS on the FS scale)
# --------------------------------------------------------------------------- #
def cmd_efficiency(a):
    rng = np.random.default_rng(a.seed)
    rows = []
    for reg, spread in regimes_from_args(a):
        for k in a.k_grid:
            acc = {m: [] for m in ("wmle", "mle", *METHODS)}
            for _ in range(a.nsim):
                day, cond, phi, th, var, raw_day, raw_cond, raw_th, tv = simulate(
                    k, k, a.sigma, rng, phi_limit=a.phi_limit,
                    spread=spread, var_df=a.var_df)
                m0 = raw_cond == 0
                acc["mle"].append(mle_nls(raw_day[m0], raw_th[m0]))
                acc["wmle"].append(wmle_nls(raw_day[m0], raw_th[m0], tv[m0]))
                for m in METHODS:
                    r = fit(day, cond, phi, th, m, a.phi_limit, var=var)
                    acc[m].append(r["kA"] if r else np.nan)
            A = pd.DataFrame(acc).dropna()
            rmse_ref = float(np.sqrt(np.mean((A["wmle"] - k) ** 2)))   # regime-correct MLE
            for m in ("wmle", "mle", *METHODS):
                r = float(np.sqrt(np.mean((A[m] - k) ** 2)))
                rows.append(dict(regime=reg, k_true=k, estimator=m,
                                 rel_bias=(A[m].mean() - k) / k, rmse=r,
                                 rmse_vs_wmle=r / rmse_ref))
    eff = pd.DataFrame(rows)
    print("EFFICIENCY vs the WEIGHTED MLE (`wmle` — the regime-correct reference; it equals")
    print("the plain `mle` when noise is homoscedastic). rmse_vs_wmle 1.0 = optimal.\n")
    for reg in eff.regime.unique():
        tag = f"  (σ-spread {a.spread})" if reg == "hetero" else ""
        print(f"--- regime: {reg}{tag} ---")
        print(eff[eff.regime == reg].drop(columns="regime")
              .to_string(index=False, float_format=lambda x: f"{x:9.3f}"), "\n")
    print("wls_fit is the delta-method linearization of the MLE — ~1.0 under homoscedastic")
    print("noise, but > 1 under heteroscedastic; wls_fit_var closes that gap.")
    return eff


# --------------------------------------------------------------------------- #
# 4b. varweight — the candidate (1−θ̂)²/Var(θ) vs (1−θ̂)²: effect size + safety
# --------------------------------------------------------------------------- #
def cmd_varweight(a):
    """Does the per-point-variance weight earn its keep?

    Sweeps {homoscedastic, heteroscedastic} × {oracle Var, noisy Var estimate} and
    compares wls_fit_var against wls_fit, ranked by RMSE vs the WEIGHTED MLE (wmle, the
    correct reference under heteroscedastic FS-scale noise). Reads off: the efficiency
    gain, and whether a NOISY variance estimate keeps it without inflating Type-I.
    """
    rng = np.random.default_rng(a.seed)
    spread_hi = a.spread
    vdf = a.var_df if np.isfinite(a.var_df) else 4.0
    regimes = [("homosced.", 0.0), ("heterosced.", spread_hi)]
    vmodes = [("Var oracle", np.inf), (f"Var est(df={vdf:g})", vdf)]
    methods = ("wls_fit", "wls_fit_var")
    rows = []
    for reg, spread in regimes:
        for vname, vd in vmodes:
            for k in a.k_grid:
                est = {m: [] for m in methods}
                cov = {m: [] for m in methods}
                rej = {m: [] for m in methods}
                mle, wmle = [], []
                for _ in range(a.nsim):
                    day, cond, phi, th, var, rday, rcond, rth, tvar = simulate(
                        k, k, a.sigma, rng, phi_limit=a.phi_limit, spread=spread, var_df=vd)
                    m0 = rcond == 0
                    mle.append(mle_nls(rday[m0], rth[m0]))
                    wmle.append(wmle_nls(rday[m0], rth[m0], tvar[m0]))
                    for m in methods:
                        r = fit(day, cond, phi, th, m, a.phi_limit, var=var)
                        if r is None:
                            continue
                        est[m].append(r["kA"])
                        cov[m].append(r["kA_lo"] <= k <= r["kA_hi"])
                        rej[m].append(r["dk_p"] < 0.05)

                def rmse(arr):
                    return float(np.sqrt(np.nanmean((np.asarray(arr, float) - k) ** 2)))

                r_wmle = rmse(wmle)
                base = dict(regime=reg, vmode=vname, k=k)
                rows.append({**base, "method": "wmle", "rmse": r_wmle, "rmse_vs_wmle": 1.0,
                             "rel_bias": (np.nanmean(wmle) - k) / k,
                             "coverage": np.nan, "type1": np.nan})
                rows.append({**base, "method": "mle", "rmse": rmse(mle),
                             "rmse_vs_wmle": rmse(mle) / r_wmle,
                             "rel_bias": (np.nanmean(mle) - k) / k,
                             "coverage": np.nan, "type1": np.nan})
                for m in methods:
                    rows.append({**base, "method": m, "rmse": rmse(est[m]),
                                 "rmse_vs_wmle": rmse(est[m]) / r_wmle,
                                 "rel_bias": (np.mean(est[m]) - k) / k,
                                 "coverage": float(np.mean(cov[m])),
                                 "type1": float(np.mean(rej[m]))})
    df = pd.DataFrame(rows)
    order = ["wmle", "mle", "wls_fit", "wls_fit_var"]
    df["method"] = pd.Categorical(df["method"], order, ordered=True)
    summ = (df.groupby(["regime", "vmode", "method"], observed=True)
              .agg(rel_bias=("rel_bias", "mean"), rmse_vs_wmle=("rmse_vs_wmle", "mean"),
                   coverage=("coverage", "mean"), type1=("type1", "mean"))
              .reset_index())
    print(f"PER-POINT Var(θ) WLS — mean over the k grid  "
          f"(nsim={a.nsim}, hetero spread={spread_hi}, days={len(DAYS)})")
    print("rmse_vs_wmle: 1.00 = matches the weighted MLE; lower is better; wmle is the ref.\n")
    print(summ.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))

    print("\nEFFECT SIZE — wls_fit_var vs wls_fit (mean over k):")
    print(f"  {'regime':12s} {'Var mode':14s} {'RMSE Δ':>8}  "
          f"{'Type-I fit→var':>16}  {'cover fit→var':>15}")
    d = df[df.method.isin(methods)]
    for (reg, vm), g in d.groupby(["regime", "vmode"], observed=True, sort=False):
        rr = g.groupby("method", observed=True)["rmse"].mean()
        t1 = g.groupby("method", observed=True)["type1"].mean()
        cv = g.groupby("method", observed=True)["coverage"].mean()
        gain = (rr["wls_fit"] - rr["wls_fit_var"]) / rr["wls_fit"] * 100.0
        print(f"  {reg:12s} {vm:14s} {gain:+7.1f}%  "
              f"{t1['wls_fit']:.3f}→{t1['wls_fit_var']:.3f}      "
              f"{cv['wls_fit']:.3f}→{cv['wls_fit_var']:.3f}")
    print("\nRead: homosced.+oracle ⇒ ~0% (identical by construction — the sanity check); "
          "heterosced.+oracle ⇒ the ceiling gain; heterosced.+noisy ⇒ the REAL question.")
    print("A noisy estimate that inflates Type-I much past 0.05 is the veto — that is what")
    print("decides whether the per-point weight earns the production default.")
    return df


# --------------------------------------------------------------------------- #
# 4c. moderate — the three arms that could make per-point Var(θ) safe to ship
# --------------------------------------------------------------------------- #
def ebayes_var(var, d, pool, d0):
    """limma/eBayes moderated variance (arm #1): shrink each estimate toward the
    pool by prior df ``d0``. ``d`` = the estimate's own df; effective df → ``d + d0``.
    ``d0 → ∞`` recovers `wls_fit` (pooled), ``d0 → 0`` recovers `wls_fit_var` (raw)."""
    var, d = np.asarray(var, float), np.asarray(d, float)
    return (d0 * pool + d * var) / (d0 + d)


def threshold_var(var, d, pool, T):
    """Hard-threshold baseline (arm #2): use the per-point estimate only where its
    df exceeds ``T``, else fall back to the pooled variance (= `wls_fit` there)."""
    return np.where(np.asarray(d, float) > T, np.asarray(var, float), pool)


def simulate_moderate(kA, kB, sigma, rng, *, spread, df_lo, df_hi, frac_lo,
                      between, phi_limit=PHI_LIMIT):
    """Two-source variance model for the moderation arms. Each cell has a
    within-peptide σ (heteroscedastic via ``spread``) whose ESTIMATE has df drawn
    from a {``df_lo`` single-peptide, ``df_hi`` multi-peptide} mix (``frac_lo`` at
    the low end), PLUS a heterogeneous **between-peptide** component (``between``,
    the source of the measured 1.25–1.72× dispersion) that the propagated PI-variance
    MISSES. Returns the two variance estimates the arms choose between:
    ``prop`` (propagated, within-only → biased low, noisy at within-df) and ``emp``
    (empirical between-peptide scatter → unbiased of the true variance but noisy at
    df = n_pep−1; falls back to ``prop`` for single-peptide cells)."""
    day = np.concatenate([DAYS, DAYS])
    cond = np.concatenate([np.zeros(len(DAYS), int), np.ones(len(DAYS), int)])
    ktrue = np.where(cond == 0, kA, kB)
    n = len(day)
    within = per_point_sigma(n, sigma, spread, rng) ** 2
    is_lo = rng.random(n) < frac_lo
    d = np.where(is_lo, df_lo, df_hi).astype(float)              # within-peptide bootstrap df
    n_pep = np.where(is_lo, 1, 3)                                # peptides/cell (df_emp = n_pep−1)
    btw = between * within * rng.gamma(2.0, 0.5, n)              # heterogeneous between comp (mean≈between·within)
    true_var = within + btw                                      # what the cell θ actually scatters by
    theta_raw = 1 - np.exp(-ktrue * day) + rng.normal(0, np.sqrt(true_var))
    prop = within * rng.chisquare(d) / d                        # misses `btw` → biased low
    df_emp = np.maximum(n_pep - 1, 1)
    emp = np.where(n_pep >= 2, true_var * rng.chisquare(df_emp) / df_emp, prop)
    phi = to_phi(theta_raw)
    keep = (phi > phi_limit) & (day > 0)
    return (day[keep], cond[keep], phi[keep], np.clip(theta_raw[keep], 0.001, 0.999),
            prop[keep], emp[keep], d[keep], day, cond, theta_raw, true_var)


def cmd_moderate(a):
    """Rank the three moderation arms against `wls_fit` / raw `wls_fit_var`, at the
    measured production regime (heteroscedastic + per-point-df mix + between-peptide
    dispersion). The winner is the arm that keeps the RMSE gain while landing
    Type-I ≈ 0.05 and coverage ≈ 0.95 — the veto the raw weight fails."""
    rng = np.random.default_rng(a.seed)
    d0_grid, T_grid = a.d0_grid, a.t_grid
    arms = (["wls_fit", "wls_fit_var(raw)"]
            + [f"ebayes(d0={x})" for x in d0_grid]
            + [f"thresh(T={x})" for x in T_grid] + ["empirical"])
    rows, base = [], {}
    for k in a.k_grid:
        acc = {m: dict(kk=[], cov=[], rej=[]) for m in arms}
        wmle = []
        for _ in range(a.nsim):
            (day, cond, phi, th, prop, emp, d,
             rday, rcond, rth, rtv) = simulate_moderate(
                k, k, a.sigma, rng, spread=a.spread, df_lo=a.df_lo, df_hi=a.df_hi,
                frac_lo=a.frac_lo, between=a.between, phi_limit=a.phi_limit)
            m0 = rcond == 0
            wmle.append(wmle_nls(rday[m0], rth[m0], rtv[m0]))
            pool = float(np.mean(prop))
            variants = {"wls_fit": None, "wls_fit_var(raw)": prop, "empirical": emp}
            variants.update({f"ebayes(d0={x})": ebayes_var(prop, d, pool, x) for x in d0_grid})
            variants.update({f"thresh(T={x})": threshold_var(prop, d, pool, x) for x in T_grid})
            for m in arms:
                v = variants[m]
                r = fit(day, cond, phi, th, "wls_fit" if v is None else "wls_fit_var",
                        a.phi_limit, var=v)
                if r is None:
                    continue
                acc[m]["kk"].append(r["kA"])
                acc[m]["cov"].append(r["kA_lo"] <= k <= r["kA_hi"])
                acc[m]["rej"].append(r["dk_p"] < 0.05)
        rw = float(np.sqrt(np.nanmean((np.asarray(wmle) - k) ** 2)))
        base[k] = float(np.sqrt(np.mean((np.asarray(acc["wls_fit"]["kk"]) - k) ** 2)))
        for m in arms:
            kk = np.asarray(acc[m]["kk"])
            rmse = float(np.sqrt(np.mean((kk - k) ** 2)))
            rows.append(dict(k=k, arm=m, rmse=rmse, rmse_vs_wmle=rmse / rw,
                             gain=(base[k] - rmse) / base[k] * 100.0,
                             coverage=float(np.mean(acc[m]["cov"])),
                             type1=float(np.mean(acc[m]["rej"]))))
    df = pd.DataFrame(rows)
    df["arm"] = pd.Categorical(df["arm"], arms, ordered=True)
    summ = (df.groupby("arm", observed=True)
              .agg(rmse_gain_vs_fit=("gain", "mean"), type1=("type1", "mean"),
                   coverage=("coverage", "mean"), rmse_vs_wmle=("rmse_vs_wmle", "mean"))
              .reset_index())
    print(f"MODERATION ARMS — spread={a.spread}, df mix {{{a.df_lo:g}@{a.frac_lo:.0%} / "
          f"{a.df_hi:g}}}, between={a.between}  (nsim={a.nsim})")
    print("target: max rmse_gain_vs_fit at Type-I ≈ 0.05, coverage ≈ 0.95 "
          "(rmse_vs_wmle 1.0 = the weighted MLE)\n")
    print(summ.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
    print("\nwls_fit = shipped default (gain 0 by definition); wls_fit_var(raw) = the")
    print("unmoderated candidate (the Type-I offender). An arm 'wins' if it keeps most of")
    print("the raw gain with Type-I back near 0.05.")
    return df


# --------------------------------------------------------------------------- #
# 5. irls — how many steps, and does dropping t=0 matter?
# --------------------------------------------------------------------------- #
def cmd_irls(a):
    rng = np.random.default_rng(a.seed)
    print("IRLS steps × the t=0 point   (nominal coverage .95, Type-I .05)\n")
    print(f"{'k':>7} {'steps':>6} {'t=0':>6} {'rel_bias':>9} {'coverage':>9} {'type1':>7}")
    rows = []
    for k in (0.075, 0.216):
        for steps in (1, 2, 5):
            for drop in (True, False):
                cov, rej, bias = [], [], []
                for _ in range(a.nsim):
                    day, cond, phi, th, *_ = simulate(k, k, a.sigma, rng,
                                                      phi_limit=a.phi_limit, drop_t0=drop)
                    X = _design(day, cond)
                    if (cond == 0).sum() < 2 or (cond == 1).sum() < 2:
                        continue
                    m = sm.OLS(phi, X).fit()
                    for _s in range(steps):
                        w = np.exp(2.0 * np.maximum(m.fittedvalues, a.phi_limit))
                        m = sm.WLS(phi, X, weights=w).fit()
                    kk = -m.params
                    ci = m.conf_int()
                    bias.append(kk[0] - k)
                    cov.append(-ci[0, 1] <= k <= -ci[0, 0])
                    rej.append(float(np.ravel(m.t_test(np.array([-1., 1.])).pvalue)[0]) < 0.05)
                b, c, r = np.mean(bias) / k, np.mean(cov), np.mean(rej)
                rows.append(dict(k=k, steps=steps, drop_t0=drop,
                                 rel_bias=b, coverage=c, type1=r))
                print(f"{k:7.3f} {steps:6d} {'drop' if drop else 'keep':>6} "
                      f"{b:9.3f} {c:9.3f} {r:7.3f}")
    print("\nONE step suffices (it is enough to make the weights exogenous); iterating further")
    print("lets them chase the noise realization. Dropping t=0 tightens coverage/Type-I.")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 6. real — impact on a real run + agreement with the nonlinear MLE
# --------------------------------------------------------------------------- #
def cmd_real(a):
    pts = pd.read_table(Path(a.run) / "riana_rollup_fractions.txt", comment="#")
    pts = pts.rename(columns={"fs": "theta"})
    out, ref, test = {}, a.reference, a.test
    # wls_fit_var is omitted here: real `riana_rollup_fractions.txt` does not yet carry a
    # per-point Var(θ) column (that plumbing is the production step this bench gates).
    for how in ("ols", "wls_obs", "wls_fit"):
        rows = []
        for (_e, prot), grp in pts.groupby(["experiment", "protein"], sort=False):
            per = {}
            for cond, cg in grp.groupby("condition", sort=True):
                th = cg["theta"].to_numpy()
                phi = to_phi(th)
                day = cg["labeling_time"].to_numpy().astype(float)
                keep = np.isfinite(phi) & (phi > a.phi_limit) & (day > 0)
                if keep.sum() >= 2:
                    per[str(cond)] = (day[keep], phi[keep],
                                      np.clip(th[keep], 0.001, 0.999),
                                      cg["labeling_time"].to_numpy().astype(float), th)
            if ref not in per or test not in per:
                continue
            conds = [ref, test]
            day = np.concatenate([per[c][0] for c in conds])
            phi = np.concatenate([per[c][1] for c in conds])
            tht = np.concatenate([per[c][2] for c in conds])
            ci = np.concatenate([[i] * len(per[c][0]) for i, c in enumerate(conds)])
            r = fit(day, ci, phi, tht, how, a.phi_limit)
            if not r:
                continue
            rows.append(dict(protein=prot, k_ref=r["kA"], k_test=r["kB"],
                             delta_k=r["dk"], p=r["dk_p"],
                             k_nl=mle_nls(per[ref][3], per[ref][4])))
        d = pd.DataFrame(rows)
        d["p_adj"] = multipletests(d["p"], method="fdr_bh")[1]
        out[how] = d

    print(f"REAL DATA — {a.run}   ({test} vs {ref})\n")
    print(pd.DataFrame([
        dict(method=m, n=len(d), sig_p_adj_05=int((d.p_adj < 0.05).sum()),
             k_ref_med=d.k_ref.median(), k_test_med=d.k_test.median(),
             dk_med=d.delta_k.median())
        for m, d in out.items()
    ]).to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\nAgreement with the NONLINEAR MLE (unbiased reference), reference arm:")
    for m, d in out.items():
        j = d.dropna(subset=["k_nl"])
        ratio = (j.k_ref / j.k_nl).replace([np.inf, -np.inf], np.nan).dropna()
        fast = j[j.k_nl > 0.15]
        rf = (fast.k_ref / fast.k_nl).replace([np.inf, -np.inf], np.nan).dropna()
        print(f"  {m:8s} ρ={j.k_ref.corr(j.k_nl):.3f}  median k/k_MLE={ratio.median():.4f}"
              f"   fast tail (k>0.15, n={len(rf)}): {rf.median():.3f}")
    print("\n  OLS reads LOW in the fast tail; wls_obs reads low throughout (endogenous")
    print("  weights); wls_fit ≈ 1.000 everywhere.")
    return out


# --------------------------------------------------------------------------- #
# 7. measure — the REAL per-point-variance regime (places production on the curves)
# --------------------------------------------------------------------------- #
def _measure_one(run):
    """Per-run per-point-variance readout + the arrays for plotting. Reports:
      1. heteroscedasticity — the spread of per-point σ (≈ the bench --spread), at both the
         raw peptide level and the COLLAPSED cell level the linear model actually weights on;
      2. effective df — how well each point's σ is pinned (bootstrap timepoints/peptide,
         peptides/collapse-cell) → the Type-I-risk axis of the varweight sweep;
      3. calibration — do peptides in a cell scatter as much as their σ claims (dispersion)?
    """
    ff = Path(run) / "riana_fit_fractions.txt"
    d = pd.read_table(ff, comment="#")
    need = {"fs", "fs_lower", "fs_upper", "labeling_time", "concat", "protein id",
            "biological_replicate", "experiment", "condition"}
    miss = need - set(d.columns)
    if miss:
        raise SystemExit(f"{ff}\n  lacks {sorted(miss)} — need a PEPTIDE-level "
                         "riana_fit_fractions.txt (not the collapsed rollup one).")
    d = d[d["labeling_time"] > 0].copy()                    # t=0 is excluded from the linear fit
    d["sigma"] = (d["fs_upper"] - d["fs_lower"]) / _PI_SPAN  # the same σ the rollup collapses on
    v = d[np.isfinite(d["sigma"]) & (d["sigma"] > 0)].copy()
    CELL = ["experiment", "condition", "protein id", "biological_replicate", "labeling_time"]

    print(f"REAL VARIANCE REGIME — {run}  (peptide-level {ff.name})")
    print(f"points with a usable PI: {len(v)}/{len(d)}   "
          f"{v['concat'].nunique()} peptidoforms, {v['protein id'].nunique()} proteins\n")

    # 1. heteroscedasticity — spread of per-point σ (maps to the bench --spread)
    def _spread(sig):
        sig = np.asarray(sig, float)
        p10, p50, p90 = np.percentile(sig, [10, 50, 90])
        return float(np.std(np.log(sig))), p50, p90 / p10
    sp_pep, med_pep, ratio_pep = _spread(v["sigma"])
    cell_var = v.groupby(CELL)["sigma"].apply(lambda s: 1.0 / np.sum(1.0 / s.to_numpy() ** 2))
    cell_sigma = np.sqrt(cell_var.to_numpy())               # σ of the inverse-variance mean
    sp_cell, med_cell, ratio_cell = _spread(cell_sigma)
    print("1. HETEROSCEDASTICITY   SD(log σ) ← compare to the bench --spread "
          "(0 = homosced.; bench default 0.7)")
    print(f"   raw peptide points : spread {sp_pep:.3f}   median σ {med_pep:.4f}   "
          f"p90/p10 {ratio_pep:.1f}×")
    print(f"   COLLAPSED cells    : spread {sp_cell:.3f}   median σ {med_cell:.4f}   "
          f"p90/p10 {ratio_cell:.1f}×   ← what the linear model weights on\n")

    # 2. effective df — reliability of each point's variance estimate
    npts = v.groupby(["experiment", "condition", "concat",
                      "biological_replicate"])["labeling_time"].nunique()
    ppc = v.groupby(CELL)["concat"].nunique()
    print("2. EFFECTIVE df   (how noisy the per-point σ estimate is → the varweight df axis)")
    print(f"   timepoints/peptide (bootstrap df) : median {int(npts.median())}   "
          f"p10 {int(npts.quantile(.1))}  p90 {int(npts.quantile(.9))}")
    print(f"   peptides/collapse-cell            : median {int(ppc.median())}   "
          f"single-peptide cells {float((ppc == 1).mean()) * 100:.0f}%\n")

    # 3. calibration — do peptides scatter as much as their σ claims? (dispersion)
    disp = []
    for _, cg in v.groupby(CELL):
        if cg["concat"].nunique() >= 3:
            stated = float(np.mean(cg["sigma"].to_numpy() ** 2))
            if stated > 0:
                disp.append(float(np.var(cg["fs"].to_numpy(), ddof=1)) / stated)
    print("3. CALIBRATION   (cells with ≥3 peptides: empirical θ-scatter / mean stated σ²)")
    if disp:
        disp = np.array(disp)
        p10, p50, p90 = np.percentile(disp, [10, 50, 90])
        het = p90 / max(p10, 1e-9)
        print(f"   dispersion ratio: median {p50:.2f}  (p10 {p10:.2f}, p90 {p90:.2f})  "
              f"n={len(disp)} cells")
        print("     1 = calibrated; >1 = σ under-stated (the formal 1/Σ(1/σ²) under-covers). "
              f"p90/p10 = {het:.1f}× ⇒")
        print("     " + ("HETEROGENEOUS across cells → the empirical / random-effects arm can help"
                         if het > 3 else
                         "fairly uniform → a uniform scale is harmless for the WLS point estimate"))
    else:
        print("   (no cells with ≥3 peptides)")
    print(f"\nPLACEMENT: run `varweight --spread {sp_cell:.2f}` (the measured cell spread) to read "
          "the real-regime effect size; the median df above says which var-df row applies.")
    return dict(sigma_pep=v["sigma"].to_numpy(), sigma_cell=cell_sigma,
                disp=np.asarray(disp, float), ppc=ppc.to_numpy(), spread_cell=sp_cell)


def cmd_measure(a):
    """Measure the real per-point-variance regime from one or more runs' peptide-level
    fractions. With ``--plot`` it overlays their variance-spread distributions for visual
    comparison. Runs come from ``--runs`` (falls back to the single ``--run``)."""
    runs = a.runs if getattr(a, "runs", None) else [a.run]
    stats = {}
    for i, run in enumerate(runs):
        if i:
            print("\n" + "-" * 70)
        stats[Path(run).name] = _measure_one(run)
    if a.plot:
        _plot_measure(stats, a.plot)
    return stats


# --------------------------------------------------------------------------- #
# plots (optional)
# --------------------------------------------------------------------------- #
def _mpl(outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return plt


def _plot_diagnose(g, sigma, outdir):
    plt = _mpl(outdir)
    fig, ax = plt.subplots(figsize=(6, 4))
    x = g["mean_theta"]
    ax.plot(x, g["resid_sd"], "o-", label="observed φ-residual SD")
    ax.plot(x, g["pred_FS_scale"], "s--", label=r"predicted, FS-scale  $\sigma_\theta/(1-\theta)$")
    ax.plot(x, g["pred_phi_scale"], "^:", label="predicted, φ-scale (flat)")
    ax.set_xlabel(r"$\theta$ (fraction new)"); ax.set_ylabel("residual SD on φ")
    ax.set_title("Is the noise FS-scale?  (rising ⇒ yes ⇒ weight the fit)")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(Path(outdir) / "diag_noise_scale.png", dpi=150); plt.close(fig)


def _plot_recover(rec, outdir):
    """Colour = estimator, linestyle = regime (solid homo, dashed hetero) — so a
    ``--regime both`` run overlays the scenarios on the same axes for comparison."""
    plt = _mpl(outdir)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    multi = rec.regime.nunique() > 1
    for reg in rec.regime.unique():
        ls = _REGIME_LS.get(reg, "-")
        for m in METHODS:
            s = rec[(rec.method == m) & (rec.regime == reg)]
            lbl = f"{m} [{reg}]" if multi else m
            axes[0].plot(s.k_true, s.rel_bias, ls, marker="o", color=_MCOLOR[m], label=lbl)
            axes[1].plot(s.k_true, s.coverage, ls, marker="o", color=_MCOLOR[m], label=lbl)
            axes[2].plot(s.k_true, s.type1, ls, marker="o", color=_MCOLOR[m], label=lbl)
    axes[0].axhline(0, ls=":", c="k", lw=.8); axes[0].set_title("relative bias in k")
    axes[1].axhline(.95, ls=":", c="k", lw=.8); axes[1].set_title("95% CI coverage")
    axes[2].axhline(.05, ls=":", c="k", lw=.8); axes[2].set_title("Type-I of the Δk test")
    for ax in axes:
        ax.set_xlabel("true k (/day)"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(Path(outdir) / "mc_recover.png", dpi=150); plt.close(fig)


def _plot_power(pw, outdir):
    """Colour = estimator, linestyle = regime (solid homo, dashed hetero)."""
    plt = _mpl(outdir)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    multi = pw.regime.nunique() > 1
    for reg in pw.regime.unique():
        ls = _REGIME_LS.get(reg, "-")
        s = pw[pw.regime == reg]
        for m in METHODS:
            lbl = f"{m} [{reg}]" if multi else m
            ax.plot(s.delta_k, s[m], ls, marker="o", color=_MCOLOR[m], label=lbl)
    ax.axhline(.05, ls=":", c="k", lw=.8)
    ax.set_xlabel("true Δk (test − reference)"); ax.set_ylabel("rejection rate")
    ax.set_title("Power of the Δk contrast (δ=0 ⇒ Type-I)")
    ax.legend(fontsize=7); fig.tight_layout()
    fig.savefig(Path(outdir) / "mc_power.png", dpi=150); plt.close(fig)


def _plot_measure(stats, outdir):
    """Visualize the REAL variance spread across run(s): (a) the per-point σ the linear
    model weights on (the heteroscedasticity), (b) the dispersion = empirical/propagated
    variance (the fixed-vs-random-effects gap), (c) peptides/cell (the effective-df axis)."""
    plt = _mpl(outdir)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    colors = [f"C{i}" for i in range(len(stats))]
    for c, (name, s) in zip(colors, stats.items()):
        cs = s["sigma_cell"]; cs = cs[np.isfinite(cs) & (cs > 0)]
        ax[0].hist(np.log10(cs), bins=60, density=True, histtype="step", color=c,
                   label=f"{name}  (SD log σ = {s['spread_cell']:.2f})")
        ax[0].axvline(np.log10(np.median(cs)), color=c, ls=":", lw=.9)
    ax[0].set_xlabel(r"$\log_{10}\,\sigma_\theta$  (collapsed cell)")
    ax[0].set_ylabel("density")
    ax[0].set_title("per-point σ — the heteroscedasticity the WLS weights on")
    ax[0].legend(fontsize=7)
    for c, (name, s) in zip(colors, stats.items()):
        dd = s["disp"]; dd = dd[np.isfinite(dd) & (dd > 0)]
        if len(dd):
            ax[1].hist(np.log10(dd), bins=60, density=True, histtype="step", color=c,
                       label=f"{name}  (median {np.median(dd):.2f}×)")
    ax[1].axvline(0.0, color="k", ls="--", lw=.9)          # dispersion = 1 (calibrated)
    ax[1].set_xlabel(r"$\log_{10}$ (empirical / propagated var)")
    ax[1].set_title("dispersion — fixed- vs random-effects gap"); ax[1].legend(fontsize=7)
    for c, (name, s) in zip(colors, stats.items()):
        ppc = np.clip(s["ppc"], 1, 8)
        vals, cnt = np.unique(ppc, return_counts=True)
        ax[2].plot(vals, cnt / cnt.sum(), "o-", color=c,
                   label=f"{name}  ({(s['ppc'] == 1).mean() * 100:.0f}% single-peptide)")
    ax[2].set_xlabel("peptides / collapse cell  (≥8 clipped)")
    ax[2].set_ylabel("fraction of cells")
    ax[2].set_title("effective-df axis — single-peptide = low df"); ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(Path(outdir) / "measure_variance.png", dpi=150); plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=["diagnose", "recover", "power", "efficiency",
                                   "varweight", "moderate", "measure", "irls", "real", "all"])
    p.add_argument("--run", default="runs/lve_atr_clean",
                   help="a rollup run dir holding riana_rollup_fractions.txt")
    p.add_argument("--reference", default="control")
    p.add_argument("--test", default="atrium")
    p.add_argument("--sigma", type=float, default=SIGMA_THETA,
                   help=f"per-point θ-scale noise SD (default {SIGMA_THETA}, from lve_atr)")
    p.add_argument("--regime", choices=["homo", "hetero", "both"], default="homo",
                   help="noise regime for recover/power/efficiency: homoscedastic, "
                        "heteroscedastic (per-point σ spread), or BOTH side-by-side in one "
                        "table/plot. (diagnose/real read real data; varweight always sweeps "
                        "both internally.)")
    p.add_argument("--spread", type=float, default=0.7,
                   help="heteroscedastic σ_θ spread — lognormal, E[σ²] held fixed so a "
                        "homo-vs-hetero comparison isolates the heteroscedasticity. Used by "
                        "the 'hetero'/'both' regimes and varweight's hetero arm (default 0.7).")
    p.add_argument("--var-df", type=float, default=float("inf"), dest="var_df",
                   help="fidelity of the per-point Var(θ) estimate fed to wls_fit_var: inf = "
                        "oracle (true σ²), finite = noisy χ²(df)/df estimate. varweight's "
                        "noisy mode uses this (default 4).")
    p.add_argument("--df-lo", type=float, default=4.0, dest="df_lo",
                   help="[moderate] variance-estimate df for single-peptide cells (default 4).")
    p.add_argument("--df-hi", type=float, default=12.0, dest="df_hi",
                   help="[moderate] variance-estimate df for multi-peptide cells (default 12).")
    p.add_argument("--frac-lo", type=float, default=0.4, dest="frac_lo",
                   help="[moderate] fraction of single-peptide (low-df) cells (default 0.4).")
    p.add_argument("--between", type=float, default=0.5,
                   help="[moderate] between-peptide variance fraction (the dispersion source; "
                        "mean dispersion ≈ 1+between, default 0.5 → ~1.5×).")
    p.add_argument("--d0-grid", type=float, nargs="+", default=[2, 4, 8, 16], dest="d0_grid",
                   help="[moderate] eBayes prior-df (shrinkage) grid to sweep (default 2 4 8 16; "
                        "d0→∞ = wls_fit, d0→0 = raw wls_fit_var).")
    p.add_argument("--t-grid", type=float, nargs="+", default=[4, 8], dest="t_grid",
                   help="[moderate] hard df-threshold grid to sweep (default 4 8).")
    p.add_argument("--runs", nargs="+", default=None,
                   help="[measure] one or more run dirs to read/overlay (falls back to --run).")
    p.add_argument("--phi-limit", type=float, default=PHI_LIMIT)
    p.add_argument("--nsim", type=int, default=2000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--k0", type=float, default=0.075, help="[power] baseline k")
    p.add_argument("--k-grid", type=float, nargs="+", default=K_GRID, dest="k_grid")
    p.add_argument("--deltas", type=float, nargs="+",
                   default=[0.0, 0.01, 0.02, 0.04, 0.08])
    p.add_argument("--plot", metavar="DIR", default=None,
                   help="write PNGs to DIR (e.g. 36_out)")
    a = p.parse_args()

    def rule(t):
        print("\n" + "=" * 88 + f"\n{t}\n" + "=" * 88)

    if a.cmd in ("diagnose", "all"):
        rule("1. DIAGNOSE — is the noise FS-scale? (the premise)"); cmd_diagnose(a)
    if a.cmd in ("recover", "all"):
        rule("2. RECOVER — bias / coverage / Type-I"); cmd_recover(a)
    if a.cmd in ("power", "all"):
        rule("3. POWER — the Δk contrast"); cmd_power(a)
    if a.cmd in ("efficiency", "all"):
        rule("4. EFFICIENCY — vs the exact MLE"); cmd_efficiency(a)
    if a.cmd == "varweight":
        rule("4b. VARWEIGHT — per-point Var(θ) effect size + safety"); cmd_varweight(a)
    if a.cmd == "moderate":
        rule("4c. MODERATE — eBayes / threshold / empirical arms"); cmd_moderate(a)
    if a.cmd == "measure":
        rule("7. MEASURE — the real per-point-variance regime"); cmd_measure(a)
    if a.cmd == "irls":
        rule("5. IRLS — steps and the t=0 point"); cmd_irls(a)
    if a.cmd == "real":
        rule("6. REAL DATA — impact"); cmd_real(a)


if __name__ == "__main__":
    main()
