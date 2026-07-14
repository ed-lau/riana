"""`linear simple` Δk — OLS vs WLS: the diagnostic + Monte-Carlo harness.

The workup behind `reports/2026-07-13_linear_model_wls.md` and the switch of the
`linear simple` rollup from unweighted OLS to WLS (`--linear-weights wls`, the default).

THE ARGUMENT, in one line. θ (fraction new) carries roughly homoscedastic *measurement*
noise on the FS scale, but the model fits φ = log(1−θ), a **log** of it. By the delta
method ``Var(φ) = σ_θ²/(1−θ)²`` — so the φ-residuals are strongly **heteroscedastic**, their
SD blowing up as θ → 1. An unweighted OLS treats that ~10× SD range as equal, which makes the
Δk test badly anti-conservative and biases k low in the fast tail.

THREE ESTIMATORS are compared throughout:
  ols      unweighted (the pre-2026-07 default)
  wls_obs  weights = (1−θ_obs)²   — the naive inverse-variance form. Fixes calibration but
           BIASES k LOW: the weight is a function of that point's OWN error, so points noise
           pushed high get down-weighted, discarding the most-negative φ and flattening the
           slope. (Endogenous weights.)
  wls_fit  weights = (1−θ̂)² = exp(2·φ̂)  — from the FITTED value (one IRLS step). The weight
           then depends only on k̂ (all n points) and t (noise-free), so it is exogenous and
           the bias vanishes: the feasible-GLS rule. THIS IS WHAT PRODUCTION NOW DOES
           (`core.linear_model.fit_linear_deltak(weights="wls")`).

All simulation is faithful to RIANA's real chain: clamp θ→[0.001, 0.999] → φ = log(1−θ) →
truncate φ > `--phi-limit` → drop t=0 → joint through-origin fit `phi ~ 0 + day:C(condition)`.

Subcommands
-----------
  diagnose    Is the noise FS-scale or φ-scale?  (real data; the PREMISE the whole thing rests on)
  recover     MC: bias / CI coverage / Type-I across a k grid, all three estimators
  power       MC: power curve for the Δk contrast (the δ=0 row is Type-I)
  efficiency  MC: how close is each estimator to the EXACT MLE (nonlinear LS on FS)?
  irls        MC: how many IRLS steps, and does dropping t=0 matter?
  real        Real-data impact on a run: k, Δk, significance, agreement with the nonlinear MLE
  all         diagnose + recover + power + efficiency

Usage
-----
    python -m tests.benchmark.bench_linear_weights all
    python -m tests.benchmark.bench_linear_weights diagnose --run runs/lve_atr_clean
    python -m tests.benchmark.bench_linear_weights recover --sigma 0.08 --nsim 500
    python -m tests.benchmark.bench_linear_weights real --run runs/lve_atr_clean \
        --reference control --test atrium
    python -m tests.benchmark.bench_linear_weights all --plot 36_out   # writes PNGs

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
K_GRID = [0.03, 0.05, 0.075, 0.11, 0.216, 0.30]   # lve_atr 5/25/50/75/95 pct + a fast one
METHODS = ("ols", "wls_obs", "wls_fit")


# --------------------------------------------------------------------------- #
# the estimators
# --------------------------------------------------------------------------- #
def _design(day, cond):
    """Through-origin `day:condition` design (one slope column per condition)."""
    X = np.zeros((len(day), 2))
    X[cond == 0, 0] = day[cond == 0]
    X[cond == 1, 1] = day[cond == 1]
    return X


def fit(day, cond, phi, theta, how, phi_limit=PHI_LIMIT):
    """Joint through-origin fit; returns per-condition k, a CI for k(A), and the Δk test."""
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
    else:
        raise ValueError(how)
    k = -m.params                                            # k = −slope
    ci = m.conf_int()
    tt = m.t_test(np.array([-1.0, 1.0]))                     # Δk = k_B − k_A
    return dict(kA=k[0], kB=k[1], kA_lo=-ci[0, 1], kA_hi=-ci[0, 0],
                dk=-float(np.ravel(tt.effect)[0]),
                dk_p=float(np.ravel(tt.pvalue)[0]))


def mle_nls(day, theta):
    """The EXACT MLE under FS-scale Gaussian noise: nonlinear LS on θ, no transform.
    (This is RIANA's nonlinear `simple` model — the unbiased reference.)"""
    try:
        return float(curve_fit(lambda x, k: 1 - np.exp(-k * x), day, theta,
                               p0=[0.05], bounds=(0, 5), maxfev=5000)[0][0])
    except Exception:
        return np.nan


def simulate(kA, kB, sigma, rng, days=DAYS, phi_limit=PHI_LIMIT, drop_t0=True):
    """One two-condition dataset through RIANA's exact chain. Returns the linear-space
    arrays plus the raw (day, θ) for the nonlinear MLE."""
    day = np.concatenate([days, days])
    cond = np.concatenate([np.zeros(len(days), int), np.ones(len(days), int)])
    ktrue = np.where(cond == 0, kA, kB)
    theta_raw = 1 - np.exp(-ktrue * day) + rng.normal(0, sigma, len(day))  # FS-scale noise
    phi = to_phi(theta_raw)                                   # RIANA clamp + log
    keep = phi > phi_limit                                    # RIANA plateau truncation
    if drop_t0:
        keep &= day > 0                                       # zero leverage + clamp artifact
    return (day[keep], cond[keep], phi[keep],
            np.clip(theta_raw[keep], 0.001, 0.999), day, cond, theta_raw)


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
    for k in a.k_grid:
        acc = {m: {"k": [], "cov": [], "rej": []} for m in METHODS}
        for _ in range(a.nsim):
            day, cond, phi, th, *_ = simulate(k, k, a.sigma, rng, phi_limit=a.phi_limit)
            for m in METHODS:
                r = fit(day, cond, phi, th, m, a.phi_limit)
                if r is None:
                    continue
                acc[m]["k"].append(r["kA"])
                acc[m]["cov"].append(r["kA_lo"] <= k <= r["kA_hi"])
                acc[m]["rej"].append(r["dk_p"] < 0.05)
        for m in METHODS:
            kk = np.array(acc[m]["k"])
            rows.append(dict(k_true=k, method=m, n=len(kk),
                             rel_bias=(kk.mean() - k) / k,
                             rmse=float(np.sqrt(np.mean((kk - k) ** 2))),
                             coverage=float(np.mean(acc[m]["cov"])),
                             type1=float(np.mean(acc[m]["rej"]))))
    rec = pd.DataFrame(rows)
    print("RECOVERY / CI COVERAGE / TYPE-I   (true Δk = 0; nominal coverage .95, Type-I .05)\n")
    for m in METHODS:
        print(f"--- {m} ---")
        print(rec[rec.method == m][["k_true", "rel_bias", "rmse", "coverage", "type1"]]
              .to_string(index=False, float_format=lambda x: f"{x:8.3f}"), "\n")
    print("SUMMARY over the k grid")
    print(rec.groupby("method").agg(
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
    for d_k in a.deltas:
        acc = {m: [] for m in METHODS}
        for _ in range(a.nsim):
            day, cond, phi, th, *_ = simulate(a.k0, a.k0 + d_k, a.sigma, rng,
                                              phi_limit=a.phi_limit)
            for m in METHODS:
                r = fit(day, cond, phi, th, m, a.phi_limit)
                if r:
                    acc[m].append(r["dk_p"] < 0.05)
        rows.append(dict(delta_k=d_k, **{m: float(np.mean(acc[m])) for m in METHODS}))
    pw = pd.DataFrame(rows)
    print(f"POWER of the Δk test at k0={a.k0} (α=0.05). The delta_k=0 row is TYPE-I.\n")
    print(pw.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
    print("\nNote: OLS's apparent sensitivity at small δ is largely its false-positive rate.")
    if a.plot:
        _plot_power(pw, a.plot)
    return pw


# --------------------------------------------------------------------------- #
# 4. efficiency — vs the EXACT MLE (nonlinear LS on the FS scale)
# --------------------------------------------------------------------------- #
def cmd_efficiency(a):
    rng = np.random.default_rng(a.seed)
    rows = []
    for k in a.k_grid:
        acc = {m: [] for m in ("mle", *METHODS)}
        for _ in range(a.nsim):
            day, cond, phi, th, raw_day, raw_cond, raw_th = simulate(
                k, k, a.sigma, rng, phi_limit=a.phi_limit)
            acc["mle"].append(mle_nls(raw_day[raw_cond == 0], raw_th[raw_cond == 0]))
            for m in METHODS:
                r = fit(day, cond, phi, th, m, a.phi_limit)
                acc[m].append(r["kA"] if r else np.nan)
        A = pd.DataFrame(acc).dropna()
        rmse_mle = float(np.sqrt(np.mean((A["mle"] - k) ** 2)))
        for m in ("mle", *METHODS):
            r = float(np.sqrt(np.mean((A[m] - k) ** 2)))
            rows.append(dict(k_true=k, estimator=m,
                             rel_bias=(A[m].mean() - k) / k, rmse=r,
                             rmse_vs_mle=r / rmse_mle))
    eff = pd.DataFrame(rows)
    print("EFFICIENCY vs the EXACT MLE (nonlinear LS on FS = RIANA's `simple` model).")
    print("The MLE is the correctly-specified estimator under FS-scale noise.\n")
    print(eff.to_string(index=False, float_format=lambda x: f"{x:9.3f}"))
    print("\nwls_fit is the delta-method LINEARIZATION of the MLE — it should sit at ~1.0.")
    return eff


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
    for how in METHODS:
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
    plt = _mpl(outdir)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    for m in METHODS:
        s = rec[rec.method == m]
        axes[0].plot(s.k_true, s.rel_bias, "o-", label=m)
        axes[1].plot(s.k_true, s.coverage, "o-", label=m)
        axes[2].plot(s.k_true, s.type1, "o-", label=m)
    axes[0].axhline(0, ls="--", c="k", lw=.8); axes[0].set_title("relative bias in k")
    axes[1].axhline(.95, ls="--", c="k", lw=.8); axes[1].set_title("95% CI coverage")
    axes[2].axhline(.05, ls="--", c="k", lw=.8); axes[2].set_title("Type-I of the Δk test")
    for ax in axes:
        ax.set_xlabel("true k (/day)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(Path(outdir) / "mc_recover.png", dpi=150); plt.close(fig)


def _plot_power(pw, outdir):
    plt = _mpl(outdir)
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for m in METHODS:
        ax.plot(pw.delta_k, pw[m], "o-", label=m)
    ax.axhline(.05, ls=":", c="k", lw=.8)
    ax.set_xlabel("true Δk (test − reference)"); ax.set_ylabel("rejection rate")
    ax.set_title("Power of the Δk contrast (δ=0 ⇒ Type-I)")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(Path(outdir) / "mc_power.png", dpi=150); plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=["diagnose", "recover", "power", "efficiency",
                                   "irls", "real", "all"])
    p.add_argument("--run", default="runs/lve_atr_clean",
                   help="a rollup run dir holding riana_rollup_fractions.txt")
    p.add_argument("--reference", default="control")
    p.add_argument("--test", default="atrium")
    p.add_argument("--sigma", type=float, default=SIGMA_THETA,
                   help=f"per-point θ-scale noise SD (default {SIGMA_THETA}, from lve_atr)")
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
    if a.cmd == "irls":
        rule("5. IRLS — steps and the t=0 point"); cmd_irls(a)
    if a.cmd == "real":
        rule("6. REAL DATA — impact"); cmd_real(a)


if __name__ == "__main__":
    main()
