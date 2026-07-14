# `linear simple` Δk — the OLS-on-φ estimator is anti-conservative; use weighted LS

**Date:** 2026-07-13 · **Data:** `runs/lve_atr_clean` (atrium vs control, 686 proteins,
12 timepoints 0–30 d) · **Status: SHIPPED** — `--linear-weights wls` is the default
(`--linear-weights ols` restores the old fit for audit); t = 0 excluded. See `CHANGELOG.md`.

> **Notation.** θ and FS are the same quantity — the fraction new (RIANA calls it `theta` in
> `core/linear_model`, `fs` in the data columns). **σ_θ** is the SD of the measurement error in
> θ; the subscript names the *variable*, matching `Var(φ)`. **"FS-scale"** is used only as an
> adjective, for *which scale* the noise is homoscedastic on (FS-scale vs φ-scale) — that
> distinction is the whole point of this report. (The reference R workup calls σ_θ `SIGMA_FS`.)

## The problem

`core/linear_model.fit_linear_deltak` fits **unweighted** OLS through the origin on the
linearized clearance φ = log(1 − θ):

```python
smf.ols("phi ~ 0 + day:C(condition)", data=fit_df).fit()
```

θ carries roughly **homoscedastic measurement noise on the FS scale**, but φ is a *log*
transform of it, so by the delta method

    Var(φ) = Var(θ)·(dφ/dθ)² = σ_θ² / (1 − θ)²

i.e. the φ-residuals are strongly **heteroscedastic** — their SD blows up as θ → 1. OLS
assumes they are equal. This was already listed under Known Limitations; this report
quantifies it and shows the fix.

## The premise holds on RIANA's own data

Fitting the current OLS per protein on `lve_atr_clean` and binning the φ-residuals by θ:

| θ bin | 0–0.2 | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–0.9 | 0.95+ |
|---|---|---|---|---|---|---|
| **observed** residual SD | 0.061 | 0.094 | 0.147 | 0.235 | 0.387 | 0.612 |
| predicted, FS-scale `σ/(1−θ)` | 0.063 | 0.080 | 0.113 | 0.186 | 0.370 | 1.59 |
| predicted, flat φ-scale | 0.060 | 0.060 | 0.060 | 0.060 | 0.060 | 0.060 |

Formal test — regress `log|resid|` on `−log(1−θ)`: **slope = 1.101, 95 % CI [1.066, 1.136]**.
The FS-scale model predicts exactly 1; the flat-φ model predicts 0. The noise is **FS-scale**,
with σ_θ ≈ **0.056**. Residual SD spans a **10× range** that OLS treats as equal.

## Consequence (Monte Carlo through RIANA's actual pipeline)

Simulated with lve_atr's design (12 days), σ_θ = 0.056, and RIANA's exact transform chain
(clamp θ→[0.001, 0.999] → φ → truncate φ > −4 → joint through-origin fit). True Δk = 0, so
every rejection is a **false positive**. 2000 sims per cell.

| estimator | mean rel-bias | mean CI coverage | **mean Type-I** | worst Type-I |
|---|---|---|---|---|
| **`ols`** (current) | −0.105 | **0.543** | **0.285** | 0.342 |
| `wls_obs` — weights `(1−θ_obs)²` | −0.081 | 0.768 | 0.080 | 0.123 |
| **`wls_fit`** — weights `(1−θ̂)²` | **+0.003** | **0.926** | **0.062** | 0.070 |

**OLS rejects a true null ~28 % of the time at α = 0.05** (nominal 5 %), and its 95 % CIs
cover 54 %. It is also biased: **+5 % in the bulk, −26 % to −38 % for k ≥ 0.2**.

**The naive fix is not enough.** `wls_obs` (weights from the *observed* θ — the obvious
inverse-variance weight) fixes calibration but introduces a **systematic downward bias**: the
weights are correlated with the response noise, so points where noise pushed θ *up* get
down-weighted, dragging k low.

**The fix is fitted-value weights.** Take one IRLS step — fit OLS, then re-fit weighting by
`(1 − θ̂)² = exp(2·φ̂)` from the *fitted* values. The weights are then a smooth, monotone
function of time (`exp(−2k̂·t)`), independent of each point's own noise.

## Confirmed on real data against an independent reference

The nonlinear `FS = 1 − exp(−k·t)` fit is the correctly-specified MLE under FS-scale noise, so
it is an unbiased reference for the point estimate. On 1154 real (protein, condition) curves:

| estimator | ρ vs nonlinear | median k/k_nl | **bulk (k≤0.15)** | **fast tail (k>0.15)** |
|---|---|---|---|---|
| `ols` | 0.804 | 1.037 | **+4.7 %** | **−14.4 %** |
| `wls_obs` | 0.829 | 0.971 | −2.4 % | −8.6 % |
| **`wls_fit`** | **0.932** | **1.0004** | **+0.1 %** | **−0.5 %** |

`wls_fit` reproduces the nonlinear k almost exactly (IQR of the ratio ±0.7 %) at every speed —
no simulation involved. OLS's **−14 % fast-tail bias is real**, on 180 real curves.

## What the estimator actually is (WLS vs IRLS vs MLE)

Worth stating precisely, because the three are easy to conflate:

- **WLS** is the *estimator*: minimize `Σ wᵢ(φᵢ − xᵢβ)²` for **given** weights. With the true
  inverse variances it is optimal. But the true variance `σ_θ²/(1−θᵢ)²` depends on the
  **true** θᵢ, which we do not have — so the weights must themselves be estimated, and *how*
  decides everything.
- **IRLS** is not a different estimator, it is the *algorithm* for that chicken-and-egg
  (fit → weights from the fit → re-fit), the same one every GLM is fitted with. "One IRLS step"
  = OLS → weights → WLS → stop.
- **The endogeneity.** `w = (1−θᵢ_obs)²` makes the weight a function of that point's **own
  error**, so down-weighting high-θ points preferentially discards the most-negative φ and
  flattens the slope (k biased low). `w = (1−θ̂ᵢ)² = exp(−2k̂·tᵢ)` depends only on `k̂` (from all
  n points) and `tᵢ` (noise-free design), so its correlation with any single εᵢ is O(1/n) and
  the bias vanishes. This is exactly the **feasible-GLS** rule: build weights from a consistent
  first-stage fit, never from the raw response. It is also why **one step suffices** — the job
  is to make the weights exogenous.

**Are we doing MLE? No — a first-order approximation to it, at ~97 % of its efficiency.** Under
homoscedastic Gaussian noise on the FS scale, the *exact* MLE is plain nonlinear least squares
on FS (`minimize Σ (θᵢ − (1−e^{−k·tᵢ}))²`) — i.e. RIANA's **nonlinear `simple` model**. `wls_fit`
is the delta-method linearization of it. Efficiency against that exact MLE (4000 sims/cell):

| k | estimator | rel-bias | RMSE | **RMSE ÷ RMSE(MLE)** |
|---|---|---|---|---|
| 0.050 | MLE (nonlinear LS) | +0.003 | 0.0032 | 1.00 |
| | `ols` | +0.018 | 0.0045 | **1.41** |
| | **`wls_fit`** | +0.012 | 0.0033 | **1.03** |
| 0.075 | MLE (nonlinear LS) | +0.001 | 0.0047 | 1.00 |
| | `ols` | +0.023 | 0.0085 | **1.80** |
| | **`wls_fit`** | +0.010 | 0.0049 | **1.03** |
| 0.150 | MLE (nonlinear LS) | +0.004 | 0.0099 | 1.00 |
| | `ols` | −0.133 | 0.0269 | **2.73** |
| | **`wls_fit`** | +0.001 | 0.0103 | **1.04** |
| 0.300 | MLE (nonlinear LS) | +0.003 | 0.0223 | 1.00 |
| | `ols` | **−0.389** | 0.1303 | **5.84** |
| | **`wls_fit`** | −0.022 | 0.0264 | **1.18** |

`wls_fit` gives up only **3–4 % of RMSE** to the exact MLE across the bulk (18 % in the fast
tail); OLS throws away **40–480 %** *and* is badly biased.

**Then why not just use the exact MLE?** Because the linear model does not exist to estimate k —
it exists to **test Δk**. The nonlinear `simple` model fits each condition *independently* and
yields no contrast. The linearization is what buys one **joint** fit with a **shared residual
variance**, a **closed-form covariance** → an exact t-test on the slope difference, CIs and BH —
all essentially free in a linear model and fiddly in a nonlinear one. `wls_fit` keeps every bit
of that machinery while recovering MLE-grade estimates, for ~4 lines.

*The purist alternative* (filed, not built): a **joint nonlinear fit** of both conditions with a
shared σ and a **Wald test** on `k_test − k_ref` from the Jacobian covariance — the exact MLE
*with* exact inference. It costs real machinery (joint NLS + delta-method contrast, plus
`curve_fit` convergence risk across ~1000 proteins) to buy back that last ~3 %.

## Two implementation details

- **One IRLS step is optimal.** 2 and 5 iterations do not improve bias/coverage and slightly
  *degrade* Type-I at high k (0.086 → 0.106) — the weights begin chasing the noise.
- **Drop the t = 0 point from the linear fit.** At t = 0 the true θ is 0, so ~half the measured
  θ go negative and are clamped to 0.001 → φ ≈ −0.001, while the through-origin model predicts
  exactly 0. Its residual is *artificially* ≈ 0, which deflates σ̂² and shrinks every SE. In a
  through-origin fit t = 0 has **zero leverage on the slope**, so dropping it costs nothing.
  Effect at k = 0.216: coverage 0.892 → **0.929**, Type-I 0.086 → **0.069**.

## Recommendation

Swap the OLS for **WLS with fitted-value weights (one IRLS step), excluding t = 0**:

```python
res0 = smf.ols("phi ~ 0 + day:C(condition)", data=fit_df).fit()
w    = np.exp(2.0 * res0.fittedvalues)          # = (1 - theta_hat)^2
res  = smf.wls("phi ~ 0 + day:C(condition)", data=fit_df, weights=w).fit()
```

`t_test` / `conf_int` work unchanged on the weighted fit, so the Δk contrast, CIs and BH
correction need no other change.

**Pros.** Type-I 0.285 → 0.062; coverage 0.543 → 0.926; the fast-tail bias (−14 % on real data)
disappears; **lowest RMSE at every k**; and because α is now honest it has *more genuine power*
(at Δk = 0.04: 0.994 vs OLS 0.847 — OLS's apparent sensitivity is largely its false-positive
rate). The linear k finally agrees with the nonlinear k (ρ 0.93, ratio 1.000), which is what a
*linearization* should do — the model exists for the joint Δk contrast, not to give a different k.

**Cons / caveats.**
- **Results-affecting default change** — needs announcing like `deberneh_2025_rss` / the Spep
  floor. On lve_atr the significant-protein count moves 329 → 288 (79 % of OLS hits survive;
  29 new ones appear). **This does not overturn the biology** — atrium really does turn over
  faster than ventricle, and the large effects survive. What changes is that the *p-values no
  longer overstate the evidence*, and the fast-tail k's stop being ~14 % low.
- Assumes θ noise is homoscedastic on the FS scale. It is *within* a protein, but per-point θ
  precision genuinely varies with peptide depth. The ideal weight is `(1−θ̂)² / Var(θ_i)`;
  `riana_rollup_fractions.txt` does not currently carry a per-point `Var(θ)`. Adding it would
  be strictly better — a clean follow-up.
- Residual coverage is 0.93–0.95, not exactly 0.95, at the fastest k — the plateau truncation
  (φ > −4) and the θ clamp still censor the tail. With proper weighting the truncation is
  arguably redundant (the weights already down-weight the saturated tail smoothly); worth a
  `phi_limit` sensitivity sweep as a follow-up.
- Does **not** rescue genuinely saturated curves (θ ≈ 1 before the last timepoints): they are
  under-identified by design, whatever the estimator.

**Not recommended:** `weights=(1-θ_obs)²` (the naive inverse-variance form) — it fixes
calibration but biases k low by ~3 % overall and ~9 % in the fast tail.
