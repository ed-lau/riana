# `linear simple` Δk — the OLS-on-φ estimator is anti-conservative; use weighted LS

**Date:** 2026-07-13 · **Data:** `runs/lve_atr_clean` (atrium vs control, 686 proteins,
12 timepoints 0–30 d) · **Status: SHIPPED** — `--linear-weights wls` is the default
(`--linear-weights ols` restores the old fit for audit); t = 0 excluded. See `CHANGELOG.md`.
**Caveat 2026-11-12:** the Monte Carlo here calibrates against envelope noise only; peptide
offsets and animal variation are treated in `2026-11-12_deltak_inference_manual.md` (§3, §6).
**Appendix (2026-11-10):** the as-shipped code, mapped line by line to the R `lm` + `emtrends`
reference, with a numerical parity check against R — see the end of this report.

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
| **observed** residual SD | 0.060 | 0.093 | 0.147 | 0.233 | 0.384 | 0.603 |
| predicted, FS-scale `σ_θ/(1−θ)` | 0.065 | 0.083 | 0.117 | 0.193 | 0.382 | 1.64 |
| predicted, flat φ-scale | 0.064 | 0.064 | 0.064 | 0.064 | 0.064 | 0.064 |

Formal test — regress `log|resid|` on `−log(1−θ)`. A pure FS-scale noise model predicts a
slope of exactly **1**; a pure φ-scale (flat) model predicts **0**:

| t = 0 | plateau truncation | slope | 95 % CI |
|---|---|---|---|
| **excluded (what the fit sees)** | **φ > −4 (production)** | **0.800** | [0.770, 0.831] |
| excluded | none | 0.881 | [0.851, 0.912] |
| *included* | φ > −4 | *1.101* | *[1.066, 1.136]* |

**Read the first row.** The two lower rows are instructive but must not be quoted: including
t = 0 inflates the slope by **+0.30**, because those points sit at θ ≈ 0 (x ≈ 0) with
*artificially tiny* residuals — the θ-floor clamp pins them (see "Two implementation details"),
which drags the regression's left end down and steepens it. That is the very artifact this
report goes on to remove from the fit, so it cannot be used to justify the fit. The plateau
truncation then costs a further ~0.08 by censoring the largest residuals at the top end.

**So the honest exponent is ≈ 0.8 (0.88 untruncated), not 1.0.** The noise is therefore
**predominantly FS-scale but not purely so** — most consistent with a *mixture*,
`Var(φ) ≈ σ_θ²/(1−θ)² + σ_bio²`, i.e. dominant FS-scale measurement noise plus a smaller
φ-scale (biological / model-misfit) floor. σ_θ ≈ **0.056**.

This does **not** weaken the case for weighting — 0.80 is ~9 SE from the φ-scale prediction of
0, and the residual SD still spans a **10× range** that OLS treats as equal. It does mean the
pure-FS Monte-Carlo below is a mild **upper bound** on OLS's real-world Type-I inflation. The
load-bearing evidence is therefore the **real-data** validation against the nonlinear MLE
(next-but-one section), which assumes nothing about the exponent and is unaffected.

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
- **The weight exponent is slightly over-aggressive.** `(1−θ̂)²` is exactly right for *pure*
  FS-scale noise (exponent 1). The measured exponent is ≈ **0.88** (untruncated), so the
  variance-optimal weight is nearer `(1−θ̂)^1.8`. In practice this costs nothing measurable —
  `wls_fit` still lands on the nonlinear MLE at ratio **1.0005** — but **note that the MLE
  reference itself assumes pure FS-scale noise**, so their agreement confirms `wls_fit`
  approximates *that* estimator, it does not independently prove the exponent. Fitting the
  exponent (weight `(1−θ̂)^{2γ}`, γ ≈ 0.88) is a clean, cheap refinement. What is *not* in doubt
  is the direction: OLS implicitly assumes exponent **0**, and the data say **0.8–0.9**.
- Assumes θ noise is homoscedastic on the FS scale *within* a protein, but per-point θ
  precision genuinely varies with peptide depth. The ideal weight is `(1−θ̂)² / Var(θ_i)`;
  `riana_rollup_fractions.txt` does not currently carry a per-point `Var(θ)`. Adding it would
  be strictly better — a clean follow-up, and it may absorb part of the apparent φ-scale floor.
- Residual coverage is 0.93–0.95, not exactly 0.95, at the fastest k — the plateau truncation
  (φ > −4) and the θ clamp still censor the tail. With proper weighting the truncation is
  arguably redundant (the weights already down-weight the saturated tail smoothly); worth a
  `phi_limit` sensitivity sweep as a follow-up.
- Does **not** rescue genuinely saturated curves (θ ≈ 1 before the last timepoints): they are
  under-identified by design, whatever the estimator.

**Not recommended:** `weights=(1-θ_obs)²` (the naive inverse-variance form) — it fixes
calibration but biases k low by ~3 % overall and ~9 % in the fast tail.

---

## Appendix (2026-11-10) — the as-shipped implementation, the R reference, and a parity check

This appendix records what the code does, so the estimator can be audited against the R
workup it was ported from. Line numbers refer to **`riana/core/linear_model.py` at commit
`af20701`** (branch `1.2.0`). The rationale in the sections above is unchanged.

### The R reference

The model was ported from `data/notebook/03_R_linearmodel_reference.Rmd` (local, not in the
repository). Its per-protein fit is:

```r
fit_interaction <- lm(protein_clearance ~ 0 + day + day:treatment, data = subset_data)
em <- emtrends(fit_interaction, pairwise ~ treatment, var = "day")
...
p_adj_ISO = p.adjust(p_ISO, method = "BH")
```

### R ↔ Python, step by step

| Step | R reference | Python (`linear_model.py` @ `af20701`) |
|---|---|---|
| Clearance | `log(1 − FS_mean)`; FS rounded to 2 dp, then clamped to [0.01, 0.99] | `to_phi`: θ clipped to [0.001, 0.999], no rounding (L292) |
| Saturated tail | kept | dropped where φ ≤ `phi_limit` (default −4, θ ≈ 0.98) (L294) |
| t = 0 | kept | dropped, except under `ols` (L306) |
| Model | `lm(φ ~ 0 + day + day:treatment)`: reference slope + difference | `ols`/`wls("phi ~ 0 + day:C(condition)")`: one slope per condition (L333) |
| Weights | none (OLS) | `wls` (default): `(1 − θ̂)² = exp(2·max(φ̂, phi_limit))` from a first OLS pass, one step (L349–350). `ols`: none. `wls-var`: also ÷ moderated Var(θ) (opt-in) |
| Per-condition k and 95% CI | `emtrends(...)$emtrends`: `day.trend`, `lower.CL`, `upper.CL` | `res.params`, `res.conf_int()`; k = −slope, so the bounds swap (L373–374) |
| Δk, SE, p | `emtrends(..., pairwise ~ treatment)$contrasts`: t-test at residual df | `res.t_test(c)` with c = +1 test, −1 reference; Δk = −effect (L384–393) |
| Multiplicity | `p.adjust(p, "BH")` over one table (one strain × compartment) | BH over one p per protein, separately in each experiment (L440–449) |

The two model formulas span the same columns, so they give identical fitted values, residual
df, σ̂², per-condition slopes and contrasts; only the raw coefficients are labelled
differently. With more than two conditions, `emtrends(pairwise ~ …)` adjusts the pairwise p by
Tukey; Python's named `--reference-condition`/`--test-condition` pair is one unadjusted
contrast.

### Code excerpt (abridged — long comments trimmed, `…` marks omitted lines)

Per-condition preparation (θ → φ, plateau truncation, t = 0):

```python
290          for cond, cg in grp.groupby("condition", sort=True):
291              day0 = cg["labeling_time"].to_numpy()
292              phi = to_phi(cg["theta"].to_numpy(),
293                           theta_floor=theta_floor, theta_ceiling=theta_ceiling)
294              t, p = truncate_plateau(day0, phi, phi_limit=phi_limit)
 …
306              keep = np.ones(len(t), dtype=bool) if weights == "ols" else (t > 0)
307              t, p = t[keep], p[keep]
308              if len(t) >= min_points_per_condition:
309                  per_cond[str(cond)] = (t, p)
```

The joint fit — OLS, then one weighted step:

```python
333          formula = "phi ~ 0 + day:C(condition)"
334          try:
335              res = smf.ols(formula, data=fit_df).fit()
336              if weights == "wls":
 …
349                  w = np.exp(2.0 * np.maximum(res.fittedvalues.to_numpy(), phi_limit))
350                  res = smf.wls(formula, data=fit_df, weights=w).fit()
```

Per-condition k with CI, and the Δk contrast (the `emtrends` equivalent):

```python
373          slope = res.params
374          ci = res.conf_int()
 …
384          def _contrast(ref: str, other: str) -> tuple[float, float, float]:
385              names = list(slope.index)
386              c_vec = np.zeros(len(names))
387              c_vec[names.index(coef[other])] = 1.0   # slope_test
388              c_vec[names.index(coef[ref])] = -1.0    # − slope_reference
389              tt = res.t_test(c_vec)
390              # k = −slope, so Δk = k_test − k_ref = −(slope_test − slope_ref).
391              return (-float(np.ravel(tt.effect)[0]),
392                      float(np.ravel(tt.sd)[0]),
393                      float(np.ravel(tt.pvalue)[0]))
 …
402          elif (len(conds) == 2 and all(coef[c] in slope.index for c in conds)
403                and (reference_condition is None or reference_condition in conds)):
 …
408              ref = reference_condition if reference_condition is not None else conds[0]
409              other = next(c for c in conds if c != ref)
410              delta_k, delta_se, delta_p = _contrast(ref, other)
```

Benjamini–Hochberg, one p per protein, within each experiment:

```python
447          for exp, grp in prot_p.groupby("experiment"):
448              prot_p.loc[grp.index, "delta_k_p_adj"] = multipletests(
449                  grp["delta_k_p"].to_numpy(), method="fdr_bh")[1]
```

The θ points come from the per-(biological replicate, timepoint) inverse-variance collapse in
`riana/core/protein.py` (`_collapse_group`, L868), which `fit_linear_deltak` receives from
`_rollup_linear` (L433).

### Parity check against R (2026-11-10)

**Method.** `runs/linear_model_parity/export.py` writes the exact points `fit_linear_deltak`
fits (after φ, plateau and t = 0 filtering) for every chamber × protein in
`runs/timeseries_atf6` with both conditions, plus the Python outputs.
`runs/linear_model_parity/r_reference.R` fits the same points with the reference
`lm(φ ~ 0 + day + day:condition)` + `emtrends(pairwise ~ condition, var = "day")`, and for
`wls` with the same `lm` given the one-step weights. BH is applied per chamber in both. Because
R receives Python's filtered points, the check covers the estimator, CIs, contrast and BH, not
the data-prep differences in the table above.

**Result.** 3,510 chamber × protein fits (7,020 condition rows, 101,778 points):

| Output | `ols`: max abs diff | `wls`: max abs diff |
|---|---|---|
| k_deg | 1.7e-14 | 6.2e-14 |
| ci_lo / ci_hi | 1.1e-12 | 2.8e-12 |
| delta_k | 2.6e-15 | 1.5e-14 |
| delta_k_se | 4.3e-15 | 1.3e-14 |
| delta_k_p | 1.5e-14 | 1.5e-14 |
| delta_k_p_adj | 1.7e-14 | 3.3e-14 |
| Proteins with BH p < 0.05 | 883 (R: 883) | 37 (R: 37) |

The Python implementation reproduces R's `lm` + `emtrends` + BH to machine precision, both
unweighted and weighted.

A side result: the unweighted fit (the original R workflow) calls **883** chamber × protein Δk
significant on ATF6, against **37** for the shipped weighted fit. ATF6 is strongly
heteroscedastic DIA data (`reports/2026-07-24_linear_wls_per_point_var.md`), the regime where
the main report shows OLS is anti-conservative.

`wls-var` (opt-in) has no counterpart in the R reference and is not covered by this check.
