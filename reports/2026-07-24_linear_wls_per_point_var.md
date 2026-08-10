# `linear simple` Δk — per-point Var(θ) WLS: a large but conditional gain; moderate the variance first

**Date:** 2026-07-24 · **Data:** synthetic Monte-Carlo + real `runs/{lve_atr_clean,
boomi_ipsc_d2o, juber_ac16_d2o}` · **Status (as written): MEASURED — not shipped.**
**UPDATE 2026-07-30: the eBayes-moderated path SHIPPED as opt-in `--linear-weights wls-var`**
(`c82cc55`; `fs_var`/`fs_df` columns) — `wls` stays the default. The effect is real
and large, but the naive estimate inflates Type-I; a variance-moderation step is required
before it can become a default. **The shipped `wls` (IRLS) remains the standing default and
recommendation.** Bench: `tests/benchmark/bench_linear_weights.py` (`varweight` / `measure`).

> **Notation.** Continues `reports/2026-07-13_linear_model_wls.md`. θ = fraction new (`fs` in
> the data), φ = log(1 − θ), the linearized clearance the Δk model fits (`φ ~ 0 + day:C(cond)`).
> **σ_θ** = SD of θ's measurement error. The shipped `wls` weights each point by `(1 − θ̂)²`
> (`exp(2·φ̂)`, from the *fitted* value); this note asks whether to add a per-point `1/Var(θ_i)`.
> **PI** = the per-point **prediction interval** (`fs_lower`/`fs_upper`, from the fit's residual
> bootstrap — the model value plus a resampled residual, so its width reads that point's
> measurement scatter). Hence `σ_θ,i ≈ (fs_upper − fs_lower)/3.29` (a 5–95 % band spans 3.29·σ)
> and **Var(θ_i) = σ_θ,i²** is the "PI-width variance" — the quantity the rollup's inverse-variance
> collapse already uses. (It is a *prediction* interval on θ, not a *confidence* interval on k̂ —
> the latter is `ci_lo`/`ci_hi`, which feeds the `k_cv` gate.)

## TL;DR

- The ideal (MLE) weight factors *exactly* into `(1 − θ̂_i)² / Var(θ_i)`: the shipped `wls`
  carries the first factor (the log-transform heteroscedasticity) and **assumes the second
  constant**. `wls_fit_var` adds it — the empirical per-point precision.
- **Production is strongly heteroscedastic**, and the collapse *amplifies* it: the per-point σ
  the linear model weights on spans **9–16× (p90/p10)** across cells, SD(log σ) ≈ **0.83–1.04** —
  well above the bench's 0.7 default. So the prize is real: under this regime the shipped
  `wls_fit` sits at **~4× the weighted-MLE RMSE**; an oracle `wls_fit_var` recovers essentially
  all of it (**+74% RMSE, ~1.0× the MLE**) with Type-I and coverage intact.
- **But the per-point variance is not known — it is estimated, and noisily.** Effective df ≈
  **6–9** (40–44 % of cells are single-peptide), and the PI-width σ's are **under-calibrated by
  1.25–1.72×** (peptides scatter more than their σ claims). Fed that estimate raw, `wls_fit_var`
  keeps most of the gain (+69–71 %) but **inflates Type-I to ~0.07 and drops coverage to ~0.91**.
- **Verdict:** big, worthwhile efficiency prize, but the naive weight is **not** a safe default —
  and the fix is now **settled by the bench**: **eBayes-moderated** per-point variance keeps ~+58 %
  of the +70 % at Type-I ≈ 0.01–0.02 (safely conservative), while a hard df-threshold and the raw
  empirical variance both **fail**. The shipped `wls` (IRLS) remains the standing default; `wls-var`
  with eBayes moderation is the validated path to evaluate for a future opt-in.

## The candidate

The Δk model fits φ = −k·t; the optimal weight is inverse variance, `w = 1/Var(φ)`, and the
delta method splits it with no approximation:

```
Var(φ_i) = (dφ/dθ)²·Var(θ_i) = Var(θ_i)/(1 − θ_i)²
   ⇒   w_i = (1 − θ̂_i)²  ·  1/Var(θ_i)
             └ transform ┘   └ per-point precision ┘
```

- **`wls_fit` (shipped)** = `(1 − θ̂_i)²` only, from the *fitted* θ̂ (one IRLS step). It corrects
  the heteroscedasticity the log **creates**, and assumes `Var(θ_i) = σ²` constant.
- **`wls_fit_var` (candidate)** = `(1 − θ̂_i)² / Var(θ_i)`. Same one-step pilot for the transform
  factor (so the *only* difference is `/Var`), times the **exogenous** per-point precision. The
  precision is safe to take from the data (it is *how noisy* the point is, independent of *which
  direction* it deviated), unlike weighting by observed θ, which biases k low.

`Var(θ_i)` is the peptide-point's PI-width variance `((fs_upper − fs_lower)/3.29)²` — the same
quantity `core.protein._weighted_theta` already inverse-variance-collapses, then discards. Under
the default `weighted` rollup the linear model's point is a `(protein, biorep, timepoint)` **cell**,
so its variance is `1/Σ(1/σ²_pep)` — the collapse's own byproduct.

## The estimators — one machinery, four weights

All four fits are the **same generalized-least-squares estimator** on the through-origin design
$\varphi = X\beta + \varepsilon$ (one slope column per condition; $\beta = -k$, so a $\Delta k$ is a
contrast $L\hat\beta$). They differ **only** in the diagonal weight matrix
$W=\operatorname{diag}(w_1,\dots,w_N)$:

$$
\hat\beta = (X^\top W X)^{-1} X^\top W \varphi,
\qquad
\hat\sigma^2 = \frac{(\varphi - X\hat\beta)^\top W (\varphi - X\hat\beta)}{N-p},
\qquad
\widehat{\operatorname{Var}}(\hat\beta) = \hat\sigma^2\,(X^\top W X)^{-1},
$$

$$
t \;=\; \frac{L\hat\beta}{\sqrt{L\,\widehat{\operatorname{Var}}(\hat\beta)\,L^\top}},
\qquad \textbf{df} = N-p \quad\text{(for every method).}
$$

$N$ = collapsed $(\text{condition},\text{biorep},t)$ points for the protein, $p$ = number of
conditions, and $\hat\sigma^2$ is the **single** residual scale estimated with $N-p$ df — so **only
the *relative* weights matter**: rescaling all $w_i$ by a constant leaves $\hat\beta$, the SE, and
$t$ unchanged (verified). Only $w_i$ changes across methods:

| method | weight $w_i$ | notes |
|---|---|---|
| **OLS** | $1$ | $W=I$; anti-conservative on φ (which is heteroscedastic) |
| **WLS (naive)** `wls_obs` | $(1-\theta_i)^2$ — *observed* $\theta$ | endogenous → biases $k$ **low** |
| **IRLS** `wls_fit` *(shipped)* | $(1-\hat\theta_i)^2 = e^{2\hat\varphi_i}$ — *fitted* $\hat\theta$ | feasible-GLS, exogenous → unbiased |
| **WLS-Var** *(candidate)* | $\dfrac{(1-\hat\theta_i)^2}{\widetilde{\operatorname{Var}}(\theta_i)}$ | + moderated per-point precision |

The transform factor is the delta-method inverse variance of the log:
$\varphi=\log(1-\theta)\Rightarrow \operatorname{Var}(\varphi)=\operatorname{Var}(\theta)/(1-\theta)^2$,
so the ideal weight is $(1-\theta)^2/\operatorname{Var}(\theta)$. **IRLS is one step:** pilot
$\tilde\beta=(X^\top X)^{-1}X^\top\varphi$, $\hat\varphi = X\tilde\beta$, then
$w_i=\exp\!\big(2\max(\hat\varphi_i,\varphi_{\lim})\big)$ — taking $\hat\theta$ from the fit (not the
observed $\theta$) is what keeps the weight exogenous.

`WLS-Var`'s per-point variance is **eBayes-moderated** (arm #1) — a precision-weighted blend of the
point's own estimate and the pooled variance:

$$
\widetilde{\operatorname{Var}}(\theta_i)
   = \frac{d_0\,s_0^2 + d_i\,\widehat{\operatorname{Var}}(\theta_i)}{d_0 + d_i},
\qquad s_0^2 = \text{pooled }\widehat{\operatorname{Var}},
$$

with $d_i$ = the *variance estimate's* own df and $d_0$ = the Smyth prior df ($d_0\to\infty$ recovers
`wls_fit`; $d_0\to 0$ the raw `wls_fit_var`). **This enters only $w_i$ — the contrast df stays
$N-p$.** The alternative that *would* move the df — moderating each protein's $\hat\sigma^2$ across
proteins with a Satterthwaite contrast df, the true limma analog — is a separate roadmap item.

## Bench design (the homo/hetero split)

`bench_linear_weights.py` runs every dataset through RIANA's exact chain (clamp θ → φ = log(1−θ)
→ truncate at `--phi-limit` → drop t=0 → joint through-origin fit). Extended here with:

- **`--regime {homo,hetero,both}` + `--spread`** — a heteroscedastic regime that draws per-point
  σ_i from a lognormal spread (`E[σ²]` held fixed, so homo-vs-hetero isolates the
  heteroscedasticity, not the noise level). First-class across `recover`/`power`/`efficiency`
  (`both` overlays the scenarios; colour = estimator, linestyle = regime).
- **`--var-df`** — the fidelity of the variance *estimate* fed to `wls_fit_var`: `inf` = oracle
  (true σ²), finite = an unbiased but noisy `σ²·χ²(df)/df`, as a PI-width from ~df observations.
- **`wmle`** — a **weighted** nonlinear MLE (`1/Var` weights) — the correctly-specified reference
  under heteroscedastic FS-scale noise (the plain `mle` is only optimal when σ_θ is constant).
- **`varweight`** — sweeps `{homo, hetero} × {oracle, noisy}` and reports the effect size.

Sanity: **homoscedastic + oracle ⇒ `wls_fit_var` is byte-identical to `wls_fit`** (the `/Var`
factor cancels), confirmed in every run. The extension is otherwise **purely additive**: in the
default homoscedastic regime the `ols` / `wls_obs` / `wls_fit` numbers are **byte-identical to the
pre-extension script** (verified by an exact element-wise comparison of the returned tables from
the same seed), so the shipped-WLS evidence in `reports/2026-07-13_linear_model_wls.md` reproduces
unchanged.

## Effect size — synthetic sweep (spread 0.7)

`varweight --nsim 500 --seed 3`, mean over k = {0.05, 0.075, 0.11, 0.216}. `wls_fit` baseline:
Type-I ≈ 0.05, coverage ≈ 0.94, and **rmse_vs_wmle ≈ 2.3×** (it leaves that much on the table).

**Heteroscedastic regime — the gain and its Type-I cost by variance fidelity:**

| Var estimate | RMSE gain vs `wls_fit` | `wls_fit_var` Type-I | coverage |
|---|---|---|---|
| oracle (∞) | **+52%** | 0.053 | 0.944 |
| df = 40 | +53% | 0.074 | 0.939 |
| df = 16 | +50% | 0.061 | 0.937 |
| df = 8 | +45% | 0.084 | 0.903 |
| df = 4 | +33% | **0.120** | 0.852 |

**Homoscedastic regime — the cost of using it when it is not needed:**

| Var estimate | RMSE Δ | Type-I |
|---|---|---|
| oracle | 0.0 % (identical — sanity ✓) | 0.052 |
| df = 16 | −3 % | 0.064 |
| df = 8 | −12 % | 0.090 |
| df = 4 | −31 % | **0.139** |

Reading: the efficiency gain is large and survives a noisy estimate, but **Type-I safety is
governed by the estimate's df** — only ≥ ~16 is clean; a df ≈ 4 estimate is a Type-I disaster in
*both* regimes. So the whole question reduces to: *how well does production estimate per-point
variance?*

## The real regime — `measure`

`bench_linear_weights.py measure` reads each run's **peptide-level** `riana_fit_fractions.txt`
(σ from the fs_lower/fs_upper PI) and reports where production actually sits:

| dataset | cell SD(log σ) | p90/p10 (cell) | timepoints/peptide (df) | single-peptide cells | dispersion |
|---|---|---|---|---|---|
| `lve_atr_clean` | **1.03** | 12× | median 9 (p10 6) | 41 % | **1.72** |
| `boomi_ipsc_d2o` | **0.83** | 9× | median 6 (p10 3) | 40 % | 1.25 |
| `juber_ac16_d2o` | **1.04** | 16× | median 6 (p10 5) | 44 % | 1.49 |

Three findings, consistent across datasets:

1. **Strongly heteroscedastic, and the collapse amplifies it.** The cell-level spread (0.83–1.04)
   *exceeds* the raw peptide spread (0.64–0.88) — the inverse-variance collapse pools a variable
   number of peptides per cell (median 2, but 1→many), so a 1-peptide cell and a 10-peptide cell
   differ enormously in precision. The linear model weights on a **9–16×** range of σ. This is
   well past the bench's 0.7 default → the prize is *bigger* than the 0.7 table shows.
2. **Moderate, risky effective df.** Each peptide's σ is bootstrapped over its ~6–9 timepoints, and
   **40–44 % of cells are single-peptide** (so their cell variance inherits that one peptide's df,
   p10 as low as 3). Effective df ≈ **6–9** — the borderline-to-risky band of the sweep above.
3. **Under-calibrated σ.** In multi-peptide cells, peptides scatter **1.25–1.72×** more than their
   stated σ predicts → the PI widths *under-state* variance, so the formal `1/Σ(1/σ²)` cell
   variance is **over-confident** (a direct cause of under-coverage, and exactly what a
   dispersion / reduced-χ² correction fixes).

## Synthesis — effect size at the measured regime

`varweight --spread 1.0` (the measured cell spread) at the measured df:

| regime · Var estimate | RMSE gain vs `wls_fit` | Type-I `fit → var` | coverage `fit → var` |
|---|---|---|---|
| hetero · **oracle** | **+74 %** | 0.034 → 0.050 | 0.945 → 0.943 |
| hetero · **df = 8** | +71 % | 0.044 → **0.069** | 0.944 → 0.912 |
| hetero · **df = 6** | +69 % | 0.036 → **0.074** | 0.947 → 0.916 |

Under production's own heteroscedasticity, the shipped `wls_fit` sits at **~4× the weighted-MLE
RMSE** (`efficiency --spread 1.0`); an oracle `wls_fit_var` recovers **all** of it at nominal
Type-I/coverage. With the *realistic* noisy estimate the efficiency prize is nearly intact
(**+70 %**), but Type-I runs **~0.07** and coverage **~0.91** — an over-rejection of ~40 %
relative, and the under-coverage the dispersion measurement predicts.

## How the θ interval is built — and why it under-states

`fs_lower` / `fs_upper` are a **within-peptide** prediction interval. The fit's residual bootstrap
(`core/fitting._fit_one_concat`) resamples each peptide's own fit residuals, refits `k*`, and takes
the 5–95 % band of `model(t_i; k*) + a resampled residual` — so the width measures how well a
peptide fits **its own** curve (its measurement scatter). The rollup collapse
(`core/protein._weighted_theta`) turns that into `σ_i = (hi − lo)/3.29` and inverse-variance-combines
the peptides in a cell, so the propagated cell variance `1/Σ(1/σ²)` is a **fixed-effects,
measurement-only** quantity.

The empirical between-peptide scatter is larger because it also carries a **peptide-heterogeneity**
component the PI structurally cannot see: peptides of one protein at one `(biorep, t)` genuinely
differ in θ — differential turnover across regions/proteoforms, sequence-dependent labeling/quant
biases, residual PTM (Met-Ox / deamidation) and MBR effects. That is a *random-effects* variance the
within-peptide bootstrap never samples. So the dispersion (empirical / propagated = 1.25–1.72×
median) is that between-peptide component made visible — a fixed-vs-random-effects gap — and `measure`
shows it is **heterogeneous** across cells (p90/p10 ≈ 15×), not a uniform rescale.

## Why moderate the variance — the eBayes logic (and where limma differs)

The machinery is limma's, but note *what* we moderate. limma (Smyth 2004) moderates each gene's own
**regression residual variance** `s²_g`, estimated from a few replicates (small df); a spuriously
tiny `s²_g` inflates that gene's t-statistic, so eBayes replaces it with a precision-weighted blend
of the gene's estimate and a pooled prior `s₀²` **and augments the t-test df** to `d_g + d₀`:

```
 s̃²_g = (d₀·s₀² + d_g·s²_g) / (d₀ + d_g)          limma contrast df → d_g + d₀
```

Our situation is one level over, and the object is **different**. The linear model's *residual*
variance `σ̂²` already has decent df (`N − p` ≈ 20 collapsed points), so there is nothing to augment
there — the Δk contrast df stays **`N − p`** (see *The estimators*). What is noisy is the **per-point
weight** `Var̂(θ_i)` (effective df ~6–9, 40 % single-peptide cells), a *separate* object. We borrow
limma's **shrinkage** to de-noise it — a precision-weighted blend of the point's own estimate and the
pooled variance:

```
 Ṽar(θ_i) = (d₀·s₀² + d_i·Var̂(θ_i)) / (d₀ + d_i)      d_i = the variance ESTIMATE's own df
```

so the raw weight stops over-trusting spuriously small variances (the Type-I inflation above). The
`d_i + d₀` here is the **variance estimate's** effective reliability — **not** the t-test df: the
moderated `Ṽar` enters only the weight `w_i = (1 − θ̂_i)²/Ṽar(θ_i)`, and the Δk contrast still uses
`N − p`. (The true-limma move — moderating `σ̂²` across proteins with a Satterthwaite contrast df —
is a **separate roadmap item**; it stabilizes sparse proteins, orthogonal to per-point weighting.)
`Ṽar` is a **continuous dial between the two estimators we already have**:

- `d₀ → ∞` (variances untrustworthy) ⇒ shrink fully to pooled ⇒ **exactly `wls_fit`**;
- `d₀ → 0` (variances trustworthy) ⇒ no shrinkage ⇒ **exactly `wls_fit_var`**.

eBayes picks `d₀` from the data, so the estimator **cannot do worse than the shipped default** and
approaches the full per-point weight only as far as the data earn it (the frontier is mapped in
*The three arms*).

**We already have the mean-variance trend.** limma-*trend* makes the prior a function of the mean
(the mean-variance relationship); our `(1 − θ̂)²` factor *is* exactly that — `Var(φ) = σ²_θ/(1 − θ)²`
is a deterministic mean-variance model, exogenous and reliable. So the clean decomposition keeps the
trend and moderates only the residual: write `Var̂(θ_i) = σ̂²_θ · r_i` (`r_i` = the point's precision
relative to typical), shrink `r_i → 1`, and weight by `(1 − θ̂_i)² / r̃_i`. Full shrink recovers
`wls_fit`, so the trend is never at risk.

**Why not a hard df threshold?** "Use per-point Var only when df > T, else `wls_fit`" is the
all-or-nothing special case of this. It fails on our data specifically because the mass sits **in**
the borderline zone: a threshold high enough to be safe (T ≈ 16) kicks most points back to `wls_fit`
(little gain captured); one low enough to keep them (T ≈ 4) re-admits the risky ones (Type-I climbs).
Shrinkage uses each borderline point *partially* instead of in-or-out — and still degrades gracefully
to `wls_fit`. (We bench the threshold as a baseline anyway, to show the smooth version wins.)

**The dispersion wrinkle.** The measured 1.25–1.72× dispersion says the σ's *under-state* true
scatter. For the WLS point estimate a uniform scale error is harmless (only relative weights matter),
but it means (a) the per-point σ's warrant *more* shrinkage, and (b) a cell has two variance
estimates — the propagated `1/Σ(1/σ²)` and the empirical between-peptide scatter — and the propagated
one under-covers. A reduced-χ² / dispersion scale (itself a shrinkage between the two) is the fix, and
a third arm to bench.

## The three arms — eBayes wins

`moderate` runs all three at the measured regime (`--spread 1.0`; a df mix of **4 @ 40 %**
single-peptide / **12** multi; between-peptide `between = 0.5` → median dispersion ≈ 1.5), against
`wls_fit` and the raw `wls_fit_var`. Mean over k = {0.05, 0.075, 0.11, 0.216}, nsim 400:

| arm | RMSE gain vs `wls_fit` | Type-I | coverage | rmse_vs_wmle |
|---|---|---|---|---|
| `wls_fit` (shipped) | 0 % | 0.037 | 0.946 | 4.13 |
| `wls_fit_var` (raw) | **+70 %** | 0.078 | 0.907 | 1.21 |
| **eBayes `d₀=2`** | **+58 %** | **0.014** | 0.973 | 1.73 |
| eBayes `d₀=4` | +52 % | 0.015 | 0.974 | 1.94 |
| eBayes `d₀=8` | +46 % | 0.018 | 0.971 | 2.18 |
| threshold `T=4 / 8` | +59 % | 0.079 | 0.908 | 1.68 |
| empirical | +51 % | **0.141** | 0.810 | 1.98 |

**eBayes (arm #1) is the clear winner.** *Any* mild shrinkage fixes the Type-I — a finer sweep at
the measured regime gives `d₀=0.5` → **+65 %** gain / Type-I **0.020**, `d₀=1` → +62 % / 0.015,
`d₀=2` → +58 % / 0.014: the raw weight's 0.078 collapses to ~0.02 by `d₀=0.5`, and the whole
`d₀ ≳ 0.5` range is safely conservative (coverage ~0.97). So `d₀` is **not** a value to hand-pick
(the grid just maps the frontier, and `d₀=2` is already conservative — more gain sits at lower `d₀`);
it should be **estimated per run à la Smyth**, from how much the per-point log-variances scatter
beyond their sampling noise. eBayes works because it moderates **every** point, including the
higher-df ones the threshold leaves raw.

**Threshold (arm #2) fails.** `T=4` and `T=8` are identical to each other and barely better than the
raw weight (Type-I 0.079): a hard cut only pools the sub-threshold points; the "passing" (df ≈ 12)
points still carry noisy variances and drive the inflation. The cliff moderates the wrong thing.

**Empirical (arm #3) fails as posed.** Swapping to the between-peptide scatter is *worse* than raw
(Type-I 0.141) — with 2–3 peptides per cell the empirical variance has df ≈ 2, far noisier than the
propagated one. Unbiased but unusable raw; it would itself need shrinkage. So the dispersion, though
real and heterogeneous, is **not** fixed by naively adopting the empirical variance.

## Recommendation

**Publish / recommend on the shipped `wls` (IRLS) now.** It is a large, validated improvement over
OLS (Type-I 0.285 → 0.062, coverage 0.543 → 0.926; `reports/2026-07-13_linear_model_wls.md`,
reproduced byte-for-byte here) and is the correct default regardless of how the per-point work lands.
Nothing below blocks it.

**The per-point path is now settled: eBayes-moderated per-point variance.** The arms above decide the
design — moderation captures ~+50–58 % of the +70 % at nominal-or-conservative Type-I, while the hard
threshold and the raw empirical variance both fail. Plan before it becomes a default:

1. Estimate `d₀` per run **à la Smyth** (from the spread of the per-point log-variances across all
   cells) so the shrinkage is data-driven, and tune toward Type-I ≈ 0.05 (`d₀ = 2` is slightly
   conservative — more gain is available nearer nominal).
2. Wire it as **`--linear-weights wls-var`** (opt-in), carrying the per-point variance through the
   collapse into `fit_linear_deltak` and moderating it there; A/B on `lve_atr`, then decide the default.
3. Optional later refinement: an eBayes-moderated **random-effects** variance (moderate the empirical
   between-peptide scatter rather than discard it) to also absorb the heterogeneous dispersion.

## Reproduce

```bash
# the effect-size sweep (homo/hetero × oracle/noisy)
python -m tests.benchmark.bench_linear_weights varweight --nsim 500 --spread 1.0 --var-df 6

# where production sits + visualize the variance spread across runs
python -m tests.benchmark.bench_linear_weights measure \
    --runs runs/lve_atr_clean runs/boomi_ipsc_d2o --plot out_dir   # → measure_variance.png

# the three moderation arms; sweep d0 to map the shrinkage frontier
python -m tests.benchmark.bench_linear_weights moderate --spread 1.0 --d0-grid 0.5 1 2 4

# play with the scenarios directly (overlaid homo vs hetero, with plots)
python -m tests.benchmark.bench_linear_weights all --regime both --spread 1.0 --plot out_dir
python -m tests.benchmark.bench_linear_weights efficiency --regime hetero --spread 1.0
```

---

## UPDATE 2026-08-09 — `runs/timeseries_atf6`: a new 5-tp DIA, 4-chamber test point for `wls` vs `wls-var`

**Data:** the full ATF6 dataset (`data/timeseries_atf6`; the complete set of which the old
`timeseries_dia` was the LV-control subset). A 4×5×2×3 factorial — **4 heart chambers**
(LA/RA/LV/RV) × **5 D₂O timepoints** (1/3/7/10/14 d) × **2 conditions** (control / ATF6-KO)
× **3 biological replicates**, 120 mouse-heart DIA runs (DIA-NN 2.5.0), RIA 4.6 %. Integrated
+ fit with current defaults (`--coefficients deberneh_2025_rss`, `simple`, `hw`), rolled up
`--model "linear simple" --reference-condition control` with `experiment = chamber`
(`integrate --experiment-column "characteristics[organism part]"`) so the Δk is control-vs-KO
**within each chamber** (33 k peptidoforms, 4 679 proteins over the collapse). Δk BH is
**per-experiment** (each chamber its own hypothesis family — the shipped default). This adds a data
point **between** the sparse 3-tp DIA (`timeseries_dia`, +2 %) and the rich 9–12-tp sets
(boomi/lauren, +10 %) that anchored the "gain scales with curve richness" thesis.

### 1. Where atf6 sits — `measure`

| quantity | atf6 (collapsed cells) | reads as |
|---|---|---|
| heteroscedasticity SD(log σ) | **0.91** (p90/p10 11×) | strongly hetero — on par with lve 1.03 / juber 1.04 (≫ bench 0.7) → the prize exists |
| effective df (timepoints/peptide) | **median 3** (p10 2, p90 5); 32 % single-peptide cells | LOW — DIA curves are gappy, so the *nominal* 5 tp become ~3 usable → the noisy-estimate regime |
| dispersion (empirical/stated σ²) | median **1.31** (p90/p10 11×) | σ under-stated + heterogeneous → the eBayes-moderation regime, not raw per-point var |

So atf6 is the "**big heteroscedastic prize, but low-df noisy variance estimate**" case the shipped
eBayes-moderated `wls-var` is built for — high spread (like the rich sets) but low df (like the
sparse ones).

### 2. Direct A/B on the production rollup (`wls` vs `wls-var`, identical inputs/gates)

- **Point estimates unchanged:** per-condition `k_deg` Pearson **0.990**, Spearman 0.998, median
  |Δ| 0.0005 (0.5 % rel). `wls-var` does not move the rates — it re-weights the **inference**.
- **Power:** significant per-chamber Δk (`delta_k_p_adj < 0.05`, per-experiment BH) rises
  **37 → 65** — a **near-superset** (1 `wls`-only, +29 `wls-var`-only). Per chamber:
  LA 18→31, LV 9→14, RA 2→5, RV 8→15.

A ~1.8× significance jump demands the reproducibility check below before it can be read as power.

### 3. Reproducibility — biorep leave-one-out (the real-replicate test)

Re-ran the rollup 3× per weight, each omitting one biological replicate (the collapse then sees
2 of 3 bioreps; peptide gates unchanged), and measured cross-fold consistency of the estimates:

| metric (cross-fold, n) | `wls` | `wls-var` | change |
|---|---|---|---|
| `k_deg` CV — median (7 219 curves) | 0.0371 | **0.0364** | **−1.8 %** (more reproducible) |
| `delta_k` SD — median (3 508 proteins) | 0.00309 | **0.00298** | −3.6 % |
| `delta_k` sign-consistent across all 3 folds | 57.4 % | **59.9 %** | **+2.5 pp** |

**Reproducibility IMPROVES under `wls-var`** — it does not degrade. Had the +29 extra hits been
Type-I inflation, `wls-var`'s estimates would scatter *more* across biorep folds and its Δk signs
would flip *more* often; instead both tighten. So the extra sensitivity is legitimate power from
the per-point-precision weighting, not over-rejection — consistent with the eBayes moderation
holding Type-I in check on real data.

### 4. Placement on the richness curve

The **+1.8 %** k-reproducibility gain sits just above the sparse-DIA point (+2 % was the *biorep
consistency* metric on `timeseries_dia`) and well below the rich 9–12-tp sets (+10 %). This tracks
**effective df, not nominal timepoints**: atf6 has 5 nominal timepoints but only ~3 usable per
peptide (gappy DIA), so it behaves like the sparse end even though its heteroscedasticity is as
strong as the rich sets. **Takeaway:** `wls-var`'s efficiency gain is governed by curve *df*, and
DIA's gappiness caps it — the strong heteroscedasticity alone is not enough. `wls` remains the safe
default; `wls-var` is a real, reproducibility-positive power gain to opt into on DIA turnover data,
with the caveat that the yield of new significant calls (~1.8×) outruns the point-reproducibility
gain (~2 %), so the extra calls are best treated as *power-recovered, still worth orthogonal
confirmation* rather than a free doubling of discoveries.

### Reproduce

```bash
# integrate → relabel experiment=chamber → fit → rollup (wls + wls-var)
python -m riana integrate data/timeseries_atf6/mzml \
    data/timeseries_atf6/quantms_results/quant_tables/diann_report.parquet \
    --sdrf data/timeseries_atf6/samplesheet_atf6.sdrf.tsv -o runs/timeseries_atf6 -W 6
python scratchpad/relabel_manifest.py runs/timeseries_atf6/riana_manifest.tsv \
    data/timeseries_atf6/samplesheet_atf6.sdrf.tsv --fix-headers
python -m riana fit --manifest runs/timeseries_atf6/riana_manifest.tsv \
    --coefficients deberneh_2025_rss --model simple --label hw -o runs/timeseries_atf6 -W 12
python -m riana rollup --manifest runs/timeseries_atf6/riana_manifest.tsv \
    --model "linear simple" --reference-condition control -W 12            # wls (default)
python -m riana rollup runs/timeseries_atf6 --model "linear simple" \
    --reference-condition control --linear-weights wls-var -o runs/timeseries_atf6/wlsvar -W 12

# regime + biorep-LOO reproducibility A/B
python -m tests.benchmark.bench_linear_weights measure --runs runs/timeseries_atf6
#   (folds: 3× leave-one-biorep-out rollups per weight; see notebooks/wls_var_investigation.ipynb)
```
