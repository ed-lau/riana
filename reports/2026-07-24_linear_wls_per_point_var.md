# `linear simple` Δk — per-point Var(θ) WLS: a large but conditional gain; moderate the variance first

**Date:** 2026-07-24 · **Data:** synthetic Monte-Carlo + real `runs/{lve_atr_clean,
boomi_ipsc_d2o, juber_ac16_d2o}` · **Status: MEASURED — not shipped.** The effect is real
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
- **Verdict:** big, worthwhile efficiency prize, but the naive weight is **not** a safe default.
  Next step (pre-registered here): **moderate the variance** — empirical-Bayes shrinkage toward
  the pooled variance + the measured dispersion scale — and re-run `varweight` to confirm Type-I
  returns to ~0.05 while keeping the bulk of the gain. Do **not** wire the raw estimator into the
  `wls` default.

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

## Why moderate the variance — the limma/eBayes logic

This is the problem limma solves for gene expression, one level up. There, each gene's variance
`s²_g` is estimated from a few replicates (small df), so some genes draw a spuriously tiny `s²_g` →
an inflated t-statistic → a false positive. eBayes (Smyth 2004) replaces each estimate with a
**moderated** variance — a precision-weighted blend of the gene's own estimate and a pooled prior
`s₀²`, with the prior df `d₀` estimated empirically from how much the `s²_g` scatter:

```
 s̃²_g = (d₀·s₀² + d_g·s²_g) / (d₀ + d_g)          effective df = d_g + d₀
```

Our cells are the "genes": `Var̂(θ_i)` is noisy (effective df ~6–9, 40 % single-peptide), so the
raw weight over-trusts spuriously small variances — the Type-I inflation above. The moderated weight
is `(1 − θ̂_i)² / Ṽar(θ_i)` with the same blend, and it is a **continuous dial between the two
estimators we already have**:

- `d₀ → ∞` (variances untrustworthy) ⇒ shrink fully to pooled ⇒ **exactly `wls_fit`**;
- `d₀ → 0` (variances trustworthy) ⇒ no shrinkage ⇒ **exactly `wls_fit_var`**.

eBayes picks `d₀` from the data, so the estimator **cannot do worse than the shipped default** and
approaches the full per-point weight only as far as the data earn it. Target here: lift effective df
from ~6–9 into the ~16+ band the sweep showed is Type-I-safe — a `d₀` ≈ 8–10.

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

## Recommendation

**Publish / recommend on the shipped `wls` (IRLS) now.** It is a large, validated improvement over
OLS (Type-I 0.285 → 0.062, coverage 0.543 → 0.926; `reports/2026-07-13_linear_model_wls.md`,
reproduced byte-for-byte here) and is the correct default regardless of how the per-point work lands.
Nothing below blocks it.

**Continue evaluating `wls_fit_var` as a refinement — do not default it yet.** The prize is real and
large (a ~4× efficiency gap under production's own heteroscedasticity), but the raw per-point weight
inflates Type-I at the measured df, so it needs the moderation above. Plan, all in the bench before
any production code:

1. Add three `varweight` arms — **eBayes-moderated** `Ṽar` (shrink `r_i → 1`, `d₀` estimated à la
   Smyth), the **hard df-threshold** baseline, and the **dispersion-scaled / empirical** variance.
2. Confirm the winner lands **Type-I ≈ 0.05, coverage ≈ 0.95** at the measured `--spread 1.0
   --var-df 6` while keeping most of the +70 %.
3. Only then wire it as **`--linear-weights wls-var`** (opt-in first, A/B on `lve_atr`, then decide
   the default). The plumbing is small — the collapse already computes the per-point variance; it
   need only be carried into `fit_linear_deltak`.

## Reproduce

```bash
# the effect-size sweep (homo/hetero × oracle/noisy)
python -m tests.benchmark.bench_linear_weights varweight --nsim 500 --spread 1.0 --var-df 6

# where production sits (per-point variance regime, real data)
python -m tests.benchmark.bench_linear_weights measure --run runs/lve_atr_clean

# play with the scenarios directly (overlaid homo vs hetero, with plots)
python -m tests.benchmark.bench_linear_weights all --regime both --spread 1.0 --plot out_dir
python -m tests.benchmark.bench_linear_weights efficiency --regime hetero --spread 1.0
```
