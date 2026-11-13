# Δk inference in Riana — technical manual: variance sources, weights and their calibration

2026-11-12 · branch `1.2.0` · simulation code `5d0ebe7` (unchanged `riana/` since) · **Status:
reference.** This is the single place that states how the `linear simple` rollup's Δk test and
the `--interaction` model handle variance, what each weighting option computes, what evidence
calibrates each test, and what is still open. Earlier reports keep their own evidence and are
cited by section; where a statement here differs from one of them, this manual is current.
Proprietary datasets are not used; the real-data checks use ATF6 and lve_atr (§7).

## 0. Summary

- **Two weighting layers, plus a random effect.** (L1) The *collapse* averages a protein's
  peptides into one point per animal × timepoint × condition, weighting each peptide by the
  inverse of its own band variance. (L2) The *regression* fits φ = log(1 − θ) against time and
  weights those points by `--linear-weights`. (L3) The *interaction model* adds a random animal
  slope. "wls" and "wls-var" name L2 only.
- **`wls` (default) is (1 − θ̂)² only** — θ̂ the fitted value from an unweighted first pass. It
  corrects the heteroscedasticity the log transform creates (the cause of OLS's 28 % false
  positives). On the rising-precursor route it is (1 − f̂)², the same formula on f.
  **`wls-var` (opt-in) is (1 − θ̂)²/Ṽ**, adding each point's measurement precision Ṽ from the
  collapse, moderated across points.
- **Neither models peptide heterogeneity or animal variation** (§3). Their *average* size reaches
  the test through the residual variance, which each protein's fit estimates, but two of their
  properties do not: a peptide's offset repeats in every point it enters, and animal variation
  grows with time. The Monte Carlo (§6) shows both matter.
- **Monte Carlo through the shipped CLI (new, §6)** — false-positive rate at α 0.05:
  - *Rising precursor, f route:* `wls` 0.042, calibrated at every turnover rate, ratio recovered;
    `wls-var` 0.255. θ on the same data: Δk biased −59 %, CI coverage 0.38.
  - *Steady precursor, the shipped θ test:* **0.163** with peptide labeling-site error and animal
    variation both on (0.073 with neither). A 2×2 knock-out separates two mechanisms:
    peptide offsets × per-condition collapse weights (**mechanism 1**, fixed by mixing the
    peptides in the same proportions in both conditions: 0.141 → 0.032) and animal variation
    growing with time (**mechanism 2**, fixed by a variance term in φ̂²: 0.119 → 0.063). Both
    fixes together: 0.053, power kept.
  - `wls-var` is no better calibrated than `wls` in any arm.
  - *Interaction:* the ratio-scale test holds (0.036); the mean-Δk and additive tests do not
    (0.145, 0.176) under animal variation, despite the random animal slope.
- **Real data (§7): mechanism 1 dominates ATF6.** Mixing the same peptides in the same
  proportions in both conditions takes ATF6 from 37 Δk calls to 0 and from 86 interaction calls
  to 2, while lve_atr keeps 244 of 250. The calls' collapsed Δk is 4× their peptides' own Δk on
  ATF6 and equal to it on lve_atr. Within-protein peptide k CV on real data is ~0.21 — larger
  than the simulation's, so §6 understates mechanism 1.
- **Decided:** the rising-precursor (f) rollup uses `wls`; `wls-var` is refused there.
- **Open (§8):** how to fix the shipped test — a head-to-head of per-point and per-peptide
  scatter (§8.1: shared proportions, a random-effects collapse, a peptide-level mixed model,
  against today's `wls` and `wls-var`), after calibrating the simulation to real data.

## 1. The chain, stage by stage

| # | Stage | What happens | Variance handled here | Code |
|---|---|---|---|---|
| 1 | Per-point solve (`fit`) | Each peptide × timepoint envelope → θ (steady precursor) or f (rising, `--kp`). Rails drop reads below −0.05 or above 1.05 (on label content when rising). | — | `core/fitting.py` (`_fit_one_concat`), `algorithms/precursor.py` |
| 2 | Peptide curve (`fit`) | One k per peptide × condition by least squares; residual-bootstrap band per point (`fs_lower`, `fs_upper`). Rising: the curve is fitted on label content and the band is mapped to f through the gain ∂A/∂f. | the peptide's own scatter → its band | `core/fitting.py` |
| 3 | Gate (`rollup`) | R² ≥ 0.8, or R² ≥ 0.6 and k_cv < 0.2; Spep ≥ 8 (D₂O); `--peptide-admission auto` = `all` for `linear simple` (a peptide enters only if it passes in both conditions). | — | `core/protein.py` (`rollup_proteins`) |
| 4 | **Collapse — L1** | Per (experiment, condition, protein, biological replicate, timepoint): θ_cell = Σ wᵢθᵢ / Σ wᵢ with wᵢ = 1/σᵢ², σᵢ = band/3.29, **from each peptide's own fit in that condition**. Ṽ = 1/Σ wᵢ, Satterthwaite df. | peptide precision (fixed-effect pooling) | `_collapse_group`, `_weighted_theta`, `_collapse_long` |
| 5 | Transform | φ = log(1 − θ), θ clamped to [0.001, 0.999]; t = 0 dropped; points with φ ≤ −4 dropped (plateau). | plateau saturation | `linear_model.to_phi`, `truncate_plateau` |
| 6 | **Regression — L2** | Per protein and experiment: φ = −k_c·t through the origin, one slope per condition. Unweighted pass → φ̂ → one weighted refit (`--linear-weights`). Δk = k_test − k_ref; SE from the weighted covariance with the residual variance estimated from this protein's points (df n − 2); t-test; BH per experiment. | log-transform heteroscedasticity (`wls`); cell precision (`wls-var`) | `linear_model.fit_linear_deltak` |
| 7 | **Interaction — L3** (`--interaction`) | Same cells and L2 weights, all experiments jointly, plus a random animal slope u·t; REML (statsmodels MixedLM, BFGS). Tests: mean Δk (between-animal containment df), ratio-scale interaction (Satterthwaite df), additive interaction (containment df). | between-animal turnover, when an animal contributes to ≥ 2 experiments | `interaction_model.fit_interaction` |

`--method pooled` skips stage 4 and fits every peptide point; it pseudo-replicates peptides and is
not used for inference (July notebook, §5).

## 2. Terms

| Symbol | Meaning |
|---|---|
| θ | fraction new from the per-point solve under a steady precursor |
| f | fraction of the protein made since labeling began (rising-precursor route); f = θ under an instant precursor |
| A | label content, the pool's mean enrichment ÷ plateau (the measured scale under a rising precursor) |
| P(t), k_p | precursor level 1 − e^{−k_p t} and its rate constant |
| g | gain ∂A/∂f; small early in a slow-precursor series, so Var(f) = Var(A)/g² is large there |
| φ | log(1 − θ) or log(1 − f); linear in t through the origin, slope −k |
| θ̂, φ̂ | fitted values from the unweighted first pass |
| σᵢ | peptide i's band SD at a point, (fs_upper − fs_lower)/3.29 |
| Ṽ | a cell's variance 1/Σσᵢ⁻² (`fs_var` in `riana_rollup_fractions.txt`); Ṽ_mod its eBayes-moderated value (d₀ s₀² + d Ṽ)/(d₀ + d), d₀ and s₀² from a robust fitFDist over all cells |
| σ² | the residual variance of one protein's weighted fit (estimated, per protein) |
| σ_animal | SD of the random animal slope in the interaction model |
| cell | one collapsed point: an animal (biological replicate) × timepoint × condition × experiment |

## 3. Variance sources and where each is handled

| | Source | Structure | Size on φ | L1 collapse | L2 `wls` | L2 `wls-var` | L3 interaction |
|---|---|---|---|---|---|---|---|
| S1 | Envelope measurement noise | independent per peptide × point; on θ roughly constant (σ_θ ≈ 0.056, lve_atr) | σ_θ/(1 − θ): grows as θ → 1 | weights peptides by their bands | **yes — this is what (1 − θ̂)² corrects** | yes | as L2 |
| S2 | Precision differences between peptides and cells | per peptide; cells with more or better peptides are more precise | varies cell to cell | yes | no (assumes equal precision; costs efficiency, not calibration) | yes, through Ṽ | as L2 |
| S3 | Rising-precursor gain (f route only) | deterministic in t, known from k_p | Var(f) = Var(A)/g²: ~10²–10³× at the first points | yes (band ÷ g) | no — early points get full weight but sit at small t, so they move a through-origin slope little and inflate σ² (conservative) | yes, but miscalibrated (§6) | as L2 |
| S4 | Peptide-specific systematic offsets: labeling-site count error (table vs truth), interference, mis-assignment | **persistent**: the same peptide reads high or low at every timepoint, in both conditions | ≈ a per-peptide slope factor; grows with t | **not modeled; interacts with L1** (below) | no | no — Ṽ is within-peptide scatter | no |
| S5 | Between-animal turnover variation | per animal (× protein); one animal per cell, so it is cell-level error shared by the cell's peptides | ≈ t·δk: variance ∝ (k·t)² ≈ φ̂² | no | no — σ² absorbs its average, but the weights assume no growth with t | no | **yes**, when an animal contributes to ≥ 2 experiments |
| S6 | Per-animal precursor enrichment error | per animal, shared by all its proteins | ≈ a per-animal slope factor | no | no | no | partly (absorbed by the animal slope) |
| S7 | Plateau saturation | late points of fast proteins | φ → −∞ | — | φ ≤ −4 dropped | same | same |

**How S4 interacts with the collapse (mechanism 1).** A peptide's weight in a cell comes from its
band in *that condition's* fit (stage 2), which is roughly constant over the condition's cells.
The band is estimated from a handful of residuals, so the same peptide gets one weight in
condition A and a different one in B. With peptide offsets, the cell value of A and the cell value
of B are then different mixtures of the same biased peptides, and the difference behaves like a
slope difference — a Δk under the null. It is a bias of the estimate, not extra noise, so no
variance model can absorb it; it is removed by mixing the peptides in the same proportions in
both conditions, or by modeling the peptide.

**How S5 breaks the weights (mechanism 2).** In a single-experiment fit with one animal per cell,
the animal deviation is part of the residual. Its average is in σ², but `wls` assumes the
residual variance is proportional to 1/(1 − θ̂)² — flat in t for slow proteins — whereas animal
variation grows as t². The late points, which carry the slope, are noisier than the pooled σ²
says, so the SE is too small. The 2026-11-10 report found the same thing from the other side: the
joint (interaction) model's SEs were 8–25 % larger than the per-chamber fit's, "mainly because
the joint model lets the animal variance grow with day² on the late, high-leverage points";
within one chamber the two variances cannot be separated.

**ATF6 evidence for S5 and S6** (2026-11-10 report §3.1, §3.5). The animal effect is significant
in 97 % of proteins, animal-to-animal CV of k ≈ 13 %. The per-animal shifts are largely
proteome-wide (split-half r = 1.00) and follow the timepoint: day-1 animals 6–38 % slower in both
genotypes (enrichment still rising), the control day-14 group 14–23 % faster (a group effect).

## 4. The weight options

| `--linear-weights` | L2 weight wₜ | Models | Status | Evidence |
|---|---|---|---|---|
| `ols` | 1 | nothing | audit only (keeps t = 0, as the pre-July fit did) | false positives 0.285, coverage 0.543 (07-13 MC) |
| **`wls`** | **(1 − θ̂)² = e^{2φ̂}**, θ̂ from the unweighted pass, one IRLS step | S1 | **default since 2026-07-13** | 0.062, coverage 0.926 (07-13 MC, S1 only); ATF6 `varExp` δ̂ = −1.10 vs the −1 this weight implies (11-10 §3.2) |
| `wls-var` | (1 − θ̂)²/Ṽ_mod | S1 + S2 | opt-in since 2026-07-30 | raw Ṽ 0.078 → moderated 0.014–0.020, +58–65 % efficiency (07-24 MC); +~10 % biorep-split consistency on rich series (07-24 notebook) |

Rules that apply to all three: the weight uses the **fitted** θ̂, never the observed θ (observed-θ
weights correlate with each point's error and bias k by −8 %; 07-13 MC); one IRLS step (more
chase the noise); the residual variance is estimated per protein (not fixed at 1), so the weights
set only the *relative* precision of a protein's points.

**On the rising-precursor route** the same formulas apply to f: `wls` = (1 − f̂)²,
`wls-var` = (1 − f̂)²/Ṽ_mod(f) with Ṽ from bands mapped through the gain. **Riana uses `wls`
there**; `wls-var` is miscalibrated on f (§6.3) and the rollup will refuse it.

**The interaction model** uses the same L2 weight (since 2026-11-12, `d59116e`) as nlme
`varFixed(~1/w)`: `wls` → (1 − θ̂)², `wls-var` → (1 − θ̂)²/Ṽ_mod, `ols` → `wls` (logged).

## 5. How the estimator got here

| Date | Change | Why | Evidence | Commit · report |
|---|---|---|---|---|
| 2026-06-13 | `linear simple`: φ-space OLS, joint fit per protein, Δk contrast + BH | a cross-condition Δk with a shared residual variance | lve_atr, 410 calls | `e558d0a` |
| 2026-07-13 | **`wls` default**; t = 0 dropped from the fit | OLS on φ ignores the log transform's heteroscedasticity | false positives 0.285 → 0.062 (S1-only MC); fast-tail k bias −14 % → ~0 | `4ae5e50` · 07-13 |
| 2026-07-30 | `wls-var` opt-in, eBayes-moderated Ṽ (robust fitFDist d₀) | recover the efficiency lost by assuming equal cell precision (S2) | raw weights 0.078 → moderated ~0.02; +70 % → +58–65 % efficiency | `c82cc55` · 07-24 |
| 2026-08-05 | per-condition gates; `--peptide-admission` (`auto` → `all` for `linear simple`) | the condition-blind gate let a peptide's pass in one condition admit it in the other; `all` makes the two conditions' peptide sets identical (paired) | — | `b26453a` |
| 2026-08-10 | per-experiment BH (multi-experiment rollups) | one BH family per tissue | — | `180e477` |
| 2026-11-10/11 | `--interaction`: random animal slope, ratio-scale primary, Satterthwaite df, `characteristics[individual]` | an animal contributes several points in a multi-tissue fit (split-plot) | ATF6: genotype 186 → 1, interaction 3 → 132 hits | `b8715f7`, `57185c0` · 11-10 |
| 2026-11-12 | the interaction follows `--linear-weights` | it had hard-coded `wls` | — | `d59116e` |
| 2026-11-12 | rising-precursor fit (`--kp`): fs holds f | an unprimed D₂O precursor rises over days | 11-11 report | `3a526c3` · 11-11 |
| 2026-11-12 | **f route uses `wls`, not `wls-var`**; envelope-level MC through the CLI | this manual | §6 | `9a2d357`, `5d0ebe7` |

**What the earlier simulations covered.** Each calibrated the test against the sources it had:

| Simulation | Level | Sources | Result for the shipped `wls` |
|---|---|---|---|
| 07-13 `bench_linear_weights.py` | cell-level θ, lve_atr 12-day design | S1 | 0.062 |
| 07-24 `bench_linear_weights.py varweight/moderate` | cell-level, per-point σ spread, noisy σ̂ | S1, S2 | 0.034–0.044 |
| 07-28 `notebooks/wls_var_sim.py` | peptide × biorep × t, collapse; peptide k-offset τ = 0.3 shared by the conditions | S1, S2, S4 | 0.047 (150 sims) |
| **11-12 `bench_rising_null.py`** (§6) | envelopes → shipped `riana fit` / `riana rollup` | S1–S6 | §6 |

The 07-28 simulation had peptide offsets but drew each peptide's weight independently in every
cell, so condition-specific mixtures averaged out over cells. The shipped collapse takes the
weight from the peptide's per-condition curve band, the same in all of that condition's cells —
the systematic difference mechanism 1 needs. None of the earlier simulations had animals.

## 6. Monte Carlo through the shipped pipeline (2026-11-12)

### 6.1 Design

`tests/benchmark/bench_rising_null.py` simulates isotope envelopes and runs the shipped
`riana fit` and `riana rollup` on them, so every stage of §1 is exercised as shipped.

| | |
|---|---|
| Design | 2 conditions × 3 tissues sampled from the same animals; one animal per condition × day; days 0, 1, 3, 7, 14, 21, 35, 56, 84 |
| Proteins | 250 per dataset, 2–6 random tryptic peptides each (9–20 residues, ≥ 9 labeling sites); true k log-uniform 0.004–0.15/day; tissue factors 1, 0.5, 2 |
| Labeling | D₂O, RIA 0.04, `deberneh_2025_rss` sites; precursor k_p 0.08/day (*slow*) or 100/day (*instant*, ≈ a priming dose) |
| Envelope noise | per channel relative 3.5 % + absolute 0.0025, × a per-peptide log-normal multiplier (SD 0.5): peptides differ in precision (S1, S2) |
| Peptide offsets (S4) | true labeling-site count = table × e^{N(0, 0.10)}, per peptide |
| Animal (S5, S6) | turnover k × e^{N(0, 0.15)} per animal × protein, enrichment RIA × (1 + N(0, 0.05)) per animal; both shared by the animal's tissues |
| Interference | 5 % of points get +20–100 % on one of m3–m5 |
| Truth | k_B = ratio × k_A in every tissue; ratio 1 (null), 1.25, 1.5, 2 (the user's R simulation grid) |
| Fits | θ (`--kp inf`) and f (k_p from the manifest); rollup `linear simple`, `wls` and `wls-var`, `--interaction` |
| Datasets | 3 per null arm, 1 per other arm; fixed seeds; settings, commit and tree state in each arm's `ratio<r>_params.json` (all `5d0ebe7`, clean) |

`eval_rising_variance.py` then re-tests the same gated peptides under the collapse and weight
variants with the shipped functions (§9). **Audit** (`evaluation_audit.tsv`): its cells match
the CLI's `riana_rollup_fractions.txt` to 4 × 10⁻⁵ and its Δk p-values the CLI's to 5 × 10⁻⁷ (the
files' rounding) over 13,633 tests; the reimplemented `wls-var` matches the shipped one to
10⁻¹⁴. Rejection rates are pooled over datasets; with ~710 tests per dataset the binomial SE at
0.05 is ~0.008, but the datasets differ by more (per-animal effects are shared by every
protein), so per-dataset ranges are given.

### 6.2 Rising precursor (k_p 0.08/day) — the f route

Per-tissue Δk, shipped (per-condition) collapse:

| Route · weights | False positives, ratio 1 (3 datasets) | … by true k: < .01 / .01–.03 / .03–.08 / > .08 | CI coverage, ratio 1 / 1.5 | Power at 1.25 / 1.5 / 2 | log₂ ratio at 1.5 / 2 (truth 0.58 / 1) | Δk bias at 1.5 |
|---|---|---|---|---|---|---|
| **f · `wls`** | **0.042** (0.031–0.050) | 0.04 / 0.04 / 0.06 / 0.04 | 0.96 / 0.95 | 0.18 / 0.50 / 0.80 | 0.60 / 1.02 | +4 % |
| f · `wls-var` | 0.255 (0.203–0.306) | 0.27 / 0.28 / 0.28 / 0.15 | 0.75 / 0.74 | 0.37 / 0.83 / 0.84 | 0.61 / 0.95 | +5 % |
| θ · `wls` | 0.050 (0.037–0.064) | **0.16** / 0.04 / 0.01 / 0.01 | 0.95 / **0.38** | 0.15 / 0.52 / 0.92 | **0.43 / 0.76** | **−59 %** |
| θ · `wls-var` | 0.045 (0.033–0.057) | 0.15 / 0.03 / 0.00 / 0.01 | 0.96 / 0.39 | 0.15 / 0.48 / 0.91 | 0.43 / 0.76 | −59 % |

- **f · `wls` is calibrated** at every turnover rate, its CIs cover, and it recovers the ratio.
  The early points' huge Var(f) is not in the `wls` weight, but those points sit at small t and
  move a through-origin slope little; they inflate σ², which makes the test slightly
  conservative.
- **f · `wls-var` is not.** Ṽ(f) is the peptide band ÷ the gain, spanning ~10³ over a series; it
  is also misspecified in time (it omits the animal variance shared by a cell's peptides). With
  weights that extreme the residual scale no longer fits the few points that carry the slope.
  Its extra power is the inflation.
- **θ under a rising precursor** reads below f, so its k is biased low (median k̂/k 0.68) and
  its Δk by −59 %; its CIs do not cover the true Δk; slow proteins are over-called.

Interaction on the same data (ratio 1; mean Δk / additive / ratio scale): f · `wls`
0.052 / 0.014 / 0.001; f · `wls-var` 0.218 / 0.281 / 0.086; θ · `wls` 0.007 / 0.053 / 0.010.
Mean-Δk power for f · `wls` at ratios 1.25 / 1.5 / 2: 0.25 / 0.71 / 0.88.

### 6.3 Steady precursor — the shipped θ test

Per-tissue Δk, θ route, instant precursor, all nuisance on:

| Collapse · weights | False positives, ratio 1 (3 datasets) | … by true k: < .01 / .01–.03 / .03–.08 / > .08 | CI coverage, ratio 1 / 1.5 | SD of z | Power at 1.5 | log₂ ratio at 1.5 |
|---|---|---|---|---|---|---|
| **per-condition · `wls` (shipped)** | **0.163** (0.121–0.188) | 0.21 / 0.18 / 0.15 / 0.07 | 0.84 / 0.82 | 1.57 | 0.88 | 0.60 |
| per-condition · `wls-var` | 0.168 (0.128–0.197) | 0.20 / 0.18 / 0.17 / 0.10 | 0.83 / 0.80 | 1.57 | 0.90 | 0.60 |
| shared · `wls` | 0.133 (0.109–0.154) | 0.19 / 0.14 / 0.12 / 0.05 | 0.87 / 0.84 | 1.43 | 0.89 | 0.60 |
| per-condition · pool2 | 0.090 (0.059–0.108) | 0.08 / 0.10 / 0.11 / 0.05 | 0.91 / 0.87 | 1.25 | 0.89 | 0.60 |
| **shared · pool2** | **0.053** (0.036–0.070) | 0.05 / 0.06 / 0.06 / 0.03 | 0.95 / 0.92 | 1.04 | 0.88 | 0.60 |

*pool2* (research only, §9): Var(φ) = c·Ṽ/(1 − θ̂)² + s·φ̂², with c and s estimated per experiment
from the leverage-corrected squared residuals of all proteins. The s·φ̂² term is variance growing
as (k·t)², the signature of animal turnover variation.

### 6.4 The two mechanisms — 2×2 knock-out

θ route, instant precursor, null; one factor family switched off at a time (interference and
per-peptide precision stay on). False-positive rate of the per-tissue Δk test:

| Peptide site error | Animal variation | per-condition · `wls` (shipped) | shared · `wls` | per-condition · pool2 | shared · pool2 |
|---|---|---|---|---|---|
| off | off | 0.073 | 0.048 | 0.062 | 0.041 |
| **on** | off | 0.141 | **0.032** | 0.134 | 0.026 |
| off | **on** | 0.119 | 0.115 | **0.063** | 0.041 |
| on | on (§6.3) | 0.163 | 0.133 | 0.090 | 0.053 |

- **Mechanism 1 — peptide offsets × per-condition collapse weights** (§3, S4). With site error
  alone the shipped test runs at 0.141; mixing the peptides in the same proportions in both
  conditions brings it to 0.032. A variance model cannot fix it (pool2: 0.134) — it is a bias of
  the estimate.
- **Mechanism 2 — animal variation growing with time** (§3, S5). With animal variation alone the
  shipped test runs at 0.119; the shared collapse does nothing (0.115); the φ̂² variance term
  brings it to 0.063.
- **Envelope noise alone: 0.073** — the level the 07-13 simulation found for `wls` (0.062–0.070);
  the shared collapse takes it to 0.048, so even per-peptide precision differences mix the
  conditions differently.

### 6.5 Interaction model

Rejection rates at ratio 1 (θ route, instant precursor; `wls`; the shipped per-condition collapse
→ shared):

| Nuisance on | Mean Δk | Additive interaction | Ratio-scale interaction |
|---|---|---|---|
| neither | 0.041 → 0.041 | 0.081 → 0.057 | 0.077 → 0.069 |
| site error only | 0.135 → **0.016** | 0.163 → 0.041 | 0.146 → 0.045 |
| animal variation only | 0.111 → 0.094 | 0.111 → 0.107 | 0.004 → 0.000 |
| both | **0.145** → 0.147 | 0.176 → 0.145 | **0.036** → 0.011 |

(CLI `wls-var`, both on: 0.186 / 0.197 / 0.030.) Mean-Δk power at ratio 1.5: 0.94.

- Mixing (mechanism 1) inflates all three interaction tests and the shared collapse removes it.
- Animal variation inflates the mean-Δk and additive tests **despite the random animal slope**.
  A likely reason (to test in §8.1): the simulated animal effect is proportional — k × e^u, the
  same *relative* shift in every tissue — while the model's random slope is additive, the same
  *absolute* shift in φ-slope in every tissue, and the tissues' k differ four-fold. ATF6's animal
  shifts were also proportional in size (6–38 % at day 1; 11-10 §3.5). The ratio-scale test,
  which works on log ratios, holds.
- Simulated σ_animal (median 0–0.006) is below ATF6’s (0.013).

### 6.6 What the simulation does not cover

- **Peptide spread is understated:** site error 10 % against a real within-protein peptide k CV
  of ~0.21–0.23 (§7); mechanism 1 is larger on real data.
- One animal per condition × day (ATF6 has three), so each peptide's per-condition band rests on
  ~9 points here vs ~14 on ATF6.
- The animal turnover effect is independent across proteins; on ATF6 it is largely proteome-wide
  (11-10 §3.5), which correlates the tests across proteins.
- No missing data (DDA gaps would make identical proportions impossible, §8.1); random sequences;
  no DIA/DDA acquisition differences; one dataset in the non-null arms.

## 7. Real data: what the per-condition collapse does to the calls

No truth here, so two comparisons that do not need one: (a) re-test the same gated peptides with
the peptides mixed in the same proportions in both conditions (the shared collapse, C in §8.1),
and (b) compare each protein's Δk with its peptides' *own* Δk — each peptide's k in one condition
minus its k in the other, from the fit's per-peptide curves, which involves no collapse at all.
Two datasets, chosen to bracket the question:

- **ATF6** (`runs/atf6_individual/fit`; DIA; ATF6-KO vs control; 4 chambers; 5 timepoints × 3
  animals per genotype). The genotype effect is small at the animal level (−1.5 % ± 3.7 %,
  11-10 §3.6), so few true per-protein Δk are expected.
- **lve_atr** (`runs/lve_atr_clean`; DDA; atrium vs left ventricle; 12 timepoints). Large real
  differences (atrium ~1.24× faster, 2026-06-28 report).

Same gates and admission as the shipped rollup; `wls` unless stated; calls = BH < 0.05 per
experiment.

| Dataset | Collapse | Weights | Δk tests | Calls | Nominal p < 0.01 | r with peptides' own Δk | Calls: median \|Δk\| (collapsed / peptides' own) |
|---|---|---|---|---|---|---|---|
| ATF6 | per-condition (shipped) | `wls` | 3,510 | **37** | 143 | 0.66 | **0.025 / 0.006** |
| ATF6 | per-condition | `wls-var` | 3,510 | 65 | 198 | 0.61 | 0.027 / 0.007 |
| ATF6 | shared | `wls` | 3,510 | **0** | 20 | 0.58 | — |
| ATF6 | shared | `wls-var` | 3,510 | 2 | 37 | 0.53 | — |
| lve_atr | per-condition (shipped) | `wls` | 389 | **250** | 225 | 0.53 | 0.021 / 0.020 |
| lve_atr | per-condition | `wls-var` | 389 | 258 | 237 | 0.49 | 0.021 / 0.020 |
| lve_atr | shared | `wls` | 389 | **244** (233 shared with the shipped set) | 219 | 0.63 | 0.021 / 0.020 |

ATF6 interaction (ratio scale, `wls`): **86 calls per-condition → 2 shared** (Cops3 O88543,
Q9JHU4); mean Δk 1 → 0; additive 137 → 3. Of the two interaction proteins the ATF6 biology
notebook highlights, Cops3 survives (BH 0.002) and Smyd1 does not (0.97). σ_animal is unchanged
(median 0.013).

**Reading.**

- Where real differences exist (lve_atr) the shared collapse keeps them: 244 of 250 calls, and the
  called proteins' Δk equals their peptides' own Δk. Where the animal-level effect is small
  (ATF6) the shipped calls disappear, and their collapsed Δk is four times their peptides' own
  Δk. ATF6's shipped Δk and interaction calls are mostly peptide mixing, not turnover.
- With the per-condition collapse ATF6 has 143 tests at nominal p < 0.01 where ~35 are expected
  with no effects at all; with the shared collapse, 20.
- `wls-var` adds calls with the per-condition collapse on both datasets, and almost none with
  the shared collapse on ATF6 — its extra ATF6 calls come through the same mixing.
- Agreement with the peptides' own Δk: the shared collapse is closer in absolute terms on both
  datasets (median gap 0.002 vs 0.003 on ATF6, 0.004 vs 0.005 on lve_atr) and better
  correlated on lve_atr (0.63 vs 0.53), slightly less on ATF6 (0.58 vs 0.66). The peptides'
  median is itself a noisy equal-weight estimate, not truth; with little true signal on ATF6,
  r mostly reflects shared noise.
- **The peptide spread is large.** Within-protein peptide k CV (gated peptides, reference
  condition): median 0.21 on ATF6 (90th percentile 0.43), 0.23 on lve_atr (0.47). The
  simulation's ±10 % site error is well below this, so §6 understates mechanism 1.

**How it happens — ATF6, LA, P97447.** Each peptide's own Δk is ≈ 0 (median −0.001, all within
±0.012), but the peptides' apparent k span 0.027 to 0.164 per day. The knockout's weights happen
to favour the two slow peptides and the control's the fast one:

| Peptide | weight share KO / control | k KO / control | own Δk |
|---|---|---|---|
| AIVAGDQNVEYK_2 | 0.09 / 0.25 | 0.164 / 0.164 | +0.000 |
| FCANTCVDCR_2 | 0.38 / 0.18 | 0.027 / 0.028 | −0.001 |
| FVFHNEQVYCPDCAK_3 | 0.14 / 0.13 | 0.130 / 0.125 | +0.004 |
| GEDFYCVTCHETK_2 | 0.07 / 0.09 | 0.128 / 0.130 | −0.001 |
| GSSVVAYEGQSWHDYCFHCK_3 | 0.05 / 0.11 | 0.152 / 0.146 | +0.006 |
| QVIGTGSFFPK_2 | 0.08 / 0.15 | 0.099 / 0.110 | −0.012 |
| YWHDNCFR_2 | 0.20 / 0.10 | 0.039 / 0.047 | −0.008 |

Collapsed Δk = −0.043 (p 3 × 10⁻⁵); shared collapse −0.007 (p 0.44).

**Why the weights differ between conditions.** A peptide's weight is 1/σ², σ from its band in
each condition's own curve fit, i.e. from the scatter of its ~14 points (ATF6) around that
curve. Two independent variance estimates on ~12 df each differ by chance: their ratio has a
5–95 % range of about 0.37–2.7 (F₁₂,₁₂). Measured ratio of each peptide's band variance
KO / control on ATF6: 5–95 % 0.29–2.54 — about what chance alone produces (the median, 0.81, is
a uniform shift: KO curves fit slightly tighter overall, which rescales every peptide alike and
does not change the mix). Nothing about the KO makes particular peptides worse. On lve_atr the
spread is wider (atrium / ventricle 5–95 % 0.22–7.92, median 1.26): chance plus real tissue
differences in peptide quality. The peptide *set* is already identical in both
conditions (`--peptide-admission all`); the proportions are not.

**Why ATF6 and lve_atr differ.** The mixing shift (per-condition Δk − shared Δk) is the same size
in both:

| | ATF6 | lve_atr |
|---|---|---|
| Δk SE, median | 0.0066 | 0.0043 |
| \|mixing shift\|, median (90th pct) | 0.0020 (0.0112) | 0.0016 (0.0092) |
| \|shift\| / SE, median | 0.35 | 0.41 |
| proteins shifted by > 2 SE | 2.7 % | 4.1 % |
| \|Δk\| / SE with the shared collapse, median | **0.80** | **3.64** |

What differs is the signal. On lve_atr the real differences dwarf the shift, so most calls
survive (the ~4 % of proteins shifted by more than 2 SE move a few calls each way: 17 lost, 11
gained). On ATF6 there is little real Δk, so the proteins the shift pushes past 2 SE are most of
what gets called. With a median of two admitted peptides per protein in both datasets, any weight
difference moves the cell value.

**Time-growing variance on real data.** The pool2 fit (§6) puts 27–33 % of the modeled variance
in the φ̂² term in every ATF6 chamber, and 25 % on lve_atr, under either collapse — consistent
with S5, though residual peptide effects also grow with t.

## 8. Open decisions

### 8.1 Plan — head-to-head of per-point and per-peptide scatter (requested 2026-11-12)

The two kinds of scatter a protein's data carry are different things, and the current options
model only one of them:

- **Per-point (within-peptide) scatter** — how far a peptide's points sit from *its own* curve:
  measurement noise (S1, S2). It is what σᵢ and Ṽ measure, and what `wls-var` weights by. It
  is independent from point to point.
- **Per-peptide (between-peptide) scatter** — how far a protein's peptides sit from *each other*:
  apparent-k differences from labeling-site error, interference, mis-assignment (S4). In
  proteomics this is the familiar peptide-level spread. It is **persistent**: a peptide that
  reads slow does so at every timepoint and in both conditions. Riana models it nowhere, and the
  per-condition collapse turns it into a Δk (mechanism 1; on ATF6, §7).

`wls-var` is not the per-peptide component: its Ṽ = 1/Σσᵢ⁻² pools the peptides as if they
measured the same value (fixed effect), so two peptides that each fit tightly but disagree by
six-fold give a small Ṽ; and it treats errors as independent between points, whereas a peptide's
offset repeats in every point it enters. A per-peptide component keeps peptide identity, so the
offset cancels from a within-peptide Δk contrast instead of entering through unequal mixing.

**Where between-peptide variance τ² can enter — three places, not exclusive:**

1. *Collapse weights* (how a cell averages its peptides): wᵢ = 1/(σᵢ² + τ²), the random-effects
   (DerSimonian–Laird) weight, τ² estimated per protein from how much its peptides disagree beyond
   their own noise. When τ² dominates σᵢ², the weights become nearly equal, so chance differences
   in σᵢ between conditions stop changing the mix. This addresses the bias.
2. *Cell variance → regression weights* (what `wls-var` would use): Ṽ_RE = 1/Σ 1/(σᵢ² + τ²),
   which counts peptide disagreement in each cell's precision — e.g. a cell where only 2 of 8
   peptides were seen is less precise than a full one. This sets how the regression trusts each
   cell; it cannot undo a bias already in the cell values (if the mixes differ, reweighting the
   cells does not unmix them).
3. *A term in the model* (E): the peptide offset b_pep is part of the mean, shared by the
   conditions, so it drops out of the Δk; the per-point variance σᵢ² then belongs to each
   peptide point, and τ² is the variance of b_pep if it is a random effect.

So the instinct to put peptide-level variance into the linear weights is right for (2), but the
larger problem the spread causes is in the means — the mixing — which only (1), a shared mix (C)
or the model term (3) removes. **Missing data decides between them:** with complete data
(ATF6 DIA) identical proportions are achievable in every cell; with gappy DDA data a peptide
missing from one cell changes that cell's mix whatever the weights, and only (3) handles it
without discarding data.

**Candidates** (all on the same gated peptides; A–C exist, D–E would be built):

| | Collapse | Variance / weights | Peptide spread is… | Status |
|---|---|---|---|---|
| A | per-condition inverse variance | `wls` | not modeled (mixing bias) | shipped default |
| B | per-condition inverse variance | `wls-var` (per-point Ṽ) | not modeled | opt-in, unused in projects; **user not convinced** |
| C | shared proportions (each peptide's weight level pooled over conditions) | `wls` | cancels through identical mixing | evaluator only (`shared_bands`) |
| D | random-effects: wᵢ = 1/(σᵢ² + τ²), shared levels as C, τ² per protein × experiment (DerSimonian–Laird) | `wls`, Ṽ_RE available for weights | a variance in the cell (more equal mixing; Ṽ includes the spread) | to build |
| E | none — every peptide point | φ = −(k_c + b_pep)·t + u_animal·t + ε, peptide slope b_pep shared by the conditions (fixed or random), animal random slope, (1 − θ̂)² weights | a persistent peptide effect; Δk is a within-peptide contrast replicated by animals | to build; generalizes `--interaction` |

**Evidence, in order** (each for the per-experiment Δk and the interaction):

1. *Calibrate the simulation first* (the deferred item; it now gates the comparison). Set the
   simulated peptide spread to the real within-protein peptide k spread (ATF6 §7, lve_atr), the
   animal variation to ATF6's σ_animal, and add proteome-wide animal shifts (11-10 §3.5) beside
   the per-protein ones. Then null rate (overall and by k), power at ratios 1.25 / 1.5 / 2,
   coverage and ratio recovery for A–E.
2. *Real-data nulls by permutation.* Permute the condition labels between animals within each
   timepoint (ATF6: 3 + 3 per day per chamber), re-run `riana fit` on the permuted manifest (the
   per-condition bands must be re-estimated), roll up under A–E. The rejection rate over
   permutations is each candidate's false-positive rate on real data — no truth needed.
3. *Agreement with the peptides' own Δk* (§7) and cross-tissue consistency.
4. *Precision* by the 07-24 biorep split — for efficiency only: it compares bioreps within one
   condition with the peptide weights fixed, so it cannot see a between-condition mixing bias.

**Decision rule.** Among the candidates that hold the null in both the calibrated simulation and
the permutations, take the most powerful; the simplest if tied. Until then `wls` stays the
default and `wls-var` stays opt-in, documented as calibrated only against S1–S2.

**Direction (user, 2026-11-13).** Try **D as the replacement for `wls-var`** (tested on real data
as above), and **E later as the thorough solution**. D and E map onto the existing
`--method weighted | pooled` split: D upgrades the collapse route, and E is the peptide-point
route done properly, so the pooled side becomes a real estimator, not a pseudoreplicated
reference, and can later become the main route. The new plan is to be written after the
non-primed D₂O dataset analysis. Points to carry into it:

- *D keeps C's shared levels* (as specified in the table). Shared proportions alone removed
  mechanism 1 in the simulation (§6.4); τ² on top makes the cell variance honest (Ṽ_RE cannot
  fall below ≈ τ²/n, which is the failure of `wls-var`'s Ṽ) and bounds the weight range, the
  problem `wls-var` has on f (§6.2), so D may also be usable on f. Test that, don't assume it.
- *τ² is not a constant.* A peptide offset is an offset in its slope (labeling-site error scales
  its apparent k), so in φ the peptides of a cell spread as ≈ c·k·t: τ grows with time, as
  `pool2`'s s·φ̂² term does. Model it as a relative spread, τ²(cell) ≈ τ_rel²·φ̂² (or the θ-scale
  equivalent), with τ_rel pooled over a protein's cells or over proteins, not one τ² per protein.
- *D cannot fix composition changes:* a peptide missing from one cell still changes that cell's
  mix (DDA gaps). Only E handles that without dropping data.
- *E with fixed peptide effects* (sum-to-zero slope offsets per protein, shared by the
  conditions) is enough for Δk, which is then a within-peptide contrast. It keeps one random
  effect (the animal slope), so the existing `--interaction` machinery (MixedLM, BFGS,
  Satterthwaite) fits it; random peptide effects are needed only to borrow strength. Longitudinal
  designs also need a sample (run) effect beside the animal slope. Check the additive-offset
  approximation (the offset is multiplicative in k) on the residuals.
- *E unifies three things:* per-experiment Δk, `--interaction`, and within-animal comparisons
  (tissue A vs B in the same animals, §8.2). One model family instead of three.
- *Naming.* E should be a new `--method` value rather than a redefinition of `pooled`: `pooled`
  also serves the exponential refit, and redefining it would silently change results.
  `pooled` can be retired once E is validated. D replaces the per-condition weights inside
  `weighted` (results-affecting when adopted), and its Ṽ_RE takes `wls-var`'s place.

**Datasets:** ATF6 (DIA, 3 animals per cell, 4 tissues), lve_atr (DDA, 2 tissues), boomi /
lauren (cell culture, rich series), the proprietary slow-precursor set (locally only), and
**EDL / SOL** (the sister comparison of lve_atr in the Study 4 manuscript) — **to do: not yet
run through current Riana**, so its peptide-level fits do not exist; run it first. The Study 4 R
workflow itself averages peptides with log₂(intensity) weights, which are near-equal and the
same in both conditions, so it does not carry the per-condition mixing of §7 to the same degree
(whether its peptide sets match between conditions has not been checked).

### 8.2 Other open items

- **Shipped results.** Any `linear simple` Δk or interaction result rests on the per-condition
  collapse. How much mechanism 1 moved it depends on the within-protein peptide spread relative
  to the Δk; §7 shows it can dominate. Results to re-examine once 8.1 decides: the ATF6 biology
  notebook (`notebooks/atf6_biology.ipynb`) and the 11-10 report's interaction hits.
- **The rising-precursor rollup** (PROJECT_REVIEW Track C, item 7): `wls` on f; refuse `wls-var`.
  Built 2026-11-13 (`9daf9f8`). It inherits whatever 8.1 decides for the collapse.
- **Within-animal comparisons** (tissue A vs B of the same animals, one experiment, condition =
  tissue). The Δk test does not read `individual`; each animal's two tissues enter as
  independent points. The animal's shared deviation (its body-water enrichment, its overall
  turnover) nearly cancels from the tissue contrast but stays in the residual variance, so the
  test is conservative with respect to animals and loses power. The fix is the paired analysis:
  φ ~ 0 + day:tissue + (0 + day | animal), tissue within animal, as a paired t-test is a
  random-intercept model; it reduces to today's fit when σ²_animal = 0. Mechanism 1 is not
  removed by pairing and may be systematic here (a peptide measured better in one tissue), so
  the test as a whole is not guaranteed conservative. Covered by E in 8.1.
- **The interaction model's animal term.** Its random slope is additive in φ (the same absolute
  slope shift in every experiment); simulated proportional animal effects inflate its mean-Δk
  and additive tests (§6.5). Check on ATF6 whether animal shifts scale with chamber k; if so, a
  proportional random slope (u·k̂_s·t) or a log-scale model. Part of E in 8.1.
- **3b ceiling drop** — dropping reads at or above the precursor ceiling selects fast proteins'
  early reads one-sidedly (k biased low there, seen while building the simulation; not
  re-measured here — fast proteins' lower power on f in §6.2, 0.24 at ratio 1.5 for k > 0.08 vs
  ~0.55 below, is consistent); cap at the ceiling instead.
- **Δk residual-variance moderation** (PROJECT_REVIEW §3, deferred) — moderates each protein's σ²
  across proteins; do it after 8.1, on the winning model.

## 9. Reproduce

```bash
# §6 — the whole grid (~2 h on 24 cores; fixed seeds), then the evaluator
bash tests/benchmark/run_deltak_mc.sh                  # OUT=runs/deltak_mc, W=22
python -m tests.benchmark.eval_rising_variance runs/deltak_mc     # re-run the evaluation alone

# §7 — real data (fits under the gitignored runs/; outputs saved beside them as
#       eval_collapse_variants.txt)
python -m tests.benchmark.eval_rising_variance runs/atf6_individual/fit --real control ATF6-KO
python -m tests.benchmark.eval_rising_variance runs/lve_atr_clean --real control atrium
```

- `tests/benchmark/bench_rising_null.py` — the simulation and the shipped-CLI runs; writes
  `<arm>/ratio<r>_params.json` (every setting, the commit, whether `riana/` or the bench differed
  from it) and `<arm>/ratio<r>_results.tsv` (rejection rates from the CLI's own output files).
- `tests/benchmark/eval_rising_variance.py` — re-tests the same gated peptides under the
  collapse × weight variants: the shipped `_collapse_long`, `fit_linear_deltak` and
  `fit_interaction`, plus the reimplemented pool2 model. Writes `evaluation.tsv` (one row per
  dataset × route × collapse × test) and `evaluation_audit.tsv` (§6.2 checks).
- The earlier simulations: `tests/benchmark/bench_linear_weights.py` (07-13, 07-24),
  `notebooks/wls_var_sim.py` + `notebooks/wls_var_investigation.ipynb` (07-28).
