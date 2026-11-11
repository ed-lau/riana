# Random animal effects (nlme-style mixed model) — feasibility on ATF6, and what it means for the statistics

**Date:** 2026-11-10 · **Data:** `runs/timeseries_atf6` (ATF6-KO vs control, 4 heart chambers,
5 timepoints, 30 animals; canonical `wls` rollup) · **Status: SHIPPED 2026-11-10** as
`rollup --interaction` (`riana.core.interaction_model`): ratio-scale interaction as the primary test
on the conservative between-animal df, δ fixed at −1, animal = condition × day × replicate. On ATF6
it finds 39 chamber-specific KO effects across 977 proteins (see `CHANGELOG.md`).
Prompted by the statistical review of the titin manuscript, which asked for an `nlme::lme`
analysis with a `varExp` variance function, a strain effect, and a random mouse slope.

## TL;DR

- **Where an animal effect exists.** Each ATF6 animal is one genotype at one timepoint, sampled
  in all 4 chambers. In the shipped **per-chamber** fit every animal is exactly one point, so the
  residual is already animal-level and a random animal effect cannot be estimated. In a
  **cross-chamber** fit each animal contributes 4 points (a split-plot): genotype varies between
  animals, chamber within them. This is the reviewer's model, with strain ↔ genotype and
  isoform ↔ chamber.
- **The animal effect is real and large.** Significant in **97 %** of the 580 proteins observed
  in all 8 chamber × genotype cells. Animal-to-animal CV of k ≈ **13 %**.
- **The reviewer's variance function agrees with the shipped weights.** `varExp(fitted)` gives
  δ = **−1.10** (median), vs the −1 implied by Riana's fixed `(1 − θ̂)²` weight. It beats
  constant variance in 97 % of proteins. A per-chamber variance (`varIdent`) adds little
  (median ΔAIC −0.9).
- **Ignoring the animal in a cross-chamber fit is badly wrong in both directions.** The genotype
  main effect is pseudo-replicated (**186 → 1** BH hits once the animal is included; SE ≈ 2×,
  df 112 → 29). The chamber × genotype interaction is under-powered (**3 → 132** hits), because
  the animal shift cancels in within-animal contrasts.
- **About a third of those interaction hits are a scale effect.** Tested on the ratio scale
  (log k_KO/k_ctrl), **82** remain: the rest are a uniform proportional KO effect that differs in
  absolute Δk only because chambers turn over at different rates.
- **The animal shifts are proteome-wide and follow the timepoint.** Day-1 animals of both
  genotypes show 6–38 % slower apparent k (body-water enrichment still rising), and the three
  control day-14 animals are 14–23 % faster. A random animal effect absorbs this; it is also
  evidence that the constant-enrichment model misfits early timepoints.
- **The notebook's "global turnover brake under ATF6-KO" does not hold at the animal level.**
  Per-animal proteome-wide apparent k, KO vs control, adjusted for day: **−1.5 % (SE 3.7 %,
  p = 0.69, n = 30 animals)**. The per-chamber medians are not independent evidence (the same 30
  animals appear in all 4 chambers).
- **Feasible in Python without R.** statsmodels `MixedLM` on weight-transformed data reproduces
  `nlme` (estimates to 1.5e-3, F to 1.3 %, identical denominator df) at **39 ms per protein**,
  provided the optimizer is BFGS (the default L-BFGS silently returns wrong fits here).

## 1. The design, and where an animal effect can exist

ATF6 has 30 animals: 2 genotypes × 5 labeling days (1, 3, 7, 10, 14) × 3 animals. Each animal is
sacrificed at **one** day and sampled in **all 4 chambers** (LA, LV, RA, RV). The animal key is
`(condition, labeling_time, biological_replicate)`; this was checked against the SDRF
`characteristics[individual]` column and the run-file prefix (e.g. `843_LA_QEHF_DIA` → animal
843).

| Fit | Points per animal | Animal effect |
|---|---|---|
| Shipped per-chamber Δk (one chamber, 2 genotypes) | exactly 1 (0 duplicates) | not identifiable; the residual is already between-animal |
| Cross-chamber model (the planned two-factor interaction) | up to 4 (one per chamber) | identifiable; required |

The same holds for most single-tissue, cross-sectional datasets: one sample per animal per fit
means no separate animal term. An animal effect is meaningful only when one animal contributes
more than one point to a fit: several tissues (ATF6), several isoforms (the titin pairs), or
repeated sampling over time.

In the cross-chamber model the layout is a **split-plot**: animal is the whole plot (genotype
and day vary between animals), chamber is the sub-plot (it varies within an animal).

## 2. Models fitted

580 of the 1,284 rollup proteins have all 8 chamber × genotype cells with ≥ 3 points (after the
shipped t > 0 and plateau filters). In that set the median protein has 120 points from all 30
animals.

Fixed effects, with explicit numeric columns so that `nlme` assigns the correct denominator df:

    φ ~ 0 + dLA + dLV + dRA + dRV + dKO + dKOc1 + dKOc2 + dKOc3

- `dLA … dRV` = day × chamber indicator: the control slope in each chamber (varies within an
  animal).
- `dKO` = day × KO: the genotype slope shift averaged over chambers (constant within an
  animal → between-animal term).
- `dKOc1…3` = `dKO` × sum-contrast code of chamber: how the KO shift differs by chamber, i.e. the
  interaction (varies within an animal).

k = −slope. Random effect: `~ 0 + day | animal` (a per-animal slope, as in the reviewer's code).
Because each animal has one day, this is a per-animal offset whose variance grows with day².

| Model | Variance | Animal effect |
|---|---|---|
| g0 | shipped weights, `varFixed(1/(1−θ̂)²)` | none (what a naive cross-chamber WLS does) |
| m1 | shipped weights | random slope |
| m2 | reviewer: `varComb(varIdent(~1 \| chamber), varExp(~fitted(.)))` | random slope |
| m3 | `varExp(~fitted(.))` only | random slope |
| m2c | constant | random slope |
| g2 | reviewer's variance function | none |

All 580 proteins fit with every model (REML; `nlme` 3.1.169 in R). 5 models × 580 proteins took
70 s on 6 cores.

## 3. Results

### 3.1 The animal effect is real and large

| | Shipped weights (m1 vs g0) | Reviewer variance (m2 vs g2) |
|---|---|---|
| Animal variance > 0 (REML LRT, boundary-corrected, BH < 0.05) | 563 / 580 (97 %) | 572 / 580 (99 %) |
| Estimate at the zero boundary | 1 | 0 |

Animal-to-animal variation in k, σ_animal / mean k: median **0.134** [IQR 0.112, 0.164].

### 3.2 The variance function: δ ≈ −1, per-chamber variance not needed

| Comparison | Result |
|---|---|
| `varExp` δ (m3) | median **−1.10** [IQR −1.58, −0.85] |
| `varExp` δ with `varIdent` (m2) | median −1.08 [IQR −1.46, −0.84]; 72 % within [−1.5, −0.5]; 100 % negative |
| `varExp` vs constant variance (m3 vs m2c) | ΔAIC median **−35.1**; better in 97 % |
| `varIdent` + `varExp` vs constant (m2 vs m2c) | ΔAIC median −37.1 [IQR −67, −20.6]; better in 98 % |
| adding `varIdent(chamber)` to `varExp` (m2 vs m3) | ΔAIC median **−0.9** [IQR −7.2, +2.7]; better in 43 % |

`varExp(~fitted(.))` models Var(φ) = σ²·exp(2δ·φ̂). With δ = −1 that is σ²/(1 − θ̂)², the
delta-method variance Riana already weights by (`reports/2026-07-13_linear_model_wls.md`). The
data estimate δ close to −1, so the shipped fixed weight is well supported. Per-chamber variance
ratios are all near 1 (medians 0.94–1.03), and `varIdent` improves AIC in fewer than half the
proteins.

### 3.3 Genotype main effect across chambers: pooling without the animal is anti-conservative

The genotype main effect (`dKO`, the mean KO shift over chambers) varies only between animals.

| | g0 (no animal) | m1 | m2 |
|---|---|---|---|
| Denominator df | 112 | 29 | 29 |
| BH < 0.05 | **186** | **1** | 0 |
| raw p < 0.05 | 262 | 32 | 30 |

SE ratio m1/g0: median 1.94 [IQR 1.66, 2.24]. Without the animal term, the 4 chambers of one
animal are counted as 4 independent replicates.

### 3.4 Chamber × genotype interaction: a large power gain, partly a scale effect

The interaction terms (`dKOc1…3`, 3 df) vary within animals, so the shared animal shift cancels.

| | g0 (no animal) | m1 | m2 |
|---|---|---|---|
| Denominator df | 112 | 84 | 84 |
| BH < 0.05 | **3** | **132** | 136 |
| raw p < 0.05 | 41 | 224 | 212 |

The 3 g0 hits are among the 132. The two interaction proteins from the earlier prototype remain
strong: Cops3 (O88543) p = 4.4e-8 and Smyd1 (P97443) p = 2.4e-10 in m1.

**Additive vs ratio scale.** Δk is an absolute difference. A KO that slows every chamber by the
same 10 % produces a larger absolute Δk in faster chambers (atria), which reads as an
interaction. Testing equality of log(k_KO/k_ctrl) across chambers instead (delta method, same
model) leaves **82** BH hits, 81 of them among the 132. So about 40 % of the additive-scale hits
are a uniform proportional effect. Cops3 (p = 5.3e-6) and Smyd1 (p = 2.8e-9) remain
chamber-specific on both scales.

### 3.5 What the animal effect is made of

The per-animal shifts are shared across proteins and grouped by timepoint:

- Split-half reliability of each animal's proteome-wide shift (median over two random halves of
  the proteins, 50 splits): **r = 1.00**.
- Share of the per-animal (BLUP) variance that lies between the 10 genotype × day groups:
  median **0.69** [IQR 0.59, 0.81], against about 0.31 expected if animals varied independently.

Per-animal apparent k relative to each protein's genotype-blind slope (median over ~3,500
protein × chamber points per animal; positive = faster):

| Day | ATF6-KO animals | control animals |
|---|---|---|
| 1 | −0.25, −0.16, −0.06 | −0.38, −0.17, −0.15 |
| 3 | −0.14, −0.13, −0.10 | −0.22, −0.16, −0.13 |
| 7 | −0.05, −0.02, +0.08 | −0.03, −0.02, +0.15 |
| 10 | −0.17, +0.03, +0.05 | −0.10, −0.10, −0.09 |
| 14 | −0.19, +0.01, +0.03 | **+0.14, +0.18, +0.23** |

Two patterns stand out:

- **Early timepoints are slow in both genotypes.** That is what a body-water enrichment that is
  still rising after labeling starts would produce under the constant-enrichment model: early θ
  falls below the line. It is a model misspecification, not animal biology (see the open
  precursor-enrichment item in `PROJECT_REVIEW.md`).
- **The control day-14 group is offset as a group** (all three animals 14–23 % faster), while
  the KO day-14 group is not. That looks like a batch or group effect. Day 14 has the most
  leverage on a through-origin slope.

A random animal slope treats all of this as animal-level noise. That is the right error
structure for inference, but it does not remove a systematic difference between genotypes in
labeling (e.g. different body-water enrichment). That cannot be separated from a real global
turnover change without measuring the enrichment per animal.

### 3.6 The "global turnover brake" at the animal level

`notebooks/atf6_biology.ipynb` reports a global turnover brake under ATF6-KO (median Δk negative
in all 4 chambers). A proteome-wide claim is replicated by animals, not proteins. Testing each
animal's proteome-wide apparent k (KO vs control, adjusted for day; 30 animals, 24 df):

| | KO − control (relative apparent k) | SE | p |
|---|---|---|---|
| All chambers | **−0.015** | 0.037 | **0.69** |
| LA | −0.025 | 0.045 | 0.59 |
| LV | −0.010 | 0.034 | 0.78 |
| RA | −0.011 | 0.039 | 0.79 |
| RV | −0.017 | 0.035 | 0.63 |

The direction agrees with the notebook, but the effect is small and not distinguishable from
zero at the animal level. The negative median in all 4 chambers is not 4 independent results:
the same 30 animals, with the shared shifts in §3.5, appear in every chamber. **The notebook's
claim should be qualified.**

### 3.7 The shipped per-chamber Δk

The joint model's per-chamber KO effects agree with the shipped per-chamber Δk (580 proteins):

| | LA | LV | RA | RV |
|---|---|---|---|---|
| r(Δk), joint model vs shipped | 0.966 | 0.940 | 0.970 | 0.949 |
| SE ratio, joint / shipped (median) | 1.08 | 1.25 | 1.13 | 1.22 |

The shipped per-chamber analysis is not pseudo-replicated (one point per animal). Its SEs are
8–25 % smaller than the joint model's, mainly because the joint model lets the animal variance
grow with day² on the late, high-leverage points. Which variance structure is right for the
per-chamber fit is **open**; it cannot be decided within a single chamber, where animal and
residual variance are not separable.

## 4. What it means for the statistics

1. **Single-tissue, cross-sectional data (most datasets): no change.** One sample per animal per
   fit means the residual already is the animal-level error.
2. **Any fit where one animal contributes several points needs the animal term.** Here that is
   the planned two-factor (chamber × genotype) model; in the titin study it is the isoform
   pairs. Without it, between-animal terms (genotype, strain) are pseudo-replicated and
   within-animal terms (chamber, isoform and their interactions) are under-powered.
3. **The correct df follow the split-plot.** Between-animal terms are tested on
   n_animals − p_between df (29 here); within-animal terms on the within-animal df (84 here).
   `nlme`'s containment rule assigns exactly this when the design columns are coded as in §2.
4. **Proteome-wide claims are animal-level claims.** Their replication is 30 animals, and the
   ATF6 global brake does not survive that test (§3.6).
5. **"Chamber-specific" needs a stated scale.** Additive Δk and ratio k_KO/k_ctrl give 132 vs
   82 interaction hits. The ratio scale matches the biological question "does KO affect this
   chamber differently" when chambers differ in baseline turnover.
6. **The fixed `(1 − θ̂)²` weight is supported** (δ̂ ≈ −1.1), and per-stratum variance is not
   needed as a default (resolves design choice 2 in the two-factor plan).

## 5. Implementation scope (Python)

**Engine.** statsmodels `MixedLM` on weight-transformed data: multiply φ, the fixed-effect
columns and the random-effect column (day) by √w, with w the shipped `(1 − θ̂)²` weight from a
first OLS pass. This is exactly `nlme::lme(..., weights = varFixed(~1/w))`. Checked against
`nlme` m1 on all 580 proteins (`runs/atf6_lme/py_lme.py`):

| Quantity | Max relative difference | Median |
|---|---|---|
| Genotype effect estimate | 1.5e-3 | 5e-7 |
| Its standard error | 3.7 % | 6e-4 |
| Interaction F | 1.3 % | 7e-5 |
| Residual σ | 3e-5 | 4e-7 |
| Animal SD σ_animal | identical except 4 proteins at the zero boundary | 2e-6 |
| Denominator df (both terms) | identical for every protein | |

Same 132 interaction hits. 39 ms per protein single-core (about 1 min for 1,300 proteins; the
rollup's `-W` workers parallelize it). No new dependency.

**Gotcha.** statsmodels' default L-BFGS optimizer stops at a wrong point on this problem
(variances around 1e-5) and still reports convergence: 65 of 580 fits failed outright and most
of the rest were wrong. BFGS (with Nelder-Mead as fallback) agrees with `nlme` on all 580.

**Denominator df.** statsmodels gives none, so the containment rule is reimplemented from
`nlme:::getFixDF`: between-animal columns get n_animals − p_between; within-animal columns get
N − n_animals − p_within, plus 1 when the model has no intercept.

**Not available in statsmodels, and whether it matters:**
- Estimating δ (`varExp`) or per-stratum variance (`varIdent`): the data support fixed δ = −1
  and no `varIdent` (§3.2). If δ must be estimated, a small custom REML that profiles over δ is
  feasible.
- Satterthwaite or Kenward-Roger df: not needed to match `nlme`, which uses containment.
- `nlme` iterates the `varExp` weights to convergence; the shipped fit uses one step, which the
  WLS report found slightly better for Type-I error.

**Data plumbing.** The model needs an animal key. For ATF6,
`(condition, labeling_time, biological_replicate)` works because the SDRF assigns each animal
one biological-replicate index, shared across its chambers. `riana_rollup_fractions.txt`
already carries all three columns. The general fix is to carry `characteristics[individual]`
from intake (already on the deferred list). Longitudinal designs need it, since one animal then
spans several days.

**Where it plugs in.** It becomes the variance model of the two-factor interaction layer
(`rollup --interaction`, PROJECT_REVIEW Track C), built on the per-experiment substrate as
planned. Per protein it would report:
- the per-stratum Δk from the joint model;
- the genotype main effect at between-animal df;
- the interaction F, p and BH-adjusted p, on the additive and/or ratio scale;
- σ_animal.

**Effort.** About 200–300 lines in core (design columns, fit, df rule, Wald tests), CLI wiring
and output columns. Tests can use stored `nlme` reference values for a few proteins. A GUI
Interaction tab would come later.

**Decisions needed:**
1. Interaction scale: additive Δk, ratio, or both.
2. Keep δ fixed at −1 (recommended) or estimate it.
3. Animal key: the composite key with a validation check now, or thread `individual` from
   intake first.
4. Revise the ATF6 notebook's global-brake statement (§3.6).
5. Whether, when the interaction model runs, the per-stratum Δk it reports should replace the
   per-chamber fits' (§3.7).

## 6. Notes on the reviewer-requested R code (titin)

The code from the review discussion maps onto this analysis. Each point below was checked by
fitting the reviewer's formulas to ATF6 (chamber in place of isoform, genotype in place of
strain; `nlme` 3.1.169):

1. **Unique animal IDs.** The mouse factor must be unique per animal, e.g.
   `interaction(mouse, strain, day)` when IDs are reused across groups. In ATF6 the equivalent
   is `(condition, day, biological_replicate)`.
2. **The strain model `ms` is fine as written.** R codes
   `phi ~ 0 + day:isoform + day:strain` with one strain-shift column, and `nlme` tests it on
   mouse-level df (29 on ATF6; the isoform slopes get the within-mouse df, 87). That is the
   correct split-plot allocation.
3. **`anova(mi, Terms = "day:isoform:strain")` does not test whether the strain effect differs
   by isoform.** With `phi ~ 0 + day:isoform + day:isoform:strain`, R codes one strain-shift
   column per isoform (4 on ATF6), so the F-test is that **all shifts are zero**: "any strain
   effect in any isoform". Each of those columns varies within a mouse, so `nlme` tests it on
   the **within-mouse df** (83 on ATF6), although the average shift is a between-mouse
   quantity. That makes the test anti-conservative for the main-effect part. The interaction
   itself is a test that the shifts are **equal**: add sum-contrast interaction columns
   alongside `day:strain` (as in §2), or compare `ms` and `mi` by an ML likelihood-ratio test.
4. **Denominator df for the response.** Coded as in §2, strain is tested on mouse-level df and
   isoform × strain on within-mouse df. A response saying both are tested "at mouse-level df"
   would be inaccurate for the interaction (and the within-mouse df is the correct one for it).
5. **δ and AIC.** Expect δ ≈ −1 (here −1.10), which is the delta-method variance of
   φ = log(1 − θ); report it with the AIC against constant variance (here −35 median). Per-group
   `varIdent` may add little; report its AIC separately rather than only the combined model.

## Reproduce

```bash
# design audit + cross-chamber table (needs runs/timeseries_atf6 and the ATF6 SDRF)
python runs/atf6_lme/prep.py
# nlme reference fits (R, nlme) and the varExp-only model
Rscript runs/atf6_lme/nlme_fit.R
Rscript runs/atf6_lme/varexp_only.R
python runs/atf6_lme/summarize.py          # §3.1–3.5
# Python MixedLM prototype: parity vs nlme, ratio-scale interaction, per-chamber Δk (§3.4, §3.7, §5)
python runs/atf6_lme/py_lme.py
# animal-level global KO effect (§3.5 table, §3.6)
python runs/atf6_lme/animal_level.py
```
