# Δk inference in Riana — technical manual (as built 2026-11-16, updated 2026-11-19)

2026-11-16 · branch `1.2.0` · code `be2d03d` · **Status: reference.** Updated 2026-11-19 for the
rollup changes of 2026-11-17/19 (one peptide per protein, admission by level, the contrast pair,
a shared peptide for Δk, dimethyl forms; §0, §7, §8.5, §8.6, §9, §10).
This manual states how the `linear` rollup estimates per-condition k and tests Δk and the
experiment × condition interaction, at both levels (`--level protein`, the default, and
`--level peptide`). It covers what each variance source does, which test holds which null, and
what is still open. It **supersedes** [the 2026-11-12 manual](2026-11-12_deltak_inference_manual.md),
which describes the pipeline before the redesign (per-condition collapse, `linear simple`,
`wls-var`) and remains the record of the Monte Carlo and real-data evidence that motivated the
redesign. The evidence for the redesign itself is in
[the 2026-11-14 evaluation](2026-11-14_wls_var_D_E_evaluation.md); this manual cites both by
section and does not repeat their derivations. Where it differs from either, this manual is
current. Proprietary datasets are not used.

## 0. Summary

- **Two levels.**
  - **`--level protein` (default):** each protein's peptides are averaged into one θ per
    biological replicate × timepoint × condition (the *collapse*), then φ = log(1 − θ) is fitted
    against time by weighted least squares. Since 2026-11-16 the collapse gives every peptide
    **one precision level across the conditions** (candidate C), so both sides of a Δk average
    the same peptides in the same proportions.
  - **`--level peptide` (opt-in, `--model linear` only):** one mixed model per protein on the
    peptides' own points (E+time): fixed peptide slope offsets shared by the conditions, random
    slopes for labeling time and animal, an animal × condition slope when the same animals carry
    both conditions. `--interaction` under it is EJ+time.
- **Both hold every Δk null** they were tested on (§8): false-positive rates at α 0.05 of
  0.028 / 0.039 (C) and 0.025 / 0.042 (E+time) on the boomi replicate split / ATF6 animal
  permutation, against 0.059 / 0.096 for the old collapse. E+time has 1.11× C's power at equal
  size. The old collapse failed the time-split check on five of six Boulder tissue pairs; C
  passes all six and E+time five (ATR–EDL reads 1.26).
- **Interaction:** at the protein level (C's points plus the random animal slope) the
  ratio-scale test sits on the bound on ATF6 (0.065) and the additive one above it (0.081);
  EJ+time holds both (0.027 / 0.041).
- **Real data:** ATF6's per-chamber Δk calls go from 35 to 0 and its ratio-scale interaction
  calls from 79 to 2 (protein level) or 0 (peptide level); both old counts sat at their
  permutation-null level. Where real differences are large (lve_atr, the Boulder tissue pairs)
  the calls stay within 6 % of the old ones; EDL–SOL, the smallest effect, has 152 calls under
  the old collapse, 136 under C and 162 under E+time.
- **Weights:** `wls`, (1 − θ̂)² from the fitted value, at both levels (and (1 − f̂)² on the
  rising-precursor route). `wls-var` was removed on 2026-11-15 (§10). `ols` is kept for audit.
- **Model set** (2026-11-16): `exponential` (was `simple`), `linear` (was `linear simple`),
  `calibration`; guan / fornasiero removed in favour of the rising-precursor solve at fit
  (`fit --kp`).
- **Admission and coverage (2026-11-17/19, §7, §8.6):**
  - one peptide that passes the gates reports a protein (`--min-peptides 1`);
  - admission `auto` is `all` for the protein-level Δk and `own` for the peptide level and the
    exponential model (protein-level `own` fails the nulls: 18 % on the ATF6 permutation);
  - a linear rollup uses only its reference and test conditions, and a Δk needs a peptide
    measured in both.
- **Open (§11):** the ATR–EDL time-split reading; longitudinal designs and the f route at the
  peptide level; whether `--level peptide` becomes the default; a peptide × condition component
  (T6 LV–RV); residual-variance moderation; flagging unmeasurable k; single-peptide tests at
  0.06 on the ATF6 permutation.

## 1. The chain, stage by stage

| # | Stage | What happens | Code |
|---|---|---|---|
| 1 | Per-point solve (`fit`) | Each peptide × timepoint envelope → θ (precursor at plateau) or fraction new f (rising precursor, `fit --kp` or the SDRF `precursor rate constant`). Reads below −0.05 or above 1.05 are dropped (on label content when rising). | `core/fitting.py` (`_fit_one_concat`), `algorithms/precursor.py` |
| 2 | Peptide curve (`fit`) | One k per peptide × condition by least squares on the exponential model; residual-bootstrap prediction band per point (`fs_lower`, `fs_upper`). | `core/fitting.py` |
| 3 | Gate (`rollup`) | Per condition: R² ≥ 0.8, or R² ≥ 0.6 and k_cv < 0.2. Spep ≥ 8 (D₂O) / 6 (¹⁸O) is applied at fit. Admission `all`: a peptide is used only if it passes in every condition it appears in. `unique` parsimony. | `core/protein.py` (`gate_peptides`) |
| 4p | **Collapse (protein level)** | Bands rescaled to the peptide's level pooled over conditions (§4.1), then per (experiment, condition, protein, biological replicate, timepoint): θ_cell = Σ wᵢθᵢ / Σ wᵢ, wᵢ = 1/σᵢ². | `_shared_bands`, `_collapse_long`, `_collapse_group` |
| 5 | Transform | φ = log(1 − θ); t = 0 dropped; points with φ ≤ `--phi-limit` (−4 ≈ θ 0.98) dropped. Protein level clamps θ to [0.001, 0.999]; peptide level caps at 0.999 and floors at −0.5. | `linear_model.to_phi`, `truncate_plateau`, `peptide_model.peptide_phi` |
| 6p | **Regression (protein level)** | φ = −k_c·t through the origin, one slope per condition; OLS pass → θ̂ → one WLS refit. Δk t-test on the protein's own residual variance; BH per experiment. | `linear_model.fit_linear_deltak` |
| 6e | **Peptide-level model** (`--level peptide`) | E+time on the gated peptide points (§5); replaces 4p and 6p. | `core/peptide_model.py` (`fit_peptide_deltak`) |
| 7 | **Interaction** (`--interaction`, ≥ 2 experiments) | Protein level: the collapsed points of all experiments jointly, with a random animal slope (§6.1). Peptide level: EJ+time (§6.2). | `interaction_model.fit_interaction`, `peptide_model.fit_peptide_interaction` |

The `exponential` rollup uses the same gate and the collapse without the band pooling, and
refits the curve per (experiment, condition, protein); it gives per-condition k, no Δk test.

## 2. Terms

| Symbol | Meaning |
|---|---|
| θ | fraction new from the per-point solve under a precursor at plateau |
| f | fraction of the protein made since labeling began (rising-precursor route); f = θ under an instant precursor |
| φ | log(1 − θ) or log(1 − f); linear in t through the origin, slope −k |
| θ̂, φ̂ | fitted values from the unweighted (OLS) first pass |
| σᵢ | a peptide point's band SD, (fs_upper − fs_lower)/3.29 |
| m_c, m | a peptide's mean σ² over its points in condition c, and over all its conditions |
| cell | one collapsed point: biological replicate × timepoint × condition × experiment |
| animal | `characteristics[individual]` when the SDRF has it; otherwise the condition × timepoint × biological replicate composite |
| b_p | peptide slope offset (fixed, sum-to-zero, shared by the conditions) |
| u_time, v_animal, w_cell, x_tissue×time | random slopes per labeling time, animal, animal × condition (or × experiment), experiment × time |
| σ_time, σ_animal, σ_cell | their SDs, written to the output |
| paired | the same animals carry both conditions (e.g. tissues of one animal as conditions) |

## 3. Variance sources and where each is handled

| | Source | Structure | Protein level (C + `wls`) | Peptide level (E+time) | Interaction: protein level / EJ+time |
|---|---|---|---|---|---|
| S1 | Envelope noise | per peptide × point, ~constant on θ, grows as 1/(1 − θ) on φ | (1 − θ̂)² weights | (1 − θ̂)² weights | (1 − θ̂)² weights |
| S2 | Precision differences between peptides | per peptide | collapse weights 1/σ² (pooled level) | not weighted (each point carries only (1 − θ̂)²; E+band tested, not chosen, eval §3) | as the level |
| S3 | Rising-precursor gain (f route) | Var(f) = Var(A)/g², large at early points | `wls` on f; calibrated (11-12 MC: 0.042) | **not evaluated** | as the level |
| S4 | Peptide offsets (labeling-site error, interference) | persistent per peptide, both conditions | same peptides in the same proportions in both conditions → the offset cancels from Δk (mechanism 1 removed) | fixed b_p; Δk is a within-peptide contrast | as the level |
| S5 | Between-animal turnover | per animal, grows with t | in the residual (the weights assume no growth with t: mechanism 2) | v_animal (and w_cell when paired) | random animal slope / v_animal + w_cell |
| S5t | Timepoint-shared deviation | every sample taken at one time departs together; ~85 % shared across ATF6 chambers | in the residual (conservative) | u_time | — / u_time + x_tissue×time |
| S6 | Per-animal enrichment error | per animal, all proteins | residual | absorbed by v_animal | absorbed by the animal slope |
| S7 | Plateau saturation | late points of fast proteins | φ ≤ −4 dropped | same | same |
| S8 | Whole-group proteome-wide shift | a set of animals turns over faster or slower as a whole | not modeled (estimand decision, §11) | partly in v_animal | partly |

The mechanisms are described in the 11-12 manual §3 and measured in its §6–§7: **mechanism 1**
(S4 × per-condition collapse weights) is a bias of the estimate that no variance model absorbs;
**mechanism 2** (S5 growing with t) makes the protein-level SE too small. C removes mechanism 1;
the peptide-level model carries mechanism 2 and S5t in its random slopes.

## 4. Protein level (default)

### 4.1 The collapse (C, `protein._shared_bands`)

Within one (experiment, protein), for each peptide (`concat`) present in at least two
conditions:

1. Each point's band variance σᵢ² = ((fs_upper − fs_lower)/3.29)².
2. m_c = the peptide's mean σ² over its points in condition c; m = its mean over all its points.
3. Every band in condition c is rescaled about fs by √(m/m_c). The time profile within a
   condition is kept; the level becomes the peptide's pooled one.

Then the inverse-variance collapse per (experiment, condition, protein, biological replicate,
timepoint): θ_cell = Σ wᵢθᵢ / Σ wᵢ, wᵢ = 1/σᵢ²; points without a usable band take the cell
median. Peptides seen in one condition, single-condition rollups and the `exponential` model
are unchanged. The collapsed points are written to `riana_rollup_fractions.txt` with their
variance and Satterthwaite df (`fs_var`, `fs_df`, informational only since `wls-var` was
removed).

**Why pooling the level is enough.** The per-condition bands were each estimated from a
peptide's ~10–14 residuals in that condition, so the same peptide drew weights that differed
two- to three-fold between conditions by chance (11-12 manual §7: band-variance ratio 5–95 %
0.29–2.54 on ATF6, about what F₁₂,₁₂ gives). With within-protein peptide k CV ~0.21, those
different mixtures produced a Δk that no peptide had. One level per peptide makes the mixture
identical (up to gaps in the data).

### 4.2 The fit (`linear_model.fit_linear_deltak`)

Per (experiment, protein), on the collapsed points with t > 0 and φ above the plateau limit:

φ = −k_c · t + ε, one slope per condition, no intercept.

- *Weights:* an OLS pass gives θ̂; one WLS refit with w = (1 − θ̂)² (`wls`, default since
  2026-07-13). Never the observed θ (it biases k by −8 %), and one IRLS step only.
- *Δk* = k_test − k_ref (`--reference-condition` / `--test-condition`, else alphabetical). Its
  SE comes from the weighted covariance with the residual variance estimated from the protein's
  own points, on points − conditions df. Two-sided t-test; BH per experiment
  (`delta_k_p_adj`).
- *More than two conditions:* naming the pair contrasts those two, but the residual variance is
  pooled over every condition in the data.
- `--linear-weights ols` restores the unweighted fit (audit only; false-positive rate 0.285 in
  the 07-13 simulation).

## 5. Peptide level (`--level peptide`, E+time)

Per (experiment, protein), on the admitted peptides' points of the compared conditions
(t > 0, φ above the plateau limit):

φ = −(k_c + b_p)·t + (u_time + v_animal [+ w_cell])·t + ε, Var(ε) ∝ 1/w

- **k_c** one slope per condition; **b_p** fixed sum-to-zero peptide offsets shared by the
  conditions, so k_c is the peptides' unweighted average k and Δk a within-peptide contrast.
- **u_time** ~ N(0, σ²_time), shared by every sample taken at that labeling time, both
  conditions; it cancels from Δk when both conditions are sampled at each time.
- **v_animal** ~ N(0, σ²_animal), nested in time. **w_cell** (animal × condition) is added
  when animals carry both conditions; the animal's shared part then cancels from Δk. Pairing
  is detected from the data: no flag.
- **Weights** w = (1 − θ̂)² from an OLS pass on the same fixed effects, normalized to mean 1.
- **Fit:** REML (statsmodels `MixedLM` on the √w-scaled data, groups = labeling time), Powell
  then Nelder–Mead (BFGS stalls with nested components and few groups).
- **Covariance:** the GLS one at the fitted components, (X′V⁻¹X)⁻¹ — not the inverse Hessian,
  which is ill-conditioned near a zero component (24 negative variances on ATF6 before this).
- **df:** animals − 1 when paired, else animals − times − 1.
- **Fallbacks** (`random_effects` column): `time` → `animal` when the design is not
  cross-sectional (an animal at several times), has fewer than three times, or the time model
  fails (per-animal slope; paired: a slope per condition with a free covariance, then one
  shared slope; BFGS then Nelder–Mead; df animals − 1 or animals − conditions) → `none`
  (WLS, no random term) when that fails too or σ²_animal = 0.
- **Per-condition k CI:** from the same covariance, at times − 1 df when the time slope is in
  the model, else the Δk df. The Δk test is what was calibrated; the per-condition interval is
  the matching containment interval, not separately validated.
- **A condition a protein lacks** (or with fewer than two animals left) leaves the other
  condition's k only, fitted by the same model on one condition. An experiment with more than
  two conditions needs the pair named.
- **Outputs** (`riana_rollup_proteins.txt`): `delta_k_df`, `n_animals`, `paired`,
  `random_effects`, `sigma_time`, `sigma_animal`; `n_points` counts peptide points and
  `R_squared` is NaN; `peptide_median_k` is written at both levels. The collapsed points are
  still written to `riana_rollup_fractions.txt` for display.
- **Estimand.** The unweighted average over peptides (the protein level weights by
  precision). E's k level sits 1–3 % above it (φ Jensen plus additive offsets on a
  multiplicative truth); Δk is unaffected within 0.15 SE; the ratio reads ~2 % low at 1.2
  (eval §8). A second-order correction, φ + σ²/(2(1 − θ)²), applies if peptide-level k values
  are reported as levels.
- **Cost:** about 1 min for ATF6's 3,507 Δk tests at `-W 14`.

## 6. Interaction (`--interaction`)

Both forms compare `--reference-condition` with `--test-condition` (or the other condition
when there are exactly two), over the experiments that have animals under both conditions,
and report per protein, each BH-adjusted across proteins:

- `interaction_*` — **ratio scale (primary):** is k_test / k_ref the same in every
  experiment? A uniform proportional effect is not an interaction on this scale.
- `interaction_additive_*` — are the Δk equal across experiments?
- `mean_delta_k_*` — Δk averaged over experiments, the condition main effect.

Per-experiment Δk tests stay in `riana_rollup_proteins.txt`.

### 6.1 Protein level (`interaction_model.fit_interaction`)

On C's collapsed points: φ = Σ_s β_s·t·[s] + β_shift·t·[test] + Σ_j γ_j·t·[test]·C[s,j] +
u_animal·t + ε, `wls` weights, REML with BFGS (L-BFGS stops early on variances ~1e-5).
Ratio-scale df: Satterthwaite (lmerTest `contestMD`, checked against lmerTest 3.2.1);
additive and mean Δk: containment (nlme). A singular fit (σ²_animal at 0, ~1 % of ATF6
proteins) is computed directly as the WLS solution. Checked against nlme on ATF6 (estimates
within 1.5e-3, F within 1.3 %). Seconds for ATF6.

### 6.2 Peptide level (EJ+time, `peptide_model.fit_peptide_interaction`)

φ = −(k_{s,c} + b_{s,p})·t + (u_time + v_animal + w_cell + x_tissue×time)·t + ε

peptide offsets within each experiment; random slopes per labeling time, animal, animal ×
experiment (cell) and experiment × time. Containment df: mean Δk animals − times − 1; additive
and ratio scale cells − animals − (S − 1)·times − (S − 1), with the tissue × condition × time
stratum pooled with the cells (S = experiments). Adds `random_effects`, `sigma_cell`,
`sigma_time`, `sigma_tissue_time`. About 18 min for ATF6 (979 proteins, `-W 14`).

## 7. Gates, admission and the model set

- **Peptide gate (rollup):** R² ≥ `--min-r2` (0.8), or the flat-curve rescue R² ≥ 0.6 and
  k_cv < 0.2; per condition. Single-timepoint experiments bypass R² and use
  `--min-fit-points` (auto 2) plus k_cv.
- **Spep floor (fit):** 8 labeling sites (D₂O), 6 (¹⁸O).
- **`--min-peptides`** (default 1 since 2026-11-17, was 2): peptides passing the gates for a
  protein to be reported in a condition.
- **`--peptide-admission`:** `any` (passed in one condition ⇒ kept in all; anti-conservative),
  `own` (kept where it passes; the conditions can rest on different peptides), `all` (kept only if
  it passes in every condition it appears in). The default `auto` (2026-11-18) is `all` for the
  protein-level linear Δk and `own` for `--level peptide` and for the exponential / calibration
  models (§8.6).
- **The contrast pair (2026-11-18):** with `--reference-condition` and `--test-condition` named, a
  linear rollup keeps only those two conditions, so admission, the shared bands and the joint fit
  do not depend on other conditions in the table; without both, an experiment may hold at most
  two.
- **A shared peptide (2026-11-17 peptide level, 2026-11-19 protein level):** a Δk needs at least
  one peptide measured in both conditions; otherwise each condition gets its k and no Δk, and the
  interaction leaves that experiment out. Peptide identity ignores the dimethyl channel labels
  (UNIMOD 36 / 199 / 330).
- **Models:** `exponential` (fit and rollup default), `linear` (rollup), `calibration` (fit,
  auto-selected on calibration manifests). `simple` and `linear simple` are accepted with a
  warning. guan / fornasiero, `--kr`, `--rp` and rollup's `--kp` are removed; a rising
  precursor is handled at fit by `--kp` (fs then holds f, and the rollup works on f unchanged).
- **`--level`** replaces `--method`: `protein` = the old `weighted`; `pooled` (every peptide
  point as an independent replicate) is removed; `--method weighted` is still accepted with a
  warning.

## 8. Calibration evidence — what each yardstick showed

Yardsticks (eval §2): **T1** boomi replicate split (20 assignments); **T2** ATF6 genotype labels
permuted between animals within timepoint (20 permutations); **T3** Δk from odd vs even
timepoints, z robust SD of the difference (band 0.9–1.15); **T6** ATF6 chamber Δk in control vs
KO animals (same band); **S** a simulation calibrated to ATF6 (false-positive rate ≤ 0.065,
coverage ≥ 0.93, power ≥ 0.9 × A). A = the old per-condition collapse with `wls`.

### 8.1 Per-experiment Δk

| Criterion | A (old) | **C (protein level)** | **E+time (peptide level, Powell + GLS)** |
|---|---|---|---|
| T1 boomi null (worst turnover quartile) | 0.059 (0.072) | 0.028 (0.041) | 0.025 (0.028) |
| T2 ATF6 null (worst quartile) | ✗ 0.096 (0.108) | 0.039 (0.046) | 0.042 (0.048) |
| T3, seven tissue-pair readings | ✗ fails 5 of 6 Boulder pairs and lve_atr (1.07–1.45) | 0.94–1.13, all pass | 0.96–1.07; ✗ ATR–EDL 1.26 |
| T6 ATF6 LV–RA (KO brake scaled out) / LV–RV | ✗ 0.85 / ✗ 0.80 | ✗ 0.67 / ✗ 0.51 (too wide) | 1.08 / ✗ 0.88 |
| S false-positive rate, coverage | ✗ 0.094, ✗ 0.898 | 0.047, 0.953 | 0.065 (on the bound), 0.936 |
| S power at ratio 1.2 vs A: nominal / equal size | 1 / 1 | 0.98 / 1.29 | 1.17 / 1.43 |

Sources: eval §3, §4b, §5. T3 on the peptide level is read against 1.04, its small-df
reference; each entry has an SE of ±0.07–0.11. The simulated null on its bound comes from the
two datasets with the largest proteome-wide group shift (S8).

### 8.2 Interaction

| Criterion | A (old) | **C + shipped interaction** | **EJ+time** |
|---|---|---|---|
| T2 ATF6 null: ratio / additive / mean Δk | ✗ 0.216 / ✗ 0.245 / 0.044 | 0.065 / ✗ 0.081 / 0.019 | 0.027 / 0.041 / 0.033 |
| S null: ratio / additive / mean Δk | ✗ 0.194 / ✗ 0.205 / 0.045 | 0.044 / 0.052 / 0.027 | 0.039 / 0.049 / ✗ 0.074 |
| S mean-Δk power at ratio 1.2: nominal / equal size | 0.741 / 0.782 | 0.736 / 0.850 | 0.919 / 0.886 |

EJ+time's mean-Δk null of 0.074 is S8 at its strongest (two of eight datasets; pooled SE
0.038); averaging tissues does not average out a shift every tissue shares.

### 8.3 Rising precursor (f route), protein level

From the 11-12 manual §6.2 (envelope simulation through the CLI, k_p 0.08/day): `wls` on f
rejects 0.042 at α 0.05, calibrated in every turnover quartile, CI coverage 0.96, ratio
recovered. θ on the same data biases Δk by −59 %. Collapse C on f and the peptide level on f
were not simulated.

### 8.4 Real calls, BH < 0.05

| Dataset · test | A (old) | C | E+time / EJ+time | Source |
|---|---|---|---|---|
| ATF6 per-chamber Δk (3,507 tests; null per permutation A 19.7) | 35 | 0 | 0 | eval §3 |
| ATF6 interaction, ratio scale (979 proteins; null per permutation A 69) | 79 | 2 | 0 | eval §4, CHANGELOG |
| lve_atr atrium vs LV (348 tests) | 234 | 235 | 229 | eval §7.3 |
| Boulder EDL–SOL (250) | 152 | 136 | 162 | eval §4b |
| Boulder LVE–EDL / LVE–SOL / ATR–EDL / ATR–SOL | 160 / 177 / 186 / 269 | 156 / 171 / 186 / 268 | 154 / 168 / 189 / 271 | eval §4b |

On ATF6, the notebook's two interaction hits: **Smyd1** is mechanism 1 (C p 0.28, EJ+time
0.78); **Cops3** survives C (p ≈ 0) but not EJ+time (0.19) — it follows an outlying LA-control
cell that EJ+time's cell variance carries and the protein-level model does not.

### 8.5 The refreshed runs (2026-11-19, current defaults)

The 2026-11-16 fits, rolled up with the 2026-11-17/19 rollup changes: one peptide per protein,
admission `auto`, the contrast pair, a shared peptide for Δk, and dimethyl forms matched without
their label (`runs/README.md`, `runs/refresh_2026-11-1{7,8,9}.sh`). Boulder's six tissue pairs
are rolled up from the one four-tissue fit. BH < 0.05 calls of Δk tests, at the defaults of each
level:

| Dataset · test | 2026-11-16 (2 peptides), protein | **Protein level** (`all`) | **Peptide level** (`own`) |
|---|---|---|---|
| Boulder LV vs atrium | 235 of 348 | 308 of 547 | 286 of 547 |
| Boulder EDL vs SOL | 136 of 250 | 196 of 424 | 210 of 424 |
| Boulder LV vs EDL / LV vs SOL | 156 of 206 / 171 of 316 | 237 of 330 / 260 of 491 | 225 of 330 / 234 of 491 |
| Boulder ATR vs EDL / ATR vs SOL | 186 of 198 / 268 of 317 | 274 of 301 / 377 of 487 | 276 of 301 / 371 of 487 |
| ATF6 per chamber, KO vs control | 0 of 3,507 | 0 of 5,513 | 0 of 5,514 |
| ATF6 interaction: ratio / additive / mean Δk | 2 / 3 / 0 of 979 | 2 / 3 / 0 of 1,515 | 0 / 0 / 0 of 1,511 |
| SCVI480 ¹⁸O, mesoderm vs iPSC | 108 of 293 | 239 of 672 | 245 of 672 |
| dimethyl duplex, rapamycin vs control (single timepoint) | 30 of 79 | 49 of 142 | 36 of 143 |

The 2-peptide Boulder column is the evaluation's two-tissue tables (eval §4b; LV vs atrium is the
lve_atr re-search). Its C column is the protein level then.

- **One peptide per protein gives 1.5–2.3× the tests** on every set, and the calls rise with
  them wherever there were calls. SCVI480's new single-peptide calls include five
  sterol-synthesis enzymes, significant at both levels and under both admissions (¹⁸O note §6.4).
- **The two levels now test the same proteins** (to within 1 %): both need a peptide measured in
  both conditions. Their calls stay within 11 % of each other, except on dimethyl's single
  timepoint (49 vs 36).
- **ATF6 has no per-chamber call at either level.** At 1 peptide under the 2026-11-17 rules, the
  protein level called one: Dazap1, right atrium, Δk −0.10, p_adj 0.04. Its control k rested
  on one peptide and its KO k on another, so the Δk compared different peptides. The
  shared-peptide rule removed that test. The interaction keeps Cops3 and Dync1h1 (§8.4).
- **Dimethyl** gained most from the label fix: before it, the light and heavy forms of a peptide
  were different peptides to the rollup, so `all` acted as `own` and the shared band level did
  not apply.
- **Pair scoping** (Boulder LV vs atrium): from the four-tissue table without it, 484 tests;
  scoped to the pair, 569, byte for byte the result of a table cut to the two tissues; with the
  shared-peptide rule, 547.
- The 2026-11-16 refresh against the archived pre-redesign rollups (June–August) changed:
  - lve_atr: 275 of 467 → 235 of 348;
  - ATF6 per chamber: 37 → 0;
  - ATF6 interaction: 86 / 137 / 1 → 2 / 3 / 0;
  - SCVI480: 225 of 652 → 108 of 293;
  - dimethyl: 33 of 77 → 30 of 79.

  Those changes came from the collapse plus the August admission and depth fixes (this
  manual's 2026-11-16 text).
- The ATF6 notebook (`notebooks/atf6_biology.ipynb`) and the ¹⁸O technical note were regenerated
  from these runs on 2026-11-19.

### 8.6 Admission, one peptide and the contrast pair (2026-11-17/19)

Evidence: `runs/deltak_eval/admission/` (local; `admission_eval.py`, `summarize.py`). The nulls
are those of §8.1 at `--min-peptides 1`, split by the protein's smaller per-condition peptide
count. Refits: T1 20 and T2 20 at the protein level, 10 and 4 at the peptide level. The first
four rows predate the shared-peptide rule; the last two rerun the protein-level default under it.

| Null, rejection at α 0.05 | Protein `all` | Protein `own` | Peptide `all` | Peptide `own` |
|---|---|---|---|---|
| T1 boomi split: every / 1 peptide / ≥ 2 | 0.036 / 0.046 / 0.028 | ✗ 0.087 / ✗ 0.103 / ✗ 0.074 | 0.035 / 0.050 / 0.023 | 0.028 / 0.045 / 0.016 |
| T2 ATF6 permutation: every / 1 peptide / ≥ 2 | 0.047 / 0.061 / 0.039 | ✗ 0.181 / ✗ 0.225 / ✗ 0.157 | 0.035 / 0.062 / 0.020 | 0.031 / 0.059 / 0.019 |
| T2 false BH calls per permutation (real ATF6 calls then) | 1.5 (1) | ✗ 370 (399) | 0 (0) | 1.0 (0) |
| T3 lve_atr time split, z robust SD (band 0.9–1.15) | 1.08 | ✗ 1.31 | — | — |
| **Shared-peptide rule:** T1 / T2, every (1 peptide) | 0.035 (0.044) / 0.046 (0.059) | — | — | — |
| **Shared-peptide rule:** T2 false BH calls; T3 | 1.25; 1.06 | — | — | — |

- **Protein-level `own` fails every null.** Each condition keeps the peptides that pass in it,
  so the two k can rest on different peptides, and the Δk carries the difference between
  peptides. On ATF6 it makes ~370 false calls per permutation, about as many as its real calls
  (399). Hence `auto` = `all` at the protein level.
- **The peptide level holds under both admissions.** Its Δk is within-peptide; a peptide that
  passes in one condition only informs that condition's k. Compared with `all`:
  - `own` gives a k for more proteins (ATF6 5,732 → 7,794; SCVI480 1,013 → 1,547) on the same Δk
    tests;
  - its Δk agree (Spearman ρ 0.97–1.00; robust SD of the difference 0.02–0.06 SE);
  - it makes 0–9 % fewer calls.

  Hence `auto` = `own` at the peptide level. The exponential model has no Δk; `own` gives each
  condition the k of its own passes (Boulder per-tissue k: 960 → 1,165 proteins).
- **One peptide is calibrated where the admission is.** The single-peptide stratum rejects
  0.044–0.062 under protein `all` and at the peptide level. That is slightly above 0.05 on the
  ATF6 permutation, and the ≥ 2 stratum compensates (0.016–0.039).
- **The protein-level Δk is paired only in part.** The gap between C and E+time, as robust SD
  of (C − E)/SE_E, grows with a protein's fraction of unshared peptides:

  | Unshared fraction | Under `all` (ATF6, Boulder LV–ATR and EDL–SOL, SCVI480) |
  |---|---|
  | < 0.3 | 0.45–0.78 |
  | 0.3–0.6 | 0.94–1.33 |
  | > 0.6 | 0.73–1.86 |

  `own` moves many more proteins into the imbalanced strata (ATF6: 765 above 0.6, against 15
  under `all`). The peptide level is the paired estimator (§9).
- **Under the shared-peptide rule `all` and `own` test the same proteins.** A shared peptide must
  pass in both conditions under either policy. Before the rule, `own` tested 3–16 % more
  (`own` vs `all` on proteins tested by both: ρ 0.93–0.99 on the Boulder pairs, 0.84 SCVI480,
  0.60 ATF6). The rule removed the cross-peptide tests (SCVI480 674 → 672, Boulder LV–ATR
  569 → 547, ATF6 5,532 → 5,513). It moved the T2 null from 0.047 to 0.046 and the false calls
  from 1.5 to 1.25 per permutation.
- **The contrast pair.** Admission `all` over a four-tissue table needed a peptide to pass in
  all four tissues. That cost each Boulder pair 8–18 % of its tests against the pair alone
  (LV–ATR 484 vs 569), so another tissue's data changed an A-vs-B result. Scoped to the named
  pair, the full-table rollup equals the two-tissue table byte for byte, at both levels.

## 9. Choosing a level

- **Protein level (default)** when the design is unpaired and the effects of interest are large
  compared with the within-protein peptide spread, or when the run must be fast. It is the
  more conservative of the two on paired contrasts (T6) and holds every Δk null with margin.
- **Peptide level** for small effects, paired designs (tissues of the same animals as
  conditions) and the interaction, where it is the only form that holds the additive null. It
  needs at least three labeling times for the time slope. It has now run on every
  multi-condition dataset here (§8.5) and stays opt-in (open item 4). Its default admission is
  `own`, which gives a k for more proteins at the same Δk tests (§8.6).
- Report the level with the results (`level` column of `riana_rollup_proteins.txt`).

## 10. How the estimator got here

| Date | Change | Why | Evidence · commit |
|---|---|---|---|
| 2026-06-13 | `linear simple`: φ-space OLS, joint per protein, Δk + BH | cross-condition Δk with a shared residual variance | lve_atr · `e558d0a` |
| 2026-07-13 | `wls` default; t = 0 dropped | log-transform heteroscedasticity | false positives 0.285 → 0.062 · `4ae5e50`, report 07-13 |
| 2026-07-30 | `wls-var` opt-in | per-point precision (S2) | report 07-24 · `c82cc55` |
| 2026-08-05 | per-condition gates; `--peptide-admission` | a pass in one condition admitted the peptide in the other | `b26453a` |
| 2026-08-10 | per-experiment BH | one family per tissue | `180e477` |
| 2026-11-10/11 | `--interaction` with a random animal slope, Satterthwaite df, `individual` | split-plot designs | report 11-10 · `b8715f7`, `57185c0` |
| 2026-11-12 | rising-precursor fit (`--kp`): fs holds f; f route uses `wls` | an unprimed D₂O precursor rises over days | report 11-11, manual 11-12 §6 · `3a526c3`, `9a2d357` |
| 2026-11-12 | known issue: per-condition collapse mixes peptides differently (mechanism 1) | ATF6 37 → 0 calls under shared proportions | manual 11-12 §7 |
| 2026-11-15 | **`wls-var` removed** | never better calibrated than `wls`; worse on the ATF6 null (0.122 vs 0.096) | eval §3, §5 Q1 · `0b275e3` |
| 2026-11-16 | **C: one band level per peptide** across conditions | removes mechanism 1; holds every Δk null | eval §3–§5 · `a31f356` |
| 2026-11-16 | models `exponential` / `linear` / `calibration`; guan / fornasiero removed; admission `all` | θ-curve fits with a lagging precursor biased k +5–42 %; one admission rule | report 11-11 · `c2b4ff1` |
| 2026-11-16 | **`--level peptide`** (E+time; EJ+time for the interaction) | the only candidate inside the G2 thresholds with more power than C | eval §4b, §5 Q3 · `be2d03d` |
| 2026-11-17 | `--min-peptides` 2 → 1; peptide level: no Δk without a shared peptide; E fallback rejects an indefinite covariance | the two-peptide rule left most proteins untested; single-peptide tests hold the nulls | §8.6 |
| 2026-11-18 | admission `auto` by level; a linear rollup scoped to its pair; dimethyl forms one peptide | protein-level `own` fails the nulls; other conditions changed an A-vs-B result | §8.6 |
| 2026-11-19 | protein-level Δk needs a shared peptide | a Δk on different peptides made ATF6's one per-chamber call | §8.5, §8.6 |

Not built (eval §5): D, a random-effects collapse (it behaved like C everywhere); E without
the time slope (conservative wherever samples share a timepoint deviation; it is E+time's
fallback); E+band (3 % more power, breaks the simulated bound); multiplicative peptide offsets.

## 11. Open items

1. **ATR–EDL time split (1.26).** The one reading where E+time's SE is too narrow; C (1.13) and
   E (1.14) are high there too. Look again when a second heart–muscle dataset exists.
2. **Longitudinal designs** (one animal at several times): E+time falls back to E; untested, and
   a run effect beside the animal slope may be needed.
3. **The f route at the peptide level:** not evaluated (needs the simulation; the proprietary
   rising-precursor set can be checked locally only).
4. **Default level:** whether `--level peptide` becomes the default after more designs. One
   argument for it: the protein-level Δk drifts from the peptide-level one as a protein's
   unshared-peptide fraction grows, even under `all` (§8.6). The rollup default model stays
   `exponential` (user, held).
5. **T6 LV–RV (0.88):** chamber-specific peptide behaviour shared by both animal sets. A
   peptide × condition variance component would carry it; not tested. It belongs in the SE if
   the estimand is the protein's peptides at large, not if it is the peptides measured.
6. **Residual-variance moderation (limma-style):** high priority, deferred; to be done on the
   peptide-level model's per-protein σ²_time and σ²_animal (PROJECT_REVIEW §3 box).
   `linear_model._fit_fdist` is kept for it.
7. **Whole-group shifts (S8):** decided 2026-11-15 to test per-protein Δk as is; a
   proteome-level test of a group's shift is a separate question, not built.
8. **k measurability:** a k whose CI reaches 0 reads as a long t½ (φ clamps θ at 0.001); a
   flag and a t½ bound in the output are a TODO (PROJECT_REVIEW §3).
9. **Additive vs proportional animal effects:** the four-tissue Boulder check found the shared
   day deviation additive in k proteome-wide, but with one mouse per day animal and time cannot
   be separated; the clean test needs several animals per day with tissues of different k.
10. **3b ceiling drop:** dropping reads above the precursor ceiling selects fast proteins' early
    reads one-sidedly on the f route; capping instead is to be revisited.
11. **Single-peptide tests on the ATF6 permutation (0.059–0.062).** Slightly above nominal at
    both levels; the ≥ 2-peptide stratum is below it, so the pooled rate holds (§8.6). Watch it
    where single-peptide proteins dominate a call list.

## 12. Reproduce

```bash
# Protein level (default) and peptide level, two conditions; add --interaction --plot for
# a multi-experiment (e.g. per-tissue) manifest with animals shared across experiments
riana rollup --manifest <fit>/riana_manifest.tsv --model linear \
    --reference-condition <ref> --test-condition <test> -W 14
riana rollup --manifest <fit>/riana_manifest.tsv --model linear --level peptide \
    --reference-condition <ref> --test-condition <test> -W 14 -o <fit>/peptide_level

# The refreshed runs this manual's §8.5 reads (local, gitignored): fits, then the rollups
bash runs/refresh_2026-11-16.sh
bash runs/refresh_2026-11-17.sh; bash runs/refresh_2026-11-17_peptide.sh
bash runs/refresh_2026-11-18.sh; bash runs/refresh_2026-11-19.sh

# §8.6: admission × level on the real sets and the nulls (local)
python runs/deltak_eval/admission/admission_eval.py null|null_peptide|real|peptide -W 22
python runs/deltak_eval/admission/admission_eval.py null --tag _shared --admissions all -W 22
python runs/deltak_eval/admission/summarize.py; python runs/deltak_eval/admission/levels_compare.py

# The evaluation (artifacts moved 2026-11-16 to runs/deltak_eval/; the reports cite the old
# runs/<name> paths — the map is in runs/README.md)
python tests/benchmark/wlsvar_eval.py null runs/deltak_eval/wlsvar_eval/atf6/t2/t2_0* \
    --ref control --test ATF6-KO --variants A C E+time-powell --cache
python tests/benchmark/wlsvar_sim.py --work runs/deltak_eval/wlsvar_sim/null_ts --ratio 1 --pooled
```

Variants that use removed code need an older Riana: `wls-var` (B, D+var) ≤ `e12980c`;
guan / fornasiero ≤ `a31f356`; the pre-C per-condition collapse ≤ `cc8a98b` (the
`eval_rising_variance` captured frames at HEAD carry C's bands). Prototypes:
`tests/benchmark/wlsvar_*.py`, `tests/benchmark/eval_rising_variance.py`,
`tests/benchmark/bench_rising_null.py`.
