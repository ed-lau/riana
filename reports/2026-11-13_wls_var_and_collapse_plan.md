# wls-var, the random-effects collapse (D) and the peptide-level model (E): scope of the evaluation

2026-11-13 · branch `1.2.0` · **Status: evaluation complete 2026-11-14; verdict,
recommendations and the user's decisions (2026-11-15) in the
[evaluation report](2026-11-14_wls_var_D_E_evaluation.md) §5.** Scope approved by the user 2026-11-13 (answers in §7). Follows
[the Δk manual](2026-11-12_deltak_inference_manual.md) §8.1, which defines candidates A–E and the
evidence order. This document fixes what is evaluated, on which data, with which tests, and the
evidence that must be in hand before anything changes in `riana/`.

## 1. What is being decided

| | Question | Possible outcomes |
|---|---|---|
| Q1 | Does `wls-var` weight cells by their real precision? | keep opt-in, with its limits documented · deprecate |
| Q2 | Should D (random-effects collapse, τ²) replace `wls-var` inside `--method weighted`? | build · do not build |
| Q3 | Should E (peptide-level mixed model) become a new `--method`? | build · defer · do not build |

**Rule (user, 2026-11-13).** A change goes into `riana/` only when the evidence shows it works on
real data. If real data cannot show that, the outcome is a written evaluation and a design note,
not code. During the evaluation, prototypes of D and E live in the evaluator
(`tests/benchmark/`), as the research `pool2` model does now. They are not CLI options.

## 2. What the existing runs can show (audit, 2026-11-13)

Real data has no ground truth, so each dataset is useful only for the tests its design allows
(§3). Gated statistics approximate the shipped rollup gate: R² ≥ 0.8, or R² ≥ 0.6 with
k_cv < 0.2; Spep ≥ 8; admission `all`.

| Dataset (run) | Design | Timepoints after 0 | Replicates per timepoint | Protein × condition groups with one admitted precursor (dropped before testing) | Multi-peptide tests: cells holding every admitted peptide (median) / fully complete | Fit made at | Tests it supports |
|---|---|---|---|---|---|---|---|
| ATF6 (`atf6_individual/fit`) | DIA, mouse heart; KO vs control × 4 chambers, all 4 from each animal | 5 (days 1–14, no day 0) | 3 animals per genotype | 40 % | 0.60 / 7 % | `33f0939` (Aug) + `individual` added | T2, T6, T4, interaction; animal-variance parameters |
| lve_atr (`lve_atr_clean`) | DDA, mouse; atrium vs ventricle, both from the same mouse | 11 (days 1–30) | 1 mouse | 44 % | 0.64 / 11 % | `2ff71be` (Jun) | T3, T4 |
| boomi (`boomi_ipsc_d2o`) | DDA, iPSC culture; one condition | 8 (1–24 h) | 2 culture plates | 43 % | 0.44 / 2 % | `e946b6f` (Jul) | **primary:** T1, T3, T3′, Q1 and τ² diagnostics |
| lauren9 (`timeseries_lauren9`) | DDA, 8-fraction 2D-LC, iPSC; one condition | 11 (0.5–24 h) | 2 (nature unconfirmed) | — | — | `e946b6f` | T1, T3 (secondary: fraction pseudoreplication, see below) |
| lauren5_7 (`timeseries_lauren5_7_…_o18`) | ¹⁸O, 8 fractions; control vs mesoderm | 11 | 2 per condition | 67 % | 0.77 / 16 % | `e946b6f` | balanced T2, T3 (secondary) |
| proprietary non-primed set | f route, 1 animal per cell | 9 | 1 | — | — | HEAD | T3 only, local only, no numbers in committed files |
| EDL / SOL | DDA, mouse; EDL vs soleus, **the same mice as lve_atr** (all four tissues at a timepoint come from one mouse) | as lve_atr | 1 mouse | — | — | not run; the user is adding the data | T3, T4; with lve_atr, the four-tissue animal diagnostic (WP1d) |

`timeseries_dia` (an ATF6 subset, 3 timepoints) and `juber_ac16_d2o` (8 % of peptides pass the
gate) are left out.

**Findings that change the plan:**

1. **No dataset has complete peptide sets, DIA included.** In ATF6, a median 40 % of cells miss at
   least one admitted peptide, and only 7 % of multi-peptide tests are complete in every cell.
   Manual §8.1 assumed that complete data, such as ATF6 DIA, allow identical proportions in every
   cell. They do not.
   - C and D therefore cannot make the two conditions' mixes identical cell by cell on any real
     dataset.
   - Gaps that occur at random add noise. Gaps that depend on the condition (a peptide less often
     seen in one tissue) can still bias C and D.
   - E is the only candidate that does not depend on complete cells.
   - Completeness depends on the depth gate (`riana fit -d`). The table's runs were fitted at
     depth 3 (ATF6, 5 timepoints) and 6 (lve_atr and boomi, 9–12 timepoints), similar
     strictness. Refitted at the default depth 3, lve_atr falls from 0.64 to 0.33. WP0 keeps
     each run's original depth.
   - Manual §8.1 to be corrected when this plan is adopted.
2. **40–67 % of protein × condition groups hold one admitted precursor.** *Corrected 2026-11-14:*
   this first read "40–67 % of tests". The rollup drops those groups (`--min-peptides 2`, counted
   on precursors, sequence + charge), so they never reach the Δk test. Every test rests on ≥ 2
   precursors, median 4. A–E differ on all of them, though a protein whose two precursors
   are one sequence in two charge states has little peptide spread to mix.
3. **Every fit predates the envelope-binning fix `d1eb480`**, and boomi and lve_atr predate the
   1.2.0 rail changes. Refit at HEAD before scoring (WP0).
4. **Which designs can give a real-data null:**
   - boomi: two culture plates per timepoint (separate plates, not reinjections; user,
     2026-11-13), so one plate can be split against the other (T1), a null that includes
     plate-to-plate biological variation. With D₂O, 8 timepoints and true replicates, **boomi is
     the primary real-data dataset.** lauren9 has the same layout, but its replicates are
     unconfirmed and it is fractionated (finding 7).
   - ATF6: 3 + 3 animals per timepoint, so genotype labels can be permuted between animals (T2).
   - lve_atr: one mouse per timepoint. Swapping tissue labels flips the sign of a large real
     effect, so it is not a null. lve_atr supports only T3 and T4.
5. **The missing design stays missing: DDA, in vivo, with ≥ 2 animals per condition per
   timepoint.** EDL / SOL comes from the same single mouse per timepoint as lve_atr. Two things
   soften the gap:
   - After gating at similar depth strictness, ATF6's DIA cells are about as incomplete as
     lve_atr's DDA cells (0.60 vs 0.64 of cells hold every admitted peptide, finding 1), so T2 on
     ATF6 already tests gaps of a similar size.
   - Gaps that follow the label or the intensity, a DDA-specific pattern, go into the simulation
     (WP3).
6. **ATF6 has 5 timepoints, too few to split in time (T3).** Its null is the permutation (T2). For
   an SE check with real effects, it splits by animals instead (T6).
7. **The lauren sets are fractionated**, and the fit still pseudoreplicates multi-file points
   (PROJECT_REVIEW, fraction-collapse gap). Use them only as a second check.
8. **Every in-vivo tissue contrast is within animal.** All four ATF6 chambers come from each
   animal, and lve_atr, EDL and soleus from one mouse per timepoint. The shipped one-factor Δk
   ignores that pairing (manual §8.2, *Within-animal comparisons*). If an animal's deviation is
   shared across its tissues, it cancels from the tissue Δk but stays in the residual variance, so
   the SE is too large: tissue contrasts should give a T3 or T6 z SD **below** 1. The lower bound
   of the z-SD threshold (§5) is there to catch that. If the deviation scales with tissue k
   (proportional, manual §8.2), it does not cancel. WP1d measures which it is.

## 3. Truth-free tests: what each one measures

| | Test | What it measures | Requirements | Datasets | Cannot see |
|---|---|---|---|---|---|
| T1 | **Replicate split.** At each timepoint, assign one replicate at random to pseudo-condition A and the other to B; refit (`riana fit` on a relabeled manifest); roll up. | False-positive rate on real peptides and real noise; the true Δk is 0 by construction. Includes mechanism 1 from chance band differences. | ≥ 2 replicates per timepoint; a refit per assignment | boomi (128 distinct assignments), lauren9 | Peptide quality that differs between real conditions (there are none); animal variation (cell culture) |
| T2 | **Permutation within timepoint.** Permute condition labels between animals; an animal's tissues move together. Refit, roll up, run the interaction. | False-positive rate of Δk and of the three interaction tests, with real animals | ≥ 2 animals per condition per timepoint; a refit per permutation | ATF6; lauren5_7 balanced (each pseudo-condition gets one sample of each real condition) | The small real genotype effect becomes noise |
| T3 | **Time split** (the user's alternate-timepoint idea). Refit the odd- and even-numbered timepoints separately; both halves span the whole range. | Δk_odd − Δk_even has mean 0 whatever the true Δk. So z = difference / √(SE²_odd + SE²_even) should have SD 1 if the SEs are honest. Var(difference) measures precision. Works with real effects present. | ≥ 8 timepoints; 2 refits per dataset | lve_atr, boomi, lauren, proprietary (local) | Bias shared by both halves, e.g. a peptide genuinely measured better in one tissue. Each half has half the points, so its bands are noisier than the full data's. |
| ↳ | *T3 on a one-condition set* (boomi) | The same z on per-protein k instead of Δk: checks the SE that Δk's SE is built from (Q1's variance model), not mechanism 1. Each half has 4 timepoints after 0, so it is fitted at depth 3. | | boomi | Mechanism 1 (needs two conditions) |
| T3′ | **Leave one timepoint out.** Predict each held-out cell from the rest. | Standardized prediction residuals per cell: does a variance model (`wls`, `wls-var` Ṽ, D's Ṽ_RE, pool2) track the actual scatter? Answers Q1 directly. | as T3, no refit for the variance check | all | Bias |
| T6 | **Animal-set split** (user, 2026-11-13). Estimate one chamber contrast (LV vs RV, LV vs RA) separately in the control animals and in the KO animals. The genotype barely changes turnover, so both sets estimate the same contrast. | z = (Δk_ctrl − Δk_KO) / √(SE²_ctrl + SE²_KO) should have SD 1. The in-vivo counterpart of T3: independent animals and all 5 timepoints in each half. Also tests whether a paired (within-animal) model gets the SE right (finding 8). | Negligible genotype × chamber interaction: drop the interaction hits and use a robust SD. Allow for the KO brake with one global scale factor, or compare on the log-ratio scale. No refit expected: the fit's curves are per chamber × genotype whichever label is the experiment (to verify in WP0). | ATF6 | Bias shared by both animal sets: mechanism 1 from chamber-specific peptide quality repeats in both halves, so T6 calibrates the SE, not that bias (T2 and T4 cover it) |
| T4 | **Peptides' own Δk** (exists, manual §7) | Bias of the collapse against a collapse-free contrast | 2 conditions | ATF6, lve_atr | Cannot rank E, which is built from the same contrast |
| T5 | **Biorep-split precision** (07-24) | Efficiency only | replicates | boomi, lauren | Mixing bias (weights fixed) |
| S | **Calibrated simulation** (`bench_rising_null.py`) | The only source of truth: power, bias, coverage, ratio recovery | calibration inputs from §2 and §4 | — | Whatever it does not simulate (manual §6.6) |

**Check the yardsticks before using them.** T1–T3 and T6 have never been tested against a known
answer. The Monte Carlo predicts how A, B and C differ: mechanism 1 inflates A (0.141 with site
error alone), C removes it (0.032), and B is no better than A (manual §6.4). So the tests are
first run on A, B and C, which already exist, with no new estimator code. T1 on boomi and T2 on
ATF6 are the two that can see mechanism 1; T3 and T6 check the SEs.

- If the real-data tests separate A from C the way the simulation predicts, they can rank D and E.
- If A and C look the same on them, real data cannot discriminate. The decision then rests on the
  calibrated simulation alone, and D and E stay design-only.

This is the "applicable to real data" gate (G1, §5).

## 4. Work packages

| WP | Content | Output | Compute | Code |
|---|---|---|---|---|
| 0 | Refit boomi, ATF6 and lve_atr at HEAD (θ route), into new run directories; integrate and fit EDL / SOL when its data are on disk. Build the relabeled manifests for T1, T2 and T3; check that T6 needs only a relabel. Measure refit time. If a T2 refit costs about an hour, decide between full refits (10–20 permutations) and a stage-2-only refit in the evaluator: curve and bootstrap on the stored per-point fs, accepted only if it reproduces the shipped k and bands on unpermuted data, as in the manual §6.1 audit. | fits; manifests; cost table | minutes each for lve_atr and boomi; about 1 h for ATF6 at `-W 16` (estimated) | scripts only |
| 1 | **Diagnostics on the refits, no new estimators.** (a) Q1: regress the log squared standardized `wls` residual on log Ṽ within protein (slope 1 = Ṽ tracks the scatter, 0 = uninformative); also by peptides per cell, since a floor rather than a 1/n fall means peptide spread that Ṽ ignores; plus the T3′ version. (b) D's τ² shape: excess spread Q − (n − 1) against θ̂, comparing θ scaling with apparent-k offset (manual §8.1 b), and the pooled τ_rel. (c) Size of τ²·g against σᵢ² at the late, slope-carrying points. If τ² is negligible there, D reduces to C plus a cell variance, and that alone decides whether D is worth more than C. (d) **The animal term, from tissues of the same animal:** correlation of each animal's residuals across its tissues (lve_atr + EDL / SOL: 4 tissues × 11 mice; ATF6: 4 chambers × 30 animals, as deviations from the timepoint mean). It gives the shared animal component, and whether it is additive or scales with tissue k (manual §8.2), which fixes the form of the animal term in E. | Q1 evidence; D's g and τ_rel; the animal term's form | minutes | evaluator functions |
| 2 | **Check the yardsticks with A, B and C:** T1 on boomi first (20–40 assignments), T2 on ATF6 (10–20 permutations, Δk and interaction), T3 on boomi and lve_atr, T6 on ATF6. | null rates, z SD, A-vs-C separation → **G1** | T2 dominates: 10–25 h of refits, or much less with the stage-2 tool; T1 and T3 minutes each; T6 no refit | evaluator + scripts |
| 3 | **Calibrate the simulation:** peptide spread to the real within-protein k CV (0.21–0.23), σ_animal to ATF6 (0.013), proteome-wide animal shifts, and missingness from the completeness in §2 (the simulation has none today). Re-run A, B and C, then confirm the calibrated simulation reproduces the A–C gaps found in WP2. | calibrated arms | about 1.5–2 h per grid on 24 cores | `bench_rising_null.py` knobs |
| 4 | **Prototype D and E in the evaluator**, scored on the WP2 yardsticks that passed G1 and on the WP3 simulation. E: fixed sum-to-zero peptide slope offsets shared by the conditions, animal random slope, (1 − θ̂)² weights, on the `interaction_model` machinery; measure its cost per protein on ATF6. Because every in-vivo tissue contrast is within animal (finding 8), E's animal term is shared across an animal's tissues (the paired model), in the form WP1d picks; T3 on lve_atr and T6 on ATF6 score it. D on f: in the simulation only, plus the proprietary set locally. | candidates A–E on every yardstick → **G2** | WP2 + WP3 again | evaluator only |
| 5 | **Write-up:** evaluation report; design note for the winner, or for why nothing changes; manual §8.1 updated. Implementation proposal only if G2 passes. | report | — | docs |

EDL / SOL comes from the same mice as lve_atr, so it gives no null and does not fill the gap in
finding 5. It is still worth running in WP0: it is the second within-animal contrast (T3, T4),
it makes the four-tissue animal diagnostic possible (WP1d), and current Riana has never fitted
it.

## 5. Decision gates

**G1, after WP1–WP2: can real data discriminate?** Two conditions together:

- T1 on boomi, or T2 on ATF6, separates A from C in the direction the simulation predicts.
- C's null rate there is near nominal.

**Q1 verdict at G1:**

- **Deprecate `wls-var`** (results-affecting; the user decides) if Ṽ is uninformative or plateaus
  (WP1a) and B is no better calibrated than A on the yardsticks.
- **Keep it opt-in, with its limits documented,** if Ṽ tracks the scatter and B holds the nulls.

**G2, after WP4: build D, E or neither.** The manual's decision rule applies (among candidates that
hold the null, the most powerful; the simplest if tied), scored as follows:

| Criterion | Threshold (to confirm) |
|---|---|
| Null rate, pooled over T1/T2 replicates | ≤ 0.07 at α 0.05, with no k bin above 0.10 |
| T3 and T6 z SD | 0.9–1.15 (below 0.9 = SE too large, e.g. unmodelled pairing) |
| Calibrated simulation | false-positive rate ≤ 0.065 overall; CI coverage ≥ 0.93 |
| Power at ratio 1.5 (simulation) | not below A's by more than 10 % (relative) |
| E's cost | acceptable for a full ATF6 rollup |

If a candidate passes only in the simulation, it is written up but not built (the rule in §1).

## 6. Out of scope

Limma-style residual-variance moderation (after the winner); the paired within-animal model as
a separate option for today's estimators (it is part of E, WP4); the default flip to `linear`; guan / fornasiero removal; the `exponential`
rename; the k_p refits; GUI work. On the f route, D is checked in the simulation, and on real data
only locally.

## 7. Answers (user, 2026-11-13)

1. **EDL / SOL:** same design as lve_atr. All four tissues (ventricle, atrium, EDL, soleus) at a
   timepoint come from one mouse. The user is adding the data.
2. **Replicates:** boomi's `_1`/`_2` are separate culture plates, not reinjections, so T1 is a
   null that includes biological variation. boomi is the preferred dataset (D₂O, enough
   timepoints). lauren's `TR1`/`TR2` are not confirmed; lauren stays secondary.
3. **lve_atr:** both tissues come from the same mouse (as in 1).
4. **ATF6:** the user suggested the chamber contrast (LV vs RV or LV vs RA) in the control and KO
   animals as two animal sets, since the KO changes turnover very little. This is T6.
5. **Prototypes** of D and E in the evaluator count as evaluation.
6. **Placement:** this scope is committed under `reports/` and linked from manual §8.1 and the
   PROJECT_REVIEW §3 box.

Not yet confirmed: the G2 thresholds in §5 are proposals, to revisit at G1.

## Reproduce the audit

`python tests/benchmark/wlsvar_audit.py <fit_dir> ...` reads `riana_fit_fractions.txt` and
`riana_fit_peptides.txt` from each run and reports:

- the design;
- the gated completeness: per protein × condition with ≥ 2 admitted peptides, the share of
  (peptide, cell) pairs observed and the share of cells holding every admitted peptide.

The relabeled manifests for T1, T2 and T3 (and ATF6's animal IDs from its SDRF) come from
`python tests/benchmark/wlsvar_manifests.py {individuals,t1,t2,t3} ...`.

## 8. Progress log (the record behind the [evaluation report](2026-11-14_wls_var_D_E_evaluation.md))

**2026-11-13, WP0.** Refits at HEAD in `runs/wlsvar_eval/` (gitignored), each at its original
depth; the stored runs are untouched.

| Run | Refit time | Gated peptides (stored → HEAD) | k of peptides gated in both |
|---|---|---|---|
| boomi, depth 6 | 6.2 min at `-W 10` | 7,546 → 7,120 | unchanged (median log₂ ratio 0.000; 0.7 % move > 0.2) |
| boomi time-split halves, depth 3 | 8.0 / 9.5 min at `-W 8` | — | — |
| lve_atr, depth 6 | 9.5 min at `-W 8` (7.9 min at depth 3, set aside) | 10,024 → 9,273 | unchanged (0.000; 2.7 % > 0.2) |
| ATF6 (with animal IDs from the SDRF), depth 3 | 65 min at `-W 16` | 70,609 → 70,691 | unchanged; gated completeness identical |

**Stage-2 refit tool (WP0 decision: built, accepted).** A T2 full refit of ATF6 costs about an hour,
and the stage-2-only refit (`tests/benchmark/wlsvar_stage2.py`) costs about as much when it
mirrors `curve_fit`: the residual bootstrap, not the fs solve, dominates both. With `--fast`, the
200 bootstrap refits of a curve are one vectorized bounded Newton solve, with the random draws
kept in the shipped order. It needs a point store (a depth-1 fit: 80 min for ATF6). It also
writes the shipped fit's rows for unfitted curves: ≥ depth raw timepoints before rail-drop, fewer
after. Under admission `all`, those rows block a peptide in the other condition. Leaving them out
admitted 1,245 extra ATF6 curves. Acceptance, against shipped fits:

| Check | Curves / unfitted rows | k, CI, bands (max \|diff\|) | Rollup Δk, p (max \|diff\|); BH calls A / B / C / C+`wls-var` | Time |
|---|---|---|---|---|
| ATF6, unpermuted, vs the depth-3 refit | identical (190,141 / 46,374) | CI ≤ 1.7e-4 and bands ≤ 1.3e-4 at R² ≥ 0.6; k identical there (differences up to 0.02 only on R² < 0.6 curves with flat objectives, never admitted) | 1.6e-5, 7e-5; 35 / 64 / 0 / 2 on both | 7.2 min at `-W 14` |
| boomi T1 assignment 1 (a relabel), vs its full refit | identical (17,791 / 7,541) | CI ≤ 3e-5, bands ≤ 1.5e-5 | 6e-7, 4e-5; 0 on both | 0.7 min |

T2 on ATF6 therefore runs through `--fast` (about 7.5 min per permutation instead of about an hour).

**First T3 reading, boomi per-protein k (one condition):** confounded by curve shape, so on
boomi T1 is the SE check. k from the odd timepoints (1, 3, 6, 12 h) is 12 % above k from the
even ones (2, 4, 8, 24 h): median log₂ ratio +0.17, largest for slow proteins. Around that
shift, z has a robust SD of 1.68 under both `wls` and `wls-var`. The cause is that −φ/t is not
constant: on the median curve it rises from 0.017 (1 h) to 0.040 (12 h), then falls to 0.034 at
24 h. A through-origin slope is dominated by its last timepoint (weight t²), so the 12 h-ended
half reads faster than the 24 h-ended half; log₂(0.040 / 0.034) ≈ 0.23. The per-protein spread
of that shape inflates z's SD whatever the SEs are. A shape shared by two conditions largely
cancels from Δk, so T3 stays planned for lve_atr's Δk.

Side observations, not pursued here (the curve shape is a modelling question outside this
scope):
- With an intercept, boomi's collapsed points show an early lag: median −0.024 in −φ, 82 % of
  proteins below 0. The through-origin k is 7 % below the with-intercept k.
- The fall at 24 h is unexplained. One candidate is slower culture growth by then; apparent k
  in culture includes dilution by growth.

**First T3 on Δk, lve_atr (atrium vs ventricle; halves 6 / 5 timepoints, depth 3; 315 proteins
with a Δk in both halves; `tests/benchmark/wlsvar_eval.py t3`):** the time split separates A from
C in the direction the simulation predicts.

| Variant | z robust SD | share \|z\| > 1.96 | robust SD of Δk_odd − Δk_even | median combined SE |
|---|---|---|---|---|
| A (shipped) | 1.45 | 0.20 | 0.0102 | 0.0072 |
| B (`wls-var`) | 1.50 | 0.24 | 0.0111 | 0.0070 |
| C (shared proportions) | 1.13 | 0.11 | 0.0093 | 0.0072 |
| C + `wls-var` | 1.05 | 0.12 | 0.0082 | 0.0071 |

- C does not widen the SE. The halves' Δk agree better, which is what removing the random part
  of mechanism 1 predicts: each half re-estimates its bands, so its mixing differs.
- `wls-var` hurts on A and helps on C. A first hint for Q1: Ṽ may carry precision information
  that the per-condition mixing currently swamps. One dataset; WP1 and T1 / T2 to confirm.
- The pairing (finding 8) predicted z SD below 1. Even A is above 1, so any pairing gain is
  smaller than the excess spread from mixing and curve shape.

**T6, ATF6 (`wlsvar_eval.py t6`; no refit).** Experiment and condition swap roles: in each genotype,
the chamber Δk (LV vs RV, LV vs RA) from the shipped linear rollup, on the two chambers alone
(admission `all` spans every condition of an experiment). z compares the control animals'
Δk with the KO animals'. The animal sets are disjoint and share the timepoints, so curve shape
cancels, unlike T3.

| Variant | LV–RV: z robust SD | \|z\| > 1.96 | robust SD of Δk_ctrl − Δk_KO | median combined SE | LV–RA: z robust SD | \|z\| > 1.96 | robust SD of Δk_ctrl − Δk_KO | median combined SE |
|---|---|---|---|---|---|---|---|---|
| A (shipped) | 0.80 | 0.047 | 0.0050 | 0.0070 | 0.90 | 0.068 | 0.0060 | 0.0083 |
| B (`wls-var`) | 0.90 | 0.077 | 0.0052 | 0.0070 | 0.94 | 0.101 | 0.0062 | 0.0082 |
| C (shared proportions) | 0.51 | 0.006 | 0.0031 | 0.0071 | 0.71 | 0.022 | 0.0048 | 0.0085 |
| C + `wls-var` | 0.54 | 0.016 | 0.0034 | 0.0072 | 0.74 | 0.040 | 0.0049 | 0.0084 |

637 proteins (LV–RV), 695 (LV–RA); per-day units.

- C cuts the disagreement between the two animal sets by 20–40 %, and the SE stays the same. Its z
  falls well below 1, as finding 8 predicts for a paired contrast whose SE still holds the
  animal variance that cancels in Δk. Under A, the mixing noise adds about as much real
  disagreement, so A's z near 1 is a coincidence of two errors.
- LV–RA has a median z of +0.25 (LV–RV: +0.02): the chamber Δk is smaller in KO animals. This
  fits a proportional global brake, which scales a Δk and shows only where the chambers
  differ (RA vs LV).
- `wls-var` adds disagreement under both collapses.

**T1, boomi (20 assignments; `wlsvar_eval.py null`).** Pseudo-conditions are the two culture
plates at each timepoint, so the true Δk is 0. About 654 proteins are tested per assignment, all
with ≥ 2 peptides. Rejections at α 0.05 are pooled over the assignments, with the range across
them in brackets.

| Variant | α 0.05 (range) | α 0.01 | By turnover quartile, slow → fast | SD of Δk / SE | BH < 0.05 calls per refit |
|---|---|---|---|---|---|
| A (shipped) | 0.059 (0.048–0.076) | 0.018 | 0.048 / 0.062 / 0.053 / 0.072 | 1.07 | 0.45 |
| B (`wls-var`) | 0.058 (0.042–0.070) | 0.019 | 0.039 / 0.057 / 0.056 / 0.081 | 1.06 | 0.40 |
| C | 0.028 (0.017–0.047) | 0.007 | 0.020 / 0.028 / 0.025 / 0.041 | 0.87 | 0.05 |
| C + `wls-var` | 0.026 (0.015–0.045) | 0.005 | 0.016 / 0.025 / 0.023 / 0.040 | 0.83 | 0.10 |

No turnover quartile goes above 0.10 under any variant. All four rise in the fastest quartile.

The simulation predicted C below A (manual §6.3: 0.133 vs 0.163), and the direction holds. A,
however, is near nominal and C is conservative. The same reading as T6: the two plates of a
timepoint share its deviation from the curve (boomi's curve shape, §8 above) and any harvest
effect. Those cancel in Δk but stay in the residual variance. G1's second condition (C's null
rate near nominal) therefore fails as written. The term C lacks is the shared (paired) one, not
a wider band. That points to C plus E's paired term (WP4) rather than C alone.

**WP1a, Q1 (`wlsvar_wp1.py q1`).** Within-protein slope of the log squared leverage-corrected
`wls` residual (θ scale) on log cell variance; 1 = the variance tracks the scatter.

| | boomi (A) | lve_atr A | lve_atr C |
|---|---|---|---|
| on log Ṽ (moderated, as `wls-var` uses it) | 0.73 ± 0.06 | 0.84 ± 0.04 | 0.88 ± 0.06 |
| on log V (raw) | 0.65 ± 0.05 | 0.78 ± 0.04 | 0.82 ± 0.05 |
| … within condition / between conditions | — | 0.87 / 0.78 | 0.84 / 0.89 |
| scatter on log peptides per cell | −0.70 | −0.89 | −0.78 |
| Ṽ on log peptides per cell | −0.82 | −1.15 | −0.97 |
| cells | 15,408 | 9,092 | 9,103 |

- Ṽ is informative: the "uninformative" branch of the Q1 verdict is out.
- Ṽ spreads more than the scatter (slope < 1), and the scatter falls more slowly with peptides
  than Ṽ (a floor). By deciles, the cells with the lowest Ṽ scatter 1.25–1.45× more than
  Ṽ says. Weighting on Ṽ therefore over-trusts the best cells. That is the form D's Ṽ + τ²
  would correct.
- Yet B is no better calibrated than A on any yardstick so far (T1, T3, T6). WP1b–c (below)
  sizes the floor: the peptide spread τ² is 2–4× Ṽ at the late points. Ṽ is a minor part of the
  real cell variance, which is why weighting on it gains nothing.

**WP1d, ATF6 (`wlsvar_wp1.py animal`; 4 chambers × 30 animals).** Residuals of the shipped `wls`
fit, centred within (protein, genotype, timepoint, chamber), so curve shape and the timepoint
mean drop out. What is left is each animal's deviation.

- **A strong animal term, shared by the chambers.** An animal's deviations correlate across
  its chambers at 0.35–0.41 per cell (φ scale; 18–26 thousand cells per pair). Its
  proteome-wide part (median relative deviation over proteins) correlates 0.89–0.95 between
  chambers. Removing that proteome-wide shift takes 16 % of the per-cell variance. The rest
  still correlates 0.19–0.26 across chambers: animal × protein effects, not only a global factor.
- **The proteome-wide shift behaves like a turnover shift, not an enrichment error.** Its SD
  across animals at days 7 / 10 / 14 is flat on the relative-k scale (0.077 / 0.076 / 0.079). An
  enrichment (RIA) error would make it grow about 1.4× over that range as θ̂ rises from 0.46
  to 0.70. On the relative-θ scale it falls 0.049 → 0.043, between the two models' predictions.
  So a per-animal turnover shift of about 8 % (SD), from 6 animals per timepoint. Days 1 and 3
  are noisier (0.10, 0.03).
- **Additive vs proportional across tissues: not decidable on ATF6.** The chambers' k differ
  too little (median ratios 0.84–1.15). Shared covariance is about constant on both scales.
  EDL / SOL vs heart, with larger k differences, is the test (as planned).
- For E: the paired term needs an animal × protein random effect shared by the tissues. A
  proteome-wide animal factor could also enter as a covariate, since it is estimated from
  every protein. It cancels in within-animal contrasts but not between animal groups. For
  WP3: the simulation's per-animal turnover term (manual §6, S5) has the right form. Its size
  should be set from the 0.077.

**WP1b–c, D's τ² (`wlsvar_wp1.py tau`).** Per collapsed cell with ≥ 2 peptides, Cochran's Q from the
peptides' band variances (the collapse's own weights), pooled by θ̄ bin (DerSimonian–Laird). The
implied relative peptide-k spread is τ_rel = τ / g(θ), with g = (1 − θ)(−ln(1 − θ)) for an
apparent-k offset per peptide.

| θ̄ bin (median) | boomi τ (τ_rel) | lve_atr τ (τ_rel) | ATF6 τ (τ_rel) |
|---|---|---|---|
| ≤ 0.05 | 0 (Q/df 0.89) | 0 (0.93) | 0 (0.82) |
| 0.075 | 0.023 (0.31) | 0.024 (0.33) | 0.014 (0.20) |
| 0.15 | 0.034 (0.26) | 0.034 (0.26) | 0.033 (0.24) |
| 0.25 | 0.038 (0.18) | 0.042 (0.20) | 0.042 (0.19) |
| 0.45 | 0.042 (0.13) | 0.051 (0.16) | 0.055 (0.17) |
| 0.65 | 0.049 (0.13) | 0.048 (0.13) | 0.062 (0.17) |
| 0.85 | 0.029 (0.10; 53 cells) | 0.050 (0.18) | 0.047 (0.16) |
| τ² / Ṽ at the late points (share of cells with τ² > Ṽ) | 1.9 (73 %) | 4.3 (86 %) | 2.1 (73 %) |

- **τ² is not negligible (WP1c).** At the late, slope-carrying points it is 2–4× the cell variance in
  all three sets, so D does not reduce to C plus a cell variance. It also explains WP1a's floor:
  most of a cell's scatter is peptide spread, which Ṽ omits.
- **Shape (WP1b): an apparent-k offset on ATF6.** τ_rel is flat at 0.16–0.19 from θ̄ 0.25 to 0.85. On
  lve_atr and boomi, τ is closer to flat in θ (about 0.04–0.05) over that range, so τ_rel sags in the
  middle (0.13). A θ-proportional offset fits none of the three (τ / θ falls 4–5× across the range).
  For D: g = the k-offset form with τ_rel ≈ 0.15–0.18 per experiment. For WP3: peptide k spread
  around 0.17, below the within-protein k CV of 0.21–0.23, which also holds noise.
- At θ̄ ≈ 0, Q/df is 0.8–0.9: the bands are slightly wide at baseline.

**T2, ATF6 (20 permutations through `wlsvar_stage2.py --fast`; `wlsvar_eval.py null` / `inter`).**
At each timepoint, genotype labels are shuffled between the six animals, and an animal's four chambers
move together. Every timepoint mixes the groups: two animals of one genotype with one of the other,
so no timepoint keeps the contrast intact. Per-chamber Δk, about 3,500 tests per permutation, all on ≥ 2
precursors. Pooled rejections at α 0.05 (± SE over the 20 permutations, which share proteins):

| Variant | α 0.05 | α 0.01 | By turnover quartile, slow → fast | SD of Δk / SE | BH < 0.05 calls per permutation | Real labels: α 0.05, BH calls |
|---|---|---|---|---|---|---|
| A (shipped) | 0.096 ± 0.013 | 0.032 | 0.108 / 0.097 / 0.083 / 0.098 | 1.07 | 19.7 | 0.121, 35 |
| B (`wls-var`) | 0.122 | 0.048 | 0.118 / 0.119 / 0.108 / 0.142 | 1.21 | 49.4 | —, 64 |
| C | 0.039 ± 0.009 | 0.007 | 0.039 / 0.041 / 0.032 / 0.046 | 0.79 | 0.4 | 0.035, 0 |
| C + `wls-var` | 0.049 | 0.013 | 0.042 / 0.047 / 0.043 / 0.066 | 0.85 | 3.2 | —, 2 |

- **Whole permutations light up together.** A rejects 0.048–0.297 per permutation, C 0.015–0.200.
  The rate tracks the permutation's proteome-wide shift (correlation 0.89 with |median z|), not leakage of
  the real KO effect (−0.36). The shift is the per-animal turnover factor of WP1d: the SD of a
  permutation's median log₂ k ratio is 0.059 (4 %).
  - In the 9 permutations with |median z| < 0.3: A 0.061 (0.048–0.083), C 0.019 (0.015–0.030), the
    same as T1 and T6.
  - In the other 11: A 0.126, C 0.056.
  - This is mechanism 2 (animal variation growing with time; manual §6.3), on real data, at proteome
    scale. C does not touch it.
- **The real KO effect is inside that spread.** Median log₂ k_KO / k_ctrl = −0.069 (KO about 4.7 %
  slower); 7 of 20 permutations shift the proteome at least as far. This agrees with the earlier
  finding that the global brake is not significant at the animal level. The real data's 35 calls under
  A compare with 19.7 per null permutation (mean; up to 0.30 of tests in one permutation).
- `wls-var` is anti-conservative on ATF6: 0.122 overall, 0.142 in the fastest quartile, 49 BH calls
  per null permutation.

Interaction (`rollup --interaction` machinery, `wls`; about 980 proteins per permutation):

| Collapse | Ratio scale (primary): α 0.05 (range) | BH calls per permutation | Additive | Mean Δk | Real labels: ratio-scale BH calls |
|---|---|---|---|---|---|
| A (shipped) | 0.216 (0.178–0.256) | 69.3 | 0.245, 95 calls | 0.044 | 79 |
| C | 0.065 (0.048–0.087) | 1.2 | 0.081, 2.4 calls | 0.019 | 2 |

- Under A, the null alone gives 69 interaction calls per permutation against 79 on the real labels.
  The shipped ATF6 interaction calls are almost all mechanism 1, as the manual's 86 → 2 suggested; this
  is the null that confirms it. (79 at HEAD vs 86 in the stored run: the HEAD gated set differs slightly.)
- Under C, the primary test is near nominal, and the real data's 2 calls match the null's 1.2.

**Status after WP1–WP2 (G1 and Q1, for the user).**

- **G1, first condition: met.** Real data separate A from C in the direction the simulation predicts,
  consistently: T1 (boomi), T2 (ATF6, Δk and interaction), T6 and T3.
- **G1, second condition (C's null rate near nominal): not met as written.** C is conservative
  wherever no global shift is present (0.019–0.028), and its z SD is below 1 on every paired or
  shared design (T1 0.87, T6 0.5–0.7, T2 0.79). Pairing, the term C does not have, explains this. On
  ATF6, the per-animal proteome-wide shift (mechanism 2) inflates whole permutations under every
  collapse, so C's pooled 0.039 averages very conservative permutations with very loose ones. Both
  are what E's animal term is for (an animal random slope, shared by the tissues). C alone is
  therefore not a finished answer. The yardsticks can grade D and E, so WP3–WP4 can proceed if the
  user agrees to read condition 2 that way.
- **Q1: the deprecation branch applies (results-affecting; the user decides).** Ṽ is informative but
  floored: peptide spread τ² is 2–4× Ṽ (WP1a–c). B is worse calibrated than A on the ATF6 null (0.122 vs
  0.096; 49 vs 20 BH calls per permutation) and no better on T1, T3 and T6.

**2026-11-14, user decisions at G1.**

1. G1's second condition is read as "C lacks the paired (animal) term": WP3 and WP4 go ahead.
2. `wls-var`: deprecate it, or replace it with C later. 1.2.0 is unreleased and no real experiment
   is known to use the flag. The change to the flag is folded into the WP5 outcome (removed, or
   repurposed for the winner), so the unreleased flag changes once, not twice.

EDL / SOL waits for one combined ATR / LVE / EDL / SOL search (user, 2026-11-14), which will
also allow further comparison designs.

**WP3, simulation calibrated to ATF6 (`tests/benchmark/wlsvar_sim.py`, `wlsvar_calib.py`).** The
manual's bench (`bench_rising_null.py`) is left as it is, so its numbers stay reproducible. The
new bench has the ATF6 design: two conditions × four tissues from every animal, three animals
per condition × day, days 1, 3, 7, 10, 14, RIA 0.046, instant precursor, fitted as ATF6 was
(θ route, depth 3, `deberneh_2025_rss`, Spep ≥ 8). New nuisance, each a knob:

- a per-peptide apparent-k offset, shared by conditions, animals and tissues (WP1b);
- a per-animal proteome-wide turnover shift shared by the animal's tissues, plus a small
  per-animal × tissue one (WP1d);
- animal × protein turnover, shared by the tissues and per tissue;
- missingness per peptide and run, heavier for low-intensity peptides (§2 finding 1).

`wlsvar_calib.py` computes the same statistics on a real and on a simulated fit, from the frame
the shipped rollup collapses. Six rounds of 300-protein datasets tuned the knobs:

| Statistic | ATF6 | Simulation (calibrated) |
|---|---|---|
| protein k, median (10th–90th) | 0.091 (0.038–0.194) | 0.087 (0.043–0.180) |
| admitted precursors per protein × condition, median / 90th | 4 / 13 | 4 / 10 |
| cells holding every admitted precursor (median) / groups fully complete | 0.58 / 4.8 % | 0.60 / 2.8 % |
| band σ at θ̄ ≤ 0.05 / 0.4–0.7 | 0.060 / 0.072 | 0.059 / 0.062 |
| Q/df at baseline | 0.81 | 0.72 |
| τ_rel at θ̄ 0.25 / 0.45 / 0.65 | 0.19 / 0.17 / 0.17 | 0.15 / 0.17 / 0.17 |
| τ² / Ṽ at the late points | 2.1 | 4.0 |
| cross-tissue correlation of the animal deviation (φ scale) | 0.39 | 0.32 |
| relative-k deviation: variance / cross-tissue covariance | 0.0235 / 0.0055 | 0.0228 / 0.0057 |
| proteome-wide shift: cross-tissue correlation / SD (last 3 days) | 0.92 / 0.077 | 0.94 / 0.066 |

Calibrated knobs: peptide k offset SD 0.17 (log); animal shift 0.085, animal × tissue shift 0.02;
animal × protein 0.06 shared by tissues and 0.05 per tissue; RIA CV 2 %; envelope noise 8.5 %
relative; missingness 8 % × intensity^−1. The proteome-wide shift statistics vary between
datasets (6 animals per day), so the shift SD was set from three rounds. Two gaps remain: the
late-point band is 15 % narrower than ATF6's (hence the larger τ² / Ṽ), and the animal deviation
correlates a little less across tissues.

**WP4 prototypes (`tests/benchmark/wlsvar_de.py`).** Evaluator code, scored with
`wlsvar_eval.py --variants` and `wlsvar_sim.py`:

- **D**: C's shared band levels, then per-peptide weights 1/(σᵢ² + τ_rel²·g(θ̄)²), with τ_rel per
  experiment from the moments of Cochran's Q, pooled over both conditions. The shipped collapse then
  runs on the inflated bands. `D` = shipped `wls` on those cells; `D+var` = shipped `wls-var` with
  Ṽ_RE.
- **E**: per (experiment, protein), on the peptides' own points: one slope per condition, fixed
  peptide slope offsets (sum-to-zero, shared by the conditions), a random slope per animal,
  weights (1 − θ̂)² from an OLS pass (`E+band`: also ÷ the peptide's shared band variance). Where
  the animals carry both conditions (tissues of one animal), the random slope is one per condition
  with a free covariance, so the animal's shared part cancels from Δk (the paired model). If that
  fit is degenerate, E falls back to one shared slope, then to no animal term. df: animals − 2
  (between), animals − 1 (paired). The estimand is the unweighted peptide-average k.
- **E, joint (`interaction_E`)**: the shipped interaction model at the peptide level, with peptide
  slope offsets per tissue and two variance components, the animal slope shared by its tissues
  and the cell (animal × tissue) slope. Tests as shipped. df: the containment rule at the cell
  level (mean Δk: animals − 1; additive and ratio: cells − animals − 2(S − 1), 84 on ATF6).
- **EJ+time**: E, joint grouped by labeling time instead, with a random slope per time (shared
  by every sample taken then, all tissues and both conditions) and three nested components:
  animal, cell and tissue × time (the part of the timepoint deviation the tissues do not share).
  Cross-sectional designs with ≥ 3 times; otherwise, or if the fit fails, as E, joint. df: mean
  Δk animals − times − 1 (24 on ATF6); additive and ratio cells − animals − (S − 1)·times −
  (S − 1) (72), the tissue × condition × time stratum pooled with the cells. BFGS stalls at its
  start with three components and five groups; Powell converges and agrees with a long
  Nelder–Mead run.

**T6 with D and E (ATF6; no refit).** Same table as above, now with D and the paired E:

| Variant | LV–RV: z robust SD | \|z\| > 1.96 | robust SD of Δk_ctrl − Δk_KO | median combined SE | LV–RA: z robust SD | \|z\| > 1.96 | robust SD of Δk_ctrl − Δk_KO | median combined SE |
|---|---|---|---|---|---|---|---|---|
| A (shipped) | 0.80 | 0.047 | 0.0050 | 0.0070 | 0.90 | 0.068 | 0.0060 | 0.0083 |
| C | 0.51 | 0.006 | 0.0031 | 0.0071 | 0.71 | 0.022 | 0.0048 | 0.0085 |
| D | 0.51 | 0.006 | 0.003 | 0.007 | 0.75 | 0.022 | 0.005 | 0.009 |
| E (paired) | 0.82 | 0.024 | 0.004 | 0.005 | 0.94 | 0.071 | 0.005 | 0.007 |
| E+band (paired) | 0.85 | 0.019 | 0.003 | 0.005 | 0.96 | 0.062 | 0.005 | 0.006 |

- **E is the first candidate whose SE fits a paired contrast.** Its animal term cancels from the
  chamber Δk, so the SE drops by 25–30 % against C, with the same disagreement between the two
  animal sets (E+band). z robust SD rises from C's 0.51–0.71 to 0.85–0.96.
- D behaves as C here: its τ² changes the weights, not the animal term that T6 is about.
- LV–RA keeps its median z of about +0.2 under every variant (the proportional KO brake, above).
  The 6–7 % of |z| > 1.96 under E includes it.
- What is still below 1 on LV–RV may be chamber-specific peptide behaviour. E puts it in the
  residual, and T6 cannot see it because it repeats in both animal sets (§3, T6 "cannot see").

**T1 with D and E (boomi, the same 20 assignments).** boomi has no animal IDs, so E's "animal" is
the culture plate at one timepoint (each pseudo-condition holds one plate per timepoint).

| Variant | α 0.05 (range) | α 0.01 | By turnover quartile, slow → fast | SD of Δk / SE | BH < 0.05 calls per refit |
|---|---|---|---|---|---|
| A (shipped) | 0.059 (0.048–0.076) | 0.018 | 0.048 / 0.062 / 0.053 / 0.072 | 1.07 | 0.45 |
| C | 0.028 (0.017–0.047) | 0.007 | 0.020 / 0.028 / 0.025 / 0.041 | 0.87 | 0.05 |
| D | 0.028 (0.012–0.044) | 0.006 | 0.018 / 0.025 / 0.028 / 0.039 | 0.87 | 0 |
| D+var | 0.014 | 0.003 | 0.008 / 0.011 / 0.012 / 0.027 | 0.74 | 0 |
| E | 0.009 (0.004–0.017) | 0.001 | 0.003 / 0.009 / 0.010 / 0.013 | 0.69 | 0.05 |
| E+band | 0.010 | 0.001 | 0.003 / 0.009 / 0.013 / 0.015 | 0.69 | 0.10 |
| E+time | 0.028 (0.015–0.043) | 0.006 | 0.033 / 0.026 / 0.027 / 0.027 | 1.14 | 0.50 |
| E+band+time | 0.033 (0.017–0.052) | 0.005 | 0.034 / 0.034 / 0.033 / 0.032 | 1.09 | 0.70 |

- **E is the most conservative candidate on T1, and the reason is the design.** The two plates of a
  timepoint share that timepoint's departure from the line (boomi's curve shape, §8 above), and E
  counts it as plate variance. That share cancels from the plate-vs-plate Δk.
- **E+time** adds a random slope per timepoint, shared by every sample taken then (both
  conditions), with the animal or plate nested in it (cross-sectional designs). On boomi the time
  term takes nearly all of that variance (σ_time ≈ 0.0055 vs σ_plate ≈ 0.0002 on assignment 1), the
  SE falls by a third, and the rejection rate rises to C's level. df is then plates − timepoints −
  1 = 7. A t₇ has SD 1.18, so an SD of Δk/SE of 1.09–1.14 is in line. It is flat across turnover
  quartiles, unlike A's fast-tail rise. On the calibrated simulation, which has no shared timepoint
  deviation, σ_time ≈ 0 and E+time gives the same answer as E (df 24 instead of 28).
- `D+var` is more conservative than D: Ṽ_RE's moderation shrinks the spread between cells.

**WP3, calibrated null: does the simulation reproduce WP2?** (`runs/wlsvar_sim/null`, 12 datasets of
300 proteins, true ratio 1, so every test is a false positive.) The calibration statistics,
averaged over the 12, match ATF6 more closely than the single tuning dataset: cross-tissue animal
correlation 0.36 (range 0.25–0.48; ATF6 0.39), relative-k covariance 0.0059 (0.0055), proteome-wide
shift SD 0.068 (0.041–0.093; ATF6 0.077).

| α 0.05 | Simulation, pooled (range over datasets) | ATF6 T2, pooled (range over permutations) |
|---|---|---|
| Δk, A | 0.175 (0.112–0.383) | 0.096 (0.048–0.297) |
| Δk, B | 0.190 | 0.122 |
| Δk, C | 0.096 (0.039–0.322) | 0.039 (0.015–0.200) |
| Δk in datasets / permutations with little shift (\|median z\| < 0.3), A / C | 0.112–0.144 / 0.039–0.050 | 0.061 / 0.019 |
| interaction, ratio scale, A / C | 0.307 / 0.084 | 0.216 / 0.065 |
| interaction, BH calls per dataset (≈ 290 proteins) / per permutation (≈ 980), A / C | 44 / 0.4 | 69 / 1.2 |

- **The ordering and the gaps reproduce:** B above A, C well below A, and the interaction false
  positives almost all mechanism 1. As in T2, whole datasets light up with their proteome-wide
  shift: dataset 4 (median z −1.6) rejects 0.38 under A and 0.32 under C.
- **The simulation runs hotter than ATF6 by about 0.05–0.08 under both A and C.** The leading
  candidate is a deviation the simulation lacks: the timepoint-shared departure from the line
  that boomi shows (T1 with E+time, above). It makes every collapse conservative on real data,
  and it is absent from the simulation. The calibration therefore ranks candidates, but its
  absolute null rates are pessimistic for ATF6.

**WP4 on the calibrated null (the same 12 datasets; `wlsvar_sim.py --score-only`, `--pooled`).**
Pooled over datasets, so each dataset's proteome-wide animal shift counts as the population null
does. Per-dataset z SD in brackets: below the pooled value wherever the shift is shared across
proteins.

| Variant | α 0.05 (range over datasets) | α 0.01 | By turnover quartile, slow → fast | 95 % CI coverage | z SD pooled (within) | BH calls per dataset | Seconds per dataset |
|---|---|---|---|---|---|---|---|
| A (shipped) | 0.175 (0.112–0.383) | 0.079 | 0.168 / 0.170 / 0.181 / 0.180 | 0.825 | 1.56 (1.36) | 54 | 35 |
| B (`wls-var`) | 0.190 (0.124–0.393) | 0.088 | 0.179 / 0.189 / 0.191 / 0.199 | 0.810 | 1.62 (1.41) | 63 | 33 |
| C | 0.096 (0.039–0.322) | 0.031 | 0.101 / 0.099 / 0.098 / 0.088 | 0.904 | 1.24 (1.01) | 9.5 | 10 |
| D | 0.091 (0.034–0.305) | 0.028 | 0.096 / 0.098 / 0.092 / 0.081 | 0.909 | 1.22 (1.00) | 7.9 | 43 |
| D+var | 0.082 (0.029–0.267) | 0.023 | 0.082 / 0.087 / 0.082 / 0.075 | 0.918 | 1.18 (0.95) | 3.8 | 21 |
| **E** | **0.059 (0.019–0.196)** | **0.014** | **0.060 / 0.061 / 0.062 / 0.054** | **0.941** | **1.08 (0.89)** | **1.1** | 16 |
| E+band | 0.067 (0.019–0.230) | 0.017 | 0.066 / 0.066 / 0.073 / 0.062 | 0.933 | 1.12 (0.90) | 3.4 | 41 |
| E+time | 0.067 (0.023–0.223) | 0.018 | 0.067 / 0.067 / 0.071 / 0.063 | 0.933 | 1.12 (0.92) | 2.3 | 37 |
| E+band+time | 0.076 (0.020–0.267) | 0.022 | 0.074 / 0.076 / 0.082 / 0.074 | 0.924 | 1.16 (0.93) | 5.3 | 34 |

Interaction (the shipped model on A's and C's cells; E's joint peptide-level model):

| α 0.05 (BH calls per dataset) | Ratio scale (primary) | Additive | Mean Δk |
|---|---|---|---|
| A | 0.307 (44) | 0.314 (45) | 0.116 (7.6) |
| C | 0.084 (0.4) | 0.089 (0.5) | 0.074 (1.3) |
| E, joint | 0.046 (0) | 0.053 (0.08) | 0.065 (0.08) |

- **E is the only candidate inside the proposed G2 simulation thresholds** (false-positive rate
  ≤ 0.065, coverage ≥ 0.93), for Δk (0.059, 0.941) and for the interaction (0.046–0.053). It holds
  even with each dataset's own proteome-wide animal shift, because its per-protein animal slope
  puts that variance in the SE. C and D remove mechanism 1 but not mechanism 2 (0.09–0.10).
- Whole datasets still light up under E when their shift is large (0.196 at median z −1.6). That
  is the right marginal behaviour for a per-protein test, but the false positives come in
  clusters, which BH does not anticipate. A shift that large is a real difference between the two
  groups of animals, not a protein-level effect (evaluation report §9, item 1).
- The simulation has no shared timepoint deviation, so E+time's extra term only adds noise there
  (0.067). The band weights (E+band) add nothing here either.

**WP3 / WP4 power at ratio 1.5 (`runs/wlsvar_sim/r15`, 4 datasets; every B protein 1.5× faster,
on every tissue).** Too strong to rank power: every variant rejects 0.994–0.998 and calls
1125–1129 of 1131 at BH < 0.05. What the run does show is how each variant estimates a real
effect, against the simulation's truth (median over tests, mean over datasets):

| Variant | log k_A − log k_A,true | log k_B − log k_B,true | Δk / Δk_true | log₂ ratio (true 0.585) | Coverage vs k0 | Coverage vs the protein's estimand | Coverage under the null |
|---|---|---|---|---|---|---|---|
| A | +0.022 | +0.032 | 1.047 | 0.601 | 0.780 | 0.827 | 0.825 |
| C | +0.029 | +0.029 | 1.031 | 0.592 | 0.855 | 0.909 | 0.903 |
| E | +0.063 | +0.039 | 0.996 | 0.552 | 0.893 | 0.928 | 0.940 |
| Em (below) | +0.053 | +0.051 | 1.051 | 0.584 | 0.872 | 0.929 | 0.927 |

- **"Coverage vs the protein's estimand" is the fair column.** The simulation sets each peptide's
  k to k0·e^koff and mislabels its sites by Spep_true / Spep_tab. A protein-level estimate
  targets the average over its peptides, not k0, so its true Δk is (r − 1)·k0·mean(e^koff ·
  Spep_true / Spep_tab). That factor is 1 under the null and differs from protein to protein
  under the alternative. Scored against k0, every variant loses 3–5 points of coverage at ratio
  1.5. Against the estimand, each is back at its null coverage. The within-dataset spread is what
  grows against k0 (E: z SD 0.89 → 1.23, back to 1.07 against the estimand), not the
  between-dataset shift. The proteins' estimands are regenerated from the simulation seed
  (`estimand_factors` in `wlsvar_sim.py`; `--pooled` reports both coverages).
- **E's k level is about 4–6 % high in both conditions**, against 1–3 % for A and C. That holds
  under the null too (+4.9 % / +4.2 %): E's estimand is the unweighted peptide average, whereas
  the collapse precision-weights the peptides. It cancels in Δk but not in the ratio when the two
  conditions' biases differ. At 1.5 they do (+6.3 % vs +3.9 %), so E's log₂ ratio is 6 % low.
- **Em (rejected):** E with multiplicative peptide offsets, slope k_c·(1 + m_p), with the offset
  columns rescaled by k̂_c / k̄ over four WLS passes, then the same MixedLM. It removes the
  ratio bias (0.584) and the asymmetry, but coverage is no better, and its null false-positive
  rate rises from 0.059 to 0.072 (coverage 0.928), with 2.8 BH calls per dataset instead of 1.1.
  The likely cause is a design that depends on the estimate, ignored by the SE. E stays
  additive. The ratio bias is small at realistic effects (below).
- **Interaction at ratio 1.5.** The ratio-scale test is still a true null (the same ratio on every
  tissue): E joint 0.040, C 0.088. The additive test is now a real effect, because atria are 15–17 %
  faster, so their Δk is larger: power C 0.285, E joint 0.161. C's additive test rejects 0.089
  under the full null, so its power figure is not comparable.

**WP3 / WP4 power at ratio 1.2 (`runs/wlsvar_sim/r12`, 4 datasets).** The size of the real tissue
effects (LA / RA 1.15–1.17 here, lve_atr's atrium 1.24). "Size-adjusted" uses each variant's own
p threshold that rejects 5 % on the 12 null datasets, pooled (so the population null, shifts
included).

| Variant | Null α 0.05 | Power at α 0.05 (range over datasets) | vs A | p threshold for a 5 % null | Size-adjusted power | vs A | BH calls per dataset (of ≈ 1120) | Coverage vs estimand | log₂ ratio (true 0.263) |
|---|---|---|---|---|---|---|---|---|---|
| A (shipped) | 0.175 | 0.760 (0.627–0.914) | 1 | 0.0038 | 0.504 | 1 | 813 | 0.823 | 0.260 |
| B (`wls-var`) | 0.190 | 0.779 | 1.03 | 0.0028 | 0.500 | 0.99 | 837 | — | 0.262 |
| C | 0.096 | 0.780 (0.663–0.923) | 1.03 | 0.020 | 0.679 | 1.35 | 829 | 0.912 | 0.254 |
| C + `wls-var` | 0.096 | 0.795 | 1.05 | 0.020 | 0.691 | 1.37 | 849 | — | 0.254 |
| D | 0.091 | 0.779 | 1.02 | 0.021 | 0.675 | 1.34 | 824 | 0.918 | 0.254 |
| D+var | 0.082 | 0.785 | 1.03 | 0.026 | 0.704 | 1.40 | 824 | — | 0.255 |
| **E** | **0.059** | **0.721 (0.539–0.891)** | **0.95** | 0.041 | **0.689** | **1.37** | 713 | **0.946** | 0.238 |
| E+band | 0.067 | 0.777 (0.615–0.921) | 1.02 | 0.035 | 0.726 | 1.44 | 802 | 0.934 | 0.240 |
| E+time | 0.067 | 0.734 | 0.97 | 0.034 | 0.678 | 1.35 | 732 | 0.936 | 0.237 |
| E+band+time | 0.076 | 0.788 | 1.04 | 0.030 | 0.723 | 1.44 | 822 | — | 0.238 |

Interaction at ratio 1.2 (every tissue's B 1.2× faster: the ratio-scale test is a true null; the
additive difference is real but small, 0.15–0.17 × a 20 % Δk):

| α 0.05 (BH calls per dataset) | Ratio scale | Additive | Mean Δk |
|---|---|---|---|
| A | 0.325 (47) | 0.354 (61) | 0.916 (255) |
| C | 0.099 (2.0) | 0.128 (3.5) | 0.920 (258) |
| E, joint | 0.058 (0) | 0.081 (0) | 0.901 (252) |

- **E costs about 5 % of A's nominal power**, inside the G2 bound of 10 %. **At equal size it has 37 %
  more** (0.689 vs 0.504). A's nominal power is bought with its false-positive rate of 0.175.
  Among the candidates that hold their size, the band weights help: E+band 0.726 at the same
  threshold, but its null rate (0.067) sits just above the 0.065 bound.
- **E keeps its coverage under a real effect** (0.946 against the protein's estimand), and its
  ratio estimate is 9 % low in log₂ (0.238 vs 0.263), the same k-level asymmetry as at 1.5.
- **The interaction holds its size under a real main effect only with E's joint model** (ratio
  scale 0.058; C 0.099). Its mean-Δk power is within 2 % of A and C.
- The spread of power between datasets (E 0.54–0.89) is the proteome-wide animal shift. It
  moves every variant together.

**E's k level (synthetic check, 100–150 proteins of 4 peptides, 3 animals per condition × 5 days,
peptide k offsets SD 0.17, constant θ noise).** E recovers its estimand, the peptide-average k,
exactly without noise. Two small effects add up to the 1–3 % seen in the simulation:

| θ noise SD | Ratio | log k_A − estimand | log k_B − estimand |
|---|---|---|---|
| 0 | 1 | 0.000 | 0.000 |
| 0 | 1.5 | +0.005 | −0.006 |
| 0.03 | 1.5 | +0.010 | −0.006 |
| 0.065 (the simulation's band σ) | 1 | +0.011 | +0.014 |
| 0.065 | 1.5 | +0.025 | +0.006 |

- **Additive offsets on a multiplicative truth** cost ±0.5 % at ratio 1.5, about 1 % of ratio.
- **Noise passed through φ = ln(1 − θ) point by point** (Jensen) adds 1–2.5 %, more for the
  slower condition. The collapse averages θ over peptides before the log, so it shows less of
  this.
- Net effect: E's ratio is 1–2 % low at strong effects, and Δk is unaffected to within 0.15 SE. A
  second-order correction, φ + σ²/(2(1 − θ)²) from the band, would remove most of the Jensen
  part if E's k values are ever reported. Not needed for the Δk test.

**T2 with D and E (ATF6, the same 20 permutations; `wlsvar_eval.py null --variants D E E+time
E+band --cache`).** Real data, no truth needed. A, B and C are repeated from the WP2 table.

| Variant | α 0.05 pooled (range over permutations) | α 0.01 | By turnover quartile, slow → fast | SD of Δk / SE | BH < 0.05 calls per permutation |
|---|---|---|---|---|---|
| A (shipped) | 0.096 (0.048–0.297) | 0.032 | 0.108 / 0.097 / 0.083 / 0.098 | 1.07 | 19.7 |
| B (`wls-var`) | 0.122 (0.073–0.322) | 0.048 | 0.118 / 0.119 / 0.108 / 0.142 | 1.21 | 49.4 |
| C | 0.039 (0.015–0.200) | 0.007 | 0.039 / 0.041 / 0.032 / 0.046 | 0.79 | 0.40 |
| D | 0.038 (0.013–0.199) | 0.007 | 0.039 / 0.038 / 0.030 / 0.044 | 0.78 | 0.30 |
| E | 0.018 (0.007–0.087) | 0.003 | 0.012 / 0.017 / 0.018 / 0.025 | 0.70 | 0.35 |
| E+band | 0.020 (0.007–0.093) | 0.003 | 0.011 / 0.020 / 0.020 / 0.027 | 0.69 | 0.25 |
| E+time | 0.041 (0.011–0.225) | 0.007 | 0.046 / 0.044 / 0.035 / 0.039 | 0.82 | 0.40 |

- **Every candidate except A and B passes the G2 real-data null** (≤ 0.07 pooled, no turnover
  quartile above 0.10).
- **E is conservative on ATF6, as on boomi (T1: 0.009),** and for the same reason. Both conditions
  sampled at a timepoint share that timepoint's departure from the line, which cancels from Δk.
  E counts it as animal variance. **E+time takes it out and lands at 0.041,** flat across
  turnover quartiles. On boomi it gave 0.028.
- **E's worst permutation is 0.087, against 0.20–0.30 for the collapse variants.** The animal
  variance it carries per protein spans the proteome-wide shift that lights up whole
  permutations. E+time gives some of that back (0.225): with the timepoint variance removed, the
  animal term is smaller and the shift-driven permutations reappear. On real data, E and E+time
  trade the same quantity: robustness to whole-group shifts against conservatism.
- The simulation lacks the timepoint-shared deviation, so it ranks E+time below E (0.067 vs
  0.059). The two real datasets where it can be checked (boomi T1, ATF6 T2) both show the
  deviation and favour E+time for size. The simulation's 0.05–0.08 excess over ATF6 under A and C
  has the same cause.

**Correction applied (E, joint).** The prototype's tissue filter let in a tissue that had
animals under only one condition: pandas' `min` skips the missing condition, which left an
unidentified contrast and degenerate p-values. Found on real ATF6, two proteins at p ~ 1e-194
and 1e-88; the shipped `_prepare` is correct (it requires both conditions). The fix also takes
the fixed-effect covariance from GLS at the fitted variance components, (Σ X′V⁻¹X)⁻¹, rather than
the joint Hessian inverse, which was ill-conditioned near a zero component (24 negative
variances on ATF6). The E, joint rows above (null, ratio 1.5, ratio 1.2) are re-scored with the
fix; no rate moved by more than 0.003.

**E, joint on real ATF6, unpermuted (fixed code; `wlsvar_eval.py inter runs/wlsvar_eval/atf6/fit
--variants A C EJ --cache`).** About 976 proteins carry ≥ 2 chambers with both genotypes.

| Interaction, α 0.05 (BH < 0.05 calls) | Ratio scale (primary) | Additive | Mean Δk |
|---|---|---|---|
| A (shipped model on per-condition cells) | 0.227 (79) | 0.294 (135) | 0.063 (1) |
| C (shipped model on shared-band cells) | 0.056 (2) | 0.114 (3) | 0.016 (0) |
| E, joint (peptide level) | 0.028 (0) | 0.050 (0) | 0.015 (0) |
| EJ+time | 0.029 (0) | 0.067 (0) | 0.051 (0) |

- **The two interaction hits named in the ATF6 biology notebook do not survive the peptide-level model.**
  Smyd1 (P97443): A p ≈ 0, C 0.28, E joint 0.77. That is the peptide-mixing artifact
  (mechanism 1), and C already removes it. Cops3 (O88543): A p ≈ 0, C p ≈ 0 (BH 0.002), E joint 0.22
  (rank 160 of 976). This fits the known outlying LA-control cell: E joint's cell (animal ×
  chamber) variance carries it, whereas the shipped model has only the animal slope. C's other
  call, Q9JHU4, is also absent under E joint.
- The real-data rejection rates under E joint sit at or below nominal on every test. That is
  consistent with no chamber-specific KO effect beyond the proportional brake (§8, T6),
  although a permutation null (T2, below) is needed to tell "none" from "conservative".
- **EJ+time** fits its time term on 965 of 979 proteins (median σ_time 0.0117, σ_animal 0.0092,
  tissue × time 0.0017, cell ≈ 0). The interaction tests barely move (ratio-scale p vs E,
  joint: Spearman 0.96; Cops3 0.19, Smyd1 0.78). The mean-Δk test comes back to nominal
  (0.015 → 0.051): the timepoint deviation no longer sits in the animal term. 75 % of proteins
  have KO slower (the brake), but no single protein passes BH.

**WP3 with the timepoint-shared deviation.** The simulation now carries the deviation that boomi
and ATF6 show (knobs in `wlsvar_sim.py`): per protein × tissue × day, shared by every animal
sampled that day in both conditions. Calibrated to ATF6's E+time fit, where σ_time by turnover
quartile (k 0.036 / 0.069 / 0.112 / 0.184) is 0.0079 / 0.0100 / 0.0128 / 0.0177, about
√(0.0072² + (0.088·k)²), with σ_animal / k ≈ 0.11:

- A turnover factor alone (`time_k_sd`) gives a σ_time / k that is flat across quartiles, at any
  size (three rounds: 0.08–0.13). ATF6's falls with k.
- A constant θ offset per peptide (a reference-envelope error, `pep_theta_sd`) adds nothing
  E+time can see.
- A factor of 0.09 plus an absolute offset per day (`time_kabs_sd` 0.0115) reproduces the
  shape. On the 8 null datasets: k 0.056 / 0.083 / 0.116 / 0.178 → σ_time 0.0096 / 0.0116 /
  0.0130 / 0.0173. ATF6's curve gives 0.0087 / 0.0103 / 0.0125 / 0.0172 at those k, so the match
  is within about 10 %, a little high for slow proteins. σ_animal / k stays at 0.10–0.11.

In this first version the deviation is drawn independently per tissue (`runs/wlsvar_sim/null_t`,
8 datasets; `r12_t`, 4; same seeds as `null` and `r12`).

| Δk, α 0.05 pooled (range) | Simulation, no timepoint term | Simulation with it | ATF6 T2 |
|---|---|---|---|
| A (shipped) | 0.175 (0.112–0.383) | **0.096** (0.052–0.213) | **0.096** (0.048–0.297) |
| B (`wls-var`) | 0.190 | 0.105 | 0.122 |
| C | 0.096 | 0.046 (0.013–0.158) | 0.039 (0.015–0.200) |
| D | 0.091 | 0.044 | 0.038 |
| D+var | 0.082 | 0.042 | — |
| E | 0.059 | 0.025 (0.008–0.078) | 0.018 (0.007–0.087) |
| E+band | 0.067 | 0.026 | 0.020 |
| E+time | 0.067 | 0.063 (0.021–0.203) | 0.041 (0.011–0.225) |
| E+band+time | 0.076 | 0.070 | — |

- **The timepoint term closes the gap between simulation and real data.** With it, the
  simulated null reproduces ATF6's permutation null for every variant measured on both: A
  exactly, the rest within 0.01, apart from E+time, which runs 0.02 hotter in the simulation.
  The worst datasets match the worst permutations too. The 0.05–0.08 excess flagged above was
  this deviation.
- **It makes E conservative (0.025), as on real data, and E+time absorbs it (0.063).** C and D sit
  at 0.044–0.046, coverage 0.95.

Power at ratio 1.2 with the timepoint term (size-adjusted on `null_t`, as before):

| Variant | Null α 0.05 | Power at α 0.05 (range) | vs A | Size-adjusted power | vs A | Coverage vs estimand | BH calls per dataset |
|---|---|---|---|---|---|---|---|
| A (shipped) | 0.096 | 0.640 (0.486–0.849) | 1 | 0.508 | 1 | 0.903 | 565 |
| B (`wls-var`) | 0.105 | 0.654 | 1.02 | 0.496 | 0.98 | 0.889 | 596 |
| C | 0.046 | 0.625 (0.447–0.853) | 0.98 | 0.636 | 1.25 | 0.954 | 529 |
| D | 0.044 | 0.623 | 0.97 | 0.637 | 1.25 | 0.954 | 515 |
| D+var | 0.042 | 0.626 | 0.98 | 0.653 | 1.29 | 0.959 | 523 |
| E | 0.025 | 0.546 (0.317–0.783) | **0.85** | 0.676 | 1.33 | 0.970 | 374 |
| E+band | 0.026 | 0.591 | 0.92 | 0.701 | 1.38 | 0.968 | 456 |
| **E+time** | **0.063** | **0.727 (0.523–0.903)** | **1.14** | 0.688 | 1.35 | **0.938** | 675 |
| E+band+time | 0.070 | 0.770 | 1.20 | 0.721 | 1.42 | 0.926 | 748 |

- **Under realistic timepoint deviation, E at its nominal α falls outside the G2 power bound**
  (0.85 × A), because it is conservative. At equal size it still beats every collapse (1.33 ×
  A against C's 1.25 ×).
- **E+time is the one variant inside every simulation threshold here.** False-positive rate
  0.063 (≤ 0.065), coverage 0.938 (≥ 0.93), and nominal power 1.14 × A, where A runs at a null of
  0.096. E+band+time has more power but breaks the size bound (0.070).
- **E+time keeps its power when the term is added (0.734 → 0.727). Everything else loses it**
  (A 0.760 → 0.640, C 0.780 → 0.625, E 0.721 → 0.546). The deviation largely cancels from the Δk
  estimate, because both conditions share it, but it inflates every SE except E+time's.

**Is the timepoint deviation shared across tissues?** The first version drew it independently
per tissue, and its interaction tests came out far colder than ATF6's (ratio scale: A 0.091,
C 0.011, E, joint 0.006). On ATF6, it is mostly shared. Per protein × chamber × day, the mean
relative-k residual over the day's animals (both genotypes), minus the animals' own
cross-chamber covariance ÷ n, correlates across chambers at **0.85** (pairs 0.80–0.92). The same
statistic reads 0.10 on the independent simulation, its floor. (`day_shared_frac` in
`wlsvar_calib.py`.) A shared deviation lands in the animal slope, shared by an
animal's tissues, which cancels from the interaction. A per-tissue one lands in the cell
variance, which does not.

`time_shared` 0.85 (`runs/wlsvar_sim/null_ts`, `r12_ts`; same seeds) reads 0.80 on the same
statistic. Δk repeats the independent version within 0.003 on every variant (A 0.094, C 0.047,
D 0.046, E 0.026, E+time 0.064, E+band+time 0.070; coverage E+time 0.936). Size-adjusted power
at ratio 1.2: A 0.486, C 0.629, D 0.629, D+var 0.657, E 0.676, E+time 0.690, E+band 0.701,
E+band+time 0.710. E+time's nominal power is 1.16 × A (0.733 vs 0.630), E's 0.87 ×. The
interaction now matches ATF6:

| Interaction, α 0.05 | Simulation, shared deviation (8 datasets) | ATF6 T2 (20 / 10 permutations) |
|---|---|---|
| Ratio scale: A / C / E, joint | 0.194 / 0.044 / 0.022 | 0.216 / 0.065 / 0.025 |
| Additive: A / C / E, joint | 0.205 / 0.052 / 0.026 | 0.245 / 0.081 / 0.035 |
| Mean Δk: A / C / E, joint | 0.045 / 0.027 / 0.018 | 0.044 / 0.019 / 0.010 |
| BH calls, ratio scale: A / C / E, joint | 17.8 / 0.1 / 0 per ≈ 280 proteins | 69 / 1.2 / 0.1 per ≈ 980 |

- **The calibrated simulation now reproduces ATF6's permutation nulls for Δk and for all three
  interaction tests**, under every variant scored on both, to within 0.01–0.03.
- **E, joint is conservative on the interaction in both** (0.022–0.035), unlike E on Δk (whose
  conservatism E+time removes). Its additive and ratio tests sit in the cell stratum, where
  only the 15 % per-tissue share of the timepoint deviation lands. Its mean-Δk test is
  conservative for E's reason: the shared deviation sits in its animal term.
- Power at ratio 1.2 with the shared deviation, mean Δk across tissues: A 0.741, C 0.736,
  E, joint 0.725, at nulls of 0.045 / 0.027 / 0.018. The ratio-scale test (a true null at 1.2):
  A 0.204, C 0.062, E, joint 0.033.

**EJ+time on the interaction** (the same datasets and permutations; `wlsvar_sim.py --score-only
--variants EJ+time`, `wlsvar_eval.py inter --variants EJ+time`; size-adjusted power from the
pooled `null_ts` thresholds, `wlsvar_sim.py --work runs/wlsvar_sim/r12_ts --ratio 1.2
--size-adjust runs/wlsvar_sim/null_ts`, which prints the Δk rows too):

| α 0.05 | Ratio scale (primary) | Additive | Mean Δk |
|---|---|---|---|
| ATF6 T2, 10 permutations: E, joint / EJ+time | 0.025 / 0.027 | 0.035 / 0.041 | 0.010 / 0.033 (worst 0.025 / 0.103) |
| Simulation `null_ts`, pooled: E, joint / EJ+time | 0.022 / 0.039 | 0.026 / 0.049 | 0.018 / 0.074 |
| … median over the 8 datasets: E, joint / EJ+time | 0.024 / 0.042 | 0.024 / 0.048 | 0.007 / 0.027 |
| Simulation `null` (no timepoint term): E, joint / EJ+time | 0.046 / 0.049 | 0.053 / 0.056 | 0.065 / 0.082 |
| `r12_ts` at α 0.05: A / C / E, joint / EJ+time | 0.204 / 0.062 / 0.033 / 0.049 (a true null) | 0.241 / 0.089 / 0.045 / 0.083 | 0.741 / 0.736 / 0.725 / **0.919** |
| `r12_ts`, size-adjusted: A / C / E, joint / EJ+time | — | 0.055 / 0.086 / 0.085 / 0.084 | 0.782 / 0.850 / 0.880 / **0.886** |

- **The ratio-scale and additive tests:** EJ+time brings E, joint's simulation rates to nominal
  (0.039–0.049), holds the ratio-scale null under a real main effect (0.049 at ratio 1.2), and
  has the same size-adjusted power as E, joint and C. On ATF6 it is as conservative as E, joint
  (0.027 / 0.041). Those tests sit in the cell stratum, which the shared deviation never reached.
- **The mean-Δk test** (protein Δk pooled over tissues) gains what E+time gains on Δk: nominal
  power 0.919 against E, joint's 0.725 (1.24 × A), and the best size-adjusted power (0.886).
  Its pooled simulation null, 0.074, is the mechanism-2 pattern at its strongest. Six of the 8
  datasets reject 0.004–0.030; the two with the largest proteome-wide shift (median z −1.18 and
  +0.90) reject 0.301 and 0.172. Averaging over four tissues shrinks the protein-level noise
  but not the animal group's shift, which every tissue shares. With 8 datasets the pooled rate
  has an SE of 0.038, so it does not separate from 0.05 or from 0.065; E+time's per-tissue Δk
  (0.064, SE 0.025) and C's (0.047, SE 0.019) are dominated by the same two datasets. ATF6's
  permutations show it more mildly (0.033, worst 0.103).
- Timing (full ATF6, 979 proteins, `-W 14`, otherwise idle): E, joint 13.6 min, EJ+time 18.2 min.
  Both reproduce the cached rejection rates exactly.

**T3 and T6 with the time term.** T6 (ATF6, chamber Δk in the control vs the KO animals):

| Variant | LV–RV: z robust SD | \|z\| > 1.96 | LV–RA: z robust SD | \|z\| > 1.96 |
|---|---|---|---|---|
| A (shipped) | 0.80 | 0.047 | 0.90 | 0.068 |
| C | 0.51 | 0.006 | 0.71 | 0.022 |
| E (paired) | 0.82 | 0.024 | 0.94 | 0.071 |
| E+band | 0.85 | 0.019 | 0.96 | 0.062 |
| E+time | 0.82 | 0.030 | 0.98 | 0.083 |
| E+band+time | 0.85 | 0.024 | 1.07 | 0.071 |

T3 (lve_atr, atrium vs ventricle Δk from the odd vs the even timepoints; one mouse per
timepoint, both tissues from it; 315 proteins):

| Variant | z robust SD | z SD | \|z\| > 1.96 | robust SD of Δk_odd − Δk_even | median combined SE |
|---|---|---|---|---|---|
| A (shipped) | 1.45 | 1.74 | 0.200 | 0.010 | 0.007 |
| C | 1.13 | 1.29 | 0.108 | 0.009 | 0.007 |
| E | 0.75 | 0.83 | 0.025 | 0.007 | 0.009 |
| E+band | 0.88 | 0.95 | 0.054 | 0.006 | 0.007 |
| **E+time** | **0.94** | **1.01** | **0.060** | 0.007 | 0.007 |
| E+band+time | 1.07 | 1.14 | 0.083 | 0.006 | 0.006 |

- **On T3, E+time and E+band+time land inside 0.9–1.15; E alone is too wide (0.75).** The
  lve_atr fit has no animal IDs, so E's "animal" is each tissue sample (condition × timepoint),
  and E runs unpaired (df 10 per half). The mouse's deviation, shared by its two tissues,
  then sits in E's animal variance although it cancels from the tissue Δk: boomi's T1
  situation. With one mouse per timepoint, E+time's time slope is the mouse slope, so it
  takes that shared part out (df 5) and the SE falls by a quarter. The halves' disagreement
  is the same under every E variant (0.006–0.007); E+time's SE matches it. *Corrected
  2026-11-14:* this bullet first called E paired here; it was not (its `paired` flag is 0 on
  every lve_atr test).
- **On T6, the time term changes nothing** (it cancels from a paired chamber contrast). LV–RA
  passes (0.98). LV–RV stays at 0.82–0.85 under every E variant, as under A (0.80), and no
  candidate reaches 0.9. LV and RV share their turnover, so this contrast carries no brake,
  and what remains is the shared-peptide gap noted at the first T6: chamber-specific peptide
  behaviour sits in E's residual and in its SE, but repeats in both animal sets, so it never
  reaches the T6 difference. Whether it belongs in the SE depends on the estimand (these
  peptides, or the protein's peptides at large); it is the same for every E variant.

**Real ATF6, unpermuted (per-chamber KO vs control Δk, 3,507 tests).**

| Variant | α 0.05 | BH < 0.05 calls | Its T2 permutation null (α 0.05 / calls per permutation) |
|---|---|---|---|
| A (shipped) | 0.121 | 35 | 0.096 / 19.7 |
| B (`wls-var`) | 0.156 | 64 | 0.122 / 49.4 |
| C | 0.035 | 0 | 0.039 / 0.40 |
| D | 0.036 | 0 | 0.038 / 0.30 |
| E | 0.014 | 0 | 0.018 / 0.35 |
| E+time | 0.043 | 0 | 0.041 / 0.40 |

- **Every candidate that holds the null sees ATF6's per-chamber KO effect at its null rate, with
  no protein passing BH.** A's 35 calls are within the range A gives on permuted labels (its
  worst permutation rejects 0.297). The proteome-wide brake (75 % of proteins slower in KO, EJ+time's cross-chamber mean Δk) is a
  shift of the whole animal group: per protein it is the mechanism-2 pattern, which a
  per-protein test should not call (evaluation report §9, item 1: a proteome-level test of the
  brake).

**T6 with the KO brake scaled out; E+time with Powell + GLS; lve_atr's real calls.**

- *Brake.* T6 LV–RA carries the KO brake, a proportional effect a median-centred SD cannot
  remove. Scaled out by one global factor, Δk_KO ≈ s·Δk_ctrl (weighted, through the origin;
  `wlsvar_eval.py t6` now reports `scale` and `z_scaled_robust_sd`), s = 0.88–0.91. LV–RA robust
  z SD, scaled (raw): A 0.85 (0.90), B 0.92 (0.94), C 0.67 (0.71), D 0.69 (0.75), E 0.88 (0.94),
  E+time 0.95 (0.98), E+band+time 1.01 (1.07). LV–RV has Δk ≈ 0, so s is not identified there;
  it stays raw.
- *E+time optimizer.* The E+time prototype used BFGS and the Hessian covariance. On the paired
  T6 contrasts that fits the time term on only 41–46 % of tests (the rest fall back to E); on
  the unpaired ones it fits 96–99 %. `WLSVAR_ET_OPT=powell` switches to Powell plus the GLS
  covariance, as in EJ+time:
  - time term on every T6 test, 2.4× faster;
  - Δk within 0.003;
  - median combined SE 10 % smaller;
  - LV–RV 0.82 → 0.88;
  - LV–RA 0.95 → 1.08 scaled (0.98 → 1.17 raw).

  T1, T2 and the simulation are not yet rerun with it.
- *lve_atr, unpermuted* (atrium vs ventricle, 348 tests at HEAD): BH < 0.05 calls A 234, C 235,
  E 218, E+time 229.

**G2 and WP5 (2026-11-14).** The G2 matrix, the decisions for Q1–Q3, the design note for E+time /
EJ+time and the open items are in the
[evaluation report](2026-11-14_wls_var_D_E_evaluation.md).
