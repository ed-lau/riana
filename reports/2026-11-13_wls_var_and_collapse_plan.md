# wls-var, the random-effects collapse (D) and the peptide-level model (E): scope of the evaluation

2026-11-13 · branch `1.2.0` · **Status: scope approved by the user 2026-11-13 (answers in §7).
Nothing is built; WP0 is next.** Follows
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

## 8. Progress log (moves into the evaluation report at WP5)

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
