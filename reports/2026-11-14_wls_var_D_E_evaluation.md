# Δk from the linear rollup: `wls-var`, the random-effects collapse (D) and the peptide-level model (E) — evaluation

2026-11-14 · branch `1.2.0` · **Status: evaluation complete (WP5 of the
[scope](2026-11-13_wls_var_and_collapse_plan.md)). User decisions 2026-11-15 (§5): remove
`wls-var` now; do not build D; build the peptide-level model after the EDL / SOL data, with E
vs E+time still to be settled (§5 Q3); test per-protein Δk as is under whole-group shifts (§9
item 1).** Prototypes: `tests/benchmark/wlsvar_*.py`. Detailed log, every table and its command:
scope §8. Background (candidates, mechanisms, the decision rule):
[Δk manual](2026-11-12_deltak_inference_manual.md) §8.1.

## 0. Summary

- **`wls-var` (B) should go.** It is never better calibrated than the shipped `wls` (A) and is
  worse on the ATF6 permutation null: 0.122 against 0.096 at α 0.05, and 49 against 20 BH calls
  per null permutation. Its per-point variance Ṽ is informative but floored: the spread
  between a protein's peptides is 2–4× Ṽ at the slope-carrying timepoints, and Ṽ leaves it out.
  1.2.0 is unreleased and no real experiment is known to use the flag. **Decided 2026-11-15:
  removed from `riana/`** while the replacement is settled.
- **D (random-effects collapse) is not worth building.** It behaves like C (shared proportions)
  on every yardstick: its τ² changes the peptide weights, not the animal term that the real-data
  tests turn on.
- **C removes the collapse-mixing bias (mechanism 1) and holds every Δk null, but its SE is too
  wide on paired contrasts.** Its T6 z SD is 0.51–0.71 (target 0.9–1.15). Under realistic
  timepoint deviation it also loses the power its size would allow.
- **E+time, the peptide-level mixed model with random slopes for labeling time and animal, is
  the only candidate that passes every G2 criterion but one.** That one is T6 LV–RV, which only
  B passes, and B fails the ATF6 and simulated nulls. E+time:
  - holds the real nulls: boomi 0.028, ATF6 0.041, every turnover quartile ≤ 0.046;
  - holds the calibrated simulation: false-positive rate 0.064, coverage 0.937;
  - has 1.16 × A's nominal power and 1.42 × at equal size;
  - gives honest SEs on lve_atr's time split (0.94) and ATF6's LV–RA animal split (0.95 with
    the KO brake scaled out);
  - fits ATF6 in about a minute.
- **For the interaction, EJ+time** (the same model across tissues) holds the ratio-scale and
  additive nulls on ATF6 and in the simulation, with the same power at equal size as E, joint
  and C. Its mean-Δk test has the most power (0.919 nominal at ratio 1.2, against A's 0.741).
  That test's null rate is inflated by whole-group proteome-wide shifts (§4, §7.1).
- **Shipped ATF6 results do not survive.**
  - The 35 per-chamber Δk calls under A are within A's permutation null.
  - The 79 interaction calls under A compare with 69 per null permutation.
  - Every candidate that holds the null makes 0 BH calls on either test.
  - The ATF6 biology notebook's Cops3 and Smyd1 interaction hits do not survive E, joint or
    EJ+time.
- **Before anything is built:** EDL / SOL (§6), which needs the combined search, and a choice
  between E and E+time (§5 Q3: simplicity against calibration and power).

## 1. What was compared

| | Collapse | Model / weights | Status before this evaluation |
|---|---|---|---|
| A | per-condition inverse variance | `wls`: (1 − θ̂)² | shipped default |
| B | as A | `wls-var`: (1 − θ̂)² / Ṽ, Ṽ moderated per point | shipped opt-in, unused |
| C | shared proportions (each peptide's weight level pooled over conditions) | `wls` | evaluator only |
| D | C, then random-effects weights 1/(σᵢ² + τ_rel²·g(θ̄)²) | `wls`; `D+var`: `wls-var` on Ṽ_RE | prototype |
| E | none: every peptide point | fixed sum-to-zero peptide slope offsets shared by the conditions; random animal slope (per-condition slopes with free covariance when animals carry both conditions); (1 − θ̂)² weights from an OLS pass | prototype |
| E+time | none | E with a random slope per labeling time (shared by every sample taken then, both conditions), animal nested in it | prototype |
| E+band(+time) | none | E (+time), weights also ÷ the peptide's shared band variance | prototype |
| E, joint (EJ) | none | the shipped `--interaction` model at peptide level: offsets per tissue, an animal slope shared by its tissues, a cell (animal × tissue) slope | prototype |
| EJ+time | none | EJ grouped by labeling time: a time slope, with animal, cell and tissue × time components | prototype |

What the candidates have to handle, as measured on real data:

- **Mechanism 1, collapse mixing.** The per-condition collapse mixes a protein's peptides in
  different proportions in each condition, so their persistent apparent-k differences become
  a Δk. C, D and E remove it.
- **Mechanism 2, proteome-wide animal shifts.** Each animal has a turnover factor of about 8 %
  (SD), shared by its tissues (ATF6, scope WP1d). A group of animals carries the mean of its
  members' factors, so whole permutations or datasets light up together. Only a per-protein
  animal term (E) carries this in the SE, and only partly (§7.1).
- **Timepoint-shared deviation.** At each labeling time, every sample shares a departure from
  the line, about √(0.0072² + (0.088 k)²) on ATF6 and 85 % shared across chambers. It cancels
  from Δk but inflates every SE that counts it as residual or animal variance. This makes C, D
  and E conservative on real data. E+time removes it.

## 2. Yardsticks

| | Test | Data | What passes |
|---|---|---|---|
| T1 | Replicate split: one culture plate per timepoint to each pseudo-condition, refit, roll up (20 assignments) | boomi (iPSC, D₂O) | ≤ 0.07 at α 0.05, no turnover quartile > 0.10 |
| T2 | Genotype labels permuted between animals within timepoint, refit, roll up (20 permutations; 10 for EJ+time) | ATF6 (4 chambers × 30 mice, DIA) | as T1 |
| T3 | Δk from odd vs even timepoints: z SD of the difference | lve_atr (atrium vs ventricle) | z robust SD 0.9–1.15 |
| T6 | Chamber Δk in the control vs the KO animals: z SD of the difference | ATF6, LV–RV and LV–RA | as T3 |
| S | Simulation calibrated to ATF6: design, gating, completeness, peptide spread, animal and proteome-wide shifts, timepoint-shared deviation (`null_ts`, 8 datasets; `r12_ts`, 4) | — | false-positive rate ≤ 0.065; coverage ≥ 0.93; power ≥ 0.9 × A |

The calibrated simulation reproduces ATF6's permutation null for every variant measured on both.
For Δk: A 0.094 vs 0.096, C 0.047 vs 0.039, E 0.026 vs 0.018, E+time 0.064 vs 0.041. For the
three interaction tests it agrees within 0.01–0.03. Before the timepoint term was added, it ran
0.05–0.08 hotter than ATF6 under every variant (scope §8). The pooled simulation rates are
dominated by the two datasets with the largest proteome-wide shift, so each has an SE of
0.01–0.04 (§4). The bound can rank the variants on the same datasets but does not fix their
absolute level more finely than that.

## 3. G2 — per-experiment Δk

| Criterion | A | B | C | D | E | **E+time** | E+band+time |
|---|---|---|---|---|---|---|---|
| T1 boomi null (worst quartile) | 0.059 (0.072) | 0.058 (0.081) | 0.028 (0.041) | 0.028 (0.039) | 0.009 (0.013) | **0.028 (0.033)** | 0.033 (0.034) |
| T2 ATF6 null (worst quartile) | ✗ 0.096 (0.108) | ✗ 0.122 (0.142) | 0.039 (0.046) | 0.038 (0.044) | 0.018 (0.025) | **0.041 (0.046)** | not run |
| T3 lve_atr z robust SD | ✗ 1.45 | ✗ 1.50 | 1.13 | not run | ✗ 0.75 | **0.94** | 1.07 |
| T6 LV–RA z robust SD, KO brake scaled out (raw) | ✗ 0.85 (0.90) | 0.92 (0.94) | ✗ 0.67 (0.71) | ✗ 0.69 (0.75) | ✗ 0.88 (0.94) | **0.95 (0.98)** | 1.01 (1.07) |
| T6 LV–RV z robust SD | ✗ 0.80 | 0.90 | ✗ 0.51 | ✗ 0.51 | ✗ 0.82 | ✗ **0.82** | ✗ 0.85 |
| S: false-positive rate | ✗ 0.094 | ✗ 0.104 | 0.047 | 0.046 | 0.026 | **0.064** | ✗ 0.070 |
| S: coverage vs the protein's estimand (ratio 1.2) | ✗ 0.898 | ✗ 0.887 | 0.953 | 0.956 | 0.975 | **0.937** | ✗ 0.926 |
| S: power at ratio 1.2 vs A, nominal | 1 | 1.03 | 0.98 | 0.98 | ✗ 0.87 | **1.16** | 1.23 |
| S: power vs A at equal size | 1 | 1.00 | 1.29 | 1.29 | 1.39 | **1.42** | 1.46 |
| Real ATF6, BH < 0.05 calls (its null: per permutation) | 35 (19.7) | 64 (49.4) | 0 (0.4) | 0 (0.3) | 0 (0.35) | **0 (0.4)** | — |

✗ = outside the G2 threshold. T6 LV–RA carries the KO brake: the chamber Δk is 9–12 % smaller in
the KO animals, a proportional effect that a median-centred SD cannot remove. As the scope
allowed, it is scaled out by one global factor (Δk_KO ≈ s·Δk_ctrl, s 0.88–0.91; `wlsvar_eval.py
t6`, column `z_scaled_robust_sd`). LV–RV carries no brake (its Δk ≈ 0 leaves s unidentified),
so it is shown raw. Not shown: D+var (as D, slightly more conservative) and E+band
(as E, conservative; 1.44 × at equal size). At ratio 1.5 every variant rejects 0.994–0.998, so
power was ranked at 1.2, the size of the real tissue effects (LA / RA 1.15–1.17 × LV, lve_atr's
atrium 1.24).

- **E+time is inside every threshold except T6 LV–RV.** No candidate that holds the nulls passes
  T6 LV–RV either. LV and RV share their turnover, so what is left is chamber-specific peptide
  behaviour that repeats in both animal sets. E puts it in its residual and SE, but it never
  reaches the T6 difference (§7.2).
- **E+time against E+band+time.** E+band+time has 3 % more power at equal size, but it breaks
  the simulated size bound (0.070) and the T3 / T6 bounds come close (1.07). On the same
  datasets it is hotter than E+time in every simulation set (by 0.006–0.009) and on T1, T3 and
  T6. It also needs C's shared bands as input. By the decision rule (the simplest if tied),
  E+time.
- **E alone is conservative wherever samples share a timepoint deviation** (boomi plates, ATF6
  animals, lve_atr's one mouse per timepoint). It counts that deviation as animal variance,
  which costs 13 % of A's power at nominal α. Its time-split SD on lve_atr is 0.75. On a
  simulation without the timepoint term, E is the better-sized of the two (0.059 vs 0.067),
  but that simulation does not reproduce ATF6's nulls.
- **Cost:** E+time takes about 1 min for the full ATF6 Δk (3,507 tests) at `-W 14`. In the
  simulation it takes 37 s per 300-protein dataset, against 35 s for A.
- **With the optimizer an implementation would use** (Powell plus the GLS covariance, §8), the
  prototype fits its time term on every paired T6 test, not on 41–46 % of them. Its T6 numbers
  move as follows:
  - LV–RV: 0.82 → 0.88;
  - LV–RA: brake scaled out 0.95 → 1.08; raw 0.98 → 1.17;
  - median combined SE: about 10 % smaller.

  Both T6 numbers stay inside or at the edge of the band. T1, T2 and the simulation have not
  been rerun this way (§9 item 4).

## 4. G2 — interaction (tissue × condition)

| Criterion | A | C | E, joint | **EJ+time** |
|---|---|---|---|---|
| T2 ATF6 null, ratio scale / additive / mean Δk | ✗ 0.216 / ✗ 0.245 / 0.044 | 0.065 / ✗ 0.081 / 0.019 | 0.025 / 0.035 / 0.010 | **0.027 / 0.041 / 0.033** |
| S null (`null_ts`), ratio / additive / mean Δk | ✗ 0.194 / ✗ 0.205 / 0.045 | 0.044 / 0.052 / 0.027 | 0.022 / 0.026 / 0.018 | **0.039 / 0.049 / ✗ 0.074** |
| … median over the 8 datasets | 0.189 / 0.202 / 0.034 | 0.044 / 0.052 / 0.011 | 0.024 / 0.024 / 0.007 | **0.042 / 0.048 / 0.027** |
| S, ratio-scale test under a real main effect (ratio 1.2: still a true null) | ✗ 0.204 | 0.062 | 0.033 | **0.049** |
| S, mean-Δk power at ratio 1.2: nominal / at equal size | 0.741 / 0.782 | 0.736 / 0.850 | 0.725 / 0.880 | **0.919 / 0.886** |
| S, additive power at ratio 1.2, at equal size (a small real effect) | 0.055 | 0.086 | 0.085 | **0.084** |
| Real ATF6, ratio-scale BH calls (its null: per permutation) | 79 (69) | 2 (1.2) | 0 (0.1) | **0 (0)** |
| Cost, full ATF6 (979 proteins, `-W 14`) | seconds | seconds | 13.6 min | 18.2 min |

- **Ratio-scale and additive tests:** both peptide-level models hold the null everywhere.
  EJ+time is nearer nominal in the simulation (0.039–0.049 against 0.022–0.026), and at equal
  size it matches E, joint and C. These tests sit in the cell stratum, which the shared
  timepoint deviation does not reach. So on ATF6 EJ+time stays about as conservative as
  E, joint (0.027 / 0.041).
- **Mean Δk** (the protein's Δk averaged over tissues): EJ+time recovers what E+time recovers
  on Δk, with 1.24 × A's nominal power. Its pooled simulated null of 0.074 is mechanism 2 at
  its strongest:
  - six of the 8 datasets reject 0.004–0.030; the two with the largest shift (median z −1.18
    and +0.90) reject 0.301 and 0.172;
  - averaging four tissues cuts protein-level noise but not the animal group's shift, which
    every tissue shares;
  - the pooled rate's SE is 0.038, so 0.074 does not separate from 0.05;
  - ATF6's permutations show the same pattern more mildly (0.033, worst 0.103).

  The fix is a question about the estimand, not about the model (§7.1).
- **Cost.** EJ+time takes 18 min for the full ATF6 interaction, against 14 min for E, joint.
  Most of it is one REML fit per protein, with peptide offsets for four tissues and three
  variance components. That is acceptable for a rollup step, but slow enough that an
  implementation should profile it, e.g. by starting each fit from the proteome-wide variance
  components instead of from zero.

## 5. Decisions

Recommendations from the evaluation; the user's decisions (2026-11-15) follow each.

**Q1, `wls-var`: remove it before 1.2.0 is released.** Under E every peptide point carries its
own (1 − θ̂)² weight, so a cell variance has no role, and D, which would have repurposed it as
Ṽ_RE, is not recommended. *Decided: remove it now, while the way forward is settled.* Removal
touches the following, with no results change for `wls` / `ols`:

- the `--linear-weights` choice (`cli.py`, `gui/tasks.py`);
- the moderation code in `core/linear_model.py` and `core/interaction_model.py`;
- the f-route refusal and the `theta_var` plumbing in `core/protein.py`;
- five test files;
- the docs pages and the CHANGELOG.

Keep or drop the `fs_var` / `fs_df` columns of `riana_rollup_fractions.txt` at implementation:
nothing in Riana would read them.

**Q2, D: do not build.** It equals C on every yardstick and addresses neither the pairing nor
mechanism 2. *Decided: not built.*

**Q3, E: build E+time as the peptide-level route, with EJ+time as its interaction form**, after
EDL / SOL (§6). Design note §8. This follows the manual's Direction (2026-11-13) for E, but
replaces D as the near-term step: C is the only collapse that holds the nulls, and its SE fails
the paired tests. Meanwhile, `wls` (A) stays the default. Its per-experiment Δk null is 0.096 on
ATF6 and its interaction null 0.216 on the ratio scale, so A's shipped calls on those designs
should be read against §7. *Decided: build the peptide-level model after EDL / SOL; E or E+time
is still to be debated (simplicity against performance).*

**E against E+time, for that debate** (numbers from §3–§4):

| | E | E+time |
|---|---|---|
| Real nulls: boomi / ATF6 | 0.009 / 0.018 (too strict) | 0.028 / 0.041 |
| Simulated null, coverage | 0.026, 0.975 | 0.064, 0.937 (at the 0.065 bound) |
| SE honesty: lve_atr time split (T3) / ATF6 LV–RA (T6) | 0.75 ✗ / 0.88 ✗ | 0.94 / 0.95 |
| Power vs A: nominal / at equal size | 0.87 / 1.39 | 1.16 / 1.42 |
| lve_atr real calls (A 234) | 218 | 229 |
| Interaction (EJ vs EJ+time): mean-Δk power, nominal; cost | 0.725; 13.6 min | 0.919; 18.2 min |
| Variance components (cross-sectional / paired) | animal / two per-condition animal slopes | time + animal / time + animal + cell |
| Optimizer, covariance | BFGS works; Hessian or GLS | Powell then Nelder–Mead; GLS required |
| Fallback | none needed | E (fewer than 3 times, or longitudinal) |

- **At equal size the two have the same power (1.39 vs 1.42).** The difference is calibration.
  E counts the shared timepoint deviation as animal variance, so its p-values are too large on
  every real design measured (boomi plates, ATF6 animals, lve_atr's one mouse per timepoint). At
  α 0.05 that costs about a quarter of E+time's calls (0.87 against 1.16 × A).
- **E+time contains E.** Its fallback for longitudinal designs and for fewer than 3 times is E, so
  building E+time means building E plus one variance component, the Powell optimizer and the GLS
  covariance. Building E first and adding the time term later costs little rework.
- **E+time is the less settled of the two.** Its simulated null sits at the bound, and the
  Powell + GLS rerun (§9 item 4) gives SEs about 10 % smaller, which could push it over.
- **The data that decides it:**
  - the Powell rerun of T1, T2 and the simulation: is E+time still ≤ 0.065?
  - EDL / SOL T3: does E's time-split SD stay near 0.75 on a second one-mouse-per-timepoint
    design?
- **Recommendation.** E+time with E as its built-in fallback, if the Powell rerun holds the bound.
  If it does not, E, with the time term deferred.

## 6. What is still missing before building

- **EDL / SOL** (the user's combined ATR / LVE / EDL / SOL search, expected in 2–3 days). It gives:
  - T3 on a second within-animal contrast;
  - the four-tissue animal diagnostic (scope WP1d), where tissue k differs enough to tell an
    additive animal slope from a proportional one (ATF6's chambers differ too little), which
    fixes the form of E's animal term;
  - a second check of E+time's time term on DDA.
- **No DDA in-vivo design with ≥ 2 animals per condition per timepoint** (scope §2 finding 5).
  ATF6's gated cells are about as incomplete as lve_atr's, which partly covers it.
- **Longitudinal designs** (the same animal at several timepoints) were not tested. There, E+time
  falls back to E, and a run effect beside the animal slope would be needed.
- **The f route** (rising precursor) was not evaluated for E. It needs the simulation, and the
  proprietary set checked locally only.

## 7. Consequences for shipped results

**7.1 The ATF6 KO brake.** 75 % of proteins are slower in KO (EJ+time's mean Δk); the median
per-chamber log₂ k ratio under A is −0.069, KO about 5 % slower. Per protein, this is the mechanism-2 pattern: a shift of the whole animal group.
It is within the permutation spread (7 of 20 permutations shift the proteome at least as far),
and no per-protein test that holds the null calls a protein. Whether a group-level shift is
itself the effect of interest is the estimand question in §9 item 1. The user decided to test
per-protein Δk as is.

**7.2 Per-chamber Δk and the interaction.** A's 35 per-chamber calls and 79 interaction calls
are at its null level (19.7 and 69 per permutation). Every candidate that holds the null makes
0 calls. Of the notebook's two interaction hits:

- **Smyd1** (A p ≈ 0, C 0.28, E, joint 0.77, EJ+time 0.78) is mechanism 1.
- **Cops3** (A p ≈ 0, C p ≈ 0, E, joint 0.22, EJ+time 0.19) follows the known outlying
  LA-control cell, which E, joint's cell variance carries and the shipped model, with only an
  animal slope, does not.

`notebooks/atf6_biology.ipynb` and the 11-10 report's interaction hits need a note, at the
user's discretion.

**7.3 lve_atr** (one mouse per timepoint, atrium vs ventricle) is a large real effect, and the
candidates keep it. On the HEAD refit (348 tests), the BH < 0.05 calls are A 234, C 235, E 218
and E+time 229 (rejection at α 0.05: 0.70, 0.69, 0.66, 0.68). E+time keeps 98 % of A's calls,
and its power is the same as A's on a design where mechanism 1 is small (manual §7).

## 8. Design note: the peptide-level model

**Per-experiment Δk (E+time).** Per (experiment, protein), on the gated peptide points with
t > 0 and φ = ln(1 − θ) above the plateau limit:

φ = −(k_c + b_p)·t + (u_time + v_animal [+ w_cell])·t + ε

- k_c is one slope per condition. b_p are peptide slope offsets, fixed, sum-to-zero, and shared
  by the conditions, so Δk is a within-peptide contrast.
- u_time ~ N(0, σ²_time) is shared by every sample taken at that labeling time.
- v_animal ~ N(0, σ²_animal) is nested in time (cross-sectional designs). When animals carry both
  conditions (tissues of one animal), a cell (animal × condition) component w is added.
- Weights are (1 − θ̂)² from an OLS pass on the same fixed effects, from the fitted θ̂, not the
  observed one, normalized to mean 1. The fit is REML.
- Δk = k_test − k_ref. Its SE comes from the GLS covariance at the fitted components,
  (Σ X′V⁻¹X)⁻¹.
- df = animals − times − 1 (between animals) or animals − 1 (paired).
- Fallbacks:
  - fewer than 3 times, or any animal at more than one time → E: the animal slope, or
    per-condition animal slopes with a free covariance when paired;
  - a failed fit → E → no animal term.

**Interaction (EJ+time).** Per protein, over the tissues that have both conditions:

φ = −(k_{s,c} + b_{s,p})·t + (u_time + v_animal + w_cell + x_tissue×time)·t + ε

The tests are those of the shipped `--interaction`. df:

- mean Δk: animals − times − 1;
- additive and ratio scale: cells − animals − (S − 1)·times − (S − 1), with the tissue ×
  condition × time stratum pooled with the cells.

**Implementation notes, all from the prototypes:**

- **Optimizer.** Powell, then Nelder–Mead. BFGS stalls at its start values with ≥ 2 nested
  components. The prototype's E+time still uses BFGS: it fits its time term on 99 % of ATF6's
  per-chamber tests and 96–97 % of lve_atr's, but only 41–46 % of the paired T6 tests, which
  fall back to E. With Powell and the GLS covariance (switch `WLSVAR_ET_OPT=powell`), it
  fits the time term on all of them, 2.4× faster. T6 stays in range (§3); T1, T2 and the
  simulation are still to rerun.
- **Covariance.** Use the GLS covariance, not the inverse Hessian, which is ill-conditioned near a
  zero component. That produced 24 negative variances on ATF6 before the fix.
- **Tissue filter.** A tissue enters the joint model only with animals under both conditions.
  The prototype's first filter did not enforce this and gave two p-values near 1e-100 on ATF6
  (fixed; the shipped `_prepare` is correct).
- **df rule.** The prototypes use containment df, which is what was validated here. The shipped
  `--interaction` uses Satterthwaite. Keep containment unless Satterthwaite is checked on the
  same yardsticks.
- **Animal key.** `characteristics[individual]`, else the condition × timepoint × biological
  replicate composite. With one animal per timepoint (lve_atr, EDL / SOL), the time slope is
  the animal slope, which is how E+time gets the pairing there without IDs.
- **Estimand.** E's estimand is the unweighted average over peptides, whereas the collapse
  weights peptides by precision.
  - E's k level sits 1–3 % above that estimand (point-wise φ Jensen plus additive offsets on a
    multiplicative truth). Δk is unaffected to within 0.15 SE.
  - Its ratio reads about 2 % low at ratio 1.2 (1.18 for 1.20).
  - A second-order correction, φ + σ²/(2(1 − θ)²), applies if E's k values are reported.
- **Multiplicative peptide offsets** (Em) were tried and rejected: no gain in coverage, and a
  false-positive rate of 0.072.
- **CLI placement** (manual §8.1 (a), options not yet decided):
  - E+time becomes a new `--method` (name open; `peptide` is a candidate), linear model only;
  - `--interaction` becomes its joint form;
  - `--linear-weights` keeps `wls` (and `ols` for audit) for the collapse route;
  - `pooled` is retired after an audit.

## 9. Open items

1. **The estimand under whole-group shifts.** A per-protein test cannot tell a group's
   proteome-wide shift (the ATF6 brake, or a chance difference between two sets of animals)
   from protein-level effects, and its false positives cluster when the shift is large.
   Options:
   - (a) a proteome-level test of the shift with animals as the units, e.g. permutation, or a
     mixed model on each animal's proteome-wide factor;
   - (b) per-protein Δk relative to the group's shift (Δk beyond the global brake);
   - (c) the per-animal proteome-wide factor (scope WP1d) as a covariate.

   (b) and (c) change what a "hit" means, which is the user's call.

   **Decided (user, 2026-11-15): test per-protein Δk as is.** The per-protein results are the
   evidence of a group shift, protein by protein. Shifts relative to the group median are a
   separate scientific question and are not built into Riana for now (recorded in
   `PROJECT_REVIEW.md`). The yardsticks stay as they are: a real group shift (the ATF6 brake)
   shows up as many proteins moving together, and a chance one still counts against the
   pooled simulated nulls (E+time 0.064; EJ+time's mean Δk 0.074).
2. **T6 LV–RV (0.82; 0.88 with Powell + GLS).** Chamber-specific peptide behaviour is shared by both animal sets. It
   belongs in the SE if the estimand is the protein's peptides at large, and not if it is the
   peptides measured. A peptide × condition variance component would carry it; not tested.
3. **Limma-style residual-variance moderation** (PROJECT_REVIEW, deferred) should be done on the
   winner: per-protein σ²_animal and σ²_time, moderated across proteins.
4. **E+time with Powell + GLS** (§8): T6 done (§3). Rerun T1, T2 and the simulation with
   `WLSVAR_ET_OPT=powell` before building. A 10 % smaller SE could move the simulated null
   (0.064) past 0.065.
5. **EDL / SOL** when the combined search lands (§6).

## 10. Reproduce

All runs are under the gitignored `runs/wlsvar_eval/` and `runs/wlsvar_sim/`. Each command and
its output table is in scope §8. `wls-var` was removed from Riana on 2026-11-15, so the variants
that use it (B, C+wls-var, D+var) need Riana at commit `e12980c` or earlier. The main commands:

```bash
# Real data: T1 / T2 (null), T3, T6, interaction; --variants picks A B C D D+var E E+band E+time …
python tests/benchmark/wlsvar_eval.py null runs/wlsvar_eval/atf6/t2/t2_0* --ref control --test ATF6-KO --variants A C E E+time --cache
python tests/benchmark/wlsvar_eval.py t3 runs/wlsvar_eval/lve_atr/t3_odd runs/wlsvar_eval/lve_atr/t3_even --ref control --test atrium --variants A C E E+time --cache
python tests/benchmark/wlsvar_eval.py t6 runs/wlsvar_eval/atf6/fit --ref LV --test RA --variants A B C D E E+time E+band+time --cache
WLSVAR_ET_OPT=powell python tests/benchmark/wlsvar_eval.py t6 runs/wlsvar_eval/atf6/fit --ref LV --test RA --variants E+time   # no --cache
python tests/benchmark/wlsvar_eval.py null runs/wlsvar_eval/lve_atr/fit --ref control --test atrium --variants A C E E+time --cache   # §7.3
python tests/benchmark/wlsvar_eval.py inter runs/wlsvar_eval/atf6/t2/t2_0* --ref control --test ATF6-KO --variants EJ EJ+time --cache -W 12
# Simulation calibrated to ATF6 (timepoint term shared across tissues), scoring, size-adjusted power
python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/null_ts --ratio 1 --datasets 8 --set time_k_sd=0.09 --set time_kabs_sd=0.0115 --set time_shared=0.85
python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/null_ts --ratio 1 --pooled
python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/r12_ts --ratio 1.2 --size-adjust runs/wlsvar_sim/null_ts
# Calibration statistics, real vs simulated (incl. the cross-tissue day-deviation share)
python tests/benchmark/wlsvar_calib.py runs/wlsvar_eval/atf6/fit --ref control --test ATF6-KO
```
