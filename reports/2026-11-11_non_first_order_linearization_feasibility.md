# Linear and mixed models for non-first-order labeling (Guan / Fornasiero) — feasibility

**Date:** 2026-11-11 (§10 added 2026-11-12) · **Status:** feasibility only, no production code ·
**Data:** simulations of an illustrative slow-precursor in-vivo design and of mouse designs, using
Riana's own IsoSpec envelope model and θ solve, and the ATF6 mouse series (`runs/timeseries_atf6`) ·
**Scripts:** `runs/nonfirstorder_linearization/` (gitignored)

**Question.** Can the `linear simple` machinery (φ-space WLS, per-experiment Δk, the random-animal
interaction model) be used when the precursor does not reach plateau instantly, i.e. where Riana
currently offers the nonlinear Guan and Fornasiero models? The proposal: fit the precursor model
first, back-calculate a precursor-corrected fraction new per point, then roll up with the simple
linear model. Precedent: Topograph (Hsieh et al., *MCP* 2012;11(11):1468–74,
[doi:10.1074/mcp.O112.017699](https://doi.org/10.1074/mcp.O112.017699)).

## TL;DR

- **Feasible, with one change to where the correction happens.** Under any shared precursor curve,
  the fraction of protein made since labeling began is still first-order, f = 1 − e^{−kt}, so
  φ = log(1 − f) = −kt is linear through the origin. Only the *label content* is not first-order. If
  each point's f can be recovered, the existing linear, mixed and interaction models apply
  unchanged.
- **The back-calculation must be made from the isotope envelope, not from Riana's θ.** θ is solved
  as a mixture of the unlabeled and the *plateau*-labeled envelope. Molecules made while the
  precursor was still rising do not fit that two-envelope picture: θ overstates their label content
  by +17 to +36 % over the first ten days of a slow-precursor design, more for peptides with more
  labeling sites. Inverting θ through the Guan formula therefore overestimates
  k by **+5 to +42 %**, peptide-dependent. **This bias also affects Riana's current
  `fit --model guan` / `fornasiero` on D₂O data.**
- **The envelope-level route works.** Solving each timepoint's envelope for f under the kinetic
  forward model (new molecules spread over synthesis times, each at the precursor enrichment of its
  time), then fitting the line with per-point information weights, recovers k within **1 %**, with
  the precision of a joint nonlinear envelope fit (spread within ~5 %).
- **The two published shortcuts do not hold for a slow precursor.** A lag intercept
  (Topograph's ln(1 − n) = −λ(t − t₀)) and correcting θ by the *current* precursor enrichment are only
  valid when k ≪ k_p. In the slow-precursor simulation they leave k biased by −18 to −84 %.
- **The precursor curve must be shared, accurate, and specific to each condition.**
  - *Per-peptide fits fail.* k_p is not identifiable per peptide: Guan's formula is symmetric in k
    and k_p, so a fast peptide's fit swaps them.
  - *A shared fit works.* Estimated across 50–800 proteins, a shared k_p is recovered within
    ~3–7 %.
  - *Errors in k_p matter.* A 25 % error in k_p biases k by 10–50 %.
  - *Conditions must not share it blindly.* If two conditions differ in precursor kinetics and are
    analysed with one k_p, nearly every protein shows a false Δk.
- **Hard limit.** When k ≫ k_p the protein's labeling just tracks the precursor, and no method
  recovers k well. With k_p = 1/14 per day, a 30 % Δk is poorly resolvable for proteins faster than
  ~0.1/day.
- **Mice (ATF6, primed):** a finite k_p (≈ 3–9/day) is detectable; correcting changes k by 1–5 %,
  differently by genotype. Small but not zero.
- **Intermittent dosing** (doses a few times a week, no priming) front-loads the precursor. A smooth
  k_p curve can then over-estimate k — by about 5 % for slow proteins and far more for fast ones in
  an idealized twice-weekly regimen — depending on the exact calendar; check it against the dosing
  record (§9).
- **Recommendation: go** — test on a real slow-precursor dataset; §6 lists checks that need no
  ground truth.

## 1. What is and is not first-order

With a shared precursor curve P(s) = p(s)/p_max (0 at the start of labeling, 1 at plateau), a protein
with turnover rate k has:

- **fraction made since t = 0:** f(t) = 1 − e^{−kt}, independent of the precursor;
- **label content** (fraction of plateau labeling): A(t) = ∫₀ᵗ k e^{−k(t−s)} P(s) ds.

Guan's model is P(s) = 1 − e^{−k_p s}, giving the closed form Riana fits today; Fornasiero adds
precursor recycling from protein breakdown (a bi-exponential P(s) in Riana's implementation). In both,
P(s) is the same for every protein. A(t; k) increases monotonically in k at fixed t, so every point
with 0 < A < P(t) maps to a unique k_i, and hence to f_i = 1 − e^{−k_i t}. That is the proposed
back-calculation; under the model it is exact.

Each sample's envelope mixes molecules made throughout labeling, each at the precursor enrichment of
its own synthesis time. Correcting by the enrichment at the sampling time alone therefore treats
earlier, less-labeled molecules as old, and under-estimates f (§4).

Three routes were compared:

| Route | What it does |
|---|---|
| **Lag intercept** (Topograph eq. 5) | φ = log(1 − θ) with a free intercept, φ = −k(t − t₀). |
| **Current-precursor correction** (Topograph's use of the per-sample RIA) | Divide by the precursor enrichment at the sampling time, φ = log(1 − θ/P(t)). Assumes all new molecules carry today's enrichment. |
| **Model inversion** (the proposal) | Map each point to k_i through A(t; k), then φ_i = −k_i t, through the origin. |

The first two are asymptotic. For Guan at large t, 1 − A ≈ (k_p/(k_p − k)) e^{−kt}, a line with
intercept log(k_p/(k_p − k)). That intercept is the same for every protein only when k ≪ k_p and
diverges as k → k_p. Topograph's own Fig. 2 shows such a residual intercept (t₀ ≈ 1 day) after its
precursor correction.

## 2. Kinetic-level simulation (label content, k_p known)

Collapsed protein points with θ-scale noise σ = 0.03. Slow-precursor design: k_p = 1/14 per day
(0.0714), labeling days 1, 3, 5, 10, 18, 31, 46, 60 and 79, one animal per condition × day. Mouse
design: 1, 3, 7, 10, 14 days × 3 animals. 2,000 simulations per cell (`sim_kinetic.py`).

**Relative bias of k, slow precursor (k_p = 1/14 per day):**

| k (/day) | production (naive) | lag intercept | current precursor | model inversion | nonlinear Guan |
|---|---|---|---|---|---|
| 0.005 | −25 % | −22 % | −22 % | +1 % | 0 % |
| 0.02 | −33 % | −27 % | −27 % | 0 % | 0 % |
| 0.05 | −51 % | −44 % | −36 % | 0 % | 0 % |
| 0.1 | −66 % | −60 % | −46 % | −2 % | 0 % |
| 0.2 | −77 % | −73 % | −55 % | −3 % | 0 % |
| 0.4 | −86 % | −84 % | −65 % | −7 % | +1 % |

- **CI coverage:** model inversion 0.88–0.94 (slow precursor) and 0.90–0.95 (mouse); the naive model falls
  to 0 from k = 0.02 on.
- **False-positive Δk rate:** model inversion 0.04–0.09, against 0.02–0.06 for the nonlinear fit.
  It is slightly high for very slow proteins, which a second reweighting step should address.
- **Power** for a 30 % Δk: close to the nonlinear fit up to k = 0.05 (0.80 vs 0.82). It drops for
  faster proteins with every method (k = 0.2: 0.13 vs 0.22) — the precursor-limited regime.

**Mouse with a mild lag (k_p = 1/day):** naive bias −5 to −36 %; the lag intercept leaves −4 to
−20 %; model inversion and the nonlinear fit are unbiased.

## 3. The precursor curve: accuracy, sharing and conditions

*Wrong k_p* (`sim_kp.py`, slow-precursor design; model inversion and nonlinear Guan behave alike):

| assumed k_p / true | 0.5 | 0.75 | 1.33 | 2 |
|---|---|---|---|---|
| k bias, k = 0.01 | +40 % | +12 % | −8 % | −15 % |
| k bias, k = 0.05 | +123 % | +30 % | −16 % | −29 % |

*Two conditions with the same k but different precursor kinetics* (k_p(B) = 0.75 × k_p(A)):

| Analysis | False-positive Δk rate |
|---|---|
| One shared k_p (either condition's, or their mean) | 16–79 % (slow precursor, k = 0.005–0.1) |
| A k_p for each condition | 0.04–0.08 |

*Estimating k_p:*

- **Shared across the proteome** (`sim_kp_estimate.py`): profiling the summed residual of
  per-protein Guan fits over a k_p grid recovers k_p within the grid resolution (±2.5 %). Its CV is
  2–7 % with 50–800 proteins, whether the tissue's proteins are slow (median k 0.01/day) or fast
  (0.1/day).
- **Per peptide** (`sim_kp_per_peptide.py`): the interquartile range of k_p spans ±30–100 % for slow
  peptides. Because A(t) is symmetric in k and k_p, a peptide with k = 0.3/day returns k_p ≈ 0.30:
  the parameters swap. Freeing k_p per peptide also widens the k spread about 4.5-fold for
  k ≤ 0.03/day and 2-fold at 0.1/day, compared with a known k_p.

## 4. The isotope envelope: Riana's θ is not the label content

Riana solves θ by least squares as a mix of the natural-abundance envelope and the envelope at the
plateau enrichment. A molecule made at a lower enrichment p has a binomial label distribution at p,
not at p_max. Its m0 depletion, 1 − (1 − p)ⁿ, is concave in p, so the solve reports more "new"
material than its share of enrichment.

The test used real peptide sequences with labeling sites from `deberneh_2025_rss`. True envelopes
were integrated over synthesis times through IsoSpec and passed to `solve_fs_d2o` (`sim_envelope.py`;
p_max = 0.04, k_p = 1/14 per day).

**θ / A(t) − 1** (Riana's θ against the Guan label content):

> **Corrected 2026-11-11.** The first version of this table reported +55 to +94 % for
> YLYEIARRHPYFYAPELLYYANK. That was a binning bug in Riana's envelope code, not physics: its
> monoisotopic mass (2952.496) sits at the edge of the ±0.5 Da windows the code used, so its isotope
> peaks landed in the wrong channels. The binning is fixed (CHANGELOG, 2026-11-11) and the table
> re-run; the other peptides were unaffected or changed by ≤ 1 point.

| Peptide (sites) | day 1 | day 5 | day 10 | day 31 | day 79 |
|---|---|---|---|---|---|
| AEFVEVTK (15) | +20 % | +22 % | +17 to 18 % | +6 to 8 % | +0.3 to 2 % |
| DDPHACYSTVFDKLKHLVDEPQNLIK (43) | +31 % | +34 % | +27 to 30 % | +8 to 14 % | +0.4 to 4 % |
| YLYEIARRHPYFYAPELLYYANK (39) | +27 % | +30 to 31 % | +24 to 26 % | +7 to 12 % | +0.4 to 3 % |

Fixes compared:

- **Correcting θ by the current precursor** (Topograph-style) instead under-reports f by −50 % on
  day 1 and −30 to −41 % on day 10.
- **A per-peptide closed-form law** replacing linear label content with m0 depletion,
  (1 − (1 − p_max P)ⁿ)/(1 − (1 − p_max)ⁿ), brings the short peptides within 2–5 % from day 3.
  It leaves −15 % (QTALVELLKHKPK) and −33 to −34 % (the 39- and 43-site peptides) on day 1, still
  −16 to −17 % on day 10 for those (`sim_envelope_fix.py`), because the solve uses m0–m5, not m0
  alone.
- **A per-point envelope solve** (`sim_envelope_noise.py`): for each timepoint, find the k that
  makes the kinetic forward model's envelope match the observed one. Then f_i = 1 − e^{−k_i t},
  φ_i = −k_i t, and the line is fitted with weights from each point's information,
  w_i = ‖∂envelope/∂k‖² / t². How it works is set out below.

**How the envelope solve works.** Today's solve has one free parameter, θ, with the enrichment
fixed at the plateau:

  E_obs(t) ≈ (1 − θ)·E(0) + θ·E(p_max),

where E(p) is IsoSpec's envelope for the peptide with n labeling sites, each deuterated with
probability p. Under a rising precursor the new molecules are not all at p_max. A molecule made at
time s carries the precursor enrichment of that time, p(s) = p_max·P(s), and under first-order
turnover at rate k it survives to t with weight k·e^{−k(t−s)}; a fraction e^{−kt} is still old. So

  E_pred(t; k) = e^{−kt}·E(0) + ∫₀ᵗ k·e^{−k(t−s)}·E(p(s)) ds.

- *One parameter per point.* The solve finds k_i, the way today's finds θ, and reports
  f_i = 1 − e^{−k_i t}. At a given t, f and k are one-to-one; parameterizing by k is what fixes the
  new molecules' mix of enrichments. A fast protein's new molecules were made mostly recently, at
  high enrichment; a slow protein's are spread over the whole period, at lower enrichment. With
  k_p = 1/14 per day, the mean enrichment of the new molecules (as a fraction of plateau) is:

  | day | k = 0.005 | k = 0.02 | k = 0.1 | k = 0.5 | precursor at that day |
  |---|---|---|---|---|---|
  | 3 | 0.10 | 0.10 | 0.10 | 0.12 | 0.19 |
  | 10 | 0.29 | 0.29 | 0.33 | 0.43 | 0.51 |
  | 31 | 0.61 | 0.64 | 0.76 | 0.87 | 0.89 |

  Correcting with the enrichment at the sampling day (last column) therefore over-states how
  labeled the new molecules are, and under-estimates f.
- *No first pass.* The k inside the integral is the unknown being solved at that point, not a value
  carried over from a peptide-level fit. Each candidate k implies its own mix of enrichments, and
  the solve keeps the k whose predicted envelope matches the observed one. Neither a Guan fit nor θ
  enters, so a biased first-pass k cannot propagate. Every point yields its own k_i; under
  first-order turnover they agree, and the line through the φ_i = −k_i t_i is what tests that. The
  only second step is the line's weights (below), evaluated at a pilot k from an unweighted line
  through the same φ_i; an error in the pilot changes the weights, not what each point estimates.
- *It reduces to today's solve.* With an instant precursor (p(s) = p_max for s > 0) the integral is
  (1 − e^{−kt})·E(p_max): today's mixture with θ = f.
- *It is cheap.* E(p) is computed once per peptide on a grid of enrichments (about 80 cached IsoSpec
  envelopes instead of 2), so the integral is a weighted sum of cached envelopes whose weights depend
  only on t, k and the precursor curve. Each candidate k costs microseconds; each point is a 1-D
  minimization, as now. Any P(s) — Guan, Fornasiero or a measured curve — changes only the weights.
- *Its weight in the line fit* follows from how strongly the envelope responds to k at that point,
  ‖∂E_pred/∂k‖², instead of (1 − θ̂)²; early points, where the precursor is still low, respond
  little and count little.

**k recovered from noisy envelopes** (3 % relative + 0.003 absolute channel noise, 120 simulations):

| Peptide (sites), k | Guan fit on θ (today) | inversion of θ | **envelope solve → line** | joint envelope fit |
|---|---|---|---|---|
| LVNELTEFAK (17), 0.01 | +5 % | +5 % | **0 %** (SD 6 %) | 0 % (SD 6 %) |
| LVNELTEFAK (17), 0.1 | +26 % | +26 % | **0 %** (SD 9 %) | 0 % (SD 8 %) |
| QTALVELLKHKPK (19), 0.1 | +22 % | +23 % | **−1 %** (SD 8 %) | −1 % (SD 7 %) |
| DDPHACYSTVFDKLKHLVDEPQNLIK (43), 0.1 | +41 % | +42 % | **0 %** (SD 5 %) | 0 % (SD 5 %) |

The weights are essential: with the simple model's (1 − θ)² weights the same route kept +2 to +15 %
bias, because the earliest points carry almost no information and must count for little.

## 5. Real data: ATF6 mice

A shared Guan k_p was profiled per chamber × genotype on the rollup's collapsed θ points
(`atf6_kp_profile.py`; θ-level, which is adequate when k_p is fast):

| | k_p (/day) | median k(instant) / k(Guan) |
|---|---|---|
| control (LA, LV, RA, RV) | 5.9, 2.9, 3.7, 3.3 | 0.95–0.98 |
| ATF6-KO (LA, LV, RA, RV) | ∞, 8.3, ∞, 9.4 | 0.98–1.00 |

The precursor rises within hours (these mice were primed), but not instantly. Correcting would raise
k by 1–5 %, more in controls. The difference is plausibly driven by the three day-1 animals in each
genotype, the animal-level shift found in the random-animal analysis
(`reports/2026-11-10_random_animal_lme_feasibility.md`). With one sample per animal, precursor lag
and animal-to-animal variation cannot be separated. The profile intervals are too narrow, because
~1000 proteins share the same few animals.

## 6. A slow-precursor dataset: what to check without ground truth

Design points that shape the analysis:

- With k_p near 1/14 per day, most proteins sit in the transition regime; the naive model would
  under-estimate k by roughly 25–86 % (§2).
- With one animal per condition × timepoint, a per-experiment Δk rests on about N − 2 residual df
  per protein (one point per animal).
- Tissues taken from the same animals make the tissue × condition comparison a split-plot; the
  interaction model's random animal effect applies, keyed on `characteristics[individual]`.

Checks that need no ground truth:

1. **Precursor:** profile k_p per condition (and per tissue). It should agree across tissues of the
   same animals if body water is the effective precursor; compare it with any value derived from the
   dosing scheme (§9).
2. **Within-protein agreement:** per-peptide k spread within each protein (geometric CV) for the
   naive model, today's Guan fit on θ, and the envelope route. The θ bias depends on the peptide's
   envelope, so under today's Guan fit, high-site peptides should read faster than their protein; the
   envelope route should remove that trend. Test: correlation of each peptide's residual
   log(k/protein k) with its number of labeling sites.
3. **Shape:** residuals against time for each route (the naive model should show early-time
   curvature), and whether the corrected line passes through the origin.
4. **Conditions and tissues:** Δk with a k_p per condition, and the tissue interaction with a random
   animal effect; report how much each changes when k_p is shared.

## 7. What it would take (if pursued)

- **SDRF / intake.** Read `characteristics[precursor rate constant]` per row, as enrichment is read
  today, and allow a per-(experiment, condition) value. Optionally estimate it: a profile of the
  shared k_p across the proteome, per condition.
- **Fit.** A precursor-aware per-timepoint solve. For each point, find the k that best matches the
  observed envelope under the kinetic forward model; report f = 1 − e^{−kt} and its variance from the
  envelope sensitivity. Per peptide, the cost is about 80 cached IsoSpec envelopes on an enrichment
  grid plus a 1-D solve per point. Guan and Fornasiero become precursor curves feeding this solve,
  not separate curve models; any measured or parametric P(s) works the same way. The peptide fit and
  admission gate run on label content (§10). The same idea applies to ¹⁸O with its envelope model.
- **Rollup.** Unchanged models on f: φ = log(1 − f), through the origin. Weights come from the
  per-point variance, through the existing `fs_var` / `wls-var` path, not the (1 − θ̂)² rule, which
  assumes an instant precursor.
- **Inference.** Treat k_p as a per-condition parameter; propagate its uncertainty (profile or
  bootstrap) or at least report sensitivity. Flag proteins with k ≳ k_p as precursor-limited.
- **Existing defect.** Today's Guan / Fornasiero fits on D₂O θ carry the §4 bias whenever the
  precursor is not at plateau. This should be noted in the roadmap whatever is decided here.

## 8. Limits of this analysis

- **Simplified noise and precursor model:**
  - θ noise was independent and normal.
  - Envelope noise used a simple relative-plus-absolute model.
  - The precursor followed Guan exactly (§9 checks intermittent dosing), the same in every animal,
    with no animal-level variation.
  - Amino-acid-specific equilibration with body water was not modeled.
- **Small test set:** the envelope test used five peptides at one enrichment (0.04).
- **Fornasiero not simulated:** its recycling curve was not tested. The envelope route only needs
  P(s), but its three precursor parameters are much harder to identify from the data than Guan's
  single k_p.
- **ATF6 is θ-level only:** the ATF6 profile uses collapsed θ, and its intervals ignore shared animals.

## 9. Deriving the precursor curve from a dosing scheme

Without a body-water measurement, the precursor has to come from the dosing scheme and the data.

**Mass balance.** With repeated D₂O doses, body water reaches a mean plateau
plateau = (D added per day, as a fraction of body water) / k_tbw, where k_tbw is the fractional
body-water turnover. The D added per dose is dose volume × D₂O fraction / body-water volume, so the
assumed body-water fraction of body mass matters directly:

- Reptile values are higher than the 50 % often assumed: green turtle 66.7 % fed / 58.7 % fasted
  ([doi:10.1242/jeb.029330](https://doi.org/10.1242/jeb.029330)); diamondback terrapin ≈ 75 %
  ([doi:10.1242/jeb.110411](https://doi.org/10.1242/jeb.110411)). Going from 50 % to 60–75 % lowers
  the k_tbw inferred from a given plateau by 17–33 %.
- Water-turnover methods and comparative values: Streit, *Comp Biochem Physiol A* 1982;72:445–54
  ([doi:10.1016/0300-9629(82)90107-4](https://doi.org/10.1016/0300-9629(82)90107-4)).
- Deuterium also exchanges into non-water hydrogen; the dilution space is about 4 % larger than body
  water in mammals, a minor correction.
- Feeding can change water turnover: green turtles' water flux was 9.6 % of body water per day fed
  against 6.1 % fasted (same reference). Where conditions may differ, estimate the precursor per
  condition and test the difference (§3: a shared k_p across differing conditions gives
  near-universal false Δk).

**Intermittent dosing front-loads the precursor.** Each dose is a step, and the first alone reaches a
large share of the plateau. For two doses a week (days 0, 3, 7, 10, …), k_tbw = 1/14 per day, scaled
to the same mean plateau (`kp_note_checks.py`):

| day | 1 | 3 | 5 | 10 | 18 | 31 | 46 | 79 |
|---|---|---|---|---|---|---|---|---|
| pulsed (fraction of plateau) | 0.23 | 0.45 | 0.39 | 0.73 | 0.83 | 1.05 | 1.04 | 0.96 |
| smooth Guan, k_p 1/14 | 0.07 | 0.19 | 0.30 | 0.51 | 0.72 | 0.89 | 0.96 | 1.00 |

**Is a smooth k_p adequate?** Fitting the smooth Guan curve to label content generated from an
idealized twice-weekly regimen (days 0, 2, 7, 9, …; `generic_dosing_check.py`):

| true half-life | 100–350 d | 23–46 d | 14 d | 7 d | ≤ 3.5 d |
|---|---|---|---|---|---|
| k_p fixed at k_tbw | +6 % | +7 to +11 % | +20 % | +54 % | fit fails |
| best-fitting smooth k_p (an effective value) | +5 % | +6 to +9 % | +16 % | +46 % | fit fails |

The best-fitting smooth k_p is an effective parameter of the approximation, not the water turnover:
under front-loaded dosing it comes out above k_tbw (a smooth curve can only mimic early label by
rising faster), so a k_p fitted to real data under a smooth model over-states k_tbw. How much it
absorbs depends on the exact calendar (when each animal was dosed relative to sampling), so check the
smooth curve against the study's own dosing record before relying on it.
Fast proteins (half-lives under about a week) are not described by a smooth curve under intermittent
dosing.

**Bounds from fast proteins.** With a smooth precursor, a protein's apparent simple-model rate can
never exceed k_p (k_p = 1/14 per day):

| true k (/day) | 0.06 | 0.1 | 0.2 | 0.3 | 0.5 | 1 | 5 |
|---|---|---|---|---|---|---|---|
| apparent rate (/day) | 0.029 | 0.037 | 0.047 | 0.053 | 0.058 | 0.064 | 0.070 |

The fastest apparent rates in a dataset therefore bound k_p from below. Under intermittent dosing the
bound loosens, because front-loaded pulses can push fast proteins' apparent rates above k_tbw.

**The plateau and slow proteins.** In the low-labeling regime θ ≈ label / plateau, so absolute k of
slow proteins scales roughly with 1/plateau. Re-analysing a 4 % plateau as 6 % or 8 % (with the k_p
the mass balance then implies) changes the absolute k of proteins with half-lives of 100–200 days by
−26 to −28 % and −38 to −41 %, while a 2-fold ratio between conditions moves only −4 % and −7 %. The
ratio-scale tests are therefore much more robust to the plateau than absolute half-lives.

**Estimating the precursor from the data.** Proteins with k ≳ 0.3/day track the precursor with a lag
of about 1/k, so each animal's fast proteins report its own precursor level, and their plateau sets
the enrichment ceiling. Caveat: the envelope shift scales with labeling sites × enrichment, so the
plateau trades off against the labeling-site table; for a species without its own table the plateau
is an effective value, consistent within the dataset but carried into absolute k.

## 10. The peptide fit and admission gate under a rising precursor

*Added 2026-11-12.* The per-point f feeds the rollup (§7). Each peptide still needs its own fit: to
report its k and, more important, to decide whether it is measured well enough to enter the protein
rollup. Today that gate is R² ≥ 0.8, or R² ≥ 0.6 with k_cv < 0.2. The gate stays as it is, for
parity with today's default; what changes is the quantity R² is computed on.

**Label content, not θ and not f.** Label content A is the peptide pool's mean enrichment divided by
the plateau enrichment. To first order it is what the envelope's mean mass shift reads, whatever
precursor model is assumed. Under an instant precursor A = f = θ. Under a rising precursor:

- θ reads the envelope as a mixture of the unlabeled and the plateau envelopes, which no longer
  gives the mean (§4).
- f is unstable early. No protein can carry more label than the precursor has reached, so the whole
  range f ∈ [0, 1] maps into A ∈ [0, P(t)]. While P(t) is a few percent of plateau, the envelope
  noise spans the whole range of f.
- An information-weighted R² in f puts nearly all the weight on the late points, where a slow
  protein's f spans little, so it fails slow proteins.

**How R² is computed.**

1. *Each point:* the envelope solve finds the k_i whose predicted envelope matches the observed one
   (equivalently f_i = 1 − e^{−k_i t_i}).
2. *Each point:* the observed label content is A_i = A(t_i; k_i), the label content of that matching
   envelope. At t = 0, A_i is the baseline read.
3. *Each peptide:* fit one parameter, k, to the A_i by unweighted nonlinear least squares, with the
   label-content curve

   A(t; k) = ∫₀ᵗ k·e^{−k(t−s)}·P(s) ds = 1 − (k_p·e^{−kt} − k·e^{−k_p·t}) / (k_p − k),

   for P(s) = 1 − e^{−k_p·s} with k_p fixed (from the SDRF). This is the Guan curve; as k_p → ∞ it
   becomes 1 − e^{−kt}.
4. R² = 1 − Σ(A_i − A(t_i; k̂))² / Σ(A_i − Ā)², unweighted, the same formula as today; k_cv comes
   from the same residual bootstrap as today.

No fit with an instant precursor comes first. The information-weighted φ line on f (§4, §7) is the
linearized form of this fit: on a real slow-precursor dataset the two gave k within a few percent
and near-identical R². The fitted k̂ also gives the model value f̂ = 1 − e^{−k̂t} at which each
point's Var(f_i) is evaluated.

**Today's `simple` fit and the rising-precursor fit, side by side.**

| | today, `simple` | rising precursor |
|---|---|---|
| per-point value | θ_i from the two-envelope solve | k_i from the precursor-aware solve → A_i (and f_i) |
| peptide curve | 1 − e^{−kt} | the Guan curve above, k_p fixed |
| fit | nonlinear least squares on θ_i, unweighted | the same, on A_i |
| R² | 1 − SS_res/SS_tot on θ | the same formula, on A |
| k_cv | residual-bootstrap 90 % CI | the same |
| passed to the rollup | θ_i | f_i with Var(f_i) |

Today's `fit --model guan` fits the same curve, but to θ, which is not label content under a rising
precursor; that is the source of its bias (§4).

**Why not gate on today's exponential fit to θ.** It would work as a "well behaved" filter, but it
scores model misfit as well as noise. Under a rising precursor the true curve is S-shaped and an
exponential cannot follow it, so well-measured peptides lose R² for the shape. On the real dataset it
gave systematically lower R² and admitted fewer peptides, while ranking peptides similarly. An R²
against the right curve measures only the scatter of the measurements.

## Reproduce

```bash
python runs/nonfirstorder_linearization/sim_kinetic.py          # §2
python runs/nonfirstorder_linearization/sim_kp.py               # §3 wrong / per-condition k_p
python runs/nonfirstorder_linearization/sim_kp_estimate.py      # §3 shared k_p (~4 min)
python runs/nonfirstorder_linearization/sim_kp_per_peptide.py   # §3 per-peptide k_p
python runs/nonfirstorder_linearization/sim_envelope.py         # §4 θ vs label content
python runs/nonfirstorder_linearization/sim_envelope_fix.py     # §4 m0-depletion law
python runs/nonfirstorder_linearization/sim_envelope_noise.py   # §4 envelope route
python runs/nonfirstorder_linearization/atf6_kp_profile.py      # §5
python runs/nonfirstorder_linearization/kp_note_checks.py       # §9 pulses, bounds, plateau
python runs/nonfirstorder_linearization/generic_dosing_check.py # §9 smooth k_p adequacy
```
