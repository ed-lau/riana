# Riana — Roadmap & TODO (v1.2.0 forward)

> **Status (2026-07-21):** `1.1.0` is **fully released** (PyPI, GitHub, Zenodo, `master`
> merged back). Active development is on the **`1.2.0` branch**. **This document is the
> forward-looking roadmap — it lists OPEN work only.** The as-shipped record lives in
> `CHANGELOG.md`; design decisions / point-in-time findings live in the `memory/` files
> and `reports/*.md`. §2 and the M1–M4 pointers are kept as collapsed anchors because
> ~36 code docstrings reference `PROJECT_REVIEW.md §2b/c/d`, `§3`, `§4.1/4.2`, and
> `Track A–E` — the anchors must survive even as their prose is trimmed.
>
> **Doc-hygiene rule (see §7):** when an item ships, move it from here to `CHANGELOG.md`
> (leave only a one-line pointer if a code anchor needs it) and update the owning
> `memory/` file's status line. Keep this file open-work-only.
>
> Maintainer: Edward Lau.

## 1. Project status

Riana is a single-author scientific Python tool for extracting and modelling isotopomer
time-series from MS1 data, used in protein-turnover research. The scientific core
(accurate-mass, kinetic models, fractional-synthesis math, IsoSpec) is lifted-unchanged
and proven; the 1.0.0 rewrite replaced the surrounding software (typed `core/` +
`algorithms/` + `io/`, streaming mzML, SDRF/manifest intake, peak detection, a PySide6
GUI). **Package layout → README; as-built feature record → `CHANGELOG.md`; open work →
§3 + §6.**

## 2. Critical findings (historical anchors — full record in `CHANGELOG.md [1.0.0]`)

Collapsed to one line each; the code refs need the anchors to resolve, and the still-open
carry-overs are flagged.

- **§2b — scientific defects:** all fixed in 1.0.0 (unreachable AA `a_0`/`label` dispatch,
  `iso0/colsums` FS drift, heuristic kinetic CI → residual bootstrap, fixed-site-count FS
  bias). FS now comes from the IsoSpec per-peptide-Spep forward/solve model.
- **§2c — algorithmic gaps:** addressed in 1.0.0 (apex/consensus peak detection, in-window
  baseline options, per-isotopomer observed-mass/drift). **Still open →** robust in-window
  **baseline subtraction** (`noise_floor` off by default) and a **cross-proportion-stable
  peak picker** (Track B).
- **§2d — architecture:** addressed in the 1.0.0 rewrite (one frozen
  `IntegrationConfig`/`FitConfig`; typed `core/` + async GUI off the main thread; PySide6
  replacing Tkinter; global-state modules restructured).

## 3. Roadmap

M1–M4 (→1.0.0) and the 1.1.0 experimental-science + QoL line are **shipped** — itemized in
`CHANGELOG.md`. This section is the OPEN work: the **v1.2.0 release scope** first, then the
standing **Track A–E** clusters (which the v1.2.0 items draw from).

### v1.2.0 — release scope

**Theme: multiplexed & single-timepoint labeling, hardened.** Finish that line cleanly;
defer the large research items so the release stays shippable.

**Shipped on `1.2.0` so far** (record in `CHANGELOG.md [Unreleased]`): masking-vectorization
perf (`searchsorted`), parallel test suite + `slow` tier, selectable two-condition linear-simple
Δk, **TMT/TMTpro + single-timepoint labeling** (pinned-isotope envelopes, isobaric SDRF collapse,
single-timepoint auto-detect at rollup + `--min-fit-points` + `n_pts<2`→NaN CI), **multi-point FS
rail-drop** (physical rails 1.05/−0.05), **dimethyl channel→sample intake** (true sample-axis
multiplexing), the quantms-1.8.0 q-value fix, and the **linear-simple WLS** default.

#### In scope — the release spine

**1. Single-timepoint hardening** *(SHIPPED 2026-07-21/22 — see `CHANGELOG.md`)*

Single-timepoint data (TMT, dimethyl, any one-labeling-time run) — the flagship 1.2.0
substrate — is now handled directly. **Shipped:** fit-step detection + `--depth` auto-relax;
the direct k solve (closed form `k̂ = −ln(1 − F̄S)/t*` for `simple`, pole-safe Brent root-find
for guan/fornasiero) replacing futile NLS in **both** `_fit_one_concat` and the rollup refit
`_fit_kdeg`, routed locally by `np.ptp(t)==0` (so multi-timepoint peptides that collapse to one
t after rail-drop benefit too); R²→NaN; the single-point `OptimizeWarning` gone. **GUI:** the
rollup tab exposes **Min fit points** and relabels **Max k_cv** (secondary rescue vs
single-timepoint primary gate) + a single-timepoint hint; the fit tab's summary and Depth
tooltip no longer mislead. All numerically identical to the old NLS. **Only remaining
(optional):** a `reports/` writeup of the derivation + ~100× A/B (low-value — the rationale is
captured in `CHANGELOG.md` and results are unchanged).

**2. A1 — dimethyl S=1 spillover: a thin shared primitive, not an elaborate gate**

Context (`reports/2026-07-06_…_spillover.md`, `reports/2026-07-13_…_intake.md`): only **S=1**
peptides (N-term dimethyl only, R-terminal, no K) are at risk — heavy iso0 lands on light
iso8 — and only **long S=1 (≥~26 res) at high FS** are badly corrupted. That is **~1–2% of
quantifiable peptides**. The end-to-end duplex run showed the heavy channel's degradation is
**mostly SNR imbalance (0.628 heavy/light, 8% dropouts), not spillover**, and that spillover
pushes Δk *against* the observed effect (self-correcting, not manufacturing signal).

- **Exact gate vs. current heuristic:** the "heuristic" is a blunt length cutoff (drop S=1
  peptides > ~20–22 res). The **exact** gate computes the per-peptide light-tail spill into
  the heavy window from the IsoSpec forward model (sequence + Spep + enrichment, worst-case
  FS=1) — so it **scales with SDRF precursor enrichment and per-peptide envelope width**,
  where a fixed length cutoff both over- and under-drops. Its real advantage, though, is that
  **the primitive is the same computation the demux needs** (below) — so it is not throwaway.
- **Recommendation (matches the maintainer's instinct):** ship a *minimal* increment — factor
  the spill primitive out of `tests/benchmark/bench_dimethyl_spillover.py`, emit a `spillover`
  reference column, and drop/flag heavy FS above a conservative default. **Low priority** (the
  corrupted corner is small, SNR-dominated, and direction-safe); the trivial length heuristic
  is an acceptable stopgap. Do **not** build elaborate gate logic — invest that effort in the
  demux instead (deferred; see Track C).

**3. Easy wins pulled into 1.2.0** *(maintainer-prioritized)*

- **Per-point `Var(θᵢ)` → depth-aware WLS weight — SHIPPED** as opt-in
  `--linear-weights wls-var` (`c82cc55`; new `fs_var` / `fs_df` columns in
  `riana_rollup_fractions`). See `CHANGELOG.md`.
- **PyInstaller standalone GUI bundle — SHIPPED** (`packaging/`, 2026-08-08).
- **`mypy --strict` on `riana/algorithms/` — SHIPPED** (CI `typecheck` job, 2026-08-08).
- **User-facing docs refresh — SHIPPED** (2026-11-11): Quarto site rewritten for 1.2, generated
  CLI Reference page (`riana_website/_gen_reference.py`, kept current by `tests/test_docs.py`),
  README slimmed to a landing page.

#### Deferred beyond 1.2.0 (with rationale)

- **A2 — SILAC proper.** Registered as *geometry* only in `riana/multiplex.py`; its
  heavy-UNIMOD constants + intake CV terms land **when SILAC-D₂O data exists**. Blocked on data.
- **A3 — deamidation** (chemical fit-merge). Its own side project: +0.984 / C13-M+1 isobaric
  overlap needs joint envelope + deamidation-proportion modelling → silent k bias otherwise.
- **A4 — adaptive-N_ISO + robust envelope. RESCOPED — see the box below.** Short version:
  ~80% shipped; the headline (adaptive *capture*) was benched and rejected as a default; the
  rule-based N_ISO we have is good enough. Only the robust matcher genuinely remains, and it is
  research-grade. **Not a 1.2.0 item.**
- **Slow precursor — per-point envelope solve + fit/rollup model redesign — DESIGNED
  2026-11-11, productionize after the report §6 real-data checks.** Answers most of the old
  *precursor-enrichment (RIA) modelling* item: under a rising precursor the envelope is rebuilt per
  timepoint by a kinetic forward model (not the plateau mixture), the precursor curve is
  `RIA_max·(1−e^{−k_p·s})` with `k_p` set per experiment from the SDRF, and `k_p` is **not** fit
  from the data by default (per-peptide it is not identifiable; proteome-wide it is only weakly
  identified on real data). **See the box below** for the agreed design, including removing
  guan/fornasiero one release after the solve ships. Still open from the old item: a measured
  (non-parametric) precursor curve, if a dataset with one appears.
- **Δk peptide collapse: per-point vs per-peptide scatter — HIGH PRIORITY, plan agreed
  2026-11-12, nothing built.** The collapse weights each peptide by its own per-condition band, so
  the two conditions average the same peptides in different proportions; with real within-protein
  peptide k spread (CV ~0.2) that is a Δk which is not turnover (*mechanism 1*). On ATF6 it
  accounts for most shipped calls (Δk 37 → 0, interaction 86 → 2 with shared proportions; lve_atr
  250 → 244). The user is also not convinced by `wls-var` even on θ (its evidence never had
  peptide offsets or animals). Plan, candidates (A shipped `wls`, B `wls-var`, C shared
  proportions, D random-effects collapse with τ², E peptide-level mixed model with an animal
  slope) and the evidence order (calibrated simulation → real-data permutation nulls → agreement
  with peptides' own Δk → precision): `reports/2026-11-12_deltak_inference_manual.md` §8.1.
  Also there: the interaction's additive animal slope vs proportional animal effects (§8.2).
  Comes **before** the residual-variance moderation below, which should moderate the winner.
  *To do first for it:* run **EDL / SOL** (the lve_atr sister comparison, Study 4) through
  current Riana — no current-version fit exists yet.
  *Order (user, 2026-11-13):* after the rising-precursor rollup stage and the non-primed D₂O
  dataset analysis and its (local) report, write a **new plan** from §8.1: does `wls-var` do what
  it is meant to (weight by measured precision) at all, and should it be replaced by a
  peptide-specific term — a peptide offset modelled in the collapse (L1), or Ṽ replaced by the
  across-peptide scatter.
- **Within-animal one-factor comparisons — not supported yet (noted 2026-11-13).** When the two
  conditions are tissues (or other strata) of the *same* animals and there is one experiment, the
  plain Δk test does not read `individual` (only `--interaction`, which needs ≥ 2 experiments × 2
  conditions, does). Each animal's two tissues are fitted as independent points, so the animal's
  shared enrichment and turnover count in the error instead of cancelling in the contrast:
  conservative, with power lost. Mechanism 1 above still applies, and may be systematic here
  (peptides measured better in one tissue). Fix: the interaction model's random animal slope with
  the tissue as the within-animal factor, φ ~ 0 + day:tissue + (0 + day | animal). Not built; to
  weigh alongside §8.1 (candidate E has the same animal slope).
- **Δk residual-variance moderation (limma-style) — HIGH PRIORITY, deferred (decided
  2026-11-11).** Shrink each protein's residual variance toward a prior estimated across proteins
  and test Δk on the augmented df, to stabilise proteins with few points. Deferred until we know
  how the shipped baseline (`wls`, `wls-var`) behaves across complex samples. Its benefit depends
  on the per-protein residual df and dispersion, which vary by design. **See the box below** for
  what to measure first, the build/skip rule, and implementation notes. Not data-blocked:
  characterization can start on the existing runs and grows as new datasets arrive.
- **Proper demultiplexing (Basisty *TurnoveR*-style).** The endgame for multiplexing: forward-model
  the light FS and **subtract** the predicted light spillover from the heavy peaks *before*
  solving heavy FS (e.g. light iso6 overlapping heavy iso9), rather than dropping the heavy
  peptide. **Rescues** long S=1 (what the A1 gate only curates out) and becomes **mandatory** for
  tighter SILAC/DIMETHYL2/4/6 spacings where even short peptides overlap. Shares the A1 primitive.

##### A4 rescope — adaptive N_ISO + robust envelope: what we are actually up against

The 2026-06 build cycle (`memory/m8_adaptive_niso_robust_envelope.md`;
`reports/2026-06-23_adaptive_niso_limited_isotopomer.md`) already shipped most of this and
**benched the core value prop to a verdict**:

- **Shipped:** B0 (`ria_max`/`adaptive_iso`/`--ria`/`--iso auto` + SDRF enrichment wiring);
  B1 (`adaptive_channel_masses` — per-peptide N_ISO from the init∪final ≥1% union, iso0 ==
  precursor m0, averaged-isotopolog accurate-mass targets); B2 (max-width NaN-padded schema);
  B3 (H4′ mix-then-normalize FS solve); B4 (`--fs`/`--fs auto` fit-time subset scoring). The
  mass-defect θ / `fs_ds` half shipped in **1.1.0**.
- **Bench verdict (ac16/cm/ipsc, decisive):** adaptive **capture** (`--iso auto`) does **not**
  earn its keep at D₂O 4.6–6% enrichment — neutral-to-negative on θ-spread, k-CV, and recovery,
  and 3–5× slower to integrate. Narrow **scoring** (`--fs 0 1 2 3`, iso0-3 sweet spot) is the
  free, fit-side keeper. **⇒ the rule-based N_ISO rules are good enough; keep `--iso auto`
  opt-in, do not default.**
- **What genuinely remains** (all deferrable): (1) the **robust observed-vs-IsoSpec matcher** —
  Huber/soft-trim/per-channel-SNR downweighting of a contaminated high channel while still
  using the clean ones; this is the real open research question of the whole effort (Track B).
  (2) **B6 calibration high-θ recovery bench** — was gated on re-searching the calibration mzML
  (verify current status). (3) A **per-peptide scoring-width optimizer** — score each peptide on
  `min(its envelope, cap)` so a flat `--fs` doesn't under-score genuinely-wide peptides (the
  "genuine per-peptide channel optimizer" Track B item; gated on the Track D animal bench).

##### Δk residual-variance moderation — characterize the baseline first

**What it is.** The `linear simple` Δk test estimates each protein's residual variance `s²_g`
from that protein's own collapsed points, with `d_g = N_g − p_g` df (`N_g` points after the
φ/t filters, `p_g` condition slopes). limma (Smyth 2004) moderates it across proteins:
`s̃²_g = (d₀·s₀² + d_g·s²_g)/(d₀ + d_g)`, and the moderated t has `d_g + d₀` df, with `(d₀, s₀²)`
from a (robust) fitFDist over the proteins. limma-trend makes `s₀²` a function of a covariate.

**How it differs from what shipped.** `wls-var` (opt-in, 2026-07-30) moderates the *per-point*
Var(θ) that sets the point weights; the contrast keeps `N − p` df, and `s²_g` stays per protein.
Residual moderation acts on `s²_g` itself and is the only lever that changes the df. Under the
shipped weights `s²_g` is a dispersion relative to the delta-method weights (it also carries
biological and animal-to-animal variation); moderating it remains valid. The two can be combined.

**Why it is deferred.** It only matters where `d_g` is small and proteins' dispersions are
similar (large `d₀`). Rich time series give `N − p` ≈ 20+ (the 2026-07-24 report: "nothing to
augment there"); single-timepoint TMT/dimethyl, sparse DIA, proteins missing timepoints, and
designs with many conditions give far fewer. We do not yet know where real projects fall.

**Step 1 — characterize the baseline (the gate).** For each dataset × design, from the shipped
`wls` (and `wls-var`) rollups, record in a `reports/` entry (one table, extended per dataset):

1. **Residual df** per tested protein: median `d_g`, and the share with `d_g ≤ 3`, `≤ 5`, `≤ 10`.
2. **Dispersion:** the distribution of `s²_g`, and a robust fitFDist `(d₀, s₀²)` per experiment
   and pooled. `d₀` against the median `d_g` predicts how far moderation would move each test.
3. **Mean–dispersion trend:** `s²_g` against `n_peptides`, mean θ, k and `N_g` (decides plain
   limma vs limma-trend).
4. **Null calibration** of the current test: within-condition splits (two pseudo-conditions from
   one condition's biological replicates) or condition labels permuted across animals within each
   timepoint; the p-value histogram and the share below 0.05 (should be ≈ α).
5. **Reproducibility:** biological-replicate leave-one-out agreement of Δk and of significance
   (the protocol of `notebooks/wls_var_investigation.ipynb`).
6. **Yield:** BH hits at 0.05, per experiment.

Existing runs to start from, spanning the designs: `runs/timeseries_atf6` (in-vivo DIA, 4 tissues
× 5 timepoints × 3 animals, 2 conditions); `runs/lve_atr_clean` (in-vivo DDA, ventricle vs
atrium); `runs/boomi_ipsc_d2o` / `_o18`, `runs/juber_ac16_d2o` / `_o18`, `runs/timeseries_lauren9`,
`runs/timeseries_lauren5_7_ipsc_mesoderm_o18` (cell culture, rich D₂O / ¹⁸O series);
`runs/dimethyl_duplex` (single timepoint, dimethyl); `runs/splatd2o_tmt16`, `runs/splatd2o_cm_tmt18`
(single timepoint, TMT; rolled up with the kinetic model only so far, and pooled plexes may not
allow a Δk); `runs/timeseries_dia` (sparse 3-timepoint DIA; not rolled up). Add complex samples as
they arrive (more tissues, more than two conditions, longitudinal sampling).

**Enabling change (small, no change to results):** write each protein's residual variance and its
df beside `delta_k_se` in `riana_rollup_proteins.txt` (e.g. `delta_k_s2`, `delta_k_df`). Items
1–3 currently have to be recomputed outside Riana. Do this first when the work resumes.

**Build / skip rule (proposed; confirm when the table exists).** Build it (opt-in first) if, in a
design Riana supports, a material share of tested proteins (≥ ~20 %) has `d_g ≤ 5` **and** the
robust `d₀` is comparable to or larger than their `d_g`, **or** the null splits show the current
test is miscalibrated in a way moderation would fix. If `d_g` is large everywhere, record the
evidence here and close the item.

**Implementation notes (for later).**

- *Where:* `riana/core/linear_model.fit_linear_deltak`, after the per-protein fits and before BH.
  Reuse the robust `_fit_fdist` already there (it estimates the `wls-var` prior). Estimate
  `(d₀, s₀²)` per experiment, matching per-experiment BH, unless item 2 shows a shared prior is
  better.
- *Surface:* an opt-in rollup flag (and a Protein-tab checkbox). Moderated SE / df / p in the Δk
  columns. `d₀` and `s₀²` go in the provenance header. Regenerate the CLI Reference
  (`tests/test_docs.py`).
- *Weights:* must work with `wls` and `wls-var` (and `ols` for audit).
- *Validation:* add a residual-moderation arm to `tests/benchmark/bench_linear_weights.py`. It has
  none today; its `moderate` subcommand moderates per-point Var(θ). Cover sparse-design
  scenarios (Type-I, power, coverage). Check R parity with limma: `lmFit` on the protein × sample
  φ matrix with weights `(1−θ̂)²` and design `0 + day:condition`, then `contrasts.fit` +
  `eBayes(robust = TRUE)`. The collapsed points share a sample axis (biological replicate ×
  timepoint), so the limma layout fits.
- *Out of scope:* the interaction model. Its variance components (animal, residual) are per
  protein; moderating them across proteins is a different problem.

##### Slow precursor — per-point envelope solve and the fit/rollup model redesign

**Status.** Design agreed with the user 2026-11-11 from the feasibility study
(`reports/2026-11-11_non_first_order_linearization_feasibility.md`). **Fit stage BUILT
2026-11-12** (items 1–3, 3a, 3b, the item-4 k_p switch and the item-6 header):
`riana/algorithms/precursor.py`, `fit --kp`, SDRF `characteristics[precursor rate constant]` →
manifest `precursor_rate_constant`. Choices made while building:
- *Var(f).* Not the envelope Jacobian of item 3: the peptide's label-content residual-bootstrap
  band (as θ's band today) divided by the gain ∂A/∂f at f̂. **But `wls-var` on f is
  miscalibrated** (Monte Carlo below) — the rollup on f uses `wls`; items 3 and 7 are revised.
- *Header:* `fs_basis theta|f|mixed` and `k_p_applied` per curve (the config's own `k_p` line is
  the deprecated guan/fornasiero parameter).
- *Scope:* D₂O only; `fs_ds` not computed on f; a curve's runs must agree on k_p; `--kp inf`
  forces the steady path; never applied to calibration runs or under guan/fornasiero.
- *Cost:* same as the θ fit (convolution envelopes per peptide × timepoint, not cached run-wide).
- On the real dataset k, R² and gate yield matched the prototype (k within ~3 %, R² ρ 0.998).

**Rollup stage (7) BUILT 2026-11-13, except the default:** rollup reads the fit's `fs_basis`
(`protein.fit_basis_header`), echoes `fs_basis`/`k_p_applied` in its headers outside the config
hash, and on f/mixed refuses `--linear-weights wls-var` and guan/fornasiero
(`protein.check_fs_basis`; CLI and GUI); `linear simple` → `linear` with the old name accepted
(warning). The Monte Carlo bench runs its f + `wls-var` arm with the refusal lifted (audit only;
reproduces the stored arm exactly). **Default flip to `linear` held:** it would make the shipped
per-condition collapse (known issue, manual §7) everyone's default Δk — flip once §8.1 decides.
**Remaining:** the default (7), deprecations (5) and the `exponential` name (4), GUI φ-axis
labels (say θ on f), validation before merge.

**Monte Carlo of the Δk / interaction tests (2026-11-12) — full write-up, numbers and
reproduction in `reports/2026-11-12_deltak_inference_manual.md` (the single source of truth for
Δk variance; do not duplicate its numbers here).** For this item: on f, `wls` is calibrated
(0.042 at α 0.05, every turnover rate; ratio recovered) and `wls-var` is not (0.255) — hence item
7. The same work found that the **shipped θ test** is miscalibrated through the per-condition
peptide collapse (mechanism 1) and through animal variation growing with time (mechanism 2), and
that mechanism 1 dominates ATF6's calls — tracked as its own item below (*Δk peptide collapse*).
The 3b ceiling drop is to be revisited (cap instead of drop).

**Why.** Without a priming dose the D₂O precursor rises over days. Today's per-point θ assumes
every new molecule carries the plateau enrichment, so under a rising precursor θ reads below the
fraction new f (the naive k is low). The guan/fornasiero curves fitted to θ are biased too,
because θ is not label content either (report §4). Under any precursor curve shared by the
proteins, f = 1 − e^{−kt} is still first-order. So the fix belongs in the per-point solve: recover
f_i with the precursor built in, and everything downstream (exponential peptide fit, R² gate,
rollup, linear Δk, interaction) runs unchanged on f.

**Design.**

1. *Where.* The per-point solve stays in `riana fit`. Integrate stays signal extraction that knows
   nothing about labeling: the solve depends on the label, the coefficient table (labeling sites),
   the plateau enrichment, `--fs` channel scoring and the rail drop (and now k_p); changing any of
   them must not mean re-reading the mzML. Fit keeps two internal steps: envelope → fs per point,
   then fs over time → k.
2. *The solve.* For each point, find the k_i whose predicted envelope
   `E_pred(t; k) = e^{−kt}·E(0) + ∫₀ᵗ k·e^{−k(t−s)}·E(RIA_max·P(s)) ds` matches the observed one
   (H4′: mix in the full-cluster basis, truncate to the scored channels, renormalize). Write
   f_i = 1 − e^{−k_i t_i} into `fs`. No first pass: k inside the integral is the unknown at that
   point. Allow f_i < 0 (k < 0 continues the formula smoothly) so noise around zero is not clipped;
   rails as today. With an instant precursor it reduces exactly to today's θ, so the default stays
   byte-identical. Build the labeled envelopes by binomial convolution on an enrichment grid
   (Track B item) — the integral mixes many enrichments.
3. *Per-point variance.* Emit Var(f_i) = σ²_env / ‖∂E_pred/∂f‖² (σ²_env from the peptide's envelope
   residuals). It feeds the rollup's inverse-variance collapse (`_weighted_theta`, in place of or
   beside the bootstrap band) and the `wls-var` weights. The φ-space information weight is
   (1 − f̂)²/Var(f_i); today's (1 − θ̂)² is its instant-precursor special case (Var constant over
   time). Under a rising precursor early points carry little information, so this weighting is
   **required, not opt-in**: plain (1 − θ̂)² weights bias k by +2–15 % (report §4). The rollup's
   per-cell `fs_var`/`fs_df` (same columns as today, now of f) carry this; `wls-var`'s eBayes
   moderation must not shrink away the precursor factor, which is known, not estimated. In the
   prototype it did not (prior df ≈ 1 against 10–16 per cell; early cells kept their ~10²–10³×
   variance within 10 %), but production should moderate the envelope-noise part and apply the
   gain after, so thin cells cannot.
   **Evaluate Var(f_i) at the model value f̂ = 1 − e^{−kt}, never at the observed f_i** (pilot k
   from the unweighted line through the origin, then one IRLS step — the same rule as the shipped
   WLS). Early in a slow-precursor series the new molecules carry only a few percent of the plateau
   enrichment, so a small baseline bias in the envelope (a few hundredths in θ) solves to f near 1.
   Variance taken at that observed f under-weights nothing and wrecks the fit; taken at f̂ those
   points get ~1/100 of a late point's weight. Treat solves at the f bound (k → ∞) as failed.
   A real-data prototype (2026-11-11) confirmed both: observed-f variance collapsed the yield, and
   f with plain `wls` weights lost almost all Δk power, while model-value variance with `wls-var`
   matched or beat θ.
   **The interaction model needs the same weights.** `interaction_model` always weights by
   (1 − θ̂)² from a first unweighted pass, whatever `--linear-weights` says. On f the low-information
   early cells then get full weight: in the prototype the mean-Δk test found almost nothing and the
   animal variance collapsed; dropping the early points restored it (but cost per-experiment Δk
   power, so the fix is the weights, not dropping points). Pass (1 − f̂)²/Var(f) into the mixed
   model (nlme `varFixed` form). **DONE 2026-11-12:** the interaction follows `--linear-weights`.
   **Superseded by the Monte Carlo (fit-stage status above): `wls-var` on f is anti-conservative
   (the interaction too); f uses `wls`.** *Confirmed beforehand (local test):* the mixed model with wls-var weights (moderated Var as
   `fit_linear_deltak` does) brought the mean-Δk test back above θ's call count, with the animal
   variance back at θ's level. With the shipped weights the earliest cells carried nearly all of the
   weighted residual sum of squares; with wls-var, a small share.
3a. *Gate: today's default, for parity (user, 2026-11-12).* Keep R² ≥ `--min-r2` OR (R² ≥
   `--rescue-r2` AND `k_cv` < `--k-cv`), with R² and `k_cv` from a peptide fit on **label content**
   (the pool's mean enrichment ÷ plateau — to first order the envelope's mean mass shift, whatever
   the precursor model):
   - per point, the envelope solve gives k_i; the observed label content is A_i = A(t_i; k_i)
     (t = 0: A_i = the baseline read);
   - per peptide, fit A(t; k) = 1 − (k_p·e^{−kt} − k·e^{−k_p·t})/(k_p − k) — the Guan curve with
     k_p fixed, i.e. ∫₀ᵗ k·e^{−k(t−s)}·P(s) ds — to the A_i by unweighted nonlinear least squares:
     today's fit code with this curve and A_i as input;
   - R² = 1 − Σ(A_i − Â_i)²/Σ(A_i − Ā)², and `k_cv` from the same residual bootstrap as today (so the
     `--rescue-r2` floor's reason, collapsed CIs at k ≈ 0, carries over unchanged).
   Under an instant precursor A = f = θ and the curve is 1 − e^{−kt}: today's gate is the special
   case. Label content is not today's θ (a two-envelope mixture at plateau, off under a rising
   precursor, report §4). Not R² in f: f is unstable near the early label ceiling (3b), and an
   information-weighted R² in f fails slow proteins. Not today's exponential fit to θ either: it
   scores the lag shape as misfit. In the prototype this fit gave k within a few percent of the
   information-weighted φ line (its linearization) with near-identical R²; the gate admitted somewhat
   more peptides than today's θ gate, with more Δk calls and θ's calls nearly all recovered. Its k̂
   also sets f̂ for the model-value Var(f_i) in 3 (in place of the pilot line). *Possible later
   simplification:* `k_cv` alone, from the analytic SE of the weighted line, was also safe on f —
   revisit on a second dataset.
3b. *Early points and the rails — DECIDED 2026-11-12 (user): rails on label content, no floor.* Early in a slow-precursor
   series no protein can carry more label than the precursor has reached, so the whole range
   f ∈ [0, 1] maps into label content [0, P(t)]. For a precursor half-life of ~10 days the first
   day's ceiling is a few percent of plateau, about the width of the envelope noise. The envelope →
   f map is then nearly flat: Var(f) is hundreds of times Var(θ) (the weights in 3 handle that), and
   **the rails, applied in f, select**. f > −0.05 there means label content above ≈ 0, so the solve
   drops the light half of the early reads and every read above the ceiling (a bound hit), and keeps
   the heavy half. In the prototype the kept first-day points read a per-point k_i more than an
   order of magnitude above the fitted k (median), still several-fold on the second sampled day. Weighted fits barely notice (an early point's pull on
   the slope scales with its gain, a few percent), but anything unweighted does: the shipped
   interaction weights, QC plots of f over time, the k_p profile. Options: apply the rails on the
   measured scale (θ or label content) rather than on f; keep early points unrailed with their tiny
   weights; censor points below an information floor (the mirror of plateau truncation); cap reads
   above the ceiling at the ceiling instead of dropping them. Independent of this, the early envelope
   reads themselves are biased (Track B "Early-timepoint findings").
   *Prototype test (same solves, three rail variants; solve bound widened so label content can read
   to −0.1).* Rails on label content (A_i > −0.05, as θ's below-natural rail; solves at the f upper
   bound = reads at or above the ceiling, dropped as failed) against rails on f, and against
   label-content rails plus an information floor (< 1 % of the peptide's best point's weight):
   - *Selection.* On the first sampled day, f rails dropped ~60 % of reads as light and ~20 % at the
     ceiling; label-content rails ~30 % and ~20 % — about symmetric, as noise around a small true
     value should be. Early median label-content residuals fell to ≈ 0 (from about +0.01–0.02).
   - *Downstream.* More proteins tested and more Δk calls (θ's calls equally recovered, same
     agreement with θ); the interaction unchanged in calls, with more proteins.
   - *Information floor:* no gain; not needed with the weights.
   - Early per-point k_i is meaningless at that gain (with label-content rails it read *negative*
     at the median: the weak-envelope light bias, seen on the measured scale) — QC views of early
     points should show label content, not f or k_i.
   **Decision (user, 2026-11-12):** rails on label content; no information floor.
4. *The precursor is the fit's one kinetic choice, and k_p is the switch (user, 2026-11-12).* With
   guan/fornasiero removed (5) the fit has one kinetic equation, first-order, so the only decision is
   whether the precursor is at plateau from the start or rising. No mode name: **k_p decides.** It
   comes per experiment from the SDRF `characteristics[precursor rate constant]` (a manifest column,
   read like RIA from `characteristics[precursor enrichment]`), with a `fit --kp` override; absent
   or `inf` = steady precursor = today's θ path, byte-identical. A finite k_p applies the precursor
   inside the per-point solve, fits the peptide on label content (3a) and writes f_i — one switch, so
   nothing in fit can count the precursor twice. One shared value across tissues and conditions
   unless the SDRF says otherwise. The fit logs which path ran and why (SDRF or CLI). `fit --model`
   keeps only `exponential` (rename of `simple`, alias kept) and `calibration` (auto from the SDRF).
5. *Remove guan/fornasiero (user, 2026-11-11, plan 2026-11-12).* Both were written for the SILAC
   side of the SILAC-vs-D₂O comparison (Hammond et al., *MCP* 2022); amino-acid labeling is no
   longer a `--label` option, so they have no native use left. On D₂O they are biased: they fit
   their curves to θ, which is not label content under a rising precursor (report §4). Fornasiero's
   label recycling has not been seen with D₂O; if it ever matters it returns as a precursor curve
   P(s) fed to the same per-point solve, not as a curve model. **Plan:**
   - *Timing.* Deprecate in the release that ships the rising-precursor solve, **not before**:
     without the solve, non-primed D₂O data would have no correction at all, and ignoring the
     precursor biases k further (low, by 25–86 % in the slow regime) than today's guan does (high,
     by 5–42 %). Deprecation = a warning on `--model guan|fornasiero` and `--kp/--kr/--rp` at
     rollup, both removed from the GUI and the docs. **Remove one release later**, with an error
     that points to k_p.
   - *Keep* `models.two_compartment_guan`: it is the label-content curve of the 3a fit (new role).
     Remove `two_compartment_fornasiero`, the guan/fornasiero single-timepoint root-find, the
     rollup's `--kp/--kr/--rp` (its `--kp 0.5` default is unrelated to any SDRF k_p) and the GUI
     spin boxes.
   - *CHANGELOG:* reproducing older guan/fornasiero analyses (e.g. Hammond 2022) needs ≤ 1.1.0.
6. *Provenance.* The fit header records what `fs` holds (`# fs_basis theta|f`) and the precursor
   (`# k_p`, with its source). With guan/fornasiero gone the double-count guard is moot (the only
   curve refit left is exponential, valid on θ and on f); rollup echoes the basis in its own
   header and log.
7. *Rollup.* Two models: `linear` (rename of `linear simple`, alias kept for a release) as the
   **default**, and the `exponential` curve refit. `--method weighted|pooled` stays a separate
   choice (collapse to one point per biological replicate × timepoint, or fit every peptide point).
   Results-affecting: the protein output becomes per-condition rows with Δk columns;
   `--peptide-admission auto` resolves to `all` (fewer peptides per protein in multi-condition
   runs); the GUI default and docs follow; CHANGELOG results-affecting entry. With f input
   (`# fs_basis f`, any finite k_p) the line and the interaction use **`wls`**; rollup errors on
   `wls-var` there (miscalibrated, Monte Carlo above).
8. *Fix in passing.* The fit header prints the CLI `ria_max` default even when the manifest's
   per-experiment enrichment is the value applied.

**Validation before merge.** Type-I of the Δk and interaction tests on f — done (above): `wls`
calibrated, `wls-var` not; the interaction under the chosen weights still to re-check. The report's §6 checks on real data (within-protein agreement,
residual shape over time, Δk and interaction); the simulation benches (k bias ≤ 1 %, coverage) —
**add an early-timepoint baseline-bias arm** (a constant +0.01–0.05 θ offset), which the original
benches lacked and which real data have; the instant-precursor path byte-identical to today; the
gate's yield (3a) rechecked with f in place of θ; the early-point rail decision (3b). **Open
from the prototype:** on identical gated peptides, within-protein k spread was ~10–15 % larger on f
than on θ (both rail variants alike) — a precision cost of the correction, a k_p mismatch, or θ's
late-point-driven k looking tighter than it is; check on the simulation benches and a second
dataset (e.g. against labeling sites). The median residual against time is also a k_p
diagnostic: a mid-series lag that the correction does not remove says the effective precursor is
slower than assumed. A k_p profile on whole envelopes is pulled toward a faster precursor by the
same early baseline bias — exclude the early, low-information points when estimating it.

#### Curation-gate reference (fit & rollup) — be explicit

Three levers are easy to confuse; document them wherever they surface (CLI help + GUI tooltips):

- **`--min-spep`** — *sequence-level admission.* Drop peptidoforms below a labelling-site floor
  before fitting. Primary gate is at `fit` (default on: hw=8, o18=6); a manifest rollup inherits
  it, and rollup re-exposes it for explicit-file inputs.
- **`--min-fit-points`** — *peptide-level biological-replicate floor* (rollup). Keep only
  peptidoforms fit on ≥ N distinct `(biorep, timepoint)` points. **Auto = 2 for a single-timepoint
  experiment, off otherwise.** Exposed in the GUI (Protein tab **"Min fit points"**, `0 = auto`).
- **`--min-points`** — *protein-level refit floor* (rollup; GUI label **"Min refit points"**,
  default 3). Min collapsed `(t, θ)` points for the protein-level refit. **Distinct from
  `--min-fit-points`** — the naming collision is a known UX wart; consider renaming to
  `--min-refit-points` for symmetry (back-compat alias).
- **R² / `k_cv` semantics:**
  - *Multi-timepoint:* R² is the gate, `k_cv` is a **secondary flat-curve rescue** — admit if
    `R² ≥ --min-r2` (0.8) **OR** (`R² ≥ --rescue-r2` (0.6) **AND** `k_cv < --k-cv` (0.2)).
  - *Single-timepoint:* R² is degenerate and **bypassed**, so `k_cv` becomes the **primary**
    gate, alongside the `--min-fit-points` replicate floor. (`k_cv = (ci_hi−ci_lo)/(2·|k|)`, a
    scale-free CV of k̂.) The GUI tooltip currently frames `k_cv` only as an R²-conditional
    rescue — fix it to state the single-timepoint primary-gate role.

#### Locked decisions (2026-06-07) — design invariants, kept for reference

1. **Labeling time is a required SDRF column** (`characteristics[labeling time]`) — the kinetic
   x-axis, deliberately a characteristic (sample-intrinsic), not `factor value[time]`, to avoid
   colliding with a drug-treatment time course. Calibration runs declare `characteristics[mixing
   proportion]` instead; fit dispatches on which column is present.
2. **Header-authoritative identity + a manifest.** Integrate freezes the full identity into each
   output's provenance header; the stage-aware `riana_manifest.tsv` is the index. Fit groups runs
   from the manifest, not by re-reading the SDRF.
3. **DIA-NN parquet intake is a fast-follow** (shipped, M6b); `PSMRecord.retention_time` exists so
   the DIA RT-prior path is designed in.
4. **Protein comparison ships as a side-by-side view first**; the linearized simple-model fit
   (`log(1−θ) = −kt`) is the bridge to the two-sample Δk stats (shipped as `linear simple`).
5. **Per-sample RIA is an SDRF column** (`characteristics[precursor enrichment]`; SDRF → global
   `--ria` default). Load-bearing at the θ-solve — it builds the labelled envelope — so it cannot
   stay a single global scalar once one mzTab spans multiple animals.

#### Track A — I/O & run-identity data model (the spine)

Shipped (CHANGELOG): `RunIdentity` + SDRF/manifest intake (M6a), DIA-NN parquet (M6b), per-run
concurrency + `--resume`, the manifest project chain, the scan↔precursor intake guard, MBR
(mzTab/DDA) with the winner-fraction policy. **Open:**
- **DIA-NN multi-fraction intake** — no data yet; guard + document.

#### Track B — integration & fit fidelity

Shipped (CHANGELOG): adaptive N_ISO at integrate (opt-in `--iso auto`), H4′ mix-then-truncate
FS solve, limited-isotopomer scoring (`--fs` / `--fs auto`), the MS1 peak precache (~8×) and the
`searchsorted` masking vectorization. **Open:**
- **Per-point `Var(θᵢ)` → depth-aware WLS weight — SHIPPED** (opt-in `--linear-weights wls-var`; see `CHANGELOG.md`).
- **Cross-proportion-stable peak picker** (Phase C v2) — the `apex_search_half_width`/`consensus`
  levers for label-invariant boundary stability are untested in production (§2c / M3 carry-over).
- **Early-timepoint findings (2026-11-11, an unpublished in-vivo D₂O DDA dataset; local notes
  only).** Three open items, none decided, all results-affecting if changed:
  - *One-sided rail drop.* Near t = 0 the true fs is ~0 and the per-point noise is wide, so a large
    share of early points falls below the −0.05 rail. Dropping only that side keeps the points that
    read high, and the surviving early fs read high (worst on noisy samples). Options: keep
    −0.1 < fs ≤ −0.05, drop only solver-bound hits, or a noise-aware rail.
  - *Low-intensity light bias.* Weak envelopes read lighter than natural — many pile up at the
    solver's −0.1 bound even with m0–m3 all present, on two Orbitrap models. Likely weak
    heavy-isotope centroids missing from some scans of the integration window. Integration-side;
    candidate fixes: per-channel scan detectability, an apex-centred window, or treating an
    undetected channel as censored rather than zero.
  - *`--fs auto` as the D₂O default?* There it reduced the early-timepoint offset and raised the
    R² ≥ 0.8 yield over full-envelope scoring, with within-protein agreement unchanged after the
    gate. Check on a second dataset before changing the default.
- **Labeled envelopes by convolution (performance; later).** For D₂O, the labeled envelope is the
  natural envelope convolved with Binomial(labeling sites, enrichment) in isotopomer units (¹⁸O: steps
  of 2); it matches correctly binned IsoSpec within 10⁻³ (2026-11-11). It would replace an IsoSpec call
  per (peptide, enrichment) with a binomial and a short convolution — what a rising-precursor envelope
  solve needs, since it mixes many enrichments
  (`reports/2026-11-11_non_first_order_linearization_feasibility.md`). **Accurate masses are not lost
  (checked 2026-11-12):** build the natural envelope E0′ with the labeling-site H fixed at ¹H (IsoSpec,
  once per peptide) and convolve intensity and intensity × mass together, each extra ²H adding
  1.00628 Da. That reproduces IsoSpec's per-channel abundances within ~10⁻⁴ and averaged-isotopolog
  masses exactly on m0–m2 and within 0.01 mDa on any channel above ~3 % abundance; the residual on
  weak tail channels (≤ 0.4 mDa at < 1 %) is IsoSpec's own 99.9 % coverage cutoff. So the mass-defect
  θ (`fs_ds`), which under a rising precursor needs the same kinetic mixture, can use it too.
- **Robust observed-vs-IsoSpec matcher** — soft-trim / per-channel-SNR / Huber downweighting of a
  contaminated high channel; the A4 remainder and the research-grade endgame. Deferred.
- **Genuine per-peptide channel optimizer** — replace the binary `--fs auto` init-width threshold
  with an outlier-aware per-peptide channel choice that balances the bias/variance tradeoff
  (more channels = more signal but more isobaric-contamination risk). Gated on the Track D animal
  benchmark; overlaps the A4 per-peptide scoring-width item.
- **`noise_floor` baseline subtraction** — implemented but off by default (detrimental in every
  test so far); needs exhaustive testing to decide keep-vs-remove, alongside the still-open robust
  in-window baseline (§2c). Re-judge on the *high-isotopomer* probe, not `m0_rmse`.
- **Integrate apex-window science decision** — the `searchsorted` perf half shipped; the open
  question is whether the apex path needs the full-concat window vs. the narrower
  `anchor ± extraction_half_width` (changes results).
- **MS2-level integration for the DIA-NN path** — explore MS2 fragment isotopomers for θ/FS;
  first check how truncated the MS2 envelope is.
- **Other rollup options** — beyond inverse-variance weighting + pooling, e.g. mixture models
  separating peptide-level vs biological-level variance.

#### Track C — fitting / modeling science

Shipped (CHANGELOG): M5 per-timepoint fraction-new, the fit-model set (simple/guan/fornasiero/
calibration), protein rollup, the `linear simple` cross-condition Δk (WLS default), the ¹⁸O
rewrite, M7 PTM-aware envelope, 2D-LC fraction collapse, TMT/TMTpro, dimethyl channel→sample
intake, the multi-point FS rail-drop, and (2026-08-09) multi-factor stratification
(`integrate --experiment-column`), `rollup --plot` / `--experiment`, and per-experiment Δk BH,
and (2026-11-10) the experiment × condition interaction with a random animal effect
(`rollup --interaction`). **Open:**
- **>2-condition Δk — full all-pairwise / Tukey.** The interim selectable pair shipped
  (`--test-condition`/`--reference-condition`, contrasted from the joint all-condition fit). The
  all-pairwise + multiplicity-correction extension builds additively on the same joint fit —
  blocked on a good ≥3-condition dataset.
- **Experiment × condition interaction — remaining.** Shipped 2026-11-10/11 as
  `rollup --interaction` and the GUI Protein tab's *Interaction* option (animal from
  `characteristics[individual]`, Satterthwaite df on the ratio-scale test, singular fits kept;
  `CHANGELOG.md`). Open: *>2 conditions* — two conditions per run today
  (`--reference-condition` / `--test-condition`). **Deferred until a dataset with three or more
  conditions exists to test on**; ties to the all-pairwise item above.
- **[HIGH PRIORITY — DEFERRED] Residual-variance moderation for the Δk t-test (limma-style).**
  A protein with few collapsed points has a noisy residual variance and an unstable t. Moderating
  it across proteins with an augmented df would stabilise them. This is distinct from the shipped
  per-point `wls-var` weights (`reports/2026-07-24_linear_wls_per_point_var.md`). Deferred
  2026-11-11 until the baseline is characterized across complex samples; the plan, the build/skip
  rule and implementation notes are in the *Δk residual-variance moderation* box under §3
  "Deferred beyond 1.2.0". Not built and not prototyped: `bench_linear_weights.py moderate`
  moderates per-point Var(θ), not the residual variance.
- **A1 — dimethyl S=1 spillover gate** (thin primitive + `spillover` column; see the scope block).
- **Deamidation** (deferred; own side project).
- **Slow precursor — envelope solve + fit/rollup model redesign** (designed 2026-11-11;
  productionize after the report §6 checks; see the box under "Deferred beyond 1.2.0").
- **Proper demultiplexing** (deferred; the multiplexing endgame; see the scope block).
- **GG-remnant (UNIMOD:121)** — the most turnover-relevant Tier-2 PTM, blocked on anti-K-ε-GG
  enriched D₂O data (none exists).
- **Cross-fraction RT-correlation MBR** — following a peptide that drifts fractions across
  timepoints; deferred (winner-fraction MBR ships the conservative policy).

#### Track D — validation (within-protein-θ, animal in-vivo)

Shipped/ongoing (CHANGELOG + `memory/`): the frozen integrator-independent within-protein-θ bench
(Hammond-style, no fractional-pool ground truth), the o18↔D₂O head-to-head, and the lve_atr /
lauren in-vivo sets. **Open:** the animal in-vivo benchmark that gates the per-peptide channel
optimizer (Track B) — and, when built, the yield-gap diagnosis follow-up
(`reports/2026-07-01_lauren_yield_gap_diagnosis.md`).

#### Track E — GUI / UX & responsiveness

Shipped (CHANGELOG): sortable tables + PNG export, GUI `-W/--workers`, advanced-knob exposure,
the Δmass/`fs_ds` overlay, isotopomer bar chart, fixed hint area, CLI progress bars, Load-results
on all three tabs, and determinate progress bars. **Open:**
- **[1.2.0 easy win] Package `riana gui` as a standalone app** — PyInstaller / py2app bundle
  (`.app`/`.exe`, real macOS Dock icon). Nearer-term than the Electron rewrite.
- **Protein (rollup) tab: show the project's experiments and conditions before a run** (low priority,
  user-suggested 2026-11-11) — e.g. "4 experiments (LA, LV, RA, RV) × 2 conditions" under the
  manifest field, read from the manifest's integrate rows, so the Reference / Test / Interaction
  choices are made knowing the design.
- **Surface per-run SDRF sample / fraction in the Integrate view.**
- **Faithful-to-smoothing chromatogram trace** — re-apply S-G at click, or an optional
  `--save-traces` `<stem>_riana_traces.parquet` sidecar (also removes the ~10 s re-read).
- **Explore other GUI frameworks** — replace PyQt with a modern Electron app (long-term).

#### Cross-cutting chores

Done: v1.0.0 release + Zenodo, the repo-hygiene audit, Snakefile retirement, manifest
schema-versioning, the parallel≡serial reproducibility test, and **`pytest-xdist` + the `slow`
tier**, `mypy --strict` on `riana/algorithms/`, and the user-facing docs refresh (all shipped on
1.2.0). **Open:**
- **`mypy --strict` beyond `algorithms/`** — extend the scoped CI check package by package.
- **Benchmark CI smoke tier** — a subsampled fast bench tier to catch science regressions
  (e.g. within-protein-θ drift) that unit tests miss.
- **SDRF as a partial config source** — auto-load `comment[modification parameters]` into the fit
  config (mass tolerance stays a separate Riana parameter).

## 4. Cross-cutting recommendations (historical anchors — all shipped in 1.0.0)

All shipped (CHANGELOG): provenance-header reproducibility (**§4.1**), the frozen
`IntegrationConfig`/`FitConfig` single source of truth (**§4.2**), `from __future__ import
annotations` throughout, the repo-hygiene audit, and the `workflow/Snakefile` retirement (Riana
stays a linear `integrate → fit → rollup` chain glued by the manifest — orchestration-agnostic).

## 5. What this project / roadmap deliberately does not include

- **AA / SILAC-as-a-fit-label** — dropped from `riana fit`; Riana fits metabolic-water labelling
  (D₂O / ¹⁸O). (SILAC/dimethyl are handled as sample-axis *multiplexing* intake, not fit labels.)
- **A bundled workflow engine** — the `Snakefile` is retired; SDRF + quantms/DIA-NN own search+ID,
  and Riana is a linear manifest-glued chain. The SDRF/manifest path is canonical; the
  per-timepoint Percolator path is kept for benchmarks/dev and gradually deprecated.
- **Complex / Bayesian kinetic models** — first-order turnover with an instant or a rising
  precursor is the whole model space (guan/fornasiero θ-curve fits to be retired; see the
  slow-precursor box); the focus is honest uncertainty, not more models.
- **Multi-omics integration, a web interface, or a plugin system** — wrong shape for a
  single-developer scientific tool.

## 6. Known limitations

- **Peak fidelity / baseline:** integration uses an apex-centred narrow window; robust in-window
  baseline subtraction and a cross-proportion-stable peak picker are open (§2c / Track B).
- **Linear Δk model — heteroscedasticity FIXED (2026-07-13), one gap remains.** `linear simple`
  now fits by **WLS** (`--linear-weights wls`, default; delta-method weights `(1−θ̂)²` from the
  fitted value); the old OLS was anti-conservative (Type-I ~28% at α=0.05) and biased −14% in the
  fast tail; t=0 excluded (`reports/2026-07-13_linear_model_wls.md`). The ideal weight
  `(1−θ̂)²/Var(θᵢ)` **SHIPPED** as opt-in `--linear-weights wls-var` (`c82cc55`; per-point
  `Var(θ)` from the emitted `fs_var` / `fs_df`). **Still open (purist alternative):** a **joint
  nonlinear fit + Wald contrast** (exact MLE + exact inference) — filed, not built, buys back
  only the last ~3% of efficiency. **Also open:** each protein's Δk test uses its own residual
  variance, which is noisy when a protein has few points (single-timepoint and sparse designs).
  Moderation across proteins is high priority but deferred (§3, the *Δk residual-variance
  moderation* box).
- **Single-timepoint fits** run futile bounded-NLS today and mislead the GUI (R²≈0/NaN); the
  fit-step detection + direct solve (§3 scope item 1) closes this.
- **Memory** is bounded — streaming/indexed mzML, one fraction at a time (`io/mzml.py`).

## 7. Doc & memory hygiene (reconciliation policy)

The 2026-07-21 audit found the planning doc and the `memory/` files had drifted from what actually
shipped (e.g. `searchsorted` and `pytest-xdist` were listed "open" though shipped; a memory index
line still called `dea65f1` "unpushed"; the "adaptive-N_ISO" scope label read as 1.2.0 in one place
and "M8" in memory). **Root cause:** three sources overlap and none is authoritative on *state*.

**Convention going forward — one source of truth per fact type:**
- **`CHANGELOG.md`** — the *shipped* record. Authoritative on what exists.
- **`PROJECT_REVIEW.md`** (this file) — *open* work only. When an item ships, delete it here (keep
  a one-line pointer only if a code anchor references it).
- **`memory/`** — decisions, rationale, and current dev state. On ship, update the owning file's
  status line the same commit.
- **`reports/*.md`** — point-in-time analyses/benchmarks. Never edited after the fact; superseded
  by a newer dated report.

**Immediate cleanups queued:** this rewrite fixes the PROJECT_REVIEW staleness; the `memory/`
files (`v1_2_0_dev_line`, `silac_dimethyl_multiplexing`, the MEMORY.md `dea65f1`-unpushed line, and
the A4 "M8 vs 1.2.0" label) still need a status refresh — do it alongside the single-timepoint work.
