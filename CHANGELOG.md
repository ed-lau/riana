# Changelog

All notable changes to Riana are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

The 1.2.0 development line (branch `1.2.0`).

### Standalone GUI bundle + strict typing for `algorithms/` — 2026-08-08

#### Added

- **`packaging/` — a PyInstaller build for a standalone Riana bundle** (a macOS `.app`,
  or a one-dir bundle elsewhere) so non-Python users can launch the GUI. The frozen
  binary is CLI+GUI (no args → GUI, args → the `riana` CLI). New `[packaging]` extra
  (`pyinstaller>=6`); build with `packaging/build_app.sh`. The build and a CLI smoke
  test are automated; GUI-window launch needs manual verification on a machine with a
  display (see `packaging/README.md`).

#### Changed

- **`mypy --strict` now gates `riana/algorithms/`** — a scoped `[tool.mypy]` config
  (widen `files` to roll out further), a CI `typecheck` job, and a tox `typecheck` env;
  `mypy` added to the `[dev]` extra. The package is annotated to zero strict errors, made
  robust across numpy stub versions (the CI/local numpy differ).

### Fit / rollup / integrate robustness guards (audit tail) — 2026-08-08

#### Fixed

- **Single-timepoint fit at a degenerate time** (`t* ≤ 0` / zero span) returns a null
  result instead of falling through to `curve_fit`, which reported the arbitrary init
  `k = 0.5` with a spuriously zero-width CI.
- **Post-rail-drop depth re-check** counts distinct surviving timepoints, not raw rows, so
  a biological-replicate peptide rail-dropped below the distinct-timepoint floor is no
  longer fit as an under-depth curve.
- **Rollup single-point protein** gets a NaN CI, not a spuriously tight zero-width bootstrap
  CI (`k_cv = 0`) — mirrors the per-peptide `n ≥ 2` guard.
- **`IntegrationConfig`** rejects `smoothing_polyorder ≥ smoothing` up front (Savitzky-Golay
  needs `polyorder < window`) instead of failing deep inside per-run integration.
- **Isobaric channel collapse** verifies a file's channels agree on labeling time /
  enrichment / replicate before merging them (labeling time is the fit x-axis).
- **GUI Protein tab**: the rollup provenance key `k_cv_max` → `k_cv` (+ `min_fit_points` /
  `test_condition`) so GUI- and CLI-produced rollups hash to the same `config_hash`; and the
  single-timepoint hint now displays (it was overwritten before it could show).

### Non-canonical residue intake guard — 2026-08-08

#### Fixed

- **Peptides with a non-canonical residue (selenocysteine `U`, pyrrolysine `O`, or the
  ambiguity codes `B`/`Z`/`J`/`X`) are now dropped at intake** (mzTab, DIA-NN, Percolator)
  instead of being silently mis-massed. Those residues had no defined atom composition —
  `aa_atoms` carried a placeholder `[0,0,0,0,0,0]` for `U`/`X`/`B`, so a `U`-containing
  peptide's neutral mass was ~150 Da too low, driving a wrong extraction window and FS/k.
  The placeholder vectors are removed and a shared `riana.utils.is_canonical_peptide` gate
  skips such PSMs (logged as a count). An audit of the reviewed human/mouse proteomes found
  non-canonical residues in only ~0.12%/0.25% of proteins (almost all selenocysteine) and
  **zero** across ~808k identified PSMs in our data, so the practical impact is nil — this
  closes a silent-wrong-mass hole rather than changing any current result.

### Rollup curation & weighting correctness fixes — 2026-08-03

A code audit of the post-1.1.0 additions surfaced one results-affecting rollup defect and
two narrower weighting/collapse bugs on opt-in paths.

#### Fixed

- **Rollup curation gates are now per-condition (RESULTS-AFFECTING for multi-condition
  rollups).** The `--min-r2` (default 0.8), single-timepoint `k_cv`, and `--min-fit-points`
  admission gates reduced their keep-set to bare `concat` (sequence+charge, identical across
  conditions) and filtered with `.isin`, so a peptide that passed the gate in **one**
  condition was admitted in **every** condition — dragging its bad-condition fit rows into
  the other conditions' protein k and into the two-condition `linear simple` Δk contrast (an
  anti-conservative OR-across-conditions union). Admission is now keyed on the full per-curve
  key `(experiment, condition, concat)`, so a peptide is dropped in exactly the conditions
  where its own fit fails. `--min-spep` is unchanged (Spep is condition-invariant by
  construction). Single-condition rollups are unaffected. New regression test
  `test_rollup_gates_are_per_condition_not_leaked_across_conditions`.

#### Added

- **`rollup --peptide-admission {auto,own,any,all}`** (default `auto`) — makes the
  cross-condition admission policy explicit and auditable. `auto` picks per-model: **`all`
  for `--model "linear simple"`** (the Δk path — a two-condition contrast resting on
  *different* peptides per condition confounds Δk with peptide identity, so the paired
  common-support basis is the right default) and **`own` for the kinetic models**. `own` is
  the per-condition fix above (each condition uses only peptides that pass the gate there);
  `any` restores the old passed-in-one⇒kept-in-all behaviour (opt-in, for back-comparison);
  `all` keeps a peptide only if it passes the gate in **every** condition it appears in
  (same peptides on both sides of the contrast, at the cost of yield). Recorded in the
  rollup provenance header. On `runs/lve_atr_clean` (control-vs-atrium, default gate), ~22%
  of both-condition peptides pass in only one condition; `all` vs `own` shifts the Δk set by
  ~18 proteins with a 0.90 Δk correlation (bulk stable, margins cleaned).
- **`--iso auto` fraction collapse no longer fabricates observed zeros.** The multi-fraction
  sum collapse aggregated iso channels with pandas `sum` (`min_count=0`), so an all-NaN
  adaptive-N_ISO channel pad (a channel a short peptidoform never extracted — "not a channel",
  distinct from an integrated `0.0`) collapsed to `0.0`, which the FS solver then scored as an
  observed-zero channel and biased FS/k low. The iso-channel sum now uses `min_count=1`
  (NaN-preserving); fixed-iso mode is unaffected (its empty channels are a real `0.0`).
- **`--linear-weights wls-var`: a point with no per-point variance now moderates to the
  prior.** A collapsed cell lacking a usable prediction-interval σ (`theta_var` NaN) got a
  constant `1.0` weight divisor while its finite-variance siblings got `1/Ṽar ≈ 1/s0²`
  (~100–1000× heavier), silently near-excluding a real θ measurement from the joint fit. The
  missing-variance fallback is now the pooled prior `s0²` (the df→0 limit of the eBayes
  moderation), putting the point on the same scale as its siblings.

### `linear simple` Δk — opt-in per-point Var(θ) weighting (`wls-var`) — 2026-07-30

#### Added

- **`riana rollup --linear-weights wls-var`** — an opt-in refinement of the WLS Δk fit that
  multiplies the delta-method transform weight `(1−θ̂)²` by the per-point precision `1/Var̂(θ)`,
  nudging the estimator toward the MLE. The weighted rollup collapse now surfaces the per-cell
  `Var(θ)` and its Satterthwaite df as two new columns — **`fs_var` / `fs_df`** — in
  `riana_rollup_fractions.txt`; `fit_linear_deltak` moderates `Var̂(θ)` eBayes-style toward a
  pooled `s0²` with a robust (Winsorized) fitFDist `d0` (fixed-2 fallback), so a noisy per-point
  variance cannot dominate. Falls back to `wls` when the points carry no `theta_var`; `wls`
  (default) and `ols` are byte-identical to before. Improves biological-replicate consistency
  ~10% on rich 9–12-timepoint curves and ~2% on sparse 3-timepoint DIA (the benefit scales with
  curve richness), so it stays opt-in — single-peptide proteins use pooled pseudoreplicates and
  see no gain. See `reports/2026-07-24_linear_wls_per_point_var.md`.

### Single-timepoint fits — direct k solve + fit-step detection — 2026-07-21

Single-timepoint data (TMT, dimethyl, any one-labeling-time run) is a flagship 1.2.0
substrate; the fit now handles it directly instead of burning futile nonlinear least
squares on it. **Numerically identical results — this is perf, correctness, and UX, not a
results change.**

#### Changed

- **`riana fit` auto-detects a single labeling timepoint and relaxes the `--depth` floor to
  1** (with an INFO log) instead of emptying the frame with the generic "no peptides survive
  --q-value / --depth" error. Curation is unchanged (it stays at rollup: k_cv +
  `--min-fit-points`, R² bypassed).
- **k is solved directly at a single timepoint, not by `curve_fit`.** When a peptide's
  surviving points all share one labeling time the 1-parameter fit collapses to
  `model(t*, k) = mean(FS)`: an elementary closed form `k = −ln(1 − F̄S)/t*` for the `simple`
  model (verified to ~1e-12 vs the bounded `curve_fit` optimum, ~100× faster) and a pole-safe
  Brent root-find for guan/fornasiero (strictly monotonic in k at a fixed `k_p`). The
  residual-bootstrap CI is computed the same way on the same RNG stream, so
  `ci_lo`/`ci_hi`/`sd` are unchanged. Applies to the per-peptide fit **and** the
  protein-rollup refit, and to any multi-timepoint peptide whose points collapse to one t
  after FS rail-drop.
- **The single point's `OptimizeWarning` ("Covariance of the parameters could not be
  estimated") is gone** — no `curve_fit` is invoked on the single-timepoint path.
- **R² is reported as `NaN` (not a misleading ≈0) for a single-timepoint fit** — there is no
  time-axis variance for it to explain; the rollup already bypasses R² there and curates on
  k_cv.
- **GUI clarity.** The rollup tab now exposes a **Min fit points** control (the peptide
  biological-replicate floor; `0 = auto` → 2 for single-timepoint, off otherwise) — previously
  tunable only from the CLI — and relabels **Max k_cv** to state its role (a *secondary*
  flat-curve rescue for a time series; the *primary* gate for single-timepoint, where R² is
  bypassed), plus shows a single-timepoint hint when every rolled protein sits at one labeling
  time. The fit tab's Depth tooltip notes the single-timepoint auto-relax, and its summary
  reads "single timepoint — R² N/A, k solved directly" instead of a spurious "0 R²≥0.9".

#### Fixed

- guan/fornasiero single-timepoint fits no longer risk a silent null from the default fit
  init `k=0.5` landing on guan's `k=k_p` pole (`config.k_p=0.5`): the direct root-find
  sidesteps NLS init-sensitivity.

Tests: `tests/test_fitting.py` (closed-form ≡ `curve_fit`, guan root-find, edge cases,
auto-relax + warning-free `fit_run`, replicate-vs-single CI); `tests/test_gui.py`
(`min_fit_points` in `build_params` + `run_rollup` forwarding).

### `linear simple` Δk — weighted least squares (RESULTS-AFFECTING DEFAULT) — 2026-07-13

#### Changed

- **The `linear simple` Δk model now fits by weighted least squares (`--linear-weights wls`,
  the new default) instead of unweighted OLS.** θ carries roughly homoscedastic *measurement*
  noise on the FS scale, but φ = log(1 − θ) is a **log** of it, so by the delta method
  `Var(φ) = σ_θ²/(1 − θ)²` — the φ-residuals are strongly **heteroscedastic**, their SD
  blowing up as θ → 1. Measured on `lve_atr`, the residual SD runs 0.060 → 0.603 across θ bins,
  and regressing `log|resid|` on `−log(1−θ)` gives slope **0.80** [0.77, 0.83] — vs 1.0 for pure
  FS-scale noise and 0.0 for φ-scale, i.e. predominantly FS-scale with a smaller φ-scale floor.
  (An earlier draft quoted 1.10; that included the t=0 points, whose clamped residuals are the
  very artifact removed below.) The old **unweighted** fit treated that 10× SD range as
  equal, which made it **anti-conservative and biased**: Monte-Carlo through RIANA's own
  pipeline puts the false-positive rate of the Δk test at **~28 % at α = 0.05** (nominal 5 %),
  95 % CI coverage at **~54 %**, and k biased **−14 %** in the fast tail (confirmed on 180 real
  `lve_atr` curves against the nonlinear MLE).

  The fit is now weighted by the delta-method inverse variance `(1 − θ)²` taken from the
  **fitted** value (one IRLS step, `w = exp(2·φ̂)`), *not* the observed θ — weighting by the
  observed θ makes each weight a function of that point's own error, which down-weights the
  points noise pushed high and biases k **low**. The weighted estimator is the delta-method
  linearization of the exact MLE (plain nonlinear LS on the FS scale) and recovers **~97 % of
  its efficiency**, while keeping the closed-form joint covariance the Δk contrast needs.
  Result: Type-I **0.285 → 0.062**, coverage **0.543 → 0.926**, bias **−0.105 → +0.003**, and
  the *lowest RMSE at every k*. Because α is now honest it also has more **genuine** power.
- **The t = 0 point is excluded from the linear fit.** In a through-origin model it has **zero
  leverage on the slope** (so k is unchanged), but its θ is pinned by the `theta_floor` clamp —
  true θ(0) = 0, so ~half the measurements go negative and clamp to the floor — which makes its
  residual artificially ≈ 0 against the model's exact 0, deflating the residual variance and
  shrinking **every** standard error. Excluding it only corrects the inference.
- `--linear-weights ols` restores the previous unweighted estimator for audit / method
  comparison. Full workup, including the Monte-Carlo and the real-data validation against the
  nonlinear MLE: `reports/2026-07-13_linear_model_wls.md`.

  **Impact to expect:** on `lve_atr` the significant-protein count moves 321 → 275 and fast-tail
  k's rise ~15 % (the old ones read low). This does **not** overturn biology — the large effects
  survive; what changes is that the p-values no longer overstate the evidence.

### True (sample-axis) multiplexing — dimethyl channel→sample intake — 2026-07-13

#### Added

- **Sample-axis multiplexing intake (dimethyl duplex).** A multiplexing label marks a
  *different sample* that co-elutes but stays separable at MS1 (the label shifts the whole
  precursor), so — unlike isobaric TMT/iTRAQ, whose channels share one MS1 cluster and are
  **merged** into an average — the channels are now **kept as distinct samples**: each is
  integrated at its own m/z, fit at its own precursor enrichment, and combined at rollup as
  a separate condition (so a duplex flows straight into the existing `linear simple`
  two-condition Δk). New **`riana/multiplex.py`** is the label registry: it maps each
  channel to its SDRF `comment[label]` CV term and to the UNIMOD mod a peptidoform carries,
  and declares the site rule + per-residue heavy shift. Dimethyl is wired (light
  `UNIMOD:36` / medium `199` / heavy `330`, +8.0444 Da per site over `S = 1 N-term + #K`
  sites); **SILAC (K+6/R+10) is registered as geometry**, proving the heterogeneous
  per-residue-shift model (K and R shift by different amounts, so the sibling-cluster offset
  is a per-residue sum, not a uniform multiple).
- **Heavy dimethyl `UNIMOD:330` and `UNIMOD:199` (DIMETHYL4) as pinned-isotope
  pseudo-elements** — the same machinery TMT/TMTpro introduced: their built-in ²H/¹³C are
  ~100 % heavy by synthesis, so they are appended as single-isotope (prob 1.0) pseudo-elements
  that shift mass without broadening the envelope. Adds `constants.D_MASS` (²H, IsoSpec's
  built-in value). 36/199/330 are whitelisted so their peptidoforms survive intake, and are
  deliberately **not** in `CHEMICAL_MODS` — a sample-axis mod must keep its channels distinct,
  not merge them onto one curve the way isobaric TMT does.
- `io/sdrf` keeps every channel row of a multiplexed sheet as its own `RunIdentity` (distinct
  sample / condition / precursor enrichment) and exposes a channel-keyed map; `io/mztab`
  routes **each PSM to its channel by the peptidoform's own label mod**; `plan_integration`
  emits **one run per channel** (shared mzML, distinct output). The manifest and fit stages
  needed no schema change — `fit_project` already resolves the RIA per
  `(experiment, condition)` group, so the channels are fit at their own enrichment.
  Validated end-to-end on the dimethyl-D₂O liver duplex; see
  `reports/2026-07-13_dimethyl_duplex_intake.md`.

#### Fixed

- **quantms ≥ 1.8.0 (OpenMS 3.6.0) mzTabs were silently unusable.** 1.8.0 stopped writing the
  *optional* `opt_global_q-value` column, so every PSM fell back to `q = 1.0` — which can never
  pass the strict `q < --q_value` gate (itself capped at 1.0) — and the entire file was
  filtered away with a misleading "relax --q_value" error. The q-value did not disappear: it is
  `search_engine_score[1]`, whose type the metadata declares (`MS:1001491 percolator:Q value` in
  1.8.0; `MS:1003115 OpenMS target-decoy q-value` in 1.7.0). `io.mztab` now prefers the explicit
  column and falls back to `search_engine_score[1]` **only when the metadata declares it to be a
  q-value**, so the FDR semantics stay exact and pre-1.8.0 files are unaffected.

### Multi-point FS rail-drop + physical-margin rails — 2026-07-05

#### Changed

- **FS rail-drop now applies to every fit, single- and multi-timepoint** (was
  single-timepoint only). A per-timepoint fraction-synthesis value that the solver
  returns beyond a small margin of the physical [0, 1] range is a failed solve, not a
  measurement, so it is dropped **before** counting fit points / depth. Re-validated at
  the production `--depth 6` curation on all five turnover sets (boomi/juber iPSC + AC16
  D₂O and ¹⁸O, lve in-vivo): admitted-peptide **yield rises 5–14 %** where rail-hits are
  common and is neutral on the cleanest set, **R² of the fitted population is universally
  cleaner**, the **matched within-protein geom-CV** (same peptides in both arms) is better
  or unchanged on every set, and **median k is unbiased** (|Δk| ≤ 0.001). The earlier
  single-timepoint-only scoping had deferred this pending re-validation; a depth-3 pass
  had shown a spurious CV regression that a depth-6 (production) pass proved was an
  artifact of the loose depth floor. Report: `reports/2026-07-05_multipoint_rail_drop.md`.
- **Rail thresholds moved from the solver clamp to a physical-margin default.** The drop
  rails are now **`FS ≥ 1.05` or `FS ≤ −0.05`** (a 0.05 margin around [0, 1]), replacing
  the old clamp-only rails (`1.199` / `−0.099`, which dropped only points the solver
  *pinned* to its bound). A drop-threshold sweep found `1.05 / −0.05` beats clamp-only on
  R², within-protein CV, and yield across the D₂O sets — it also removes
  solved-but-implausible points — while staying wide enough not to eat a genuine near-1.0
  plateau or a near-0 t₀ anchor (the exact physical bound `1.0 / 0.0` did, cutting AC16
  yield 16 %). Confirmed non-harmful on the ¹⁸O and single-timepoint TMT regimes.

#### Added

- **`FitConfig.fs_rail_drop` (CLI `--fs-rail-drop / --no-fs-rail-drop`, default on)**
  gates the whole rail-drop; `--no-fs-rail-drop` reproduces the pre-1.2.0 multi-timepoint
  numbers (rail-hits then fall to the R² gate downstream). **`FitConfig.fs_rail_hi` /
  `fs_rail_lo`** override the drop thresholds per fit (e.g. the pre-sweep clamp `1.199 /
  −0.099`, or the aggressive physical bound `1.0 / 0.0`). All flow into the provenance
  header via the full-config record.

### TMT / TMTpro + single-timepoint labeling — 2026-07-05

Process an isobaric TMTpro D₂O turnover experiment end-to-end (integrate → fit →
rollup), validated on a 16-plex AC16 single-timepoint dataset.

#### Added

- **TMT / TMTpro peptides are modelled and processed.** Their built-in ¹³C/¹⁵N
  (≈100% heavy by synthesis) are injected as **pinned single-isotope pseudo-elements**
  in the IsoSpec envelope and `unimod_mass`, so the precursor mass and the envelope
  stay in exact lockstep (`unimod_mass(2016)` = 304.2071; the naive light-atoms +
  mass-override that desynced `get_envelope` is avoided). TMTpro (`UNIMOD:2016`) and
  TMT 6/10/11-plex (`737`) tokenise and are **chemical fit-merge mods**: the
  multiplexed samples share one MS1 cluster, so RIANA reports their
  intensity-weighted-average turnover and merges the N-term-labelled vs unlabelled
  peptidoforms onto one curve. New `constants.mod_fixed_isotopes`.
- **Isobaric SDRF intake.** A TMT/iTRAQ SDRF (detected from `comment[label]`)
  collapses its per-channel rows to **one run per data file**; the condition is the
  `|`-joined set of the channels' treatments — so a split-batch design (all-control
  files vs all-treatment files) flows into the linear-simple 2-sample Δk, while a
  pooled plex carries a combined label. Label-free intake is unchanged.
- **Single-timepoint fitting + curation.** A one-labeling-timepoint experiment is
  auto-detected from the data and curated **without R²** (degenerate at a single x):
  the new `rollup --min-fit-points N` biological-replicate gate (auto **2**, tunable)
  plus the `k_cv` relative-uncertainty gate, with the R²/rescue machinery bypassed and
  logged. `fit --depth 1` recovers k analytically from the replicate points.

#### Fixed

- **Single-point fits report NaN uncertainty, not a spurious zero-width CI.** A
  one-point residual bootstrap collapses to `k_cv = 0` (undefined uncertainty read as
  perfect certainty), which sailed through the `--k-cv` gate; now `n_pts < 2` → NaN
  CI/k_cv, so a single-replicate peptide is correctly excluded. General fix, surfaced
  by the single-timepoint data.
- **FS rail-hits are dropped from the single-timepoint fit** before counting fit
  points / depth. The FS solver clamps an unphysical fit to its ±bound (−0.1 / 1.2) as
  a diagnostic; such a point is not a real measurement, and railing to the *same*
  bound in every replicate had manufactured fake "perfectly replicated" peptidoforms
  (identical FS → R²=NaN, k_cv=0). Single-timepoint only for now; the same drop for
  multi-point fits (same criterion) is a re-validation follow-up.

### Rollup — selectable two-condition pair for the linear-simple Δk — 2026-07-03

#### Added

- **`rollup --test-condition` names the Δk comparison condition** (alongside the
  existing `--reference-condition` baseline), so `--model "linear simple"` can
  contrast any **chosen pair** of conditions — even in a project with **more than
  two** conditions, which previously produced a per-condition k for each but no Δk.
  An interim for multi-group projects ahead of full all-pairwise/Tukey: the pair is
  contrasted from the **joint (all-condition) fit** (`core.linear_model`), so it
  reuses the existing model and the eventual multi-group extension is additive rather
  than a rewrite. `--test-condition` requires `--reference-condition` and must name a
  *different* condition (a self-contrast is rejected — it would degenerate into a
  spurious slope-vs-zero test, `delta_k = −k` with `p ≈ 0`, not a Δk); both are
  validated against the conditions present in the data (a typo fails with the
  available choices). The GUI Protein tab's Reference / Test
  inputs are now **dropdowns auto-populated from the manifest's conditions**.
  **Honesty caveat** (surfaced in the CLI `--help` and the GUI tooltip): the joint
  fit still pools the residual variance over *all* conditions in the project, so
  scope the SDRF/project to the conditions you actually mean to compare. Auto mode is
  unchanged — with no `--test-condition`, a Δk is emitted only for a protein that has
  exactly two conditions.

### Integrate — binary-searched isotopomer m/z window — 2026-07-02

#### Performance

- **`integrate` binary-searches the isotopomer m/z window on dense scans.** The
  per-(scan, isotopomer) extraction matched centroids with an O(n)
  `np.abs(mz − target) ≤ delta` mask over every peak in the spectrum — the residual
  per-PSM cost after the MS1-decode precache, and the dominant one on dense
  fractionated Orbitrap runs where `use_range` spans the whole concat scan span. On
  sorted centroid m/z the ±delta window is contiguous, so `_window_sum` now brackets
  it with two `np.searchsorted` calls (O(log n)) once a scan clears a measured
  ~7k-centroid crossover, falling back to the mask on sparse spectra where the mask
  is faster (~1.3–2.7× faster window evaluation at 10–40k centroids; neutral below).
  It is **byte-identical** to the mask either way — the search only narrows the
  candidates (padded one index each side for ULP safety) and the exact predicate
  makes the final selection over a contiguous ascending slice — verified bit-for-bit
  on the sample1 integrate across both `use_range` paths and all output columns.

#### Fixed

- **`IndexedMzML.preload_peaks` verifies each MS1's m/z is non-decreasing** (the
  binary-search extractor's precondition) and raises `DataError` on a pathological
  file, rather than silently under-summing.

### Testing — parallel suite + a `slow` tier — 2026-07-02

#### Changed

- **The test suite runs under `pytest-xdist`.** `pytest-xdist` is added to the
  `[dev]` extra, and the CI test job + `tox` now pass `-n auto` to spread the
  CPU-bound suite (IsoSpec, integration, `curve_fit`, large-mzML decode) across
  cores; `pytest-cov` combines the per-worker coverage. It is opt-in per invocation,
  not in `addopts`, so single-file / `-x` / pdb runs stay serial.
- **The heavy big-mzML / full-pipeline tests are marked `@pytest.mark.slow`** (17
  items — the sample1 / ac16 integrate gates, the CLI / GUI end-to-end runs, and the
  `IndexedMzML` decode checks). `pytest -m "not slow"` skips them for a fast
  inner-loop run — **~24 s (parallel) vs ~5.5 min for the full serial suite** — while
  CI still runs the full tier. The `slow` marker is registered in `pyproject.toml`.

## [1.1.0] — 2026-07-01

The experimental-science line, opened after the 1.0.0 N_ISO finish, plus a
quality-of-life / curation pass. The science: the ¹⁸O (H₂¹⁸O) rewrite (reverse-model
coefficients, production fit, kinetic validation), the mass-defect → θ second turnover
estimate with Δmass/Δspacing QC, modern D₂O labelling-site tables, and the pyteomics 5.x
upgrade. The curation + correctness: a label-aware **Spep gate**, LC-fraction collapse with
winner-fraction MBR and a fraction-aware mass merge, and an internal label-taxonomy
cleanup. The UX: CLI progress bars for `fit`/`rollup`, an isotopomer bar chart and a fixed
hint area in the GUI, the GUI narrowed to the SDRF/manifest path, and a results-display
performance pass so the GUI stays responsive on the large result frames the recent
multi-file series produce. See `PROJECT_REVIEW.md` §3. Entries are grouped by the work that
produced them.

> **Results-affecting defaults (read before upgrading a pipeline):** the default D₂O
> coefficient table changed to `deberneh_2025_rss` (was the 1983 tritium values), a
> Spep curation floor is now applied by default (8 for `hw`/D₂O, 6 for `o18`), and the
> rollup peptide **R² gate now defaults to `--min-r2 0.8`** (was off). Pass the prior
> table explicitly, `--min-spep 0`, and `--min-r2 0` to reproduce 1.0.0 numbers.

### GUI — results responsiveness on large frames — 2026-07-01

GUI results-display items promoted from the 1.1.1 plan into 1.1.0 (worked one at
a time): make the results tables and the chromatogram view usable on the large
frames the multi-file iPSC/cardiac series now produce.

#### Performance

- **The Integrate chromatogram no longer stutters per row-click on a large
  frame.** Selecting a peptide derived its scan span with a
  `df[df["concat"] == x]["scan"]` scan of the *whole* results frame — ~8 ms per
  click at 248k rows, ~160 ms on a millions-row fractionated concat, on the UI
  thread. The `concat → (min_scan, max_scan)` map is now computed **once** per
  result (one vectorised groupby, ~75 ms at 248k rows, off the interactive path)
  and the click is an **O(1) dict lookup**. The span is identical to the old
  scan (a test pins parity, incl. an MBR `scan == -1` row); a peptide absent from
  the map falls back to the extractor's default window.
- **The results tables (Integrate / Model / Protein) no longer lag on large
  frames.** `DataFrameTableModel.data()` — Qt's per-cell, per-repaint hot path —
  read every cell through `DataFrame.iat`, whose pandas scalar-lookup overhead
  (~7 µs/cell) cost **4–13 ms per repaint** for a screenful, so scrolling,
  hovering, and row-selection were visibly janky once the frame was large (a
  multi-file `integrate` concat is 10⁵–10⁶ rows — the 24-file LVE/ATR set alone
  is ~248k rows × 40 cols, and a fractionated series is millions). The model now
  caches each column as a **native-dtype numpy array** once per (re)assignment
  and indexes those in `data()`: **~12× faster repaints** (measured 12.4 → 1.0 ms
  on a 247k-row × 40-col frame; the residual is unavoidable Qt marshalling), with
  the rendered string **byte-identical** to the old path — a test pins this
  across float / int / str / bool, before and after a header sort. The cache is
  per-column so dtypes are preserved (floats stay float for the `:.4g` format);
  there is deliberately **no whole-frame stringify** (that would be ~2.5 s on a
  247k-row frame, seconds-to-minutes on a fractionated one). Building it is
  O(columns) — each `to_numpy()` is a cheap view — so it rides the existing model
  reset for free.

#### Added

- **The Integrate results table has a filter box and a display cap.** A
  multi-file `integrate` concat is 10⁵–10⁶ rows (the 24-file LVE/ATR set is
  ~248k; a 384-file fractionated iPSC series is millions), and handing all of
  them to the view makes sort / selection / memory the bottleneck no matter how
  fast the model is — and nobody scrolls millions of rows to find a peptide. The
  table now shows at most **5,000 rows** at once, with a filter box that narrows
  by **sequence / protein id / concat** (case-insensitive substring) and a
  "showing N of M" note. Only the *view* is capped: the full result is kept in
  memory (and on disk in each run's `_riana.txt`), so sorting, selection, and the
  chromatogram/isotopomer views on the shown rows stay bounded and instant
  regardless of run count. The filter is debounced (250 ms) and matches a single
  precomputed `str.contains` key, so it stays responsive on the millions-row
  frame (~0.6 s/search there, sub-100 ms on the 248k set).

#### Changed

- **The Protein (rollup) tab now takes a `riana_manifest.tsv`, not a fit-output
  directory.** The GUI runs the SDRF / manifest project path only — matching the
  Integrate *SDRF* and Model *Manifest* fields — and the rollup tab was the last
  holdout still asking for a folder. It now locates the fit outputs from the
  manifest's `stage="fit"` rows (`fit_outputs_from_manifest`, the same resolver
  `rollup --manifest` uses), writes `riana_rollup_proteins.txt` /
  `riana_rollup_fractions.txt` next to the manifest, and records the
  `stage="rollup"` rows — so one manifest drives the whole `integrate → fit →
  rollup` chain in the GUI exactly as on the CLI. The fit-directory rollup stays
  CLI-only (`riana rollup <fit_dir>`).

### GUI — display saved results from a manifest / project without recompute — 2026-07-01

#### Added

- **A "Load results" button on the Model and Protein tabs** displays the fit /
  rollup results already saved next to a manifest, without re-running the
  (minutes-long) fit or rollup. When a manifest carrying `stage="fit"` (Model) /
  `stage="rollup"` (Protein) rows is entered, a hint appears ("✓ saved results
  found — Load to view, or Run to recompute") and the button enables. The load is
  **full-fidelity**: it reconstructs the exact in-memory shape the workers return
  — the concat-indexed fit frame with its per-timepoint
  `t`/`fs`/`evidence`/`metox`/`fs_ds`/`dmass`/`dspacing` list-cells, and the
  rollup table plus the per-protein refit `points` — from the `_peptides` /
  `_proteins` scalar summary combined with the `_fractions` substrate, so the
  tables **and** the fitted-curve / Δspacing / φ-space views render identically to
  a fresh run. The curve's kinetic model is read from the output's provenance
  header (only `k_deg` is in the table; the curve shape needs the model), via a
  new `read_provenance_header`. New off-the-Qt-loop workers `load_fit_results` /
  `load_rollup_results` (`riana.gui.tasks`).
- **The Integrate tab loads a prior project from its Output dir.** Because
  `integrate` *creates* the manifest (a manifest *input* would invert the data
  flow), the Output dir doubles as the project locator: point it at a folder that
  already holds a `riana_manifest.tsv` with `stage="integrate"` rows and the same
  Load button displays those runs' concatenated `_riana.txt` outputs — through the
  same filtered + row-capped view a fresh run uses, so a 384-file project stays
  responsive. The isotopomer bars work off the row immediately; the
  chromatogram-on-select works when the mzML folder is also set (`file_idx` → mzML
  is re-resolved from it) and degrades to table + bars otherwise. New worker
  `load_integrate_results`.

### GUI — determinate fit / rollup progress bars — 2026-07-01

#### Added

- **The Model and Protein tabs show a determinate progress bar during a fit /
  rollup**, replacing the indeterminate busy-bar. The work runs in a pool process
  (or a `-W` thread), so it can't touch Qt; it now reports throttled
  `(done, total)` through a `multiprocessing.Manager` queue (picklable across the
  spawn boundary, unlike a raw `mp.Queue`), and a main-thread `QTimer`
  (`riana.gui.progress.ProgressPump`) drains the queue to its latest value and
  updates the bar — so only the timer callback touches Qt. The bar stays busy
  until the first update lands, then goes determinate; a manifest fit over several
  condition curves refills per curve (the core's per-phase counter). It reuses the
  same `progress_callback` the CLI bars already drive — the GUI just wraps it to
  enqueue — so the two surfaces share one progress path. Closes the last GUI
  results-display item promoted into 1.1.0.

### GUI + CLI — a different `-o` on the manifest path forks a derived project — 2026-07-01

#### Changed

- **`fit --manifest` / `rollup --manifest` (and the GUI Model / Protein tabs) now
  *fork* into a different Output dir instead of silently ignoring it.** Before,
  the Output dir was ignored on the manifest path — outputs were always written
  next to the input manifest and that manifest was updated in place, so pointing
  `-o` at a new folder did nothing *and* re-running clobbered the original
  project's fit/rollup. Now a **default / same-folder** `-o` still updates the
  project in place (unchanged), but a **different** `-o` writes the outputs there
  and seeds a **new, self-contained manifest** in that folder with the upstream
  stages' rows (`integrate` for a fit; `integrate` + `fit` for a rollup) carried
  over as **absolute paths** — so the forked folder is a valid project (a later
  `rollup` / GUI "Load results" works on it) that reuses the already-computed
  upstream stages, and the **input manifest is left untouched**. Running a variant
  (different coefficients / model / curation) into a separate `-o` is now safe: a
  pristine input run is never mutated just by choosing a different output folder.
  New `core.pipeline.resolve_manifest_write`; both surfaces call it so they can't
  diverge.
- **A fork refuses to write into a folder that already holds a *different*
  project's manifest**, rather than merging into it. Appending a fork's rows to a
  foreign `riana_manifest.tsv` produced a hybrid manifest (rows from two projects)
  that then failed to load; the fork now errors with a clear message (choose an
  empty output folder) unless the folder is empty or a prior fork of the same
  source. An empty / same-source folder forks as before.

### Rollup — R² gate on by default (`--min-r2 0.8`) — 2026-07-01

#### Changed

- **`rollup --min-r2` now defaults to `0.8`** (was off). The peptide R² admission
  gate is **results-affecting** and needed for good within-protein geometric CV —
  inverse-variance weighting alone under-curates (shown in `reports/`). The
  flat-curve rescue still admits a well-measured low-R² peptide (`R² ≥ --rescue-r2
  0.6` **and** `k_cv < --k-cv 0.2`), so slow/flat curves aren't lost. Pass
  **`--min-r2 0`** (any value ≤ 0) to disable the gate entirely (the prior
  default). The GUI Protein tab's *Min R²* spin defaults to `0.8` to match; `0`
  there still means off.

### Output — full settings in the provenance header — 2026-07-01

#### Changed

- **Every output's provenance header now records the full settings**, one
  `# key value` line per config option, instead of a hand-picked few. The header
  was already stamped with `# riana` / `# git` / `# config_hash` / `# id_source`
  plus a small `extra` set (model, label, …); it now also writes the complete
  config dict that `config_hash` is computed from — so `riana_fit_peptides.txt`,
  `riana_rollup_proteins.txt`, and each `_riana.txt` are fully self-documenting and
  a run is reconstructable by reading the header, not just by matching the hash.
  Values are flattened to one line each; `extra` still overrides a shared key
  (e.g. the resolved coefficients path). Readers already skip the header
  (`comment="#"`), so parsers and the golden/parity tests are unaffected.

### Rollup — scale-free relative-uncertainty curation gate — 2026-07-01

#### Changed

- **The rollup curation admit is now a scale-free relative-uncertainty rescue,
  replacing the absolute `--alt-k` / `--alt-se` thresholds.** The `--min-r2` gate
  admits a peptide if `R² ≥ min_r2` OR (flat-curve rescue) `R² ≥ --rescue-r2 AND
  k_cv < --k-cv`, where `k_cv = (ci_hi − ci_lo) / (2·|k|)` is the rate constant's
  **relative uncertainty** — a scale-free coefficient of variation of k̂. This
  rescues well-measured but *flat*-curve peptides (slow turnover / low dynamic
  range) whose R² is pathologically low even for a good fit — the issue Lau
  *Nat Commun* 2018 (and Sadygov's d2ome) address by gating on the rate
  constant's confidence interval. Unlike the old absolute `k`/SE thresholds,
  `k_cv` needs **no retuning** across time-series ranges or k units (/day vs /h),
  which is why `--alt-k` / `--alt-se` were confusing to set. Defaults: `--k-cv 0.2`
  (set ≤ 0 to disable the rescue → R²-only gate), `--rescue-r2 0.6`.
  - The `--rescue-r2` floor is **load-bearing, not cosmetic**: without it the
    relative-uncertainty gate admits degenerate k≈0 rail-hits whose bootstrap CI
    collapses to a spuriously tight `k_cv ≈ 0` (with deeply negative R²). On
    noisier / short-window data that population is large and it *wrecks*
    protein-level ranking (ac16 D₂O↔¹⁸O ρ 0.61→0.43, lauren 0.61→0.53); any floor
    above the negative-R² band restores it, and 0.6 preserves ranking best while
    keeping within-protein geom-CV clean (reports 2026-06-26 / 2026-07-01).
  - **API change** (removes two options added in 1.0.0). The gate is niche — it
    only fires under `--min-r2`, which is **off by default** — so the default
    pipeline output is unchanged; only runs that set `--min-r2` *and* relied on
    `--alt-k` / `--alt-se` are affected, and there only the rescued subset shifts.

#### Added

- **`k_cv` reference column in `riana_fit_peptides.txt`** — the per-peptide
  relative uncertainty of k̂, emitted alongside `ci_lo` / `ci_hi` so it can be
  curated on directly (it is what the rollup `--k-cv` gate computes internally).
- **GUI Protein tab: `Max k_cv` and `Rescue R² floor` spinboxes** surface the new
  rescue next to `Min R²` (the misleading "k ≤ 0.025 & SE ≤ 0.05" hint is gone).

### Integrate — mass-based intake guard + large-run robustness/throughput — 2026-06-30

#### Changed

- **The intake guard is now mass-based (scan↔precursor), replacing the scan↔RT
  guard.** Each run verifies that a sample of mzTab `spectra_ref` scans point at the
  matching **precursor m/z** in this mzML, instead of reconciling the mzTab-reported
  retention time. The RT guard false-positived on legitimately **OpenMS-aligned**
  runs (median offsets of a few minutes while the scans were correct), blocking real
  data; the precursor check is immune to RT alignment (it compares mass) and still
  catches a wrong mzML↔mzTab pairing / quantms filename-prefix scramble. Config
  `check_scan_rt` → `check_scan_id`, `scan_rt_tol_min` (min) → `scan_precursor_tol_ppm`
  (default **10 ppm**); CLI `--no-rt-check`/`--scan-rt-tol` → `--no-id-check`/`--precursor-tol-ppm`.

#### Fixed

- **One bad run no longer freezes the whole batch.** A per-run exception used to
  escape the `ProcessPoolExecutor` `with` block, whose `shutdown(wait=True)` then
  blocked on *every* other submitted task before surfacing — an indefinite hang (no
  output) on a many-file run from a single failing file. Failures are now logged and
  skipped, so the rest of the batch completes.
- **Provenance git-SHA computed once, not per output file.** `make_provenance` forked
  `git` for every file's header; on macOS, forking from the multi-threaded pool main
  process intermittently deadlocked in the child's `pthread_atfork` handlers and froze
  the run. Now cached and warmed single-threaded before the pool spawns.
- **Empty scans no longer crash a run.** A zero-peak spectrum (`defaultArrayLength=0`
  — rare but real in crash-recovered / some ProteomeXchange files) made peak access
  raise `KeyError: 'm/z array'` (surfaced by the precache, which decodes every MS1);
  empty scans now yield empty arrays — zero signal — instead.

#### Performance

- **Per-run MS1 peak precache** — each MS1 spectrum is decoded once per run rather than
  re-decoded for every overlapping per-PSM RT window (~8× faster `integrate` on dense
  runs; output is numerically identical).

### GUI — fixed hint area + tooltip audit — 2026-06-28

#### Added

- **A fixed hint area in the status bar** mirrors the tooltip of whatever control
  the mouse or keyboard focus is on, so help is visible immediately instead of only
  after a hover-hold. An app-wide event filter reads each widget's existing
  `toolTip()` (walking to the nearest ancestor that has one), so there is nothing to
  keep in sync.
- **Tooltip coverage audit** — added help text to 24 previously-bare form controls
  across the Integrate / Model / Protein tabs (mzML & search-ID & SDRF & output
  paths, q-value, window anchor, integration ½-width, manifest, coefficients, label,
  depth, RIA, parsimony, min-peptides/points, and the MBR / intake-guard dials), so
  the main path and the common knobs all explain themselves in the hint area.

### GUI integrate — isotopomer abundance bar chart — 2026-06-28

#### Added

- **The Integrate tab shows a relative isotopomer (m0..mN) bar chart** beside the
  chromatogram. Selecting a peptide row now renders both the RT-domain trace and
  the abundance-domain envelope (the integrated `isoN` areas, normalised to sum 1,
  palette-matched so m{i} is the same colour in both views). It reads straight off
  the results row — no mzML round-trip — so it updates instantly, and has its own
  PNG export.

### Progress bars for `fit` / `rollup` — 2026-06-28

#### Added

- **Live progress for the long CLI loops.** `riana fit` and `riana rollup` now show a
  progress bar over peptidoforms / protein groups (previously they ran silent for
  minutes with only a start/end line). A dependency-free renderer
  (`riana.progress.ProgressReporter`): a single in-place carriage-return bar on a
  TTY, or a clean log line every 10 % when piped / redirected (so logfiles and CI
  output stay readable). Covers both fit paths (manifest curves and explicit files)
  and both rollup models (weighted/ODE refit and `linear simple` Δk). The core
  functions take an optional `progress_callback(done, total)`; `0`-cost when unset.
  (`integrate` already reports per-run; the GUI Integrate tab already has a bar.)

#### Changed

- **The GUI now runs the SDRF / manifest path only.** The Integrate tab requires an
  SDRF (the search-ID file is the quantms mzTab / DIA-NN parquet) and the Model tab
  requires a `riana_manifest.tsv` — the bare-Percolator integrate (psms.txt + sample
  name) and the explicit-timepoint-file fit are removed from the UI. Those paths
  reach all of Riana's current features (per-run identity, manifest chaining,
  fraction collapse, per-experiment RIA) only via the SDRF, and the explicit-file
  fit pseudo-replicates fractions; GUI users are all on the SDRF path. Both legacy
  intakes remain **CLI-only** (`riana integrate <mzml> <psms>`, `riana fit a.txt …`)
  for dev/testing and bench scripts. Removed the Integrate tab's *Sample* field and
  the Model tab's *Timepoint files* list.

### Spep curation gate (`--min-spep`, label-aware default) — 2026-06-28

#### Added

- **`fit --min-spep N` / `rollup --min-spep N`** — a curation floor on a
  peptidoform's labelling-site count (Spep). Peptidoforms below it are dropped
  **before fitting**, so under-powered curves — too few sites for the isotopomer
  envelope to shift measurably as FS goes 0 → 1 — never reach the fit results,
  metrics, or rollup. It complements the R² gate, which can't catch a noise-driven
  high-R² low-site fit. Exposed on the GUI Model tab (`Min Spep`, "auto" = label
  default). `0` disables it.

#### Changed

- **A default Spep floor is now applied** (it was off): **8 for `hw`/D₂O, 6 for
  `o18`** (`FitConfig.min_spep` resolves `None` → the label-aware default). The
  ¹⁸O floor is lower because its +2 Da-per-site shift makes fewer sites
  measurable. This is **results-affecting** — a fit with the default now omits
  the lowest-Spep peptidoforms (typically a small tail; the gate is a weak lever
  once `--fs`/R² curation is applied). Pass `--min-spep 0` for the prior
  unfiltered behaviour, or tune per sample (the floor rises with shorter time
  series and slower turnover). Defaults derived in
  `reports/2026-06-28_spep_curation_gate.md`.

### LC-fraction collapse policy + fraction-aware mass merge — 2026-06-27

#### Added

- **`riana fit --fraction-collapse sum|anchor`.** Selects how LC fractions /
  technical replicates of the same `(peptidoform, charge, biological replicate,
  labeling time)` — identified by the SDRF `comment[fraction identifier]` — are
  combined into one kinetic point before fitting (manifest path). `sum` (default)
  sums each `isoN` channel across fractions; `anchor` keeps only the single
  highest-total-intensity fraction (legacy parity). In both cases the intensities
  are combined **before a single FS is solved** — fractions are never fit as
  independent points and their FS values are never averaged.

#### Fixed

- **The fraction merge now intensity-weights the per-channel mass / QC columns.**
  `_merge_fractions` previously summed the `isoN` intensities but took the *first*
  fraction's `iso{N}_obs_mz` / `iso{N}_ppm_error` / `apex_snr`, so the summed
  envelope carried one arbitrary fraction's masses — which would corrupt the
  mass-defect (`fs_ds`) estimate on fractionated data. Mass/error columns are now
  intensity-weighted by their channel, `apex_snr` by the row's total intensity, and
  `n_scans` takes the max. (Latent since the M6a pre-wiring; only became a real bug
  once `fs_ds` shipped this line.)
- **The explicit-files fit path warns on pseudo-replication.** `riana fit a.txt
  b.txt …` (no `--manifest`) does not collapse fractions; it now warns when rows
  share a `(peptide, sample)` and points to the SDRF `--manifest` path.

#### Changed

- **MBR is restricted to each precursor's winner fraction.** Across the runs of a
  curve, match-between-runs now fills a peptide's holes only in the **single LC
  fraction where it has the most identifications** (ties broken by best q-value,
  then fraction number), instead of independently in every fraction it touches.
  This is the conservative minimal policy — recover few bona-fide signals where a
  peptide reliably elutes, rather than maximise transfers — and it deliberately does
  not model cross-fraction RT drift (a deferred refinement). A structural no-op for
  single-fraction data.

> Validated against the real fractionated iPSC D₂O SDRF (192 files = 12 timepoints
> × 2 bioreps × 8 fractions); end-to-end fit validation awaits the mzML/quantms
> search. Cross-fraction RT-correlation MBR and DIA-NN multi-fraction intake remain
> follow-ups (no data yet).

### Internal label taxonomy → string labels (`"D2O"` / `"O18"`) — 2026-06-27

#### Changed

- **The internal labeling-chemistry selector is now a string.** `get_peptide_distribution`,
  the Δspacing solver `solve_fs_d2o_ds`, the fit's `_fs_ds_points`, and `core/fsynthesis`
  take `label="D2O"` / `"O18"` (and the legacy `"AA"`) in place of the opaque integers
  `{1 = ²H in-vivo, 2 = ²H in-vitro, 3 = ¹⁸O, 4 = AA}`. The old in-vivo/in-vitro 1-vs-2
  split is **retired** — D₂O cell-specificity lives in the fit's coefficient table, not the
  label. The user-facing CLI flag (`--label hw|o18`) is unchanged; `hw` maps to `"D2O"` and
  `o18` to `"O18"` at the fit boundary.

### ¹⁸O kinetic fit + RIA-from-manifest + coefficient unification — 2026-06-26

#### Added

- **`solve_fs_o18_ds` — the ¹⁸O Δspacing `fs_ds`.** The ¹⁸O analog of the D₂O mass-defect
  estimate, delegating to the shared nonlinear spacing core (`label="O18"`). Because ¹⁸O is
  a +2 Da label, only iso1 is flat; it scores iso2–4 (iso3 = ¹⁸O + ¹³C carries signal). This
  is a **completeness/symmetry estimate, not a reliable second estimate** — the ¹⁸O-vs-¹³C
  mass-defect difference (~2.5 mDa) is ≈½ of D₂O's per-mass-unit, so per-peptide IQR spans
  the bounds. The intensity `solve_fs_o18` stays primary. The GUI `fs_ds` overlay now uses
  the correct ¹⁸O model on ¹⁸O fits.

#### Fixed

- **Per-curve RIA is resolved from the manifest.** `fit_project` now reads each curve's
  `precursor_enrichment` from `riana_manifest.tsv` instead of falling back to the 0.06
  default, so series at a different enrichment (e.g. the iPSC ¹⁸O series at RIA 0.0897) fit
  at the right precursor RIA.

#### Changed

- **¹⁸O kinetic fit validated on AC16 + iPSC D₂O-vs-¹⁸O time series — no reverse-model change
  required.** ¹⁸O ≈ D₂O on the curated set (iPSC Spearman(k) 0.66 at peptide and protein
  level; within-protein robust geomCV ~0.15, identical between labels) and out-curates D₂O at
  iPSC. Recommended curation: R² ≥ 0.8, Spep ≥ 5, depth 6. `--fs auto` is **wrong for ¹⁸O**
  (it drops iso4) — use the default full envelope. (Report `reports/2026-06-26_o18_kinetic_fit.md`.)
- **D₂O and ¹⁸O coefficient tables unified on a bootstrap-OOB freeze.** Both table families
  now report coefficient = bootstrap mean and R² = out-of-bag (CSV gains `oob_r2, ci_lo,
  ci_hi, boot_frac_nonzero`); the ¹⁸O table moved off the single 80/20 split (numbers move
  negligibly — OOB R² 0.892). The in-vitro D₂O table was retrained on the current default
  integration.
- **¹⁸O coefficient presets renamed for provenance:** `o18_ac16` → `juber_2026_o18_ac16`
  (in-vitro), `o18_previs` → `rachdaoui_2009_o18` (in-vivo mouse).

### D₂O labelling-site coefficient tables — 2026-06-25

#### Added

- **Two modern LC-MS-derived per-AA D₂O tables.** `ilchenko_2019` (Ilchenko/Sadygov 2019,
  Table 2 N_aa) and `deberneh_2025_rss` (Deberneh 2025, Table S1 RSS — the authors'
  recommended method). Both tighten within-protein k agreement ~18–21% vs the 1983 tritium
  values (median robust geometric k-CV 0.199 → 0.158 / 0.162 on the LVE manifest) — strong
  evidence the 1983 values are suboptimal for LC-MS D₂O. (Report
  `reports/2026-06-25_d2o_coefficient_tables.md`.)

#### Changed

- **`deberneh_2025_rss` is the new default coefficient table** (was the 1983 tritium values).
  The per-cell-line presets are renamed `ac16`/`ipsc`/`cm` → `alamillo_2025_ac16` /
  `alamillo_2025_ipsc` / `alamillo_2025_cm`, and the legacy `commerford` → `commerford_1983`,
  for provenance clarity. The new tables do **not** move the fs↔fs_ds residual — that residual
  is not the coefficient table.

### Mass-defect θ / Δspacing QC + calibration model — 2026-06-25

#### Added

- **`fs_ds` — a drift-robust second turnover estimate from the mass-defect (Δspacing).** Each
  neutromer's accurate-mass shift (the DeuteRater signal, Naylor/Price 2017) is inverted back
  to fraction-new by the **nonlinear init↔final mixture-spacing curve** (`solve_fs_d2o_ds`),
  replacing the linear `ΔSₓ/ΔSₓmax` ratio that over-read mid-range. Per-timepoint, M0-internal
  (global m/z drift cancels), empirically t0/f0-anchored, weighted median + MAD over iso0–3.
  It is a **cross-check, never a replacement** for the intensity FS (it is ~2.4× noisier per
  point); its *disagreement* with the intensity FS is the useful product. New `fs_ds` column.
  (Report `reports/2026-06-25_mass_defect_theta.md`.)
- **GUI Δmass / Δspacing QC.** The Model tab gains a **Fit / Δspacing / Δmass** toggle and an
  anchor checkbox, overlaying `fs_ds` on the Fit graph with agreement statistics.
- **A `calibration` fit model.** A through-origin FS-vs-mixing-proportion recovery line
  (slope → k_deg, R² = recovery quality) for mixing-series calibration runs, auto-dispatched
  in `fit_project` when the SDRF declares `characteristics[mixing proportion]`; GUI 1:1 plot.

### ¹⁸O (H₂¹⁸O) reverse model + production fit — 2026-06-24

#### Added

- **`riana fit --label o18` works end-to-end.** The ¹⁸O envelope is a 3-isotope-O shift
  (¹⁶O/¹⁷O/¹⁸O; the enrichment dilutes the natural pool and lifts ¹⁸O by the labeling
  fraction), matching the NB90c reverse-model oracle to 5e-5. The D₂O path is untouched.
- **A length-model ¹⁸O coefficient structure** — `Spep = b·(L−1) + c_D·D + c_E·E + c_N·N +
  c_Q·Q + c_S·S` (6-param). Serine added to NB90c's DENQ (+3.4 pt held-out R², 21σ; T/Y null);
  backbone `L−1` confirmed by the data and by Previs (MCP 2009 — one ¹⁸O per peptide bond, the
  terminal residue's O back-exchanges in tryptic digest). Trained on AC16 (2083 peptides).
  Two presets bundled: the in-vitro AC16 table and the in-vivo mouse reference (reproduces
  Previs' worked example LGEYGFQNAILVR = 16). (Report `reports/2026-06-24_o18_reverse_model.md`.)

### Pyteomics 5.x upgrade — 2026-06-24

#### Changed

- **Pyteomics unpinned to `>=5,<6`** (the 1.0.0 release pinned `<5`), with explicit `psims`
  and `lxml` dependencies. Verified 233-test parity vs 4.7.5 (sample1 golden to 1e-3). 5.0's
  `map()` multithreading was evaluated and not adopted — it is a full-sweep primitive that
  doesn't fit integrate's targeted random access, and file-level `-W` already covers the
  GIL-bound work.

## [1.0.0] — 2026-06-24

The breaking 1.0 release (GitHub release + Zenodo DOI): a new package structure,
peak detection, baseline subtraction, mzTab + DIA-NN intake, protein rollup, and
a Qt GUI, closing with the Track B N_ISO line — the standing calibration
benchmark driver, init-width-keyed `--fs auto` widening, adaptive `--iso`, and
GUI parity — plus DIA-NN phospho proteoforms and the `lxml` / `pyteomics<5`
dependency fixes. `v1.0.0` is tagged at the N_ISO finish. See `PROJECT_REVIEW.md`
§3. Entries below are grouped by the work that produced them.

### DIA-NN phospho proteoform sites (M7 Stage B, DIA path) — 2026-06-23

#### Added

- **DIA-NN phosphopeptidoforms roll up as distinct proteoforms.** `io.diann` now
  maps DIA-NN's localized `Protein.Sites` to the biological-mod proteoform suffix
  (e.g. `P35486_pS293`, `_pS1332_pS1333` for two sites) — **byte-identical to the
  DDA/mzTab key format** — gated on `PTM.Site.Confidence` (default ≥ 0.75, the
  class-I localized cutoff; sub-threshold folds into the bare protein). DIA-NN
  lists *every* modified site, so the constitutive fixed Carbamidomethyl (C) and
  Met-Ox (M, a chemical mod merged at fit) are excluded — only
  `constants.BIOLOGICAL_MODS` define a proteoform. The PTM columns are read only
  when present, so a no-mod / older report still loads (empty `mod_sites`). New
  `read_diann(min_site_confidence=…)`. Validated on the cardiac variable-phospho
  DIA series (577 sited proteoforms from 789 phospho precursors).

### GUI exposure of `--iso auto` + `--fs` scoring (Track B / Track E) — 2026-06-23

#### Added

- **GUI parity for the Track B knobs.** The Integrate tab's *Isotopomers* field
  now accepts `auto` (adaptive N_ISO) and a single index `N`, with a new
  *Precursor enrichment* (RIA) spin that shapes the adaptive final envelope; the
  Model tab gains an *FS scoring* dropdown (`full envelope` / `auto` / `iso0-N`)
  wiring `FitConfig.score_channels` / `fs_auto`. Both build the same frozen config
  the CLI does (GUI↔CLI parity before the 1.0.0 release).

### Per-peptide `--fs auto` widening + CLI single-int `--fs`/`--iso` (Track B) — 2026-06-23

#### Added

- **`riana fit --fs auto` — per-peptide limited-isotopomer widening.** Instead of
  a flat channel count, the fit picks per peptidoform: iso0-3 for typical peptides,
  but **widens to all captured channels for peptides whose natural-abundance (θ=0)
  envelope is broad** (`init_envelope_width ≥ 6`), where iso4-5 carry clean,
  model-predicted signal that flat iso0-3 would truncate. The criterion is the
  **natural envelope width — purely compositional, RIA-invariant** — so the
  threshold is a named constant (`core.fitting.FS_AUTO_*`), not a user dial.
  Derived on the calibration mixing series (crossover = init width 6, identical on
  ac16/cm/ipsc) as a strict Pareto win over flat iso0-3, and ties an RIA-dependent
  N_ISO-keyed alternative on the headline metric while needing no per-experiment
  retuning (report `2026-06-23_adaptive_niso_limited_isotopomer.md`;
  `bench_niso_crossover.py`). Free at fit — reuses the cached init envelope
  `solve_fs_d2o` already builds. New `algorithms.isotope_dist.init_envelope_width`
  + `FitConfig.fs_auto`.

#### Changed

- **`--fs` and `--iso` accept a single index `N`** (= the leading channels
  iso0..isoN), the easy form replacing the leading-contiguous list: `--fs 3` =
  iso0-iso3, `--iso 5` = the m0-m5 envelope (now the `--iso` default). An explicit
  list is still accepted for a non-contiguous set (the o18 `0 6` pair), and `--fs`
  also takes `auto`. Fixes a latent bug where the `--fs` *list* form
  (`--fs 0 1 2 3`) was always rejected (tuple-vs-list comparison).

#### Benchmarks / tooling

- **`tests/benchmark/bench_niso_crossover.py` — crossover derivation.** Per-integer
  N_ISO *and* init-width MAE strata + the heuristic A/B (flat / niso / init-width
  binary / graded) on the fixed v1.0.0 integrate; the tool behind the `--fs auto`
  decision.

### Standing calibration benchmark driver (Track B) — 2026-06-23

#### Benchmarks / tooling

- **`tests/benchmark/run_calibration_benchmark.py` — one-command calibration
  A/B harness.** Thin orchestrator over `run_integrate_v1_0_0.run_line` +
  `bench_fs_method_compare`: per cell line it integrates (or reuses a complete
  `integrate_outputs/<label>/`), scores |θ−f| recovery vs the ground-truth
  mixing proportion, writes the standing `runs/calib_<line>/<label>/`
  (`recovery/`, `config.json`) layout, and appends a frozen line to
  `BASELINE.md`. The default config + `--fs 0 1 2 3` reproduces the committed
  `benchmark_results/<line>/v1.0.0_fs0123/` anchor, so integration-knob / MBR
  tuning is an A/B against a recorded baseline. Design:
  `2026-06-23_calibration_benchmark_harness.md`.

### Adaptive N_ISO, H4′ FS solve, and `--fs` limited-isotopomer scoring (Track B) — 2026-06-23

#### Added

- **`riana fit --fs 0 1 2 3` — limited-isotopomer scoring.** Fit the per-timepoint
  fractional synthesis on a *leading subset* of isotopomer channels (e.g. iso0-iso3)
  rather than the full integrated envelope, to reduce sensitivity to co-eluting
  contaminants in the high channels — "integrate wide, fit narrow" (Sadygov & Currie,
  JPR 2025). On the D₂O calibration mixing series (ac16/cm/ipsc) this **tightens FS
  recovery** (within-±0.05 +1.8–2.8 pp at every mixing proportion, lower IQR and bias)
  and is a **pure improvement over the 0.9.0 baseline**, which it leaves byte-identical
  when unused. `FitConfig.score_channels`. Two guards: a **run-level** check that the
  requested channels were actually integrated (so `integrate --iso 0 1 2 3` +
  `fit --fs 0 1 2 3` is valid, and a too-narrow integrate errors clearly), and a
  **per-peptide clamp** (a peptidoform whose envelope ends before the subset is scored
  on the channels it has; a 1-channel peptidoform returns NaN rather than crashing).
- **`riana integrate --iso auto` — adaptive N_ISO (opt-in).** Runs the IsoSpec forward
  model per peptidoform at integrate time and sets the isotopomer channel set from the
  envelope (init∪final ≥1 % abundance, conservative Commerford upper-bound labelling
  sites), replacing the fixed `--iso` tuple; `iso0` is the precursor m0 and the
  extraction target is the averaged-isotopolog accurate mass. **Off by default** —
  evaluated on LVE and the calibration series and found neutral-to-slightly-negative on
  recovery/θ-spread and ~3–5× slower to integrate, so it stays opt-in (the science win
  is the fit-side `--fs` lever, not wide capture). `--ria` sets the precursor enrichment
  (also read from the SDRF `characteristics[precursor enrichment]`).

#### Changed

- **`solve_fs_d2o` mixes the init/final envelopes in the full-cluster basis before
  truncating to the scoring channels (the "H4′" order).** The previous
  normalize-each-then-mix order is correct only when the natural and labelled envelopes
  carry the same in-window mass fraction; the new order is exact under any truncation,
  which `--fs` limited-isotopomer scoring requires. A near-no-op at low labelling (LVE
  RIA 4.6 %: within ~1e-3), load-bearing for narrow scoring and high θ.

#### Benchmarks / tooling

- `bench_fs_method_compare.py` and `bench_within_protein_theta.py` gained
  `--score-channels`; `run_integrate_v1_0_0.py` gained `--adaptive` / `--ria`. New
  committed anchor `benchmark_results/<line>/v1.0.0_fs0123/` (recovery at the recommended
  `--fs 0 1 2 3`). Reports: `2026-06-23_adaptive_niso_limited_isotopomer.md` (the
  capture-vs-scoring investigation + no-regression-vs-0.9.0), `2026-06-23_robust_turnover_cv.md`
  (derivation of the `1.4826·MAD(ln k)` within-protein k-CV),
  `2026-06-23_calibration_benchmark_harness.md` (standing-benchmark design).

### Sortable result tables + graph export (Track E easy wins) — 2026-06-21

#### Added

- **Sortable result tables on all three GUI tabs** (Integrate / Model / Protein).
  `DataFrameTableModel.sort()` reorders the backing frame in pandas (vectorised,
  C-level) on a header click, and the tables set `setSortingEnabled(True)`. The
  sort is done in the model over the *raw* column, not via a
  `QSortFilterProxyModel` over the formatted cells, for two reasons: numeric
  columns sort numerically (a proxy would compare the `:.4g` display strings
  lexicographically — `"100" < "9"`), and the row-click handlers keep mapping a
  view row straight through `dataframe.iloc[row]` with no `mapToSource`
  translation. The pandas sort also scales to large result frames where a proxy
  doing per-cell Python comparisons would not — the "large-table responsiveness"
  half of the easy-win pair.
- **Save graph… (PNG) on the chromatogram and fitted-curve views.** A button under
  each embedded plot exports the live `PlotItem` via pyqtgraph's `ImageExporter`
  (`riana/gui/export.py`), replacing the removed `--plotcurves`. The button is
  disabled in the placeholder state and enabled once a peptide/protein is plotted.
  Raster only: pyqtgraph's `SVGExporter` throws on plots carrying scatter symbols
  (the observed-point / fold series these views always draw), so an SVG option
  would crash on the common case.
- **Protein + condition in the Model-tab curve header.** `CurveView.plot_fit` now
  takes optional `protein`/`condition` and joins them with the peptide `concat` in
  the plot title (`LSLIIR_2 • sp|… • control`), so an exported figure self-identifies
  which protein and — for a multi-condition manifest fit — which condition group it
  shows. The export filename picks up the condition too.
- **Legend moved outside the data area.** The embedded curve + chromatogram views
  now host the pyqtgraph legend in its own right-hand column
  (`riana/gui/plotting.plot_with_external_legend`) instead of anchored inside the
  ViewBox, where its MBR ▲ / fold ◇ sample glyphs were easily mistaken for plotted
  data points. The PlotItem's own legend is reparented, so `plot(name=…)` still
  auto-populates it; graph export now writes the whole **scene** so the external
  legend is still captured.
- **Application icon.** `riana gui` now sets a window/app icon from a packaged
  512 px PNG (`riana/gui/resources/riana.png`, loaded via `importlib.resources`).
  On macOS the Dock icon of an unbundled `python` process is owned by the launcher,
  so this drives the window icon; a true Dock icon still needs an `.app` bundle.

#### Fixed

- **Model tab plotted the wrong condition group for multi-condition (manifest)
  fits.** Pre-rollup, a manifest fit produces one result row per `(experiment,
  condition)` for each peptide; the Model-tab table showed them all but
  `_on_row_changed` looked the row back up by `concat` and took `.iloc[0]` — always
  the first (alphabetically-first) group — so the fitted-curve view, its point
  census (observed / MBR / folded), and the `k_deg`/CI could describe a *different*
  condition than the `n_points`/`n_mbr`/`n_metox`/`n_clean` shown on the selected
  row. The selection now matches on the row's `(experiment, condition)` keys
  (`_match_group`), and the `experiment`/`condition` columns were added to the
  Model-tab table so the per-condition rows are no longer indistinguishable. The
  text outputs (`riana_fit_peptides.txt`, `riana_fit_fractions.txt`) were always
  correct — both carry the right per-condition counts; this was a display-only
  mismapping in the GUI.

### M7 Tier-1 — K-acetyl proteoform key + GUI fold-point display (Track C/E) — 2026-06-21

#### Added

- **Lysine acetylation as a biological proteoform key.** `UNIMOD:1` joins
  `BIOLOGICAL_MODS` with prefix `ac`, so an internal **K-acetyl** peptidoform rolls
  up as its own unit (`P12345_acK106`). **Protein N-terminal** acetyl — the *same*
  `UNIMOD:1` — still folds into the bare protein: `io.mztab._proteoform_sites`
  skips the mzTab `pos 0` (N-terminal) occurrence, so only the side-chain (`pos≥1`)
  acetyl keys. The fuller per-experiment chemical-vs-biological mod policy is a
  recorded far-future item.
- **GUI fold-point display.** The Model-tab fitted-curve view (`CurveView.plot_fit`)
  now draws **chemical-mod-folded** points (consolidated into the curve by
  `core.fitting._fit_key`) as a distinct purple ◇ series, alongside the existing
  MBR orange △ and direct-ID blue ○. The per-point `metox` flag is generic over
  `CHEMICAL_MODS`, so it covers Met-Ox today and **TMT** the moment it lands — no
  further GUI change. `fitting._build_output_df` now carries the per-point `metox`
  list on the wide frame (dropped from the slim `riana_fit_peptides.txt`, kept
  in-memory for the curve).

### Advanced integration knobs surfaced on both CLI and GUI (Track E) — 2026-06-20

#### Added

- **All `IntegrationConfig` tuning dials are now reachable on both surfaces** — six
  were frozen at their dataclass default with no flag/widget anywhere, and several
  more existed on only one surface. New CLI options (in a `rich_help_panel`
  **"Advanced integration"** group): `--ppm-alert` (previously GUI-only),
  `--prominence-k`, `--width-rel-height`, `--apex-n-consensus`, `--scan-rt-tol`,
  `--smoothing-polyorder`. The five MBR options moved into a **"Match-between-runs
  (MBR)"** help panel.
- **GUI: a collapsed "Advanced…" group** on the Integrate tab (the form is now in a
  `QScrollArea`) carrying every previously CLI-only knob — smoothing window +
  poly-order, mass difference, apex-search ½-width, prominence-k, width-rel-height,
  apex-N-consensus, an optional extraction-½-width override, write-intensities, the
  scan↔RT guard toggle + tolerance — plus a checkable **MBR** sub-group (the 5 MBR
  knobs). `build_config()` now passes the full field set, so the GUI can build any
  config the CLI can (the no-drift contract). Baseline / apex-selection / drift-alert
  moved from the everyday form into Advanced to match the CLI's panel split.
- Tests `test_build_config_advanced_widgets_flow_through` +
  extended `test_build_config_defaults_round_trip` (GUI ↔ config parity).

### Match-between-runs (MBR) for the mzTab/DDA path (Track A) — 2026-06-18, feature-complete 2026-06-20

#### Added

- **`integrate --mbr` — gated pure RT-transfer MBR for DDA.** Donor = a precursor
  identified at q≤0.01 in ≥2 runs of its `(experiment, condition)` curve group; a
  robust per-run RT offset (median over shared IDs) places a synthetic
  `scan=-1`+RT `PSMRecord` flagged `evidence="mbr"`, which flows through the
  existing M6b `resolve_rt_anchored_scans` + apex re-detect. A **two-part
  MBR-only quality gate** admits the transfer: `--mbr-min-snr` (apex
  prominence/local-noise; `inf`=fail, a sparse MAD=0 trace) and `--mbr-min-scans`
  (nonzero scans in the window), defaults **4 / 3**, uncapped. Assembled in
  `plan_integration` (mzTab is whole-experiment, so cross-run donor assembly is
  free); the scan↔RT guard runs only on the directly-scanned subset.
- **`fit --exclude-mbr` / `rollup --exclude-mbr`** opt-outs. New output columns:
  `apex_snr` / `n_scans` (integrate); `n_points` / `n_mbr` / `n_metox` / `n_clean`
  (fit, rollup) + per-point `evidence` (fractions). GUI: orange-triangle MBR points
  + table census + per-run gate-drop count in the log.
- **Benches:** `bench_missingness.py`, `bench_rt_alignment.py`,
  `bench_mbr_quality.py`, `bench_mbr_ab.py`; design `reports/2026-06-17_mbr_v1_design.md`.

#### Validation

- **Fit A/B verdict:** *ungated* MBR is harmful (R²>0.95 −30%, pollutes clean
  curves); *gated* MBR is neutral at strict R²>0.95 and **net-positive at the
  in-vivo gates (+180 at R²>0.8) with no pollution** → shipped uncapped. Eval-only
  within-protein-θ / yield / k_deg-CV show a reasonable yield-for-consistency
  tradeoff (+37 proteins at R²>0.8 for ~+3.5% scatter).
- **RETRACTED 2026-06-20** the earlier calibration "MBR mis-quantifies at high
  label" result — root-caused as an mzTab↔mzML RT-axis mismatch (the calibration
  mzTab is `.raw`-searched, RT axis 2–4 min off the local `.mzML`, and
  `--no-rt-check` bypassed the guard). Not a labeling effect. Maintainer TODO:
  re-search quantms on the exact `.mzML`, re-run the sweep without `--no-rt-check`.

### M7 — PTM-aware envelope (Track C) — 2026-06-13 → 2026-06-17

#### Added

- **Atom-vector extension `[C,H,O,N,S]` → `[C,H,O,N,S,P]`** across `count_atoms` /
  `_calc_atom_mass` / `constants` (`aa_atoms`, `iso_abundances`, mass vector) and
  the IsoSpec `IsoParamsFromDict` formula, so phospho's P (monoisotopic) + its
  3 O shape the envelope. Verified byte-identical on unmodified peptides (mass Δ 0,
  envelope Δ 3e-18) and against the frozen 5-element M2 oracle.
- **Variable mods threaded IO→integrate→fit.** Parsed bracketed UniMod masses now
  enter the IsoSpec `formula` via a curated UNIMOD-id-keyed `mod_atoms` table
  (`{id: [C,H,O,N,S,P]}`), and modified peptidoforms integrate at their **own**
  m/z + envelope instead of being dropped. Each `[UNIMOD:N]` form is a distinct
  `concat`; mod hydrogens stay out of `num_labeling_sites` (mod enrichment unknown).
- **Proteoform-aware rollup keys.** `PSMRecord.mod_sites` carries the biological-mod
  site in protein coordinates (`pS34476` from mzTab `start`+`pos`); a
  `BIOLOGICAL_MODS={21}` (phospho) set drives a split where regulated mods get a
  distinct key (`A2ASS6_pS34476`) and constitutive/chemical mods (N-term Ac,
  Met-Ox, CAM) fold into the bare accession.
- **Chemical-mod integrate-separate / fit-merge (Met-Ox).** `CHEMICAL_MODS={35}`;
  `core.fitting._fit_key` strips chemical-mod tokens so the oxidized and unoxidized
  forms of a peptide land on one turnover curve (they share the same D₂O clock),
  with per-row FS on each form's own envelope.

#### Removed

- **Retired `-X/--ignored_mods` and `-F/--forced_mods`** (SILAC-era dual-channel
  machinery) from CLI, `IntegrationConfig`, and the `mod{offset}` channel path in
  `integration.py` — `integrate_run` now emits `iso{N}` directly at the PSM's own
  m/z. Removed the `ignored_mods` plumbing through `io/mztab`, `io/diann`,
  `io/percolator`, `pipeline`, `gui`. The `iaa` special-case for Carbamidomethyl
  is gone — CAM routes through `mod_atoms[4]` (`UNIMOD:4`) like any other mod.

#### Validation

- New baseline `runs/lve_atr_m7/` (git `ef8e3ae`): 24 runs, 20,955 converged
  peptidoforms, 1,986 proteins; `[UNIMOD:35]` in all 24 integrate outputs but 0
  fit-key rows (merged); phospho proteoform keys (`Q02566_pT2`, `Q9JJW5_pT107`)
  landed. 178-test suite + the 415s ac16 byte-identical gate pass.

### Parallelism — `-W/--workers` everywhere; `-t/--thread` removed

#### Removed

- **`-t/--thread` removed from `integrate`, `fit`, and `rollup`** (and the
  `threads` field on `IntegrationConfig`/`FitConfig`, and the GUI "Threads"
  spins). Measured 2026-06-13: thread-level parallelism gave **no** speedup and
  was often *slower* than serial — the hot loops hold the GIL (IsoSpec FS, scipy
  `curve_fit` + bootstrap, and pyteomics' mzML XML parse). rollup `-t 8` ran 551s
  vs ~435s serial / 60s at `-W 8`; integrate degraded monotonically (719 → 806 →
  1156 → 1548s for threads 1→2→4→8). **Use `-W/--workers`** (process-level), the
  real lever: `integrate -W` across files, `fit -W` across peptides, `rollup -W`
  across proteins — all deterministic regardless of N.

### Track C — `linear simple` model + cross-sample Δk

#### Added

- **`riana rollup --model "linear simple"` — the linearized cross-sample Δk
  model** (`core/linear_model.py`). A model choice mutually exclusive with the
  nonlinear ODE models: transforms fraction-new θ → clearance φ = `log(1 − θ)`
  (linear through the origin for the simple model) and fits a protein's conditions
  **jointly** in one no-intercept OLS with a `day:condition` interaction
  (statsmodels), yielding a per-condition `k_deg` (= −slope) **and** a pairwise
  **Δk** test with a shared-variance p-value, then Benjamini-Hochberg across
  proteins. Output schema `PROTEIN_LINEAR_COLUMNS` adds `delta_k` / `delta_k_se` /
  `delta_k_p` / `delta_k_p_adj`. CIs are **analytic** (statsmodels `conf_int` /
  `t_test`), not the residual bootstrap the nonlinear models use. New flags
  `--phi-limit` (plateau truncation, default −4 ≈ θ 0.98 — linear-only, since the
  saturated tail is floor-noise that flattens the through-origin slope) and
  `--reference-condition`. Adds a `statsmodels>=0.14` dependency. Validated on
  `data/timeseries_lve_atr` (atrium vs ventricle: 780 proteins in both chambers,
  atrium 1.24× faster, 410 significant at BH p_adj<0.05).
- **`riana rollup -W/--workers`** — process-level parallelism for the GIL-bound
  per-protein refit (mirrors `fit -W`); byte-identical to serial, ~7.2× cores.

#### Fixed

- **Two-condition `fit_project` crash** — the final `pd.concat` of per-curve
  frames raised `Can only compare identically-labeled` because each frame carried
  a `fractions_long` DataFrame on `.attrs` (pandas reconciles attrs by equality
  across frames once there is >1 curve). Now dropped before the concat.

### Track A / E — variable-mod drop + GUI φ-space visualization

#### Added

- **`io/mztab` variable-mod drop** (`drop_variable_mods`, default on) — mirrors
  the DIA-NN path: PSM rows carrying a non-fixed UniMod (Oxidation, Phospho,
  N-term Acetyl, …) in the `modifications` column are dropped rather than
  integrated at the *unmodified* m/z until M7. Carbamidomethyl (`UNIMOD:4`) is
  kept. On the `timeseries_lve_atr` PTM search this drops ~2.3% of PSMs.
- **GUI φ-space plotting + CI ribbons** — the Protein tab's Model combo offers
  `linear simple` (with φ-limit + reference-condition knobs shown only for it);
  selecting a protein overlays both conditions in φ-space (`CurveView.plot_linear`)
  with the through-origin k lines, truncated points drawn hollow, and the Δk +
  p_adj in the title. Both the θ-space fit curves and φ-space lines now shade a
  **confidence ribbon** from each k CI. The fit and rollup tabs gained a
  **Workers** spin (dispatched on a main-process thread when >1 to avoid nested
  pools).

### M6b — DIA-NN parquet intake (Track A)

#### Added

- **`riana/io/diann.py` — DIA-NN `report.parquet` intake.** The DIA counterpart
  of `io/mztab.py`: reads a quantms-diann DIA-NN report (≥ 2.2.0; validated on
  **2.5.0** with the `diann` parquet output), disaggregates it by the `Run`
  column, and attaches each run's `RunIdentity` from the **same `io/sdrf.py`**
  (DIA is auto-detected from `comment[proteomics data acquisition method]`), so
  every `PSMRecord` carries its full identity exactly like the DDA path. The
  peptide mass is recomputed via `mass_calc.calculate_ion_mz` (Carbamidomethyl(C)
  counted) so the m/z target derivation is identical across DDA/DIA;
  `Precursor.Mz` is used only as a self-check (logged median |ppm|). `Protein.Ids`
  (`;`-joined) is normalized to Riana's `,` convention. Decoys are dropped
  (the report is already FDR-filtered, so this is defensive).
- **RT-anchored extraction (`core/integration.resolve_rt_anchored_scans`).** DIA
  has no precursor-specific MS2 scan to anchor on, so `io/diann` emits `scan = -1`
  and carries DIA-NN's inferred apex `RT` in `PSMRecord.retention_time`. At
  integrate time each PSM's apex RT is resolved to the nearest MS1 scan in *this*
  mzML, and the **existing scan-based extraction runs unchanged** — DIA-NN is
  "just another ID + RT source"; Riana still pulls the MS1 isotopologues from the
  mzML, and the apex finder re-centres on the true MS1 apex (the ≤1-cycle anchor
  offset is well inside `apex_search_half_width`). The integrator detects the DIA
  case from `scan < 0` (a no-op on the DDA path).
- **`dia` install extra** (`pyarrow`), lazily imported in `io/diann` so a
  DDA-only install stays lightweight (with a clear `pip install riana[dia]` hint
  if the parquet path is used without it).

#### Changed

- **`plan_integration` dispatches the ID reader on `SdrfTable.acquisition`** —
  DIA → `read_diann` (parquet), else `read_mztab` (mzTab). Both CLI and GUI build
  the same per-run `RunTask`s, so the DIA path comes for free on both surfaces.
- **The intake scan↔RT scramble guard is skipped for DIA** (it is circular there:
  the scan is *derived* from the reported RT, so it always reconciles).
  `resolve_rt_anchored_scans` does an **RT-in-bounds** check instead — raising
  when most of a report's apex RTs fall outside the paired mzML's MS1 RT span,
  the DIA analog of a wrong mzML↔report pairing.

#### Caveat (M7)

- DIA peptidoforms carrying a **variable modification** (e.g. Oxidation(M),
  `UniMod:35`) are **dropped by default** (`drop_variable_mods`) — like the DDA
  mzTab path, the envelope is keyed on the *stripped* sequence, so a variable-mod
  form would integrate at the unmodified m/z until M7 threads mods into the
  envelope. On the validated cardiac DIA set this is ~1.8% of rows and costs only
  the few precursors seen *exclusively* as a modified form.

#### Validation

- Cardiac in-vivo D₂O DIA set (mouse left ventricle, 9 runs = 3 timepoints
  {3, 7, 14 days} × 3 biological replicates; DIA-NN 2.5.0, 10 ppm from the SDRF).
  Ran the full `integrate → fit → rollup` chain end-to-end.
  - **Intake/extraction is clean:** 196,310 PSMs across the 9 runs (~20–24k
    precursors/run after the ~1.8% variable-mod drop); **98.5% of rows extract a
    non-zero m0** (<0.1% all-zero) and the median |m0 ppm error| is **0.8–2.1 ppm
    per run** — i.e. the RT→MS1 anchor lands the right precursor on-target. (The
    ~1.5% MS1 miss is expected: DIA-NN IDs off MS2, so some precursors have
    weak/absent MS1.)
  - **Fit R² is lower than DDA, as anticipated for MS1-on-DIA:** 24,827 peptides
    converged; R²med **≈0.42**, ≥0.8 **≈22%** (vs the DDA LVE turnover baseline
    R²med ~0.87 / ≥0.8 ~61% at the same 10 ppm) — DIA's wide-window MS1 imports
    co-eluting interference into the m1–m5 channels, and the curve is only 3
    timepoints at RIA 4.6%. **The central estimate is sound:** k_deg median
    **0.080/day** (IQR 0.036–0.176; 94% in a plausible range), biologically
    reasonable for cardiac turnover. There is no DDA counterpart of these exact
    samples, so this is a "roughly against other datasets" sanity check, not a
    paired benchmark.
  - **The pooled R² is dominated by replicate spread, not extraction** (internal
    A/B): the 3 bioreps are *independent animals* (terminal sampling) pooled into
    one 9-point curve, so it absorbs inter-animal variance. Re-fitting each
    replicate alone (a clean 3-point curve) ~doubles the headline: R²med
    0.42→**0.68–0.72**, ≥0.8 21.5%→**38–42%**, while k_deg median is unchanged
    (0.068–0.075/day). So the central estimate is stable and the pooled-curve
    scatter is a study-design property (also the source of within-protein curve
    spread), not a DIA intake defect.

### Run-identity keys — formalize the fit/rollup grouping (Track A follow-up)

#### Changed

- **`RunIdentity.curve_key` → `RunIdentity.group_key`**, fixed to
  `(experiment, condition)` — the *actual* run-grouping `recombine_for_fit` uses.
  The old `curve_key = (experiment, sample, biological_replicate)` was **unused
  and wrong**: `sample` is per-run (it would have split every run into its own
  curve) and biological replicates are *pooled* as independent points, not
  separated.
- **The fit/rollup grouping keys are defined once** in `riana.records` —
  `GROUP_KEY_COLUMNS = (experiment, condition)`, `CURVE_KEY_COLUMNS = (…, concat)`
  (one fitted peptide curve — "per concat per condition"), `PROTEIN_KEY_COLUMNS =
  (…, protein)`. `recombine_for_fit` / `fit_project` (group), `core.protein`
  (`_GROUP_KEYS`, output schema), and `fit_project`'s curve-uniqueness guard all
  reference these instead of re-listing column names, so the stages can't drift.
  The typed `RunIdentity` is still unpacked into columns at fit recombination by
  design — a curve/protein spans many runs, so it is not one identity object.

### Fit — process-level parallelism (`riana fit -W/--workers`)

#### Added

- **`riana fit -W/--workers N` — process-pool fit.** The per-peptide fit
  (IsoSpec forward-model FS + residual-bootstrap CIs) is **GIL-bound**, so the
  existing `-t/--thread` path barely scales (measured ~1.25 effective cores at
  `-t 8` on the cardiac DIA set, 28,916 peptides). `-W` dispatches `fit_run`'s
  per-concat map over a `ProcessPoolExecutor` (shared frame/coefficients pickled
  once via a pool initializer; only `concat` strings cross per task), reaching
  ~8× on 8 workers. `FitConfig.workers`; threaded/serial stays the default.
- The GUI Model tab is unaffected (it builds `FitConfig` with the default
  `workers=1`).

#### Changed / Fixed

- **Bootstrap seeding is now content-stable** — the per-peptide RNG seeds from
  `hashlib.blake2b(concat)` (`core.fitting._concat_seed`) instead of the
  process-salted built-in `hash()`. This makes the bootstrap CIs **identical
  regardless of `-W` (and reproducible run-to-run)**, mirroring the rollup's
  `_group_rng`; a test pins `workers=1 ≡ workers=2`. (The old `hash()` seed was
  silently non-reproducible across processes/runs.)

### Output hygiene (fit / rollup)

#### Changed

- **Rollup outputs/stage renamed** for symmetry with fit: `riana_protein.txt` →
  **`riana_rollup_proteins.txt`** (+ a new **`riana_rollup_fractions.txt`** with
  the collapsed inverse-variance-weighted `(t, θ)` the GUI curve plots), and the
  manifest stage `protein` → **`rollup`**. Pre-rename manifests still load
  (`protein` aliased to `rollup`).
- **Manifest rows gain a `created_at` timestamp** (ISO-8601 UTC). The manifest is
  an *index*, not a history — a re-run overwrites the output and **replaces** its
  row (`config_hash` = settings, `git_sha` = code or `"unknown"` without git,
  `created_at` = when).

- **`riana_fit_peptides.txt` is now a scalar summary** — the per-timepoint
  `t`/`fs`/`fs_lower`/`fs_upper` list-cells are dropped from the *written* file
  (the per-timepoint detail, *with* biorep labels, lives in
  `riana_fit_fractions.txt`). The in-memory result keeps them so the GUI curve is
  unaffected. `core/fitting.peptide_summary`.
- **Estimate outputs are rounded to ~6 significant figures** (`%.6g`) — the fit
  peptides / fractions and protein files. Readable and well beyond the precision
  of k / θ / R² (reproducibility is guaranteed by the provenance header, not
  bit-exact floats). The integrate `_riana.txt` is left full-precision (its
  `isoN_obs_mz` needs ppm-level digits, and it is parity-gated).

### Track E — GUI rewiring onto core/pipeline

#### Changed

- **`integrate_project` split into plan / dispatch / finalize** — `plan_integration`
  (SDRF+mzTab → per-run `RunTask`s), `_integrate_results` (yields each frame as it
  finishes), and `finalize_run` (write `<stem>_riana.txt` + manifest row). The
  split lets the GUI drive the *same* per-run unit over its own pool.
- **Integrate is now crash-resilient.** The main process **writes each run's
  `_riana.txt` + appends its manifest row as that run completes** (workers only
  integrate), so an interruption keeps the runs already finished — the old
  buffer-all-then-write design lost everything on a kill. New **`integrate
  --resume`** skips runs already in the manifest at the current `config_hash`
  (same SDRF/mzTab inputs assumed), so a crashed run continues where it stopped.
- **GUI Integrate tab** now goes through `core/pipeline`: an **SDRF** field routes
  the search-ID file as the quantms mzTab (identity-stamped per-run outputs +
  `riana_manifest.tsv`), and both paths are **file-parallel** via a *Workers*
  spinbox (`asyncio.as_completed` over the shared pool, bounded one-mzML-per-run —
  no nested process pools), on top of the per-run *Threads*.
- **GUI Model tab** gained a **Manifest** field (`fit_project` via the new
  `tasks.run_fit_manifest` worker) and now also writes the M5
  `riana_fit_fractions.txt`.
- **GUI Integrate tab reflects the SDRF mass tolerance** — picking an SDRF fills
  the Mass tolerance spinbox from its `comment[precursor mass tolerance]` (the
  user can still override), the GUI equivalent of the CLI's SDRF resolution; the
  stale `50 ppm` spinbox default now matches the config's `10 ppm`.

#### Added

- **`riana rollup --thread`** (+ a Protein-tab Threads spinbox) — the per-protein
  refit (`curve_fit` × bootstrap) now runs over a thread pool. Per-protein RNG
  streams (`_group_rng`) keep the result identical regardless of thread count.
- **Manifest project chain (`integrate → fit → rollup`).** `fit --manifest` now
  writes its outputs **next to the manifest** (the project dir; `-o` is ignored
  there) and records `stage="fit"` rows; **`rollup --manifest`** (the fit dir
  argument is now optional) reads those rows to find the fit outputs, writes
  `riana_rollup_proteins.txt` (+ `riana_rollup_fractions.txt`) next to the
  manifest, and records `stage="rollup"` rows.
  So one `--manifest` drives the whole chain and the manifest indexes every
  stage. The GUI Model/Protein tabs follow the same rooting. (Previously `fit`
  never touched the manifest — it was integrate-only despite the documented plan.)

### Track C — protein rollup

#### Added

- **`riana rollup`** — rolls the `riana fit` per-peptide outputs up to one
  turnover estimate per `(experiment, condition, protein)`, writing
  `riana_rollup_proteins.txt` (a `method` tag + `k_deg`/`ci_lo`/`ci_hi`/`R_squared`,
  the data-structure counts `n_peptides`/`n_replicates`/`n_timepoints`/`n_points`,
  and a comparison `peptide_median_k`) + `riana_rollup_fractions.txt` (the
  long/tidy collapsed `(t, θ)` behind each refit — the GUI curve substrate).
  `core/protein.py`. (The manifest stage is `rollup`; older `protein`-stage
  manifests still read, aliased.)
- **`--method {weighted, pooled}`** (default `weighted`) — `weighted` is the
  **biorep-aware per-timepoint refit**: within each `(protein,
  biological_replicate, labeling_time)` the peptides' fraction-new θ are collapsed
  by an inverse-variance weighted average (σ from the M5 prediction interval),
  then one `k_deg` is fit to the collapsed `(t, θ)` across timepoints *and*
  bioreps (independent points → honest dof). `pooled` fits all peptide×timepoint
  points directly (pseudoreplication; for comparison). The median/harmonic point
  estimators were dropped as standalone methods (trivially DIY on the peptide
  file; the median is still carried as `peptide_median_k`). The linearized fit +
  a cross-sample Δk test are a deferred follow-up.
- **`--parsimony {unique, isoform}`** (default `unique`) — protein attribution is
  a **summarize-time** decision: a shared peptide's envelope blends both
  proteins' turnover, so it can't be attributed. `unique` keeps only
  single-accession peptides; `isoform` additionally folds isoform-only-shared
  peptides into the canonical entry **unless** an isoform in the group carries
  its own unique peptide (a dataset-wide pass over the peptide↔protein map,
  adapted from `02_R_parsimony_reference.Rmd` — Riana-`,` separator, base-accession
  fallback, UniProt `-N` *and* JCAST `-J1`/`-J2` isoform suffixes; θ weighted by
  the M5 inverse variance, not `log2(Int)`).
- **Optional peptide R² admission gate** (`--min-r2`, off by default) — hard-excludes
  peptides whose kinetic fit doesn't follow the model, *complementing* the
  inverse-variance weighting, with a JCI-style slow-turnover admit (`--alt-k` /
  `--alt-se`) so legitimately slow peptides (low R² only because θ barely moves)
  survive. Off by default keeps the inverse-variance-only result as a clean A/B baseline.
- **GUI Protein tab** — a 3rd tab over the same `rollup_proteins` core
  (`gui/protein_tab.py` + the Qt-free `tasks.run_rollup` pool worker): pick a fit
  output dir, choose parsimony/model/R²-gate, run, write the rollup outputs, and
  click a protein row to see its collapsed `(t, θ)` points + refit curve.

#### Removed

- **`integrate --unique`** (and `IntegrationConfig.unique_only` + the GUI
  checkbox). Integrate now extracts **all** peptides — shared peptides are valid
  per-peptide measurements; the unique/isoform filtering moved to `riana rollup
  --parsimony`, where protein attribution belongs.

### M5 — per-timepoint fraction-new (Track C)

#### Added

- **`riana_fit_fractions.txt`** — `riana fit` now also writes a long-format,
  one-row-per-`(concat, biological_replicate, labeling_time)` table with the
  per-timepoint fraction-new `fs` and prediction-interval bounds `fs_lower` /
  `fs_upper`. This is the substrate the protein rollup consumes. (The wide
  `riana_fit_peptides.txt` is a scalar summary — see "Output hygiene" above.)
- **`build_fractions_long`** (`core/fitting.py`) + `out.attrs["fractions_long"]`
  carry the long table through `fit_run` → `fit_project` (tagged with
  `experiment` / `condition`) → CLI without re-deriving θ.

#### Changed

- **Bootstrap is now a unified residual bootstrap** (fixed t-design: resample
  the kinetic-fit residuals, refit) feeding *both* the `k_deg` CI and the new
  per-timepoint band. `fs_lower` / `fs_upper` are a **prediction interval**
  (`model(t_i; k*) + resampled-residual`), so they reflect each peptide's
  measurement scatter, not just curve uncertainty. More robust than the prior
  pairs bootstrap on sparse 3–5-point curves. `k_deg` / `R²` are unchanged (they
  come from the main fit); `k_deg` CI bounds shift slightly (no test pins them).

### Track A — intake scan↔RT guard

#### Added

- **Per-run scan↔RT reconciliation guard** in `integrate_run`: each PSM's mzTab
  `spectra_ref` scan is mapped to this mzML's MS1 retention time and compared to
  the mzTab-reported `retention_time`; if the per-run **median** offset exceeds
  `scan_rt_tol_min` (default **2.0 min**) it raises `DataError`. This catches the
  quantms filename-prefix scan-scramble (mzML basenames that are prefixes of one
  another) and wrong mzML↔mzTab pairings — previously silent. The median gate is
  robust to the run-dependent ProteomicsLFQ alignment offset (≤~0.9 min measured
  on real output); a scrambled run sits tens of minutes off. No-ops on the
  Percolator path (no `retention_time`).
- **`--no-rt-check`** flag (config `check_scan_rt`) to override the guard for a
  run known to be correctly paired; **`IndexedMzML.rt_for_scans`** vectorized
  scan→RT lookup.

### Pre-1.0.0 chores

#### Changed

- **Version → `1.0.0`** (dropped the `.dev1` marker).
- **Docs** are no longer committed as rendered HTML. The Quarto source under
  `riana_website/` now renders to a gitignored `_site/` and is published to
  GitHub Pages by a new `.github/workflows/docs.yml`. *Manual one-time step:* set
  the repo's Pages source to "GitHub Actions" before merging to `master`.

#### Removed

- **The bundled `workflow/Snakefile`** and its `config_template.yaml`. Riana is
  orchestration-agnostic: quantms / DIA-NN own search + ID upstream, and Riana is
  a linear `integrate → fit` chain glued by the manifest (below). Compose the
  subcommands into whatever workflow already runs them.

### M6a — run-identity data model & intake refactor (Track A)

#### Added

- **`RunIdentity`** (`records.py`): the SDRF-sourced identity
  (experiment, sample, data_file, bio/tech replicate, fraction,
  labeling_time | mixing_proportion, condition, acquisition,
  precursor_enrichment) attached to every PSM at intake and carried
  header-authoritatively to fit. The experiment type is *declared by which
  independent-variable column is present*. `PSMRecord` also gains
  `retention_time` (the DIA RT-apex prior) and `identity`.
- **`io/sdrf.py`** — `read_sdrf()` reads Riana's documented SDRF column subset
  (handling duplicate `comment[modification parameters]` columns), validates it,
  and exposes a `sample_map` keyed by mzML stem. It deliberately does **not**
  read `comment[precursor mass tolerance]` (the tight search window).
- **`io/manifest.py`** — a schema-versioned, stage-aware `riana_manifest.tsv`
  (`stage = integrate | fit | protein`) that indexes every stage output with its
  identity; the project glue between subcommands.
- **`core/pipeline.py`** — the shared orchestration extracted from
  `cli.integrate` / `gui.tasks`: `integrate_project()` (SDRF → mzTab → one
  `<mzml_stem>_riana.txt` per run, identity in the provenance header, manifest
  rows) and the fit-time recombination (`recombine_for_fit` / `fit_project`)
  that groups runs into kinetic curves by `(experiment, condition)`, merges
  fractions of the same `(biorep, timepoint)` at peptide level, and takes the
  curve x-axis from the identity — not from the `sample` string.

#### Changed

- **`io/mztab.py`** — `read_mztab(path, sample_map=…)` is now the primary path:
  it joins per-`ms_run` identity by location basename and carries the mzTab
  retention time. The bare `sample=` call is kept for legacy/single-run tests.
- **`riana integrate --sdrf SDRF`** drives the identity-keyed mzTab intake (one
  output per run + manifest). Without `--sdrf`, the single-mzML Percolator path
  is unchanged (the demoted testing/legacy tier).
- **`riana fit --manifest MANIFEST`** groups runs into curves from the manifest
  identity. The positional `_riana.txt` path stays for legacy/single-curve fits.
  `fit_run` gained a `time_column` so the manifest path reads the timepoint from
  the identity; a parity test pins the manifest path's `k_deg` to the legacy
  path's.

### M4 Phase 2 — PySide6 GUI

#### Added

- **`riana gui` subcommand** launches a PySide6 GUI (`riana/gui/`). PySide6 /
  qasync / pyqtgraph are imported lazily, so they never load on the core CLI
  path; install them with the new optional extra: `pip install 'riana[gui]'`.
- **Model tab** — a form mirroring `riana fit` builds the *same* frozen
  `FitConfig` (shared `__post_init__` validation). The fit runs on the
  `ProcessPoolExecutor` via the Qt-free `riana.gui.tasks.run_fit` worker (reads
  the per-timepoint `_riana.txt` files, loads the `--coefficients` table, calls
  the shared `core.fitting.fit_run`), so output matches the CLI. Per-peptide
  results land in a table; selecting a peptide draws its `(t, fraction-new)`
  points and the fitted kinetic-model curve in an embedded pyqtgraph view
  (`riana.core.models` functions evaluated directly — no worker round-trip).
- **Integrate tab** — a form whose fields mirror `riana integrate` and build the
  *same* frozen `IntegrationConfig`, so its `__post_init__` is the single shared
  validator (CLI and GUI cannot drift). Integration runs off the Qt event loop:
  CPU work is awaited on a `ProcessPoolExecutor` one fraction at a time (via
  qasync), keeping the UI responsive with honest per-fraction progress. Results
  land in a table, the per-fraction mass-accuracy **drift summary** (median/MAD
  ppm, suggested shift, `--ppm-alert` flag) is shown inline (no separate
  Calibration tab), and selecting a peptide draws its isotopomer XICs in an
  embedded **pyqtgraph chromatogram** with the integrated window shaded. The
  `_riana.txt` (+ `.drift.json`) output is identical to the CLI's.
- `riana.core.integration.extract_peptide_trace` + `PeptideTrace` — extract one
  peptide-charge's per-isotopomer XICs and integration window (reuses the
  integrator's own extraction/boundary code), giving the previously orphaned
  `records.Chromatogram` dataclass a producer.
- `riana.io.mzml.list_mzml_files` / `mzml_stem` — the shared mzML
  directory-layout helpers the CLI and the GUI worker both use to assign the
  same fraction order.
- `tests/test_gui.py` — Qt-free worker tests (assert the GUI integrate path
  matches the committed `sample1` golden within 1e-3, and that `run_fit` fits a
  synthetic D₂O series) + headless `pytest-qt` smoke tests for both tabs, the
  form→config validation, and the fitted-curve plot.

#### Changed (BREAKING)

- **`riana integrate --iso` default is now `0 1 2 3 4 5`** (was `0 6`). The D₂O
  fit matches the observed envelope against an IsoSpec forward model over the
  contiguous m0-m5 channels, so the legacy `0 6` pair is not fittable by the new
  engine. `--iso` is still free-form for other workflows (e.g. SILAC cluster
  extraction via `-F`). The reframed `--fs` help calls out the post-M4 plan.
- **`riana fit` now errors clearly when the integrate output lacks the
  isotopomers the D₂O solver needs** (m0-m5). Previously `--iso 0 6` data was
  silently misaligned against the model, producing garbage `fs`/`k_deg`; it now
  raises with a "re-run integrate with `--iso '0 1 2 3 4 5'`" message.

#### Removed

- **`riana fit --plotcurves`** — it was a no-op in the new fit engine (the
  legacy static-PNG path was never wired in). Inspect fitted curves
  interactively in the GUI Model tab instead. `--fs` is **kept but currently
  ignored** (reserved for a post-M4 channel-subset envelope SSE; a warning is
  logged if it is set).

### M4 Phase 1 — Typer CLI + `--engine legacy` removal

#### Changed (BREAKING)

- **The CLI is now a Typer app** (`riana/cli.py`); the argparse `riana/main.py`
  is gone. The console entry point is `riana = riana.cli:main`. The typed
  pipeline (`core` + `io`, driven by the frozen `IntegrationConfig` /
  `FitConfig`) is the **only** engine — the `--engine legacy|new` flag is
  **removed**. Reproduce 0.9.0 integration with `--peak-rt ms2
  --integration-half-width 1.0` (numerically within 1e-3; pinned by a committed
  golden, `tests/data/sample1/sample1_riana.v0_9_0.txt`).
- **List flags take a single comma/space-separated token**, not argparse
  `nargs='+'`: `-i "0 1 2 3"` (was `-i 0 1 2 3`); same for `-X` / `-F`.
- **`riana fit --coefficients` is now REQUIRED** for `--label hw`. Pass a
  bundled preset (`commerford` | `ac16` | `ipsc` | `cm`) or a path to a
  `(amino_acid, coefficient)` CSV. The old fixed-site-count default is gone.
- **`riana fit --label` is now `{hw, o18}`** (was `1|2|3|4`). The deuterium
  in-vivo/in-vitro split collapses into one `hw` (heavy water / D₂O) mode —
  cell/tissue specificity comes from the coefficient table, not the label.
  **Amino-acid / SILAC fitting (old `--label 4`, `-a/--aa`) is dropped**:
  integrate still extracts SILAC peaks via `-F/--forced_mods`; do the
  L/(H+L) curve fit downstream. **`o18` is recognized but errors** ("being
  reimplemented post-M4"); o18 *integration* is unaffected.

#### Added

- Bundled per-AA D₂O coefficient presets under `riana/data/coefficients/`
  (`commerford` literature + `ac16`/`ipsc`/`cm` calibration tables);
  `riana fit --coefficients <name>` resolves a preset, else a filesystem path
  (`core.fitting.load_aa_coefficients` / `available_coefficient_presets`).
- `tests/test_cli.py` (Typer `CliRunner` smoke tests).

#### Removed

- Legacy modules `riana/{main,riana_integrate,riana_fit,spectra,peptides,project}.py`
  and the re-export shims `riana/{accmass,fsynthesis,models}.py`. The broken
  Tkinter GUI (`riana_ui/`) is deleted (PySide6 GUI lands in M4 Phase 2).
- Tests that A/B-tested against the live legacy modules now compare against
  committed golden fixtures captured from 0.9.0 before deletion
  (`sample1_riana.v0_9_0.txt`, `percolator_parity_v0_9_0.csv`,
  `mzml_ms1_v0_9_0.npz`).

### M3 Week 1 — skeleton + lifts

#### Added

- New package layout: `riana/core/`, `riana/algorithms/`, and `riana/io/`
  subpackages. The numerically-sensitive science modules are *lifted* into
  their new homes unchanged rather than rewritten — `accmass` →
  `algorithms/mass_calc`, `utils.get_peptide_distribution` →
  `algorithms/isotope_dist`, `models` → `core/models`, `fsynthesis` →
  `core/fsynthesis`.
- `riana/records.py` — typed, frozen data records (`PSMRecord`,
  `Chromatogram`, `IsotopomerPeak`) replacing the untyped `pandas` rows and
  positional lists the 0.9.0 pipeline passes around. Field names are anchored
  to the existing `*_riana.txt` column schema.
- `riana/config.py` — frozen `IntegrationConfig` / `FitConfig` dataclasses,
  the single typed source of truth shared by the CLI and (M4) GUI.

#### Changed

- Version is now `1.0.0.dev0` (PEP 440) for the duration of the M3 rewrite.
- During the rewrite the legacy flat modules (`accmass`, `models`,
  `fsynthesis`, `utils`) remain as thin re-export shims, so the 0.9.0
  `riana integrate` / `riana fit` CLI keeps working as a regression gate.
  The shims are removed with the legacy pipeline in M3 Week 4.

#### Fixed

- `fsynthesis.calculate_a0` tested `label == 'aa'` (a string) while callers
  dispatch with an integer label, so the amino-acid `a_0` branch was
  unreachable and `--label 4` experiments silently used the natural-abundance
  baseline (`PROJECT_REVIEW.md` §2b). Corrected to `label == 4`. This changes
  fit output for amino-acid-labeling experiments.

### M3 Week 2 — I/O layer (both ID paths)

#### Added

- `io/percolator.py` — typed Percolator parser (no exception-as-control-flow);
  emits `PSMRecord`s.
- `io/mztab.py` — quantms/OpenMS mzTab intake (`read_mztab` → `PSMRecord`s + an
  `ms_run → filename` map), the second ID pipeline.
- `io/mzml.py` — `IndexedMzML`, indexed/streaming reads (one fraction in memory,
  not the whole run).
- `io/writers.py` — TSV writer with a provenance header (riana version, git SHA,
  config hash, id source). Readers skip it with `comment='#'`.

### M3 Week 3 — core integration pipeline + peak-detection toolkit

#### Added

- `core/integration.py` — rewritten against the typed records, streaming per-PSM.
- `algorithms/peaks.py` — `find_apex`, `consensus_apex`, `detect_peak`, co-elution
  check, SNR/symmetry/quality scoring. `algorithms/baseline.py` — `noise_floor`
  (+ SNIP/AsLS) options.
- Per-isotopomer mass-accuracy columns (`iso{N}_obs_mz`, `iso{N}_ppm_error`) and a
  per-fraction drift sidecar (Phase D); `--ppm-alert`.
- `riana integrate --engine new` — dispatcher onto the typed pipeline.

### M3 Week 4 — fitting rewrite + §2b science fixes

#### Added

- `core/fitting.py` (consumes `IntegrationResult`); `riana fit --engine new`.
- Provenance headers on `_riana.txt` / `riana_fit_peptides.txt` (Phase F3).

#### Changed

- **Fractional synthesis is now computed by the IsoSpec forward/solve model**
  (`algorithms/isotope_dist.solve_fs_d2o`: per-peptide Spep + full-envelope
  least-squares) instead of the fixed-site-count `m0`-analytic relation. This
  closes the ≈ −0.5 pseudo-time `k_deg` recovery bias the M3 Week 0 baseline
  documented (median `k_rel_err` → ~0 with `--engine new`).

#### Fixed

- AA `a_max` dispatch (§2b); FS-denominator drift; kinetic-fit confidence via
  bootstrap.

### M3 peak-detection spike (pre-M4, 2026-06) — **default integration changed**

A time-boxed, benchmark-gated investigation of chromatographic peak detection and
baseline subtraction, validated on the D₂O calibration series (ac16/ipsc/cm at
0–100%) and an independent in-vivo mouse set.

#### Changed (BREAKING — default integration output differs from 0.9.0 by design)

- **Default integration is now an apex-centred narrow window**: `peak_rt="apex"`,
  `integration_half_width=0.15`, `apex_selection="tallest"`, `baseline="none"`.
  It beat the 0.9.0 fixed ±`r_time` rectangle on every test set (envelope RMSE
  vs IsoSpec; cross-proportion mixing linearity), generalizing across cell line,
  organism, and ID pipeline. **0.9.0 behaviour is reproducible with
  `--peak-rt ms2 --integration-half-width 1.0`** (the parity tests pin this).
- `integration_half_width` tracks the chromatographic peak width — dial it to
  your gradient (≈0.1 sharp UPLC … ≈0.2–0.33 broad).
- Renamed config field / CLI flag `r_time` → `extraction_half_width` (it
  conflated the *extraction* trace with the *integration* window; now distinct).
  `--r_time` kept as a legacy alias.

#### Added

- `peak_rt` options `apex` (default) and `consensus` (median apex over m0..m3 —
  labelling-independent, contamination-robust; prefer at high D₂O). Tunable
  knobs: `integration_half_width` (float|`auto`), `extraction_half_width`,
  `apex_selection` {tallest,nearest}, `apex_search_half_width`, `apex_n_consensus`,
  `prominence_k`, `width_rel_height`. CLI: `--peak-rt`, `--integration-half-width`,
  `--baseline`, `--apex-selection`.
- Benchmarks: `bench_mixing_linearity.py` (model-free mᵢ:mA-vs-proportion R²),
  `bench_zero_sweep.py` (fast single-proportion RMSE-vs-IsoSpec sweep; percolator
  + mzTab intake), `run_mixconfirm.py`, `run_integrate_invivo.py`.

#### Removed

- `baseline="linear"` (Skyline-style local-linear) — it over-subtracts on the
  narrow on-peak boundaries and failed on every line; rationale retained in
  `PROJECT_REVIEW.md` §2c.

---

## [0.9.0] — Unreleased

This is a stabilization release. It fixes correctness defects on the existing
algorithm, modernizes packaging, and adds CI — but does not yet change the
underlying integration or fitting algorithms. A larger, breaking 1.0 release
will introduce a new package structure, peak detection, baseline subtraction,
Qt GUI, and mzTab intake; see `PROJECT_REVIEW.md` for the roadmap.

### Changed (BREAKING)

- **Mass tolerance now means ±N ppm half-width.** `-m N` / `--mass_tol N`
  previously divided the requested ppm by 2 before applying as a half-width,
  so `-m 50` actually integrated only `±12.5` ppm around the theoretical m/z.
  As of 0.9.0, `-m 50` integrates `±50` ppm (a 100 ppm wide window). This
  doubles the effective mass window vs. 0.8.x. To reproduce 0.8.x integrated
  values, halve your `-m` argument.
- Empty results, missing mzML files, unreadable Percolator files, and
  out-of-range fraction indices now raise typed `RianaError` subclasses
  (`DataError`, `IntegrationError`, `ModelingError`) instead of calling
  `sys.exit` or raising bare `AssertionError`/`Exception`. Scripts that
  caught `SystemExit` will need to catch `riana.exceptions.RianaError`.
- The `riana preprocess` subcommand is removed. The implementation was
  inoperative (referenced removed NumPy APIs); MBR will return in a later
  release with a working implementation.

### Fixed

- `tqdm` progress bars now show the correct total (was `max(range(N)) = N-1`,
  off by one).
- `riana fit --plotcurves` no longer crashes — `plot_model()` was being called
  with the wrong keyword argument (`model=` vs. `model_to_use=`).
- `plot_model()` no longer leaks plot commands to the active matplotlib axes
  (it now draws onto its own `Figure` consistently).
- `except ValueError or IndexError:` clauses (which only caught `ValueError`)
  replaced with `except (ValueError, IndexError):` in both `riana_integrate.py`
  and `riana_fit.py`.
- `flanking aa` column in standalone-Percolator parsing previously stored the
  same scalar in every row due to a `[1]` index typo; now stores the per-row
  flanking residues.
- Error message in `get_isotopomer_intensity` no longer raises
  `UnboundLocalError` when integration fails in the `use_range=False` branch.
- `np.trapz` (deprecated in NumPy 2.0) replaced with `np.trapezoid`.
- `np.int` (removed in NumPy 1.20) usage eliminated alongside dead-code
  deletion.
- Logger is now cached by `(name, out_path)` instead of `name` alone, so
  successive runs writing to different output directories no longer share a
  file handler pointing at the first run's directory.
- Crux vs. standalone Percolator format is now detected by inspecting the
  file header explicitly, not by catching a `KeyError` while parsing.
- Standalone-Percolator `PSMId` parsing uses an anchored regex against the
  MSFragger format `filename.scan.scan.charge_index` instead of repeated
  string `.split()` indexing.
- Raw-string regex literals throughout to silence `SyntaxWarning:
  invalid escape sequence` introduced in Python 3.12.

### Added

- `riana.exceptions` module with `RianaError`, `DataError`,
  `IntegrationError`, `ModelingError`, `ValidationError`.
- GitHub Actions CI running pytest + coverage on Python 3.10–3.12.
- `[project.optional-dependencies] dev` extra installs pytest + pytest-cov
  (`pip install -e ".[dev]"`).

### Removed

- `riana/riana_preprocess.py` (mostly commented-out scaffolding; the
  `riana preprocess` subcommand was non-functional).
- `documentation/PROJECT_EVALUATION.md`, `documentation/ROADMAP.md`,
  `documentation/MASS_ACCURACY_SPECIFICATION.md`, `documentation/README.md`,
  and root `TODO.md`. Their content is consolidated into `PROJECT_REVIEW.md`.
- `setup.py` (project metadata fully moved to `pyproject.toml`, PEP 621).
- `rx` and `snakemake` as declared runtime dependencies (rx was only used by
  the broken Tk GUI; snakemake is a workflow-runner users install themselves
  if they choose to use `workflow/Snakefile`).

### Repackaged

- `pyproject.toml` (PEP 621) is the single source of project metadata and
  version, dynamically read from `riana.__version__`.
- Minimum Python is now 3.10 (was 3.9). Tested on 3.10, 3.11, 3.12.

---

Prior to 0.9.0 the project used `CHANGES.txt`. Entries below are copied
verbatim from that file for posterity.

## [0.8.2]

* Added the `-X`, `--ignored_mods` argument to `riana integrate` for
  modification(s) to ignore in the search result for calculating true
  peptide mass. Must match the exact string in the search-engine output.
* Added the `-F`, `--forced_mods` argument to `riana integrate` for
  modification(s) to always look for in the search result. Useful for
  amino-acid labeling experiments. `riana fit` only performs curve-fitting
  for unmodified peptides.
* The mass-defect argument has been replaced by `-D` `--mass_difference`,
  a float value specifying the exact mass difference between each
  isotopomer. Defaults to 1.003354835 Da (C12 vs. C13).

## [0.8.1]

* Minor bug fixes.

## [0.8.0]

* Added support for automatic isotopomer selection (`--fs Auto`).
* Added experimental support for in-vitro heavy-water labeling
  (`--label hw_cell`).

## [0.7.3]

* Added graphical user interface (GUI). Launch with `riana/riana_ui.py`.

## [0.7.2]

* Included additional Snakemake options in the config files.
* The isotope argument now takes multiple numerical values separated by
  spaces, e.g. `2 4 6 8`.

## [0.7.1]

* Added the `-D` mass-defect parameter to `riana integrate`.
* Added the `-p` flag to `riana fit` to toggle plotting fitted curves.
* Added the `-w` flag to `riana integrate` to write pre-integration
  intensities.

## [0.7.0]

* Riana now supports a Snakemake pipeline performing protein database
  search, filtering, integration and fitting.
* Native support for curve-fitting and best-fit-curve plots.
* Added one-pool, two-compartment (Guan et al.) and three-exponent
  (Fornasiero et al.) models.
* Match-between-runs temporarily disabled to support the Snakemake
  workflow; will return in a later version.

## [0.6.4]

* Added support for standalone Percolator results for v.3.0.5.
* Project layout: one `psms.txt` for target PSMs per percolator folder.
* When using standalone Percolator, peptide masses are calculated de novo.
* Toggle match-between-runs with `--mbr` / `-b`.
* No longer writes results of individual fractions separately.

## [0.6.3]

* Fixed an issue where MBR all-NaN slice caused an error.
* Restructured project directory: mzML files must be in an `mzml/`
  subfolder, Percolator files in a `percolator/` subfolder.

## [0.6.0]

* Started implementation of match-between-runs; restructured project for
  distribution via pip.

## [0.5.0]

* Updated to use pymzml 2.2.
* Multi-threading via `--thread`.
* Riana now loads all spectra needed for integration into memory.
* Substantial speed gain vs. 0.4.0; can finish a sizeable fraction
  (~1 GB raw) in 10 min.
* User-definable MS1 mass tolerance via `--masstolerance` (default 100 ppm).

## [0.4.0]

* Python 3.5+, up-to-date scipy/numpy/pymzml.
* Multi-fraction runs supported — mzML files in the directory must be in
  the same order as the Percolator indices.
* Riana now takes Percolator tab-delimited files instead of mzid (mzid
  support will return in a future version).
