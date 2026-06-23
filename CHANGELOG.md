# Changelog

All notable changes to Riana are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [1.1.0] — Unreleased

Post-1.0 development line. Track B–E follow-ons (per-peptide `--fs` keyed on the
adaptive envelope, MBR tuning on the standing calibration benchmark, …).

## [1.0.0] — 2026-06-23

The breaking 1.0 release: a new package structure, peak detection, baseline
subtraction, mzTab intake, and a Qt GUI. See `PROJECT_REVIEW.md` §3 for the
roadmap. Entries below are grouped by the work that produced them. (The git tag
and Zenodo code DOI follow at release.)

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
