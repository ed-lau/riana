# Riana — Project Review & Roadmap

> Status: M4 Phase 2 done (1.0.0.dev1 — PySide6 GUI: `riana gui` Integrate +
> Model tabs run end-to-end via qasync + a ProcessPoolExecutor; Integrate has an
> inline drift summary and a pyqtgraph chromatogram, Model has fitted-curve
> plots; deps ship as a `[gui]` extra). M4 Phase 1 done before it (Typer CLI
> replaces argparse; the typed pipeline is the only engine, `--engine legacy`
> removed; legacy modules + Tk GUI deleted; `riana fit` requires
> `--coefficients`, with bundled presets). M3 done before that (restructure,
> streaming I/O, dual ID intake, apex-default peak integration, IsoSpec
> forward-model FS). The post-M4 planning round is **done** (2026-06-07): the
> record-only M5–M8 numbering is retired and reorganized into five tracks (A
> I/O & run-identity model, B integration fidelity, C fitting science, D
> validation infra, E GUI/UX) — see §3 "Post-M4 roadmap". **Pre-1.0.0 chores +
> M6a (the run-identity model) are now done on branch `m3-rewrite`** (version →
> 1.0.0, Snakefile retired, CI Quarto docs deploy, `RunIdentity` + `io/sdrf.py` +
> `io/manifest.py` + `core/pipeline.py`, `integrate --sdrf` / `fit --manifest`).
> Track D shipped (within-protein-θ benchmark + frozen LVE bench set), which
> surfaced + fixed two integration bugs: the quantms **filename-prefix**
> scan-scramble (de-prefix mzML names) and **mass_tol 50→10 ppm** (now read from
> the SDRF; centroid search-tolerance, not a profile-width window). Engine also
> gained a best-q apex anchor, `apex_search_half_width=0.25`, bounded file-
> parallelism (`--workers`), and a profile-mzML intake warning. The **intake
> scan↔RT guard is shipped** (2026-06-10): per run, `integrate_run` errors if the
> median |spectra_ref-scan→mzML-RT − mzTab-RT| exceeds 2 min (the prefix-scramble
> signature; ≤~0.9 min on clean ProteomicsLFQ output), with `--no-rt-check` to
> override. **M5 (per-timepoint fraction-new) is shipped** (2026-06-10): `riana
> fit` now also writes `riana_fit_fractions.txt` — long-format `(concat, protein,
> biorep, t, fs, fs_lower, fs_upper)` with prediction-interval bounds from a
> unified residual bootstrap; this is the substrate the protein rollup
> inverse-variance-weights. **Track C protein rollup — `unique` first cut shipped**
> (2026-06-10): `riana rollup` writes `riana_rollup_proteins.txt` (median-of-peptide-k +
> the biorep-aware per-timepoint weighted refit), parsimony is now a
> summarize-time decision (`integrate --unique` removed so integrate extracts all
> peptides). **`--parsimony isoform` + the GUI Protein tab are shipped too**
> (2026-06-10): `isoform` folds isoform-only-shared peptides into the canonical
> entry unless an isoform carries its own unique peptide (a dataset-wide pass,
> adapted from `02_R_parsimony_reference.Rmd` with Riana-`,` / NA-fallback / `-N`+`-J`-suffix
> fixes + an optional `--min-r2` gate); the GUI grew a 3rd **Protein** tab over
> the same `rollup_proteins` core, with a per-protein refit curve on row select.
> **M6b (DIA-NN parquet intake) is SHIPPED** (2026-06-11): `io/diann.py` reads the
> DIA-NN `report.parquet` (validated on 2.5.0), `plan_integration` dispatches the
> reader on `SdrfTable.acquisition`, and DIA's apex RT is resolved to the nearest
> MS1 scan (`resolve_rt_anchored_scans`) so the scan-based extraction runs
> unchanged. Validated end-to-end on a cardiac in-vivo DIA set (98.5% m0
> coverage, ~1–2 ppm mass accuracy; fit R²med ~0.42 — lower than DDA's ~0.87, as
> expected for MS1-on-DIA, but k_deg median 0.08/day is biologically sound).
> Also shipped alongside: **`riana fit -W/--workers`** (process-pool fit — the
> GIL-bound IsoSpec/bootstrap fit went ~1.25→~7.8 effective cores; deterministic
> per-concat seed).
> **Track E GUI rewiring shipped** (2026-06-10): `integrate_project` split into
> plan/dispatch/finalize so the GUI Integrate tab drives the same per-run unit
> over its own pool (Workers spinbox, no nested pools) + an SDRF path; GUI Model
> gained a manifest path (`fit_project`) — integrate/fit now go through
> `core/pipeline` on both surfaces. **Also shipped (2026-06-11):** the **manifest
> project chain** (`fit`/`rollup --manifest` root at the manifest dir + write
> `stage="fit"`/`"rollup"` rows, each with a `created_at`); rollup **`--method
> {weighted, pooled}`** (point estimators dropped; the median stays as a
> `peptide_median_k` column; `n_replicates`/`n_timepoints` added); the GUI **SDRF
> mass-tolerance reflect**; **output hygiene** (slim `riana_fit_peptides.txt`,
> ~6-sig-fig estimate outputs, `riana_protein.txt`→`riana_rollup_proteins.txt` + a
> tidy `riana_rollup_fractions.txt`); and **crash-resilient integrate** (per-run
> incremental writes + `integrate --resume`). **Shipped 2026-06-12/13 (this
> session):** the **`io/mztab` variable-mod drop** (mirrors `io/diann`, so the
> PTM-search `timeseries_lve_atr` set integrates its unmodified peptidoforms
> cleanly); **rollup `-W/--workers`** (process pool, ~7.2× cores, byte-identical);
> the **`linear simple` rollup model** — the cross-sample Δk milestone: φ=log(1−θ)
> through-origin OLS fit jointly across a protein's conditions → per-condition k +
> a **Δk test** (statsmodels, analytic CIs) + Benjamini-Hochberg, with linear-only
> plateau truncation (`--phi-limit`); validated on `data/timeseries_lve_atr`
> (atrium vs ventricle: 780 proteins both chambers, atrium 1.24× faster, 410
> significant at BH p_adj<0.05); the **GUI φ-space plot + CI ribbons** (and `-W` in
> the fit/rollup tabs); and the **removal of `-t/--thread`** everywhere (GIL-bound,
> no speedup — `-W` is the sole parallelism lever). **Next:** **M7 PTM-aware
> envelope** (now the headline gap — thread parsed mod masses through the IsoSpec
> envelope so modified peptidoforms integrate at the *correct* m/z, plus
> proteoform-aware rollup keys; design in §3 Track C + `m7_ptm_envelope_design`),
> then the parked/short-term items (see §3 "Handoff — ordered priorities"). Track
> B fidelity gateable by Track D.
> Maintainer: Edward Lau. Last reviewed: 2026-06-13.

This document consolidates and supersedes the prior `documentation/` folder
(`PROJECT_EVALUATION.md`, `ROADMAP.md`, `MASS_ACCURACY_SPECIFICATION.md`). The
old contents are folded in below where still load-bearing; the rest is cut.

## 1. Project status

> **Update (2026-06, 1.0.0.dev1):** the three-step plan below is complete
> through M4 — 0.9.0 stabilization, the calibration dataset, and the 1.0.0
> rewrite (typed pipeline, streaming I/O, peak detection, mzTab intake, and a
> PySide6 GUI) have all shipped. The narrative below is the original pre-rewrite
> evaluation that motivated the plan, kept for context; the as-built record is
> CHANGELOG `[1.0.0]` + §3. The post-M4 planning round is done — the next
> batch of work is organized into five tracks in §3 "Post-M4 roadmap".

Riana is a single-author scientific Python tool for extracting and modeling
isotopomer time-series from MS1 data, used in protein-turnover research. The
core scientific functions — accurate-mass calculation, kinetic models,
fractional-synthesis math, IsoSpec wrapping — are correct and proven. The
surrounding software (data parsing, integration pipeline, Tkinter GUI,
packaging, tests) shows typical scientist-coder accumulation: tight coupling,
weak tests, dead code, brittle string parsing.

Rather than incrementally patch the structure, the plan is:

1. **0.9.0 — Stabilization.** Fix correctness defects in place. Modernize
   packaging and add CI. No restructure.
2. **Calibration dataset** (parallel work). A 9-level D₂O mixing series
   becomes the regression check for every future change.
3. **1.0.0 — Aggressive rewrite.** New package structure, typed records,
   streaming mzML, async pipeline, Skyline-style peak detection, mzTab
   intake, PySide6 GUI. The numerically-sensitive scientific functions are
   **lifted unchanged** rather than rewritten.

The 0.9.0 → 1.0.0 boundary is a deliberate breaking change. The user base is
small and scientific software is allowed to evolve.

## 2. Critical findings (post-evaluation)

### 2a. Correctness defects fixed in 0.9.0

See `CHANGELOG.md` for the full enumerated list (B1–B13 plus the
mass-tolerance semantic change). Highlights:

- Mass tolerance was applied at half the requested ppm (`/2` in the
  `delta_mass` computation). Now `-m N` means `±N` ppm.
- `riana fit --plotcurves` crashed because `plot_model()` was being called
  with the wrong keyword argument.
- The `flanking aa` column in standalone-Percolator output stored the same
  scalar for every row.
- `except ValueError or IndexError:` only ever caught `ValueError`.
- Logger was keyed by name alone, so subsequent runs wrote to the first
  run's output directory.

### 2b. Scientific defects deferred to 1.0.0

- **Amino-acid-labelling `a_max` path is unreachable.** `fsynthesis.py`
  checks `label == 'aa'` (string), but `riana_fit.py` dispatches with
  `label = int`, so AA experiments silently fall into the natural-abundance
  branch.
- **Fractional-synthesis denominator drift.** `iso0 / colsums` silently
  changes meaning depending on the `-i` choice passed upstream to
  `integrate`. Same data + different `-i` ⇒ different `mi`.
- **Kinetic-fit uncertainty.** Reported `sd` is `sqrt(diag(pcov))` of
  `k_deg` only; the plotted confidence band uses a heuristic
  (`k_deg ** 2 / (k_deg + sd)` as "lower bound") that is not a defined
  statistical CI. Use proper bootstrap CIs in 1.0.
- **Fixed labelling-site model biases fractional synthesis.** `riana fit`
  derives FS from the monoisotopic peak alone — `calculate_fs_m0` takes
  `mi = iso0/colsums`, then inverts the analytic relation
  `a_max = a_0·(1−ria_max)^n` with a site count `n` from `calculate_label_n`,
  a fixed per-peptide model. The M3 Week 0 `bench_fit_recovery.py` baseline
  (2026-05-20) shows this yields FS ≈ 0.6·f on the calibration mixing
  series: pseudo-time `k` recovery is biased low by ≈ −0.5 (median fitted
  `k` 0.24 AC16 / 0.28 iPSC vs target `k_deg₀` 0.5; stable across both cell
  lines and the R²≥0.9 subset). M2's `bench_fs_recovery` confirms FS ≈ f is
  recoverable (bias −0.017) once the per-peptide Spep is fitted with the
  IsoSpec forward model — so the gap is the fixed site-count model, not
  integration. Fix in M3 Week 4 (below).

These are deferred to 1.0.0 because each interacts with the planned
data-model rewrite (typed records, dataclass-based config) and is cleaner
to fix there than to patch in place.

### 2c. Algorithmic feature gaps — not bugs, but quantifiable shortcomings

These are not defects in what the algorithm does; they are limits in how
sophisticated the algorithm is. Whether they actually matter for Riana's
outputs is what the calibration dataset (Milestone 2) is designed to
measure.

1. **No chromatographic peak detection.** Integration is over a fixed RT
   window: every MS1 scan within `±r_time` of the PSM span is summed, then
   `np.trapezoid` over RT. There is no boundary detection. Outcome:
   integrated value includes anything in window — co-eluting peptides,
   baseline, tail of neighboring isotopologue.
2. **No background / baseline subtraction.** Same root cause: by integrating
   a rectangle, baseline is included proportionally to RT width.
3. **`polyorder=1` Savitzky-Golay smoothing** is mathematically equivalent
   to a moving average and distorts peak heights. The smoothed trace then
   directly feeds the area integral, so toggling `-S` changes the integrated
   value.
4. **No mass-domain refinement.** Centroid-summing in window with no
   observed-mass tracking throws away mass-accuracy information that would
   diagnose calibration drift and confirm correct peak assignment.

#### Recommended approach to (2c)

The reference target is Skyline ([peak picking documentation](https://skyline.ms/wiki/home/software/Skyline/page.view?name=tip_peak_calc)).
For 1.0.0, an `algorithms/peaks.py` module will implement:

- **Peak detection:** `scipy.signal.find_peaks` with prominence threshold
  on each isotopomer's XIC; the PSM scan provides a strong RT prior.
- **Co-elution grouping:** isotopomers of the same peptide should co-elute;
  if iso0's apex is more than e.g. 2× MS1 cycle time away from iso1's apex,
  flag the peptide and fall back to fixed-window integration.
- **Boundary determination:** `scipy.signal.peak_widths` at `rel_height=0.05`
  (95% of apex below peak).
- **Baseline:** local linear between detected boundaries (Skyline default),
  with SNIP and AsLS as benchmark alternatives via `pybaselines`.
- **Quality scoring:** S/N from baseline residual + peak symmetry.

Whether each of these improvements actually wins is decided by the
calibration dataset, not by intuition.

**Spike outcome (2026-06-05) — these were tested; the calibration data
decided** (full record in the `m3-peak-detection-spike` memory; pending the
full-mixing confirm):

- **Width is the primary lever.** A narrow integration window
  (`integration_half_width` ≈ 0.1–0.2 min) beats the wide 0.9.0 ±0.33
  rectangle — on the model-free mixing-linearity metric and on 0%-envelope
  RMSE vs IsoSpec, across all three lines (ac16/ipsc/cm). It removes
  background by *exclusion* rather than subtraction.
- **Apex-centring (`peak_rt="apex"`) is a real but second-order win** — it
  only matters once the window is narrow (a tight window must sit on the apex
  or it clips); at the wide width it's a wash. The apex finder uses a
  prominence gate + nearest-to-MS2-RT selection (not "tallest").
- **FWHM (`width_rel_height=0.5`) beats 5%** for detected boundaries, but a
  well-chosen fixed narrow width beats FWHM-auto (which over-clips long
  peptides). The optimum is mildly line-dependent (0.1–0.2) → length-adaptive
  width is the open refinement (M8).
- **`linear` (Skyline local-linear between boundaries) is discarded.** It is
  the correct Skyline algorithm but assumes boundaries at the chromatographic
  *floor*; our narrow on-peak boundaries make it over-subtract real signal,
  catastrophically with tight windows (worst row on every line). Removed from
  `algorithms/baseline.py` 2026-06-05. `none` is the default; `noise_floor`
  is competitive-but-no-op at low labelling, kept for future tuning (its
  quality leaks the extraction width — off-peak-estimate TODO in baseline.py).
- The config field `r_time` was renamed `extraction_half_width` (it conflated
  extraction with integration; those are distinct knobs now).

### 2d. Architecture findings (addressed in 1.0.0 rewrite)

- CLI (`main.py`) and GUI (`riana_ui/`) duplicate validation/config logic
  with non-identical types — the surfaces drift.
- `integrate_all(args)` is 230 lines doing six things; the GUI calls it
  directly and blocks the Tk main loop. The async dispatch with `rx` does
  not actually move work off the main thread.
- The Tkinter GUI references a `console` module that is not in the tree,
  so the GUI is broken on a fresh clone. Replaced by PySide6 in 1.0.
- `logger.py` global dict + `__init__.py` glob-import + `params.py`
  module-level globals together prevent isolated test runs.

## 3. Roadmap

Sequencing: M1 (0.9.0), M2 (calibration dataset + benchmarks), M3 (restructure
+ peak-detection spike), and M4 **Phase 1** (Typer CLI + `--engine legacy`
removal) have all shipped — see CHANGELOG `[1.0.0]`. M4 **Phase 2** (PySide6 +
async Qt GUI) is shipped too: the `riana gui` Integrate and Model tabs both run
end-to-end. The post-M4 planning round is done — the next batch of work is the
five-track "Post-M4 roadmap" at the end of this section (it supersedes the old
record-only M5–M8).

### M1 — Stabilization → 0.9.0 — DONE

Smallest viable bugfix release on the existing layout. Fixes B1–B13,
modernizes packaging, adds CI. Behavior-preserving except for the
documented mass-tolerance semantic correction. See `CHANGELOG.md`.

### M2 — Calibration test dataset (acquisition + dual-ID processing) — DONE

**What shipped (2026-05-17).** Benchmark scaffolding lives in
`tests/benchmark/`; calibration artifacts in
`tests/data/calibration_d2o_mixing/{ac16,ipsc,cm}/` (`cm/` was added during M3
— see the M2 addendum below). The approach diverged from
the original plan below — there is no external ground truth for per-peptide
isotope envelopes (no animal calibration curve), so "predicted vs observed
m0/mA" is not directly scorable. Instead the benchmark ports NB87a
(`data/notebook/87a_…IsoSpec_AC16.ipynb`): integrate output → per-peptide
Spep via an IsoSpec forward model → per-AA non-negative regression →
fractional-synthesis recovery vs the nominal mixing proportion. The escape
from circularity is coefficient *stability* across N_ISO, smoothing, and cell
line.

Scripts: `build_ground_truth.py`, `bench_aa_coefficients.py`,
`bench_fs_recovery.py`, `bench_n_iso_sweep.py`, `bench_smoothing.py`,
`run_integrate_v0_9_0.py`, plus `_helpers/forward_model.py`. The
aa-coefficient and N_ISO-sweep ports reproduce NB87a bit-exact (Δ ≤ 6e-7).

v0.9.0 baseline (committed under `benchmark_results/v0.9.0/`):

- AC16: 1512 peptides, coeff train/test R² 0.86/0.82, FS bias −0.017.
- iPSC: 2964 peptides, coeff train/test R² 0.77/0.76, FS bias −0.018.
- N_ISO sweep: test-R² peaks at N_ISO=4 for both lines.
- Smoothing sweep: per-AA coefficients shift with window size (≤0.12 AC16,
  ≤0.22 iPSC at S=9) but FS-recovery bias is insensitive to smoothing
  (±0.002) — confirms §2c point 3 (SG distorts areas) while showing the
  calibration verdict is robust.
- Cross-line: AC16 vs iPSC coefficients differ by up to 0.54 (Met) — *not*
  cleanly transferable; revisit before adopting a single frozen table.

Deferred from M2: the mzTab→RIANA adapter (M3); a Zenodo deposit (raw `.raw`
files are already citable on JPOST — `JPST002443` AC16, `JPST003556` iPSC —
so no separate deposit is needed); `pseudotime_map.csv` and the fit-module
benchmark (`riana fit` correctness work interacts with the M3 rewrite and is
deferred to post-M3/M4). `integrate_outputs/` (~450 MB of `_riana.txt`) is
gitignored — regenerable via `run_integrate_v0_9_0.py`.

Other notes from the baseline run: v0.9.0 was run with `-m 15` (±15 ppm under
the 0.9.0 semantic); the snakemake-era reference effectively used ±7.5 ppm due
to the pre-0.9.0 `/2` bug, so the two are not strictly mass-window-matched.
The bundled `config_template.yaml` passes `-D D` (a label) which no longer
parses — the current CLI takes a float isotopomer mass-step; `run_integrate`
omits `-D` and uses the default 1.003354835.

**M2 findings → M3 implications.**

1. *Median FS bias is a blunt metric for peak-detection work.* It moved only
   ±0.002 across the entire smoothing sweep. An M3 benchmark that judges peak
   detection by median FS bias will see almost nothing. The metrics that
   actually move with integration quality are: the **spread** of FS recovery
   around each nominal proportion (RMSE/IQR, not the median — `fs_recovery.csv`
   already carries every per-(peptide, fraction) row); per-AA coefficient
   std errors and R²; the count of peptides passing the R²>0.95 curation gate;
   and the *shape* of the N_ISO sweep (if peak detection cleans iso5/iso6, the
   post-N_ISO=4 R² degradation should flatten).
2. *The curation filter hides peak-detection's main win.* The R²>0.95 gate
   discards exactly the co-eluting / low-SNR peptides where peak detection
   helps most. M3 benchmarks (`bench_peak_boundary.py`, `bench_baseline.py`)
   **must also report metrics on the uncurated population**, or they will
   systematically understate the improvement.
3. *Per-cell-line frozen coefficient tables are viable* — the original
   "predicted vs observed m0/mA" target is recoverable, just per-line not
   universal. Within-line coefficient drift is small (≤0.09 AC16 across
   integrate versions, ≤0.12 across the smoothing sweep). A frozen table is a
   **constant**, so its bias cancels in any method-vs-method comparison: scoring
   integration A and B against the same frozen predicted m0/mA preserves the
   relative ranking, which is what regression-gating needs. Plan: bootstrap one
   table per cell line from the best available integration, freeze + version it
   (`d2o_aa_coefficients_<line>.csv`), and add `bench_m0_ma_recovery.py` scoring
   observed-vs-predicted m0/mA RMSE — the sensitive per-peptide metric (1) calls
   for. Caveats: the iPSC table is noisier (R² 0.77 vs 0.86), so its absolute
   numbers are less trustworthy though still usable as a constant reference; a
   new cell type needs its own re-derived table. Bonus diagnostic: if M3's
   better integration makes the AC16 and iPSC tables *converge*, that is
   evidence the 0.54 cross-line divergence was an integration artifact rather
   than real biology.

**M2 addendum — third calibration line (2026-05-22).** A third D₂O mixing
series was acquired and wired in under `data/calibration_cm/` and
`tests/data/calibration_d2o_mixing/cm/`: contractile human iPSC-derived
cardiomyocytes (iPSC-CM), 9 proportions, JPOST `JPST003582`. It is processed
identically to ac16 and ipsc — same Crux+Percolator IDs, same pinned `-m 15`
integrate config — and re-uses every benchmark script via `--line cm`.

Rationale: ac16 and ipsc are both **proliferative**, and a dividing cell
dilutes isotopic label through division independently of protein turnover;
iPSC-CM is **post-mitotic**, isolating turnover from division, and is the more
physiologically relevant model. The third line turns M2 finding 3's two-point
cross-line comparison into a three-point one with a biology axis.

*v0.9.0 baseline, all 9 proportions:*

| line | curated n | OOB R² | m0_rmse curated | m0_rmse uncurated |
|------|-----------|--------|-----------------|-------------------|
| ac16 | 1512      | 0.848  | 0.018           | 0.062             |
| ipsc | 2964      | 0.766  | 0.024           | 0.072             |
| cm   |  564      | 0.754  | 0.025           | 0.074             |

cm came back the noisiest line: its `time50` fraction is a weak acquisition
(16,090 vs ~20k target PSMs; 5,671 vs ~13–16k integrated peptides at `-q 0.01`).
Because NB87a curation requires a peptide observed at all 9 proportions,
`time50` alone bottlenecks cm to 564 curated peptides — a third of ipsc's — and
its frozen table inflates near-zero-labeling residues (Lys 0.235, Tyr 0.258,
Phe 0.40) where a starved 20-parameter fit absorbs integration noise.

*Relaxed-coverage A/B (drop `time50`).* The "observed at all 9" rule is a
**coverage** requirement, separate from the R²>0.95 quality gate; `time50` is a
bad *run*, not a bad set of peptides. Re-curating cm on the 8 surviving
proportions (`--drop-proportion 50`, new flag on `bench_aa_coefficients.py` /
`bench_m0_ma_recovery.py`) gives:

| cm variant | curated n | OOB R² | mean coef boot-std | mean \|Δ\| vs ac16 / ipsc |
|------------|-----------|--------|--------------------|---------------------------|
| 9/9        |  564      | 0.754  | 0.147              | 0.172 / 0.194             |
| 8/9 (−t50) | 1817      | 0.782  | 0.094              | 0.165 / 0.157             |

Dropping one weak fraction recovers 3.2× the peptides, cuts coefficient
bootstrap noise by 36% (0.147→0.094, toward ipsc's 0.064), and collapses the
unphysical low-labeling coefficients (Lys 0.235→0.022, Tyr 0.258→0.065). The
table also moves *toward* both proliferative lines — most toward ipsc, cm's
parental line. The cm self-shift 9/9→8/9 is 0.11 mean (0.26 max), comparable to
the cross-line deltas themselves: **a large part of cm-9/9's apparent
cross-line divergence was small-N noise, not biology.** Both tables are kept
(`d2o_aa_coefficients_cm.csv` = 9/9; `d2o_aa_coefficients_cm_drop50.csv` = 8/9,
the recommended cm reference); the A/B is the record.

Implications: (a) M2 finding 3's "do the frozen tables converge under better
integration" diagnostic now has three legs, and cm-8/9 already sits closest to
ipsc; (b) this is concrete input for the M3 curation-gate revisit — the
coverage rule should tolerate a known-bad fraction rather than discard every
peptide missing from it; (c) cm remains the noisiest line even at 8/9 (boot-std
still ~1.5× ipsc), so it is the sharpest stress test for whether M3 peak
detection tightens the per-AA fit.

The original M2 plan is retained below for reference.


Goal: a ground-truth benchmark dataset that survives every later
milestone. MS data already exists (parallel project); the outstanding
work is re-searching with two ID pipelines and building benchmark
scaffolding.

**Experimental design (already acquired).** Cells cultured in 6% D₂O for
≥10 doublings → effectively complete proteome labelling. Lysate from
labelled cells mixed with lysate from unlabelled cells at 9 nominal
heavy fractions: `0%, 12.5%, 25%, 37.5%, 50%, 62.5%, 75%, 87.5%, 100%`.
Multiple technical replicates per level.

**Why it's a uniquely good Riana benchmark:**

1. **Tests `integrate` directly.** For every peptide and every mixing
   fraction `f`, the expected isotopomer envelope is computable from the
   peptide sequence + 6% D₂O enrichment via `get_peptide_distribution()`.
   RMSE between observed and predicted m0/mA across the curve is the
   headline metric.
2. **Tests `fit` independently of biology.** Pseudo-time trick: choose
   `k_deg₀`, compute `t_i = -ln(1-fᵢ)/k_deg₀` per level, relabel samples
   `time<t_i>`, run `riana fit`. A correct fit module recovers `k_deg₀`
   for every peptide. One source of variability (integration) instead of
   the usual biological/kinetic cocktail.
3. **Edge cases come for free.** `0%` ⇒ false-positive m1+ floor.
   `100%` ⇒ direct check on `a_max = a_0·(1-ria_max)^n`.

**Outstanding tasks:**

1. Re-search the raw mzMLs through two pipelines, same FASTA and FDR:
   - Percolator path → `*.target.psms.txt` (matches 0.9.x).
   - mzTab path via quantms → `*.mzTab` (matches 1.0).
2. Build `ground_truth.csv` from the sample-prep spreadsheet.
3. Compute `predicted_distributions.csv` offline using the existing
   `get_peptide_distribution()` machinery.
4. Build benchmark scripts in `tests/benchmark/` (see below).

**Dataset layout:**

```
tests/data/calibration_d2o_mixing/
├── README.md                       # protocol, instrument, FASTA, search params
├── mzml/                           # raw mzML (.gz), HOSTED ON ZENODO
├── ids/
│   ├── percolator/                 # *.target.psms.txt
│   └── mztab/                      # *.mzTab from quantms
├── ground_truth.csv                # nominal heavy fraction per (peptide, replicate)
├── predicted_distributions.csv     # theoretical m0/mA per (peptide, fraction)
├── pseudotime_map.csv              # f → t for fit-module test
└── benchmark_results/
    └── <git_sha>/                  # frozen diffs across releases
```

Raw mzMLs are deposited to Zenodo under a dataset DOI separate from the
code DOI; repo carries CSVs + scripts + a `make calibration-data` target.

**Benchmark scripts:**

- `bench_integrate_recovery.py` — runs `riana integrate` on the series,
  compares observed vs. predicted m0/mA, reports RMSE, bias, R²
  per peptide and aggregated by length/labelling-sites/intensity bin.
- `bench_fit_recovery.py` — applies pseudo-time mapping, runs
  `riana fit`, reports `(k_deg - k_deg₀) / k_deg₀` per peptide. Repeat
  for Guan/Fornasiero by inserting synthetic precursor lag into the
  pseudo-time map.
- `bench_smoothing.py` — sweeps `--smoothing ∈ {None,3,5,7,9,11}`,
  reports RMSE shift. Quantifies whether SG is helping or hurting.
- `bench_peak_boundary.py` (post-M3) — fixed-window vs. detected-boundary
  vs. Skyline-style.
- `bench_baseline.py` (post-M3) — no baseline / linear / SNIP / AsLS.
- `bench_id_path.py` (post-M3) — Percolator-ID vs. mzTab-ID on identical
  raw data; per-peptide agreement.

**M2 deliverables:**

- Re-searched ID files (Percolator + mzTab) archived in
  `tests/data/calibration_d2o_mixing/ids/`.
- Zenodo-deposited raw data with DOI.
- `ground_truth.csv`, `predicted_distributions.csv`, `pseudotime_map.csv`.
- `bench_integrate_recovery.py` + `bench_fit_recovery.py` running on 0.9.0,
  baseline numbers committed to `benchmark_results/v0.9.0/`.

### M3 — Aggressive restructure → 1.0.0 — DONE

Delivered over Weeks 0–4 plus a pre-M4 peak-detection spike (2026-06). The
itemized record is in `CHANGELOG.md` (`[1.0.0]`); the high-level outcome:

- **New package layout + typed records.** `core/` (integration, fitting,
  models, fsynthesis), `algorithms/` (mass_calc, isotope_dist, peaks, baseline,
  calibration), `io/` (mzml, percolator, mztab, writers); frozen
  `IntegrationConfig`/`FitConfig` + `PSMRecord` etc.; science modules lifted
  unchanged.
- **Streaming I/O + dual ID intake.** Indexed/streaming mzML (one fraction in
  memory), typed Percolator parser, and quantms **mzTab** intake; provenance
  headers on every output.
- **Integration rewrite + peak detection** (then behind `riana integrate
  --engine new`, now the default and only engine): apex/consensus detection,
  baseline options, per-isotopomer mass-accuracy + drift. The typer/click CLI
  rewrite was deferred to M4 and shipped in M4 Phase 1 (`--engine` removed).
- **Fitting rewrite + §2b fixes** (`riana fit --engine new`): FS via the IsoSpec
  forward/solve model (per-peptide Spep + full-envelope least-squares), closing
  the ≈ −0.5 pseudo-time `k_deg` recovery bias; AA `a_max` dispatch,
  FS-denominator, and bootstrap CIs fixed.
- **Peak-detection spike → default integration changed.** Benchmark-gated on the
  D₂O calibration series (ac16/ipsc/cm, 0–100%) + an in-vivo mouse set: an
  **apex-centred narrow window** (`peak_rt="apex"`, `integration_half_width=0.15`,
  `apex_selection="tallest"`, `baseline="none"`) robustly beats the 0.9.0 fixed
  rectangle and generalizes across cell line, organism, and ID pipeline. Now the
  **default**; 0.9.0 reproducible via `--peak-rt ms2 --integration-half-width
  1.0`. `consensus` (median apex over m0..m3) is the high-D₂O alternative; all
  window/apex knobs are tunable. Rationale (incl. discarded `linear` baseline)
  in §2c; full record in the `m3-peak-detection-spike` memory.

**Regression gates met:** explicit-`ms2` `sample1` integration within 1e-3 of
0.9.0; `bench_fit_recovery` `k_deg₀` recovery; Percolator/mzTab agreement;
one-fraction memory ceiling. **Open refinements (deferred):** tight
`apex_search_half_width` (one untested lever); typer/click CLI (M4); adaptive
width / N_ISO (M8); Phase C v2 cross-proportion-stable picker.

<details><summary>Original M3 Week 0–4 plan, target layout, and verification
(as-planned; superseded by the summary above and CHANGELOG)</summary>

0. **Week 0 — benchmark infrastructure (no rewrite code).** M2 finding 3:
   bootstrap one frozen per-cell-line coefficient table from the best v0.9.0
   integration, freeze + version it as
   `d2o_aa_coefficients_<line>.csv`, and add `bench_m0_ma_recovery.py`
   scoring observed-vs-predicted m0/mA RMSE on **both the curated and the
   uncurated** (pre-R²>0.95-gate) peptide populations — M2 finding 2. Also
   land `bench_peak_boundary.py` and `bench_baseline.py` as *runnable
   stubs against v0.9.0* (fixed-window only for now) so Week 3 has a
   regression gate the moment peak detection is added; both report curated
   and uncurated metrics. Build `pseudotime_map.csv` + `bench_fit_recovery.py`
   here too (M2 deferred them, but Week 4 now needs them — see below). None
   of this depends on the rewrite, so it lands first.
1. **Week 1 — skeleton + lifts.** New layout (below). Lift `accmass`,
   `models`, `fsynthesis`, `constants`, `utils.get_peptide_distribution`
   into their new homes. Apply two science-layer bug fixes:
   - `plot_model` `model_to_use=` kwarg (already done in 0.9.0).
   - `fsynthesis` `label == 'aa'` → `label == 4`.
   Define typed records: `PSMRecord`, `Chromatogram`, `IsotopomerPeak`,
   `IntegrationConfig`, `FitConfig`.
2. **Week 2 — I/O layer (both ID paths).** `io/percolator.py` as a typed
   parser (no exception-as-control-flow; M2 already produces both formats
   so the mzTab adapter lands here too). `io/mztab.py` via
   `pyteomics.mztab`. `io/mzml.py` with indexed/streaming read via
   `pyteomics.mzml` — never hold more than one fraction in memory.
3. **Week 3 — core integration pipeline.** Rewrite `core/integration.py`
   against the new types. Add peak detection, baseline subtraction,
   mass-accuracy outputs inline. Validate each addition against the Week 0
   benchmarks before committing — judge by the *sensitive* metrics (M2
   finding 1: FS-recovery spread, per-AA R²/std errors, R²>0.95 gate count,
   N_ISO-sweep shape), not median FS bias. Re-run `bench_n_iso_sweep.py`
   after peak detection: if it cleans iso5/iso6, the post-N_ISO=4 R²
   degradation should flatten.
4. **Week 4 — fitting rewrite + §2b science fixes** *(re-scoped 2026-06-04
   after Week 3 shipped)*. Original spec also included a typer/click CLI
   rewrite; that is **deferred to M4** so it lands alongside the GUI's
   shared config-driven API surface, and so Week 4 stays focused on the
   load-bearing science change (closing the ≈ −0.5 k-recovery bias the
   Week 0 baseline documented). `--engine new` keeps the existing
   argparse surface through Week 4; the CLI rewrite happens once. The
   peak-detection engine revisit (Phase C v2 — cross-proportion-stable
   boundaries; the M3 boundary-stability finding) is its own
   planning effort post-fit-rewrite.

   Rewrite `core/fitting.py` to consume `IntegrationResult` records,
   **and apply the §2b scientific fixes**:
   - AA `a_max` dispatch (`label == 4`), FS-denominator drift, bootstrap
     kinetic-fit CIs.
   - **Replace the m0/mA-analytic FS calculation with the IsoSpec
     forward/solve model.** Drop `calculate_fs_m0` + `calculate_label_n`'s
     fixed site-count model (§2b) in favour of the per-peptide Spep fit +
     full-envelope `solve_fs` validated in M2 — the approach in
     `tests/benchmark/_helpers/forward_model.py` used by
     `bench_fs_recovery.py` (IsoSpec forward envelope at natural vs. Spep-
     labelled enrichment, FS by least-squares against the observed
     envelope rather than from `iso0` alone). The M2 forward model lifts
     into `algorithms/isotope_dist.py`; IsoSpecPy becomes a runtime
     dependency of `riana fit`, not just a benchmark one. This is what
     closes the ≈ −0.5 k-recovery bias the Week 0 `bench_fit_recovery.py`
     baseline records.
   Regression-gated by `bench_fit_recovery.py` (`k_deg₀` recovery). Output
   provenance header (git SHA, riana version, config hash) via
   `io/writers.py`. Concrete target: close the Week 0 baseline's median
   k_rel_err of −0.51 (ac16) / −0.44 (ipsc) toward 0 by replacing the
   fixed-site-count analytic FS with the per-peptide Spep + IsoSpec
   forward FS.

   Note: the fit fixes deliberately change fit output, so the
   `tests/data/sample1/` smoke test (below) can no longer demand bit-near
   identity for the fit stage — it gates *integration* m0/m6 only; fitting
   is gated by `bench_fit_recovery.py` recovering `k_deg₀`.

**Target layout:**

```
riana/
├── __init__.py
├── __main__.py
├── cli.py                   # typer/click dispatcher
├── exceptions.py            # already exists in 0.9.0
├── config.py                # frozen dataclasses for both CLI + GUI
├── records.py               # PSMRecord, Chromatogram, IsotopomerPeak, ...
├── pipeline.py              # async stage composition
├── core/
│   ├── integration.py       # NEW
│   ├── fitting.py           # NEW
│   ├── models.py            # LIFTED
│   └── fsynthesis.py        # LIFTED (fix label==int bug)
├── algorithms/
│   ├── mass_calc.py         # LIFTED from accmass.py
│   ├── isotope_dist.py      # LIFTED from utils.get_peptide_distribution
│   ├── peaks.py             # NEW: detection, boundaries, SNR, symmetry
│   ├── smoothing.py         # NEW: SG polyorder ≥ 2, AsLS, SNIP
│   └── calibration.py       # NEW: per-peak ppm error, drift summary
├── io/
│   ├── mzml.py              # NEW: indexed/streaming
│   ├── percolator.py        # NEW: typed parser
│   ├── mztab.py             # NEW: quantms intake
│   └── writers.py           # NEW: TSV/JSON with provenance
├── constants.py             # LIFTED
└── gui/                     # NEW (Milestone 4): PySide6
```

Deleted: `riana_ui/`, `riana/spectra.py` (replaced by `io/mzml.py`),
`riana/peptides.py` (replaced by `io/percolator.py`),
`riana/riana_integrate.py` (replaced by `core/integration.py` +
`pipeline.py`), `riana/riana_fit.py` (replaced by `core/fitting.py`),
`riana/project.py`.

**Mass-accuracy output (folded into M3, was a separate spec doc):**

Per-isotopomer columns: `iso{N}_obs_mz`, `iso{N}_ppm_error`,
`iso{N}_snr`, `iso{N}_quality`. Per-sample summary footer: median ppm,
MAD ppm, suggested calibration shift. CLI flag `--ppm-alert <ppm>`
(default 20) emits warnings via the logger when systematic drift exceeds
threshold. Optional `--json-out` for downstream tooling.

Skip: real-time monitoring class, calibration dashboard, automated
correction. Keep it simple.

**Verification (regression-gated by the Week 0 benchmarks):**

- Calibration benchmark holds or improves vs the committed v0.9.0 baseline.
  Judge by the *sensitive* metrics, not median FS bias (see M2 findings):
  FS-recovery spread, per-AA coefficient R²/std errors, peptide count through
  the R²>0.95 gate, N_ISO-sweep shape, and observed-vs-predicted m0/mA RMSE
  (via `bench_m0_ma_recovery.py` against the Week 0 frozen tables) on **both
  the curated and uncurated** peptide populations.
- M3 integration benchmarks must be run with the same `-m 15` mass window as
  the committed v0.9.0 baseline, or the comparison is invalid (the baseline
  is *not* mass-window-matched to the snakemake-era reference — see M2 notes).
- `bench_fit_recovery.py` recovers `k_deg₀` per peptide within tolerance on
  the pseudo-time-mapped series — this is the gate for the Week 4 fit rewrite
  and its §2b fixes.
- Percolator-ID and mzTab-ID paths produce m0/mA values that agree within
  tolerance (cross-format A/B is itself a validation of the mzTab adapter).
- `tests/data/sample1/` end-to-end smoke test produces a `_riana.txt` whose
  per-peptide *integration* m0/m6 agree with 0.9.0 within 1e-3 relative
  tolerance (near-identical — same numerical core). The fit stage is exempt:
  the §2b fixes change fit output by design (see Week 4 note).
- Memory peak on a 2 GB mzML drops from "all of it" to "one fraction
  worth."
- Bonus diagnostic (M2 finding 3, now three lines): if the better M3
  integration makes the ac16 / ipsc / cm frozen tables *converge* — and cm's
  residual low-labeling-residue inflation collapses — the cross-line divergence
  was an integration artifact, not biology. The M2 addendum already shows a
  large part of cm's divergence was small-N noise; M3 peak detection is the
  test of the rest. Score cm against `d2o_aa_coefficients_cm_drop50.csv`.

</details>

### M4 — Qt + CLI rewrite, legacy removal

**Phase 1 — Typer CLI + `--engine legacy` removal — DONE (2026-06-06).** The
typer/click CLI rewrite deferred from M3 Week 4 landed as `riana/cli.py` (Typer);
the argparse `main.py` and the whole legacy pipeline
(`riana_integrate`/`riana_fit` + the `accmass`/`fsynthesis`/`models` shims +
`spectra`/`peptides`/`project`) are deleted, along with the broken Tkinter
`riana_ui/`. The typed pipeline is the only engine — `--engine` is gone;
0.9.0 integration is reproduced with `--peak-rt ms2 --integration-half-width
1.0` (pinned to a committed golden within 1e-3). `riana fit` now **requires
`--coefficients`** (bundled presets `commerford`/`ac16`/`ipsc`/`cm` under
`riana/data/coefficients/`, or a path); `--label` collapsed to `{hw, o18}`
(cell specificity is the coefficient table, not the label), amino-acid/SILAC
fitting dropped (the SILAC `-X/-F` extraction knobs were later retired in M7
Stage A3), and `o18` is recognized
but errors pending its post-M4 rewrite. The A/B-against-legacy tests were
converted to committed-golden comparisons. See CHANGELOG `[1.0.0]` M4 Phase 1.

**Phase 2 — PySide6 + async GUI — DONE (2026-06-07).** PySide6 (LGPL) + `qasync`
(bridges asyncio with the Qt event loop) under `riana/gui/`. Long-running CPU
work via `ProcessPoolExecutor` driven from async tasks. `pyqtgraph` for fast
embedded chromatogram / fitted-curve inspection; matplotlib only for static
export. Entry point: a lazy `riana gui` subcommand in `cli.py`
(PySide6/qasync/pyqtgraph are imported only inside it, so they never load on the
core CLI path). GUI deps ship as a `[gui]` extra so the core CLI stays
lightweight. Both tabs share the Qt-free worker layer (`riana/gui/tasks.py`,
unit-tested without a display) and a `DataFrameTableModel` (`riana/gui/models.py`).

*Shipped (the vertical slice that proves the architecture):* the **Integrate**
tab runs end-to-end. Its form builds the *same* frozen `IntegrationConfig` from
widget values, so `__post_init__` is the single shared validator (CLI and GUI
cannot drift). Integration is awaited on the process pool one fraction at a time
(responsive UI, honest per-fraction progress) through the Qt-free workers in
`riana/gui/tasks.py` (which call the identical `core.integration.integrate_run`
→ numerics provably match the CLI; gated in `tests/test_gui.py` against the
`sample1` golden within 1e-3). Calibration is **folded into Integrate**: the
per-fraction `DriftSummary` (median/MAD ppm, suggested shift, `--ppm-alert`
flag) shows inline in the results panel — no separate Calibration tab.
Peptide-row selection draws the isotopomer XICs in a pyqtgraph
`ChromatogramView` with the integrated window shaded, backed by the new
`core.integration.extract_peptide_trace` / `PeptideTrace` (which also gives the
orphaned `records.Chromatogram` its first producer). The per-fraction
orchestration is *mirrored* from `cli.integrate` (not refactored) to keep the
tested CLI path untouched; a shared `core/pipeline.py` extraction is a noted
future cleanup.

The **Model** tab is the fit counterpart: its form builds the same frozen
`FitConfig`, and the fit runs as a single batched job on the pool via
`tasks.run_fit` (reads the per-timepoint `_riana.txt` files, loads the
`--coefficients` table, calls the shared `core.fitting.fit_run`), so output
matches `riana fit`. Selecting a result row plots that peptide's
`(t, fraction-new)` points + the fitted kinetic curve in a pyqtgraph `CurveView`
(`core.models` functions evaluated on the GUI thread — pure math, no worker
round-trip). `plot_curves` / `fs_formula` are not surfaced (no-ops in the new
fit engine; the interactive curve replaces the legacy `-p` static plots).
(`rx`, `sv_ttk`, `pandastable`, the missing `console` shim, and Tkinter are
already gone after Phase 1.)

*Fit/integrate consistency fix (2026-06-07, from the Phase 2 review):* the
`riana integrate --iso` default changed from the legacy `0 6` pair to the
contiguous **m0-m5** envelope, because `solve_fs_d2o` matches the observed
envelope against the IsoSpec forward model over contiguous channels from m0 — the
`0 6` pair silently misaligned (observed m6 vs predicted m1) → garbage `k_deg`.
`fit_run` now **guards** on the canonical set (`_REQUIRED_D2O_ISOTOPOMERS =
0..5`) and errors clearly if absent. `--plotcurves` (a no-op) was removed;
`--fs` is kept but ignored, reserved for a post-M4 *channel-subset* envelope SSE
(use fewer high isotopomers to dodge co-eluting contaminants — literature-backed;
deferred alongside the M8 adaptive-N_ISO work).

### Post-M4 roadmap (planning round, 2026-06-07)

M1–M4 shipped the rewrite (typed pipeline, streaming I/O, peak detection,
mzTab intake, Typer CLI, PySide6 GUI). This round retires the record-only
M5–M8 numbering and the scattered deferrals (Phase C v2, the animal
within-protein benchmark, the `--fs` subset, the o18 rewrite) and reorganizes
everything — plus substantial new feature work surfaced in the session — into
**five tracks** with named milestones and a sequence.

**The through-line.** Nearly every I/O pain point is downstream of one missing
abstraction: a **run-identity data model**

```
RunIdentity = (experiment, source/sample, condition/group,
               biological_replicate, technical_replicate, fraction,
               labeling_time, acquisition, precursor_enrichment)
```

attached to every PSM at intake and carried — *header-authoritatively* — to
fit and protein rollup. Today identity is positional and conventional
(`file_idx` as a sort index, timepoint encoded in the filename,
one-Percolator-file-per-sample by habit); that is the fragility. The SDRF is the
carrier for both acquisition modes — DDA via quantms and DIA via quantms-diann,
which emits its own SDRF (a near-identical variant, a few columns different) —
so the identity model is shared; only the PSM/quant file format differs (mzTab
vs DIA-NN parquet).

#### Handoff — ordered priorities (next sessions, as of 2026-06-13)

The Track C cross-sample Δk milestone (`linear simple`) is **done**; the headline
gap is now M7. Recommended order:

1. **M7 — PTM-aware envelope (next; the headline).** Load-bearing: variable-mod
   peptidoforms are currently *dropped* (the `io/mztab`/`io/diann`
   `drop_variable_mods` interim) because they would integrate at the unmodified
   m/z. Full design + next-tier roadmap in Track C below +
   `[[m7_ptm_envelope_design]]`; scope **locked 2026-06-13**. Stage it: (a)
   **forward-model atom accounting** — add each parsed mod's composition to the
   IsoSpec envelope *and* the integrate-side target m/z (two separate fixes),
   extending the atom vector `[C,H,O,N,S]` → `[C,H,O,N,S,P]` and routing the
   hardcoded Carbamidomethyl through a curated `UNIMOD:id → [C,H,O,N,S,P]` table;
   **starter set = phospho-STY + protein N-term Acetyl** (the two
   reliably-identifiable-in-*un*enriched-data mods); (b) **proteoform-aware rollup
   keys** (`P12345_pS235` for phospho; N-term Ac / Met-Ox / deamidation / CAM fold
   into the bare accession) — both mzTab (`start`+pos) and DIA-NN (`Protein.Sites`)
   already carry the protein-coordinate site, **no FASTA needed**. Also **retires
   `-X/--ignored_mods` + `-F/--forced_mods`** (SILAC-era). Build order A1 → A3 → A2
   → B. Substrate: `data/timeseries_lve_atr` (PTM search already in hand).
2. **Expose hidden integrate knobs as clearly-marked *advanced* options
   (short-term, cheap, independent).** Audit the full `IntegrationConfig` knob set
   (`--peak-rt`, `--apex-selection`, `--integration-half-width` vs
   `--extraction-half-width`, `--baseline`, `--smoothing`, …); surface them in CLI
   help + a collapsible "Advanced" group in the GUI forms (pairs with the
   smoothing-in-GUI Track E item). Good warm-up; improves the fidelity work.
3. **Cut a real `1.0.0` tag + repo hygiene (deliberate hygiene chunk).** A clean
   line in the sand before more features pile on (move 1 GB+ of personal outputs
   out of `data/`, drop stray root outputs — see Cross-cutting chores). Can follow
   M7 but shouldn't slip indefinitely.
4. **Match-between-runs re-explore for the mzTab/DDA path (mid-term).** Start with
   the *measurement* — compare per-peptide missingness DIA vs DDA (`runs/lve_dia`
   vs `runs/lve`) — before building anything; DIA-NN already propagates, so DDA is
   where MBR would pay. Track A item below.
5. **User-facing docs refresh (large; deliberate, not a feature side-effect).**
   Stale post-M3; docstrings are the interim source of truth.

**Parked / blocked (do when unblocked, not ahead of the above):**
- **>2-condition pairwise Δk contrasts** — `fit_linear_deltak` handles exactly two
  qualifying conditions today; the pairwise/Tukey extension is **blocked on a good
  ≥3-condition test dataset** (user to provide). Forward-compatible (per-condition
  k for any N; Δk left NaN when ≠2).
- **Track B integration fidelity** + the `calibration` 4th model + guan/fornasiero
  validation — demand-driven, gateable by the Track D within-protein-θ bench.
- **o18 rewrite** — needs the NB90b frozen coefficient table; `fit --label o18`
  errors until then.
- **GUI framework feasibility (Electron/Tauri/web)** — long-term decision doc, only
  if Qt becomes a concrete graphing ceiling.

#### Locked decisions (2026-06-07)

1. **Labeling time is a required SDRF column** — `characteristics[labeling
   time]`, not parsed from `source name` or the filename. Deliberately a
   *characteristic* (sample-intrinsic), not `factor value[time]`: it is the
   kinetic-curve x-axis, and keying on `characteristics[labeling time]` avoids a
   real collision with a drug-treatment time course (where the *treatment* time
   is the genuine `factor value[time]` while the metabolic-labeling duration is a
   separate property). The current `data/timeseries_lve/samplesheet_lve_sdrf.tsv`
   does *not* carry it (time is only in `LVE_d0` / `…time0`), so adopting Riana's
   SDRF convention means adding this column. Open sub-decision: the value/unit
   format (bare number with a documented unit vs. a unit token — LVE is days).
   Riana reads its own documented subset, so it does not depend on generic SDRF
   tooling treating the studied variable as a factor value. The calibration
   mixing series declares its type with the sibling `characteristics[mixing
   proportion]` instead; fit dispatches on which column is present (kinetic
   models vs. the `calibration` 1:1 recovery model — see Track C).
2. **Header-authoritative identity + a manifest.** Integrate freezes the full
   identity into each output's provenance header (a frozen SDRF snapshot →
   reproducible, SDRF-independent at fit time); a stage-aware
   `riana_manifest.tsv` is the index. Fit groups runs from the manifest rather
   than re-reading the SDRF.
3. **DIA-NN parquet intake is a fast-follow** (its own milestone), but
   `PSMRecord` grows `retention_time` now so the DIA RT-prior path is designed
   in, not bolted on.
4. **Protein comparison ships as a side-by-side view first**; in-app
   cross-sample statistics are deferred. The bridge is the linearized
   simple-model fit (`log(1−θ) = −kt`): once that supports a shared-variance
   two-sample linear model (marginal-mean k per sample, Δk test), the stats
   layer becomes reachable.
5. **Per-sample RIA is an SDRF column** (optional `characteristics[precursor
   enrichment]`; resolution order SDRF → global `--ria` default). It is
   load-bearing at the per-timepoint theta-solve — it builds `_get_final_env`
   in the IsoSpec forward model, which is what recovers θ in the first place —
   so it cannot stay a single global scalar once one mzTab spans multiple
   animals.

#### Track A — I/O & run-identity data model (the spine)

- **M6a — identity model + intake refactor. DONE (2026-06-07, branch
  `m3-rewrite`).** Shipped as `RunIdentity` (`records.py`) + `io/sdrf.py` +
  `io/manifest.py` (schema-versioned, stage-aware) + `core/pipeline.py`
  (`integrate_project` one-file-per-run with identity-in-header + manifest;
  `recombine_for_fit`/`fit_project` grouping curves by `(experiment,
  condition)`), wired as `integrate --sdrf` / `fit --manifest`. Curve grouping
  was user-confirmed: one curve per condition, `characteristics[biological
  replicate]` are independent replicate points (not separate curves). Percolator
  is the demoted single-mzML tier. Variable mods are parsed/exposed but not yet
  threaded into the envelope (that stays M7). A parity test pins the manifest
  path to the legacy path. **Deferred to follow-ups:** DIA-NN intake (M6b),
  bounded file-parallelism (below), GUI rewiring onto `core/pipeline.py`
  (Track E), threading SDRF mods into integration. Original spec below.
- **M6a (original spec).** `io/sdrf.py` + a documented
  Riana-read column subset (`comment[data file]` → mzML join key; `source
  name`; `characteristics[biological replicate]`; `comment[technical
  replicate]`; `comment[fraction identifier]`; `characteristics[labeling time]`
  (required for turnover, the kinetic-curve x-axis) OR `characteristics[mixing
  proportion]` (calibration fixtures) — the independent-variable column declares
  experiment type and fit dispatches on it; `factor value[...]` → condition/group
  for side-by-side + future cross-group stats (key on whichever `factor
  value[...]` column(s) exist — `factor value[condition]` the canonical case,
  so `factor value[disease]`/`[genotype]` work too);
  `comment[proteomics data acquisition method]` → DDA/DIA;
  `characteristics[precursor enrichment]` → RIA (optional)); validation with
  clear errors; everything else quantms fills is ignored. Identity +
  `retention_time` on `PSMRecord`. mzTab becomes the primary path
  (`read_mztab(path, sample_map)` where `sample_map` comes from the SDRF);
  Percolator is demoted to a single-mzML testing/legacy tier (don't rescue
  multi-fraction Percolator). **Output: one `<mzml_stem>_riana.txt` per run**
  (no timepoint-name collisions), full identity + RIA frozen into the header.
  **Stage-aware `riana_manifest.tsv`** (`stage = integrate | fit | protein`) as
  the project index. **Fit-time recombination** — group runs by `(experiment,
  sample, bio_rep)`, merge fractions of the same `(sample, rep, timepoint)` at
  peptide level, order by time — lands in a shared **`core/pipeline.py`** (the
  extraction the GUI "mirror, don't refactor" note flagged), consumed by both
  CLI and GUI. Read `comment[modification parameters]` from the SDRF so variable
  mods need not be set on the CLI; **and read `comment[precursor mass tolerance]`
  as the integration window** (CLI `--mass_tol` overrides; default 10 ppm).
  *(Reverses the original rec4: it had said NOT to inherit the search tolerance,
  using a deliberately wider ~50 ppm "signal-capture" window. That was a
  profile-vs-centroid category error — on **centroid** mzML, ~50 ppm imports
  co-eluting interference into the heavy isotopomer channels. Verified 2026-06-09
  on LVE: tightening 50→10 ppm lifted R²med 0.69→0.87, ≥0.8 41→61%, within-protein
  θ 0.16→0.10, out-of-range θ 27%→8%; all current mzMLs are MS1 centroid. The
  search tolerance IS the right integration window for centroid data. A
  profile-mode mzML triggers an intake warning.)*
- **M6b — DIA-NN parquet intake. SHIPPED (2026-06-11, branch `m3-rewrite`).**
  `io/diann.py` reads the DIA-NN `report.parquet` (validated on **2.5.0** with the
  `diann` parquet output, ≥ 2.2.0 as emitted by quantms-diann); `plan_integration`
  dispatches the reader on `SdrfTable.acquisition` (DIA → parquet, else mzTab), so
  the same `io/sdrf.py` spine and both CLI + GUI surfaces get DIA for free — only
  the parquet reader is new. DIA has **no MS2 anchor**, so the reader emits
  `scan = -1` + carries DIA-NN's apex `RT`, and `core/integration.resolve_rt_anchored_scans`
  maps it to the nearest MS1 scan in the mzML at integrate time (+ an RT-in-bounds
  wrong-mzML check); the scan↔RT scramble guard is skipped for DIA (circular
  there). Riana still extracts the MS1 isotopologues itself — DIA-NN is "just
  another ID + RT source." **Variable-mod peptidoforms are dropped** (the M7
  caveat — they'd integrate at the unmodified m/z; ~1.8% on the cardiac set).
  Needs `pyarrow` (the `[dia]` extra, lazily imported). Validated end-to-end
  (`integrate → fit → rollup`) on the cardiac in-vivo DIA set under
  `data/timeseries_dia` → `runs/lve_dia`: **98.5% m0 coverage, ~1–2 ppm mass
  accuracy** (extraction is on-target), fit **R²med ~0.42 / ≥0.8 ~22%** (lower
  than the DDA LVE ~0.87 / ~61% — DIA's wide-window MS1 imports co-eluting
  interference into m1–m5, and the curve is 3 timepoints at RIA 4.6%), but
  **k_deg median 0.08/day** is biologically sound. *Open follow-ups:* a stronger
  DIA-specific QC than RT-in-bounds; if DIA turnover matters, an MS1-interference-
  robust envelope matcher (overlaps the Track B adaptive-N_ISO matcher).
- **Integrate concurrency + GUI rewiring (Track E, SHIPPED 2026-06-10).**
  `integrate_project` was split into `plan_integration` (SDRF+mzTab → per-run
  `RunTask`s) / `_integrate_results` (yields each frame as it finishes) /
  `finalize_run`. **Crash-resilient (2026-06-11):** the main process writes each
  run's `_riana.txt` + appends its manifest row *as that run completes* (workers
  only integrate), so an interruption keeps finished runs; `integrate --resume`
  skips runs already in the manifest at the current `config_hash`. The GUI
  Integrate tab now drives the *same* plan + per-run unit over its **own** shared
  pool via `asyncio.as_completed` bounded by a *Workers* spinbox (one mzML per
  concurrent run) — cross-file parallelism with **no nested process pools** —
  and gained an **SDRF** field that routes the search-ID file through
  `core/pipeline` (identity-stamped `<stem>_riana.txt` per run + manifest), on
  top of the per-run *Threads*. The GUI Model tab gained a **Manifest** field
  (`fit_project` via `tasks.run_fit_manifest`) and now also writes the M5
  `riana_fit_fractions.txt`. So integrate/fit go through `core/pipeline` on both
  surfaces. **Rollup threading shipped** (2026-06-10): the per-protein refit is
  dispatched over `--thread` / the Protein-tab spinbox with per-protein RNG
  streams (`_group_rng`), so the result is identical regardless of thread count.
  *Still TODO:* the SDRF mass-tolerance resolution the CLI does (the GUI uses the
  spinbox value explicitly).
- **Manifest project chain (SHIPPED 2026-06-10).** `fit --manifest` /
  `rollup --manifest` now root outputs at the **manifest's folder** (the project
  dir; `-o` ignored there with a warning) and write back their `stage="fit"` /
  `stage="rollup"` rows (`record_stage_rows` + `fit_outputs_from_manifest` in
  `core/pipeline.py`), so one `--manifest` drives `integrate → fit → rollup` and
  the manifest indexes every stage. `rollup`'s `FIT_DIR` argument is now optional
  (XOR `--manifest`). The aggregate fit/protein rows carry a coarse
  experiment-level identity (`aggregate_identity`) — indexing/provenance only.
  (Closes the gap where `fit` never updated the manifest despite the M6a plan.)

> **⚠️ quantms input gotcha — mzML filenames MUST NOT be prefixes of one another
> (verified 2026-06-08, quantms/OpenMS ~1.7.0).** ConsensusID/ProteomicsLFQ
> matches mzML/spectra by a **filename-prefix** rule, so when one basename is a
> prefix of another (`…_time1` vs `…_time10`/`…_time15`; `…_time3` vs `…_time30`)
> the mzTab `spectra_ref` **scan is pulled from the sibling file** — the scan↔file
> association silently scrambles for every entangled run. On the LVE series this
> wrecked integration on 8 of 12 runs (their `spectra_ref` scans indexed a
> sibling's mzML; e.g. `time30`'s scans → `time3`'s mzML to 0.13 min). Proof was
> single-variable: **renaming to zero-padded non-prefix names (`time00…time30`),
> same 2023 conversion, restored scan-match to ~90% on all 12.** Also previously
> seen as the reason Sage couldn't be added to ConsensusID. **Mitigations:**
> (1) zero-pad / otherwise de-prefix mzML basenames before a quantms run;
> (2) **anchor mzTab integration on `spectra_ref` scan, NOT the reported
> `retention_time`** — RT carries an alignment-frame offset (~0.1 min on clean
> runs, but the mzTab RT is OpenMS-aligned, not raw) so scan is the correct,
> confound-free key *once filenames are clean*; (3) **Track A intake guard
> (SHIPPED 2026-06-10):** `integrate_run` reconciles, per run, each PSM's
> `spectra_ref scan → mzML RT` against the mzTab `retention_time` and raises
> `DataError` when the **per-run median** offset exceeds `scan_rt_tol_min`
> (default **2.0 min**). The median (not per-PSM) gate is robust to the ~10% of
> PSMs that legitimately mismatch and to the run-dependent ProteomicsLFQ
> alignment offset — measured at ≤~0.9 min on the de-prefixed LVE runs (time00
> median 0.92, time03 0.22, time30 0.36; all 100% within 2 min), vs tens of
> minutes on a scrambled run (~25× separation). Errors loudly with the
> de-prefix/zero-pad remedy in the message; `--no-rt-check` (config
> `check_scan_rt=False`) overrides. No-ops on the Percolator path (no
> `retention_time`). **DIA caveat:** this is scan-anchored, so it does not
> transfer to M6b — DIA-NN reports an inferred per-run apex RT with no MS2-scan
> anchor to reconcile against; the DIA analog (reported RT within mzML bounds +
> run-column matches the paired mzML) is a separate, weaker check for M6b.

- **Re-explore match-between-runs (MBR) for the mzTab/DDA path (spiked 2026-06-12,
  mid-term).** DDA misses peptides stochastically across the time series, so a
  turnover curve can lose timepoints to identification gaps rather than real
  absence. MBR (transfer an ID/RT to runs where the peptide was observed but not
  picked for MS2) would fill those. **First step is measurement, not code:
  compare the per-peptide missingness of the DIA vs DDA sets** (`runs/lve_dia` vs
  `runs/lve`/`lve_atr`) — DIA-NN already does cross-run propagation, so DIA is
  likely complete enough to *not* need MBR; quantify the gap before building
  anything. If warranted, scope MBR at the intake layer (it is an ID-assembly
  concern, hence Track A) — there is a stale `data/mbr_test` fixture from the
  pre-rewrite era to revisit. Pairs with the Track D missingness metric.

#### Track B — integration fidelity (research cluster)

Plan these together — they are complementary integration-fidelity levers, and
the high-isotopomer channels are exactly where peak detection, baseline, and
N_ISO all interact:

- **Phase C v2 — cross-proportion-stable peak picker.** The one untested lever
  from the M3 spike: tight `apex_search_half_width` (~0.15) and/or
  `consensus_apex` for boundaries stable across labelling proportions (the
  boundary-stability finding).
- **Adaptive N_ISO + robust matcher.** Set the per-peptide channel count from
  the IsoSpec plateau envelope instead of a fixed integer, paired with a
  contaminant-robust observed-vs-IsoSpec matcher (Huber / soft-trim /
  per-channel-SNR weighting) so the extra channels don't import co-eluting
  isobars. The open research question is that matcher.
  - **⚠️ HARD CONSTRAINT — normalize the FULL IsoSpec envelope BEFORE truncating
    to N_ISO** (the "H4′" finding, verified 2026-06-08). `solve_fs_d2o` today
    normalizes init/final *separately over the truncated channels* then mixes:
    `(1−f)·norm(init_trunc) + f·norm(final_trunc)`. That equals the physical
    `norm_trunc[(1−f)·init_full + f·final_full]` **only when init and final carry
    the same fraction of mass inside the window.** The labeled envelope spills
    past the window (more, the higher the RIA / Spep), so normalize-then-mix drops
    the θ-dependent denominator and FS↔θ goes **non-linear** — the math won't work.
    Negligible at LVE's RIA 4.6% (in-window mass fraction ΔS ≤ 0.04 ⇒
    current == correct to 3 dp), which is why it is *deferred*; but it becomes
    load-bearing the moment adaptive N_ISO widens/varies the window, or at high
    RIA / high Spep / ¹⁸O. **So build the mix-then-normalize FS solve as part of
    this item** (mix full-envelope *abundances*, THEN truncate + normalize), not
    after. The same trap kills a naive iso0/iso1-ratio interpolation — you must
    mix the abundances, then take the ratio.
- **`--fs` channel-subset SSE.** Wire the reserved flag: solve FS over a chosen
  isotopomer subset to dodge contaminated channels. (Subject to the same
  full-envelope-normalization constraint above.)
- **Dual-mode FS (abundance + mass-defect shift).** D₂O labelling shifts the
  intensity-weighted accurate mass of each isotopomer (the 2H−1H mass defect
  differs from 13C−12C), so the per-channel mass shift over the init (natural)
  envelope is an *orthogonal* FS estimator. Let the model solve FS both ways and
  combine/cross-validate; for low-abundance peptides mass accuracy can beat
  spectral accuracy, so this is a robustness win that pairs with the adaptive
  N_ISO matcher. (The near-term QC half of this idea is in Track E.)

Gated by both the mixing-series benchmarks and the new animal benchmark
(Track D).

#### Track C — fitting / modeling science

- **M5 — persist per-timepoint fraction-new (SHIPPED 2026-06-10).** `riana fit`
  writes `riana_fit_fractions.txt` — long-format one row per `(concat,
  biological_replicate, labeling_time)` with `fs` and prediction-interval bounds
  `fs_lower` / `fs_upper`. The CIs come from a **unified residual bootstrap**
  (fixed t-design, resample residuals, refit) that now feeds *both* the `k_deg`
  CI and the per-timepoint band: `fs_lower`/`fs_upper` are a **prediction
  interval** — `model(t_i; k*) + resampled-residual` — so they reflect each
  peptide's measurement scatter (the quantity the protein rollup inverse-variance
  weights), not just curve uncertainty. On noise-free data the band collapses to
  ~0; on a sparse 3-5 point curve the fixed design is more robust than the old
  pairs bootstrap. `build_fractions_long` + `out.attrs["fractions_long"]` carry
  it through `fit_run` → `fit_project` (tagged with experiment/condition) → CLI;
  the wide `riana_fit_peptides.txt` also gains `fs_lower`/`fs_upper` list-cells
  for the GUI curve view. **This is the substrate the protein rollup consumes.**
- **Fit-model set + the calibration model.** The kinetic models `{simple, guan,
  fornasiero}` are *all wired end-to-end* already (models math lifted unchanged;
  `_MODELS` dispatch → `curve_fit`; CLI `--model`; `FitConfig.model` validation;
  GUI combo + k_p/k_r/r_p spinboxes) — they fit `k_deg` with `k_p/k_r/r_p` held
  *fixed* at config values. **Caveat (verified 2026-06-07):** only `simple` is
  validated — guan/fornasiero have **no test or benchmark** (every gate, incl.
  the −0.5 bias close, ran on `simple`) and their `k_p/k_r/r_p` defaults
  (0.5/0.05/10.0) are placeholders. Add **`calibration` as a 4th model** — a 1:1
  line fit of observed FS vs. known mixing proportion (slope ≈ 1, bias, R² = the
  recovery metric; the productized "option b" recovery mode), dispatched
  identically and reusing the GUI curve view. Its x-axis is `characteristics[mixing
  proportion]`, so **fit dispatches on the experiment-type SDRF column**:
  `characteristics[labeling time]` → kinetic models; `characteristics[mixing
  proportion]` → `calibration`. **Future science (as-needed):** validate
  guan/fornasiero against the in-vivo animal labeling data (hard to model) and
  decide whether the precursor parameters (k_p / k_r / r_p) are *fitted* or
  *supplied* — both are fixed today, so meaningful two-compartment use needs real
  precursor priors. Unlikely near-term.
- **Protein rollup** (new milestone). **Shipped 2026-06-10** — `riana rollup`
  (`core/protein.py`) writes `riana_rollup_proteins.txt` + `riana_rollup_fractions.txt`
  (a `method` tag + one `k_deg` /
  CI / R², the structure counts `n_peptides` / `n_replicates` / `n_timepoints` /
  `n_points`, and a comparison `peptide_median_k`), grouped by `(experiment,
  condition, protein)`; `--parsimony {unique, isoform}`;
  an optional peptide R² admission gate (`--min-r2`, off by default, with a JCI
  slow-turnover admit); threaded refit (`--thread`); and a **GUI Protein tab**
  (`gui/protein_tab.py` + `tasks.run_rollup`) with a per-protein refit **curve**.
  **`--method` (decided 2026-06-10, after evaluating the menu below):**
  `weighted` (default, the biorep-aware inverse-variance per-timepoint collapse)
  or `pooled` (all peptide×timepoint points, pseudoreplication-naive, for
  comparison). The **point estimators were dropped as methods** — median/harmonic
  over peptide k are a trivial `groupby` the user does on the peptide file (and
  harmonic is outlier-sensitive on the low-k tail); the median is still carried
  as a `peptide_median_k` column. **Next milestone:** the **linearized
  `log(1−θ) = −kt` fit + cross-sample Δk test** — *not* weaker (strong basis;
  estimates marginal means + a two-group Δk in one linear model), it's the
  genuinely new statistical capability. (Manifest stage-row wiring is now done —
  see the manifest project chain above.) The menu that was evaluated:
  - **Parsimony is a summarize-time decision, not integrate-time** (locked
    2026-06-10). A shared peptide's envelope blends both proteins' turnover, so
    its isotope signature can't be attributed — the rollup is unique-by-default.
    `integrate --unique` is therefore **removed**: integrate extracts *all*
    peptides (shared included — still valid per-peptide measurements), and
    `riana rollup --parsimony {unique, isoform}` (default `unique`) does the
    attribution. **`unique`** (shipped): drop any peptide whose `protein id` is
    multi-accession; attribute the rest to the bare UniProt accession.
    **`isoform`** (shipped 2026-06-10; adapted from `02_R_parsimony_reference.Rmd`
    by Juber/Lau): like MSFragger gene/protein parsimony but with a specific
    isoform rule — if the protein group spans **multiple genes**, reject
    (cross-gene turnover is unattributable); if it is **multiple UniProt entries
    of one protein (isoforms)**, then for the isoform set this peptide maps to,
    ask whether **any** of those isoforms has its own *unique* peptide (mappable
    to it alone): **no** → assign this shared-among-isoforms peptide to the
    **canonical** entry; **yes** → exclude it. Rationale: physically-present
    isoforms (evidenced by isoform-unique peptides) make the shared signal
    genuinely ambiguous, while isoforms that don't actually exist / are
    low-abundance (no unique evidence) shouldn't decimate "unique" peptide counts
    just because an isoform-bearing database makes most peptides nominally shared
    — assumes most isoforms are low-abundance or absent. It is a **dataset-wide
    pass** over the peptide↔protein map (`_resolve_parsimony`), built once from
    the peptide table and applied to both the peptide and fraction frames.
    **Deviations from the R reference** (scrutinized, not copied): Riana joins
    accessions with `,` not `;`; an all-isoform group with no canonical present
    falls back to the shared **base** accession (the R left `collapsed_uniprot`
    NA); the isoform suffix is UniProt `-N` *and* JCAST `-J1`/`-J2` (`-J?\d+`,
    requiring digits so a bare trailing `-` is not a suffix); and the θ collapse
    weights by the **M5 inverse variance**, not the R's `log2(Int)`. An optional
    `--min-r2` peptide gate (off by default, with the R's JCI slow-turnover
    admit) is available as a complementary A/B to the inverse-variance weighting.
  - **Point estimators over fitted peptide k:** median (robust default);
    harmonic mean (rate-correct — mean of half-lives ↔ harmonic mean of k — but
    outlier-sensitive on the low-k tail).
  - **Protein-level refit:** pool the qualifying peptides' `(t, θ)` and fit one
    k. Two flavours — all peptide points (pseudoreplication risk: peptides are
    not independent replicates) vs. the **per-timepoint weighted-average
    collapse** (the pseudoreplication-safe default; open sub-decision = weight by
    inverse variance from the peptide CI vs. by intensity). **Biorep-aware:** the
    collapse is *within* `(protein, labeling time, biological replicate)` — only
    peptides of the same protein in the *same animal* are pseudoreplicates;
    *different* `characteristics[biological replicate]` are genuine replicates and
    stay as independent points, giving the refit honest degrees of freedom.
  - **`linear simple` — a new rollup model (SHIPPED 2026-06-12).** A *model choice*
    alongside `simple`/`guan`/`fornasiero` and **mutually exclusive with them by
    nature**: those are nonlinear ODE fits via scipy `curve_fit`; `linear simple`
    is OLS in φ-space. **Implemented:** `core/linear_model.py` (`fit_linear_deltak`
    + `to_phi` / `truncate_plateau`) driven from `core/protein._rollup_linear`
    (reuses the weighted/pooled collapse via the extracted `_collapse_group` /
    `_collapse_long`); CLI `--model "linear simple"` + `--phi-limit` (default −4) +
    `--reference-condition`; output schema `PROTEIN_LINEAR_COLUMNS` (per-condition
    k + `delta_k`/`delta_k_se`/`delta_k_p`/`delta_k_p_adj`). **CI provenance:**
    analytic (statsmodels `conf_int` / `t_test` on the OLS coefficient
    covariance), **not** the residual bootstrap the nonlinear models use — the
    linearized model has a closed-form covariance, so `n_boot`/`random_state` are
    ignored on this path. **Validated on
    `runs/lve_atr`:** 780 proteins in both chambers, atrium 1.24× faster, 71%
    faster in atrium, **410 significant at BH p_adj<0.05** (347 of them
    atrium-faster) — consistent with the descriptive paired Δk. Tests:
    `tests/test_linear_model.py` (6) + `test_rollup_linear_simple_delta_k`.
    **GUI φ-space plotting — DONE 2026-06-13:** "linear simple" is selectable in
    the Protein tab Model combo (φ-limit + reference-condition knobs shown only for
    it); selecting a protein row overlays both conditions in φ-space via
    `CurveView.plot_linear` (clearance points with truncated points hollow, the
    through-origin k line per condition, the φ-limit threshold line, and the
    Δk + p_adj in the title). **Parked — >2-condition pairwise contrasts:** the Δk
    needs exactly two qualifying conditions today; the >2-condition pairwise/Tukey
    case is deferred until there is a good ≥3-condition test dataset (user,
    2026-06-13). `fit_linear_deltak` already emits per-condition k for any N and
    leaves `delta_k`/p NaN when ≠2, so the path is forward-compatible. The design
    below records the rationale.
    - *Transform.* φ = `log(1 − θ)`; clamp θ to ~[0.01, 0.99] and drop non-finite
      (the R reference does this). φ = 0 at t = 0, so the line is **through the
      origin**.
    - *Model.* Per protein, the no-intercept interaction OLS
      `φ ~ 0 + day + day:condition` (R: `lm(protein_clearance ~ 0 + day + day:treatment)`):
      the per-condition slope is `−k`, and the `day:condition` term is the **Δslope
      = −Δk**. Two conditions → marginal-mean k per condition (with CI) + a pairwise
      Δk contrast p-value, then **Benjamini-Hochberg** across all proteins, filter
      `p_adj < 0.05`. Reference (do **not** follow — study-specific compartment /
      strain / JCAST aspects no longer apply): `data/notebook/03_R_linearmodel_reference.Rmd`.
    - *⚠️ Plateau truncation (required, the key deviation from the R reference;
      user 2026-06-12).* Once φ saturates (θ near its measurement ceiling, φ ≈ −4),
      later timepoints are noise around the floor, not slope — including them
      **breaks linearity and drags the through-origin slope flat**. So per curve,
      **truncate at/after the first timepoint that reaches saturation**. **Knobs
      (user wants these exposed):** a configurable **φ limit** (the truncation
      threshold) — `−3…−4` is sensible (≈ θ 0.95–0.99; beyond that is measurement
      noise) — plus the θ clamp. The nonlinear models **need** the plateau (their
      asymptote parameter fits it), so truncation is **linear-only**. This is
      **inseparable from `linear simple`**; build them together, not as a
      follow-on.
    - *EMMEANS — use statsmodels (user decision 2026-06-12).* numpy/scipy have
      **no** emmeans/emtrends/marginal-means machinery (scipy.stats stops at
      `linregress` / `ttest` / `f_oneway`). Although this simple fixed-effects
      design *could* be hand-rolled (β=`lstsq`, cov=`σ²(XᵀX)⁻¹`, contrast
      `cᵀβ ± t·√(cᵀ cov c)`, BH ~5 lines), the user prefers **adding statsmodels
      as a dependency** and building on it (`OLS.t_test`/`f_test` for the contrasts
      + `stats.multitest.multipletests` for BH): we will want it later anyway (the
      R notes defer mixed / `limma`-style models — "limma later"), and a
      well-established package is easier to maintain than bespoke contrast math.
      → add `statsmodels` to deps when `linear simple` lands.
    - *Why linear is not weaker* (user, prior): strong mathematical basis; it
      estimates marginal-mean k **and** a two-group Δk in one model — the genuinely
      new statistical capability. Guan/Fornasiero add a precursor-rise term and are
      not linearizable, so their cross-sample equivalent needs nonlinear
      mixed-effects (deferred).
  Reads both the fit and integrate outputs via the stage-aware manifest; groups
  samples by `factor value[...]` (condition/group) and supports opening multiple
  groups side-by-side for visual comparison — the substrate for the deferred
  cross-group stats.
- **o18 rewrite.** Same SSE-vs-IsoSpec → fit machinery as D₂O, but FS is
  computed with o18 isotope mass/proportion **and a different Spep model**: the
  labelling sites are only the oxygen-bearing residues, so the coefficient
  schema is a regression on `length, #D, #E, #N, #Q, #S` (possibly #T/#Y) — *not*
  a per-AA-over-20 dict. This needs a distinct `--label o18` dispatch (different
  coefficient format) and a frozen o18 coefficient table built via
  `data/notebook/90b_O18_LearnAALabelingSites_IsoSpec_AC16.ipynb` (the analog of
  87a for D₂O). Until then `riana fit --label o18` errors; o18 *integration* is
  unaffected.
- **M7 — PTM-aware envelope** (demand-driven). Thread the parsed bracketed mod
  masses through `get_peptide_distribution` / `solve_fs_d2o` so the IsoSpec
  envelope shape is correct for modified peptides (today the FS solve uses the
  unmodified-backbone envelope; the per-cell-line coefficient tables are fit on
  the backbone and don't change). `_fit_one_concat` already separates
  `seq_with_mods` from `seq`; M7 threads the mods alongside instead of dropping
  them. Prioritize only for PTM-focused datasets (e.g. labelled
  phosphoproteomics). Three concrete pieces (user context, 2026-06-11):
  1. **Atom accounting in the forward model.** A mod's extra atoms must be added
     to the IsoSpec `formula` *even though the mod is not D₂O-labeled* (we don't
     know its enrichment, and the forward model mixes initial + fully-labeled
     envelopes — the extra C/H/O/N/S still shape the envelope and shift the
     channel masses). The envelope formula comes from `count_atoms` →
     `IsoParamsFromDict({"C":…,"H":…,…})` in `algorithms/isotope_dist.py`; M7 adds
     each parsed mod's composition there. **Starter set (locked 2026-06-13):
     phospho (S/T/Y, `UNIMOD:21` = `[0,1,3,0,0,1]`) + protein N-term Acetyl
     (`UNIMOD:1` = `[2,2,1,0,0,0]`)** — trimmed to the two highest-value,
     reliably-identifiable mods (the next tier + its binding constraint are in the
     roadmap below). **Composition source = a small curated `mod_atoms` table**
     keyed by UNIMOD id (`{id: [C,H,O,N,S,P]}`), **not** the UniMod XML ontology
     (too heavy for ~5 entries); today `constants.mod_atoms` holds only `IAA`
     (which becomes `UNIMOD:4`).
  2. **Generalize Carbamidomethyl(C) to a normal UniMod.** Today it is special-
     cased: `count_atoms(iaa=True)` adds `mod_atoms['IAA']=[2,3,1,1,0]` per
     cysteine, so C is effectively `[5,8,2,2,1]` (actual C `[3,5,1,1,1]` + IAA) in
     the forward model. M7 should route CAM through the same per-mod composition
     machinery as any other UniMod (`UNIMOD:4`) rather than the hardcoded `iaa`
     flag — one code path for all fixed/variable mods.
  3. **Proteoform-aware rollup.** A phosphopeptide must **not** collapse into its
     unmodified protein (that hides the PTM's effect on turnover). At rollup,
     count the PTM-site residues and append them to the accession as a distinct
     proteoform key, e.g. `P12345_pS235`, so each modified form rolls up as its
     own unit. User has R scripts to reference for the site-naming convention.
     This extends `core/protein.py` parsimony, which today keys on bare
     accession.

  Non-obvious structural implications (noticed while grounding in the code,
  2026-06-11):
  - **Two separate fixes, don't conflate.** (a) *Integrate-side target m/z* —
    today `mass_calc.calculate_ion_mz` recomputes from the **bare** sequence, so a
    modified PSM is extracted at the unmodified m/z (the reason we currently
    *drop* it). M7 must add the mod mass to the extraction target. (b) *Fit-side
    envelope shape* — the mod atoms must enter the IsoSpec `formula`. The drop
    sidesteps both; M7 needs both.
  - **Phosphorus is not in the `[C,H,O,N,S]` atom vector.** `count_atoms` returns
    5 elements and `constants.iso_abundances` is 5-long; phospho (HPO₃) adds a
    **P** (and 3 O, 1 H). P is monoisotopic so it does not broaden the envelope,
    but supporting it means extending the atom vector + abundance list, not just
    the mod table. The 3 O *do* affect envelope shape.
  - **Mod hydrogens default to non-labelable.** The forward solver subtracts
    `num_labeling_sites` (from `label_deuterium_de`/`label_oxygens`) from H before
    handing to IsoSpec. Conservative default: mod-contributed atoms are static
    (not added to the labelable count) since the mod's D₂O enrichment is unknown
    — matches the user's "account for the atoms even if not necessarily enriched."
  - **Site localization gates the proteoform key — and both primary formats give
    the protein-coordinate site directly (no FASTA needed).** Checked 2026-06-11:
    - *mzTab* carries `start`/`end` (the peptide's protein-coordinate start) plus
      the peptide-relative mod position in the `modifications` column
      (`pos-UNIMOD:id`, 1-based, 0 = N-term), so the protein site =
      `start + pos − 1` (e.g. MYL4 `start=93` + Ox at pos 10 → M102; TITIN
      `start=34473` + Phospho at pos 4 → S34476). Shared peptides carry a
      comma-separated `start` per accession.
    - *DIA-NN parquet* carries `Protein.Sites` **pre-formatted** as `[acc:res+pos]`
      (e.g. `[Q9Z1P6:M87]`, multi-site `[P55264:C139,C142]`) plus
      `PTM.Site.Confidence` (1.0 = confident) and `Site.Occupancy.Probabilities`.
      It also lists CAM cysteine sites, so the proteoform-key builder must filter
      to the variable mods of interest.
    - *Percolator* (demoted/testing tier) has no protein-coordinate site — **drop
      PTM support there** (user decision 2026-06-11) rather than require a FASTA.
    Only confidently-localized sites become distinct `_pS###` keys; decide a
    bucket/drop policy for ambiguous ones (gate on DIA-NN `PTM.Site.Confidence`).

  **Locked planning decisions (2026-06-13):**
  - **Phospho is in v1** → extend the atom vector `[C,H,O,N,S]` → `[C,H,O,N,S,P]`
    across `count_atoms` / `_calc_atom_mass` / `constants` (`aa_atoms`,
    `iso_abundances`, the mass vector) / the IsoSpec `IsoParamsFromDict` formula.
    P is monoisotopic (no envelope broadening) but phospho's 3 O *do* shape it.
  - **Curated UNIMOD-id-keyed `mod_atoms` table** (not the UniMod ontology); CAM
    routed as `UNIMOD:4`, retiring the `iaa` special-case flag → one code path for
    fixed and variable mods. Unmodified-peptide masses/envelopes must stay
    byte-identical; the frozen M2 oracle (`_helpers/forward_model.py`) stays
    5-element as the independence check and is **not** touched — instead assert
    production still matches it on unmodified peptides.
  - **Proteoform rollup = "differential-turnover vs constitutive/artifactual"
    split** (the bio-vs-artifactual decision, refined). Mods with regulated,
    site-specific turnover — **phospho** now, K-acetyl / GG later — get distinct
    `_pS###` proteoform keys. Constitutive or artifactual mods — **protein N-term
    Acetyl** (co-translational, ~constitutive per N-terminus), **Met-Ox**,
    **deamidation**, **CAM** — **fold into the bare accession** (same turnover
    unit), but still get atom accounting so their peptidoforms integrate at the
    right m/z + envelope instead of being dropped. So N-term Ac is in v1 for
    *envelope fidelity / peptide retention*, phospho for *proteoform turnover*.
    A `BIOLOGICAL_MODS` (gets-its-own-key) set in `constants` drives the split.
  - **Retire BOTH `-X/--ignored_mods` and `-F/--forced_mods`** — SILAC-era
    dual-channel machinery (SILAC fitting already removed in M4). Drop the CLI
    options (`cli.py`), the `IntegrationConfig` fields (`config.py`), and the
    `mod{offset}` channel machinery in `integration.py` (hardcode the single
    `mod0` path — the current default), plus the `ignored_mods` plumbing through
    `io/mztab`, `io/diann`, `io/percolator`. Percolator keeps working for
    unmodified peptides only (no protein-coordinate site → no PTM support there).
  - **Build order:** **A1 — DONE 2026-06-13** (atom-vector `[C,H,O,N,S]` →
    `[C,H,O,N,S,P]`; `mod_atoms` now a curated UNIMOD-id table; CAM via
    `mod_atoms[4]`; `count_atoms(mods=…)` + `get_peptide_distribution(mods=…)`
    thread variable-mod composition; verified byte-identical on unmodified peptides
    — mass Δ 0, envelope Δ 3e-18 — and against the frozen 5-element oracle; 178
    tests pass) → **A3 — DONE 2026-06-13** (retired `-X/--ignored_mods` +
    `-F/--forced_mods` from CLI/`IntegrationConfig`; dropped the `mod{offset}`
    channel machinery in `integration.py` — `integrate_run`/`_extract_per_psm`
    now emit `iso{N}` directly at the PSM's own m/z; removed `ignored_mods` from
    `calculate_ion_mz` + the mztab/diann/percolator/pipeline/gui plumbing;
    byte-identical — the 415s ac16 baseline gate + 177 others pass)
    (independent cleanup; removes dead surface before threading) → **A2** thread
    mods IO→fit (load-bearing; needs a **peptidoform-distinct `concat`** so
    phospho ≠ the unmodified form of the same sequence; mod H stays out of
    `num_labeling_sites`) → **B — DONE 2026-06-13** (proteoform rollup keys:
    `PSMRecord.mod_sites` = biological-mod site in protein coords — `pS34476` from
    mzTab `start`+`pos`, BIOLOGICAL_MODS={21} only — threaded IO→integrate(`mod
    sites` col)→fit→`riana.core.protein._resolve_parsimony`, which appends it so a
    phosphopeptidoform rolls up as `A2ASS6_pS34476`; N-term Ac / unmodified → bare
    accession; two peptides covering one site share the key. DIA-NN deferred
    (no phospho fixture yet). Validated on `data/timeseries_lve_atr` — 259 phospho
    PSMs keyed; full suite + ac16 byte gate pass).

  **Roadmap — which mods come next, and the binding constraint.** The hard gate is
  **identifiability in a search over *un*enriched data**: no PTM-enrichment
  D₂O-labeling dataset exists yet, so we can only measure turnover of PTM forms
  detectable in the ordinary global-proteome runs. That, not envelope difficulty,
  is what sequences the list.
  - *Tier 0 (v1):* **phospho-STY**, **protein N-term Acetyl** — high value, reliably
    found unenriched (N-term Ac is near-universal/high-stoichiometry; abundant
    phosphosites do show up without enrichment, just fewer).
  - *Tier 1 — biological, own key, acceptable unenriched yield:* **Lysine acetylation
    (K-ac, `UNIMOD:1` — same `[2,2,1,0,0,0]` composition as N-term Ac, side-chain
    site)**; low stoichiometry unenriched, but abundant metabolic enzymes / histones
    give real sites → gets its own `_acKxxx` key. **K/R methylation** (mono/di/tri,
    `UNIMOD:34/36/37` = `[1,2,0,0,0,0]` / `[2,4,…]` / `[3,6,…]`). Cheap once the
    machinery exists — composition only, no new atom-vector work.
  - *Tier 1b — artifactual / chemical, fold-into-bare, envelope fidelity:* **Met-Ox
    (`UNIMOD:35` = `[0,0,1,0,0,0]`)** and **deamidation N/Q (`UNIMOD:7` =
    `[0,-1,1,-1,0,0]`)** — ubiquitous, trivially identified, not turnover units, but
    accounting their atoms **recovers the abundant modified peptidoforms we
    currently drop**, raising bare-protein peptide counts. **Chemical mods need a
    peptide-level merge, not just envelope accounting** — see the boxed design
    below; deamidation is the harder of the two (isobaric-overlap regime).

  **Chemical-mod handling — integrate-separate, fit-merge (design, user 2026-06-13).**
  A purely chemical mod (Met-Ox, and the chemical part of deamidation) happens
  *post-synthesis*, so it does **not** reset the D₂O clock: the oxidized and
  unoxidized forms of a peptide share the *same* FS-vs-time signature. They should
  therefore fold into the same proteoform **and the same peptidoform** (one
  turnover curve), not be fit as two separate underpowered curves.
  - *Mechanism:* **integrate the forms separately** (each at its own clean m/z +
    envelope — Met-Ox is +15.995, well resolved; A2 already makes each `[UNIMOD:N]`
    form a distinct `concat`), then **merge at the fit level** via a
    *chemical-mod-stripped fit key*: strip chemical mods (Ox) from the grouping key
    so their per-timepoint FS points land on one curve; keep biological mods
    (phospho) distinct. Solves the use_range/ms2 "which m/z to integrate?" problem
    (you never integrate a blended m/z). The fit key is the layer *between* the
    integrate `concat` and the Stage-B protein/proteoform key.
  - *Empirical motivation (LVE_ATR PTM search, depth≥4):* of 545 Met peptidoforms
    seen oxidized, **451** appear as *both* Ox and non-Ox (same cond/charge) — i.e.
    currently double-fit as two curves; merging gives **4** genuine depth-rescues
    (neither form alone ≥4) and **51** power-gains (already fittable, +points). The
    raw "fittable-series count" *drops* under merging (deduplication), so it is the
    wrong metric — the gain is consolidation/power, not count.
  - *Deamidation is the hard case — conditionally tractable, mass-gated.* +0.98401 Da
    (N→D / Q→E) sits only **0.0193 Da below the C13 M+1** (`1.00335`), i.e. a
    separation of **≈ 19340 / M_neutral ppm** (~19 ppm at 1000 Da → ~6 ppm at
    3000 Da; charge-independent in ppm, but higher z lands the peak at lower m/z
    where the Orbitrap resolves better). Because deamidation is **non-stoichiometric
    and time-correlated** (old/unlabeled protein is most deamidated), an unresolved
    deam-M0 bleeds into the non-deam M+1 *in lockstep with the labeling state* — a
    smooth, systematic bias that can warp `k` and pass the R²>0.95 curation silently.
    Two gates decide separability: (1) MS1 resolution actually centroided them apart
    (R ≈ M/0.0193 — ~52k at 1000 Da, ~104k at 2000 Da, lost by ~3000 Da or at fast-DIA
    15–30k MS1), and (2) a tight extraction window (±3–5 ppm, not the default ±10).
    **Where both hold** (small/mid peptides, high-res MS1, tight window) → it reduces
    to the Met-Ox case (integrate-separate + fit-merge). **Where they don't** (large
    peptides, low-res MS1, wide window) → merged centroid, irrecoverable → **flag /
    exclude the peptide (both forms)**. The *rigorous* alternative for the
    unresolvable regime — **jointly modelling both forms' FS/D₂O IsoSpec envelopes
    AND the (unknown, time-varying) deamidation proportion** — is a **standalone
    side project**, deferred (user 2026-06-13). So the near-term policy is a
    per-peptide **mass/resolution gate** (`19340/M_neutral` vs the window ppm), not
    blanket merge or blanket exclude. The bio-vs-chemical distinction is **moot** — both are
    the identical +0.984 mass shift, unactionable; FS-merge is valid for either, so
    the only question is separability. When deamidation lands, first measure the MS1
    resolution + deam-peptide mass distribution (e.g. on LVE) to size the separable
    fraction. This is a **fit-stage** concern; none of it is in v1 (Ox/deam are
    dropped), so it does not affect Stage B.
  - *Tier 2 — highest biological value, BLOCKED on data:* **Ubiquitin/ISG15
    GG-remnant (`UNIMOD:121`, GlyGly = `[4,6,2,2,0,0]`)** — it literally *is* the
    degradation tag, so it is the most turnover-relevant PTM imaginable, **but**
    diGly-remnant peptides are essentially undetectable without anti-K-ε-GG
    enrichment, and no enriched D₂O dataset exists. Deferred until such data lands
    (same blocker for the acyl class — succinyl / malonyl / crotonyl — that needs
    enrichment). Worth naming now precisely because the science case is strongest:
    the moment an enriched D₂O-diGly dataset appears, GG jumps the queue.

#### Track D — validation infrastructure (unblocked by M6a)

- **Animal within-protein-θ benchmark — BUILT (2026-06-10).** The 4th dataset
  (`data/timeseries_lve`, mouse in-vivo D₂O, mzTab + SDRF) has no fractional-pool
  ground truth, so it is scored by *within-protein θ/FS variance per timepoint*
  (Hammond 2022): a better integrator minimizes the spread among peptides of the
  same protein. Shipped as `build_lve_bench_set.py` (frozen, integrator-
  independent set: proteotypic peptides of ≥3-peptide proteins at q<0.01,
  abundance tertiles) + `bench_within_protein_theta.py` (per-timepoint θ via the
  production `solve_fs_d2o`; overall / per-stratum / curated-vs-uncurated;
  `compare_methods` so a Track B variant slots in on the identical set). Baseline
  at the adopted defaults (apex@0.15 + best-q anchor + 10 ppm): within-protein-θ
  robust-SD median **0.10** (curated 0.044), correct abundance ordering. Building
  it surfaced and fixed **two integration bugs**: (1) the quantms **filename-prefix
  scan-scramble** (Track A gotcha above) and (2) **mass_tol 50→10** (centroid
  search-tolerance, not profile width — LVE R²med 0.69→0.87). Now the orthogonal
  regression gate for Track B.
- **Curation-gate revisit.** Treat both knobs as tunable — the R²>0.95 threshold
  *and* the "observed at every proportion" coverage rule (it should tolerate a
  known-bad fraction rather than discard peptides wholesale). Report the R²
  *distribution* shift before/after a fidelity change, not just curated counts.
- **DIA test dataset (TODO, enables M6b).** A small quantms-diann run lands under
  `data/timeseries_dia` (DIA-NN `report.parquet` ≥ 2.2.0 + the DIA-variant SDRF),
  the fixture the `io/diann.py` adapter and DDA/DIA auto-detect are built and
  tested against.

#### Track E — GUI/UX & responsiveness (surfaces the engines above)

- **Sortable tables** (`QSortFilterProxyModel` over `DataFrameTableModel`) — an
  easy win.
- **Save/export graphs** — replaces the removed `--plotcurves`; pyqtgraph export
  for interactive views, matplotlib for static; all tabs.
- **Δmass-over-time QC display** — the near-term half of dual-mode FS: surface
  the per-isotopomer accurate-mass shift over the init envelope across the time
  series, a high-level QC for advanced users.
- **Surface smoothing in the GUI + make the displayed trace faithful.** The CLI
  `-S/--smoothing` (Savitzky–Golay, `polyorder=2`) is missing from the Integrate
  form, so GUI runs can't enable it. Worse, the chromatogram re-extracts the
  *raw* XIC from the mzML at row-click, so the plotted trace does not reflect the
  S-G smoothing integration actually applied. Fix both: add the form control, and
  either re-apply the run's `IntegrationConfig.smoothing` at click or — better —
  read the persisted smoothed trace from the `--save-traces` sidecar (which also
  removes the ~10 s re-read).
- **Fit progress + real parallelism** (engine is Track C / `core/pipeline.py`,
  surfaced here). Today the GUI submits one opaque `run_fit` to a single worker
  with an indeterminate "busy" bar, and the inner peptide parallelism is a
  GIL-bound `ThreadPoolExecutor`. Make `fit_run` chunk-aware (peptide chunks),
  dispatch chunks across the **ProcessPool**, emit a progress callback the bar
  maps to. The CLI gets the same speedup. **Reproducibility caveat:** chunking
  across processes requires a deterministic per-peptide bootstrap seed (e.g. from
  `concat`), or output becomes scheduling-dependent — add a parallel == serial
  test (rec below).
- **Modernize look** — low priority, ready-made only (a QSS theme or
  `qt-material` / `qtmodern`).
- **`-W/--workers` in the GUI fit + rollup forms — DONE 2026-06-13.** Integrate
  already had it (the "Workers (files)" spin + semaphore dispatch over the shared
  pool). Added a **Workers** spin to the Model (fit) and Protein (rollup) tabs.
  **Nesting guard (verified the hard way — a `ProcessPoolExecutor` inside a pool
  worker raises `BrokenProcessPool` here):** when `workers > 1` the tab dispatches
  the run on a **main-process thread** (`run_in_executor(None, …)`) so `fit_run` /
  `rollup_proteins` create their pool from the *main* process, not nested inside a
  shared-pool worker; `workers == 1` keeps the shared pool (one CPU task off the
  main process). Per-peptide/-protein deterministic seeds keep the result
  worker-count-independent.
- **Expose the hidden integrate knobs in CLI + GUI as clearly-marked *advanced*
  options (spiked 2026-06-12, short-term).** Several `IntegrationConfig` knobs the
  CLI already takes (`--peak-rt`, `--apex-selection`, `--integration-half-width`
  vs `--extraction-half-width`, `--baseline`, `--smoothing`, `--mass-difference`)
  are power-user dials that are easy to mis-set; the GUI exposes only a subset.
  Audit the full knob set, surface them all (CLI help + a collapsible "Advanced"
  group in the GUI forms), and **label them advanced** so the default path stays
  simple. Pairs with the smoothing-in-GUI item above.
- **GUI-framework feasibility report — is Qt a ceiling? (spiked 2026-06-12,
  long-term.)** If pyqtgraph/Qt limits graphing or interactivity as the protein /
  Δk / φ-space views grow, write a feasibility note comparing alternatives
  (Electron + JS charting, Tauri, a web/server split, Dash/Streamlit for the
  analysis surfaces) against the cost of leaving the native PySide6 app. Decision
  doc only — no migration implied; revisit only if a concrete Qt limit blocks a
  needed view.

#### Persist the integration signal (optional sidecar)

The GUI re-extracts a peptide's XIC from the mzML on every row click (~10 s). Add
an optional `--save-traces` → `<stem>_riana_traces.parquet` (keyed `(concat,
isotopomer) → (rt[], intensity[])`), never in the `.txt`; default off for batch
(size), on for GUI/diagnostic workflows. This is purely a GUI-speed/diagnostics
feature — the protein refit consumes per-(peptide, t) θ (the M5 intermediate),
not the raw XIC, so the two are decoupled.

#### Suggested sequence

1. **Pre-1.0.0 chores** — repo hygiene + cut a real `1.0.0` tag (below). A clean
   line in the sand before piling on features.
2. **M6a** — the spine; unblocks Track D and the protein layer; pulls in
   `core/pipeline.py` and integrate concurrency.
3. **Track D** — stand up the orthogonal animal gate while M6a is fresh.
4. **M5** — small, but gates the protein layer.
5. **Track C protein rollup + Track E polish** — the protein tab needs M5 + the
   M6a manifest; GUI polish anytime.
6. **Track B** — the fidelity research, now with two gates ready.
7. **o18 / M7 / M6b / dual-mode FS** — slot in by demand and bandwidth.

#### Cross-cutting / pre-1.0.0 chores

- **Cut a real `1.0.0`** (currently `1.0.0.dev1`, last tag `v0.9.0`) + a Zenodo
  *code* DOI. M1–M4 is a complete rewrite and deserves the tag.
- **Repo hygiene** (also §4.5): remove the stray root outputs
  (`riana_fit_peptides.txt` etc.), `riana_website/` (mode 0700), the committed
  `docs/` Quarto HTML, and move 1 GB+ of personal outputs out of `data/`.
- **User-facing docs refresh — post-M3, still stale (large TODO).** The narrative
  docs were not redone after the M3/M4 rewrite, so they lag the current CLI/GUI,
  the SDRF/manifest model, the rollup, and now `linear simple`. Until that
  happens, the **docstrings are the source of truth** and should carry the
  non-obvious contracts explicitly — e.g. `linear simple` CIs are analytic, not
  bootstrap (now noted in `core/linear_model.py`). Track a full docs pass as its
  own chunk of work, not a side-effect of feature PRs.
- **Snakemake — retire the bundled `workflow/Snakefile`** (supersedes §4.4). With
  Percolator demoted and quantms / DIA-NN handling search + ID end-to-end
  upstream, Riana's own pipeline collapses to a linear `integrate → fit →
  protein` chain whose glue is the manifest (each step reads the prior's manifest
  rows). Snakemake's DAG/partial-rerun value largely evaporates; ship the three
  CLI subcommands (optionally a thin `riana run` convenience wrapper) and stay
  orchestration-agnostic so Riana composes into whatever workflow already runs
  quantms.
- **mypy --strict** rollout (§4.3); `py.typed` already ships.
- **Test-suite runtime (~931s; profiled 2026-06-13).** ~845s lives in ~7 real-mzML
  tests, all single-fraction `integrate_run` calls (`test_ac16_time0_matches_
  committed_baseline` alone is **415s**; six `sample1`/cli integrate tests are
  51–114s). **`-W` does not help these** — `max_workers` is an `integrate_project`
  cross-*run* lever, but these integrate one fraction each and `integrate_run`'s
  per-PSM loop has no internal pool, so `n_parallel = min(workers, 1) = 1`. The
  effective levers (deferred — recorded, not done): (1) a **`slow` pytest marker**
  on the heavy real-data regression gates so local dev runs `pytest -m "not slow"`
  (the 415s ac16 test is already `skipif`-guarded on gitignored heavy inputs, so it
  only runs where those exist); (2) **`pytest-xdist` (`-n`)** for GIL-free
  cross-test process parallelism (~2× — the 415s ac16 test becomes the wall-clock
  pole). Single-fraction `integrate_run` is intentionally left serial (cross-run
  `-W` covers the multi-file production case).

#### Additional recommendations (added this round)

1. **Benchmark CI smoke tier.** The `bench_*` suite runs manually; add a
   subsampled fast tier to CI so a science regression (e.g. a fidelity change that
   worsens within-protein θ) fails the build.
2. **Reproducibility as an explicit invariant.** Add a parallel-fit == serial-fit
   test (forces the deterministic per-peptide seeding); protects the
   header-authoritative frozen-output story.
3. **Schema-version the manifest/header** so the breaking output-granularity
   change (timepoint-named → mzml-named + manifest) and future format changes are
   detectable downstream.
4. **SDRF as a partial config source** (scoped): read mods from the SDRF; keep
   integration mass tolerance a separate, wider Riana parameter (not the search
   tolerance).
5. **Sequence the breaking-change migration with the Snakefile retirement** so
   there is never a broken intermediate workflow state.


## 4. Cross-cutting recommendations

These apply during and after the rewrite:

1. **Reproducibility.** Stamp output files with git SHA, riana version,
   and a hash of input CLI args. Goes in the first line of `*_riana.txt`
   as a comment.
2. **Configuration.** Frozen `dataclass` `IntegrationConfig` consumed by
   both CLI and GUI. Single source of truth, type-checked.
3. **Type safety.** `from __future__ import annotations` everywhere; add
   `py.typed` marker; run `mypy --strict` on `riana/algorithms/` first
   (smallest blast radius).
4. **Snakemake — superseded by the post-M4 round: retire the bundled
   `workflow/Snakefile`.** With Percolator demoted and quantms / DIA-NN owning
   search + ID end-to-end, Riana's own pipeline is a linear `integrate → fit →
   protein` chain glued by the manifest; ship the three CLI subcommands
   (optionally a thin `riana run` wrapper) and stay orchestration-agnostic. See
   §3 "Post-M4 roadmap" → cross-cutting chores.
5. **Repo hygiene.** Move 1+ GB of personal experiment outputs out of
   `data/` (Zenodo for the M2 dataset bits; gitignore the rest). Remove
   `riana_website/` (mode 0700 directory; either commit cleanly or
   delete). Move `docs/` (generated Quarto HTML) to a `gh-pages` branch
   or rebuild via CI rather than committing to `master`.

## 5. What this plan deliberately does not include

- No "AI/ML peak detection," no cloud, no plugin architecture, no
  database backend, no web UI. Earlier evaluations suggested all of
  these; none are right for a single-author academic tool with this user
  base.
- No new kinetic models beyond the existing three. Focus on getting the
  existing three to produce reliable numbers with honest uncertainty.
- No multi-omics integration, no Kubernetes, no "plugin ecosystem with
  >10 community plugins." Same reason.

## 6. Known limitations (honest list)

- Memory: fixed in 1.0 via streaming/indexed reads (`io/mzml.py`,
  `IndexedMzML` — one fraction in memory). The 0.9.0 list-of-arrays
  loader (`spectra.py`) is removed.
- Peak fidelity: see §2c. The calibration dataset quantified the impact;
  the M3 spike adopted an apex-centred narrow window as the default.
- AA / SILAC fitting: dropped from `riana fit` in M4 (the buggy
  `--label 4` `a_0` path is gone). The SILAC dual-channel extraction knobs
  (`-F/--forced_mods`, `-X/--ignored_mods`) were retired in M7 Stage A3 — mods
  are now handled uniformly via UniMod; integrate extracts one channel per
  peptidoform at its own m/z.
- o18 fitting: recognized (`--label o18`) but errors pending a post-M4
  rewrite (length + selected-residue coefficients, not a per-AA dict).
  o18 *integration* is unaffected.
- GUI: the broken Tkinter `riana_ui/` is deleted (M4 Phase 1); the PySide6
  replacement shipped in M4 Phase 2 (`riana gui` — Integrate + Model tabs).
- The bundled `workflow/Snakefile` is an example, not a tested pipeline.
  Treat as a starting point.
