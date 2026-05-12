# Riana — Project Review & Roadmap

> Status: 0.9.0 (stabilization). Next: 1.0.0 (aggressive rewrite, lift the
> science). Maintainer: Edward Lau. Last reviewed: 2026-05-12.

This document consolidates and supersedes the prior `documentation/` folder
(`PROJECT_EVALUATION.md`, `ROADMAP.md`, `MASS_ACCURACY_SPECIFICATION.md`). The
old contents are folded in below where still load-bearing; the rest is cut.

## 1. Project status

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

Sequencing: M1 has shipped (you're reading the 0.9.0 review). M2 runs in
parallel with M1 (no code dependency between them). M3 is gated on M2's
ground-truth tables existing.

### M1 — Stabilization → 0.9.0 — DONE

Smallest viable bugfix release on the existing layout. Fixes B1–B13,
modernizes packaging, adds CI. Behavior-preserving except for the
documented mass-tolerance semantic correction. See `CHANGELOG.md`.

### M2 — Calibration test dataset (acquisition + dual-ID processing)

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

### M3 — Aggressive restructure → 1.0.0

Branch from the `v0.9.0` tag. Pre-req: M2 baseline numbers recorded so
the rewrite is regression-gated.

**Approach: walking skeleton, then fill in.** ≈3–4 weeks of focused work.

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
   mass-accuracy outputs inline. Validate each addition against M2
   benchmarks before committing.
4. **Week 4 — CLI + pipeline glue + fitting.** Build `cli.py`
   (typer or click). Rewrite `core/fitting.py` to consume
   `IntegrationResult` records. Output provenance header
   (git SHA, riana version, config hash) via `io/writers.py`.

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

**Verification (regression-gated by M2):**

- Calibration benchmark recovers nominal labelling fractions at least as
  accurately as 0.9.0 baseline, peptide by peptide.
- Percolator-ID and mzTab-ID paths produce m0/mA values that agree within
  tolerance (cross-format A/B is itself a validation of the mzTab adapter).
- `tests/data/sample1/` end-to-end smoke test produces a `_riana.txt`
  whose per-peptide m0/m6 agree with 0.9.0 within 1e-3 relative tolerance
  (it should be near-identical — same numerical core).
- Memory peak on a 2 GB mzML drops from "all of it" to "one fraction
  worth."

### M4 — Qt + async GUI rewrite

PySide6 (LGPL) + `qasync` (bridges asyncio with the Qt event loop) under
`riana/gui/`. Long-running CPU work via `ProcessPoolExecutor` driven from
async tasks. `pyqtgraph` for fast embedded chromatogram inspection;
matplotlib only for static export. Tabs: Integrate, Model, Calibration
(after M3). Entry point: `riana gui` subcommand.

Drops: `rx`, `sv_ttk`, `pandastable`, the missing `console` shim,
Tkinter dependencies entirely.

### M5 — Flexible D₂O reporting (record only; defer design)

User-requested: per-timepoint **fraction-new** as a first-class output
(today, `riana fit` reports per-peptide `k_deg` only; per-timepoint `fs`
is intermediate and not persisted). Likely surface: `riana fit
--emit-fraction-new` producing a long-format table
`(concat, sample, t, fs, fs_lower, fs_upper)`, with CIs via bootstrap
residuals.

The math is already mostly in `fsynthesis.py` — the gap is in serializing
intermediate `fs` per-timepoint instead of only the fitted constant.
This stays a planning placeholder until scoped further.

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
4. **Snakemake.** Keep `workflow/Snakefile`; update rules to use whichever
   CLI surface lands in 1.0. Add a parallel `workflow/Snakefile.quantms`
   after M3.
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

- Memory: full mzML loaded into a Python list of numpy arrays
  (`spectra.py`). 1.0 fixes this with streaming/indexed reads; 0.9.0
  carries the limit.
- Peak fidelity: see §2c. The calibration dataset will quantify the
  impact; until then, treat integrated areas as comparable across
  samples but not as absolute peak areas.
- AA labelling (`--label 4`): `a_0` dispatch bug (§2b) means the fit
  uses an incorrect baseline. Fix lands in 1.0.
- GUI: broken on a fresh clone today (`console` import missing).
  Replaced wholesale in 1.0.
- The bundled `workflow/Snakefile` is an example, not a tested pipeline.
  Treat as a starting point.
