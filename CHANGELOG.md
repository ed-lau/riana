# Changelog

All notable changes to Riana are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.0] — Unreleased

The breaking 1.0 release: a new package structure, peak detection, baseline
subtraction, mzTab intake, and a Qt GUI. See `PROJECT_REVIEW.md` §3 for the
roadmap. This section is built up iteratively as the milestones progress;
entries below are grouped by the work that produced them.

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
