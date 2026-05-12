# Changelog

All notable changes to Riana are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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
