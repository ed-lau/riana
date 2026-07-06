# Riana calibration benchmarks

Benchmark + analysis scripts, scored against the D₂O / ¹⁸O calibration mixing
series (`../data/calibration_d2o_mixing/`, `runs/calib_*`) and the in-vivo animal
sets. The current calibration baselines + dataset provenance live in
`reports/2026-06-23_calibration_benchmark_harness.md`; the shipped-feature record
is `CHANGELOG.md`. **Each script's module docstring is the authoritative detail —
run with `--help`;** the catalogue at the bottom is the index.

## Peptide populations ("curated" vs "uncurated")

The m0/mA recovery benchmarks (`bench_m0_ma_recovery.py`,
`bench_peak_boundary.py`, `bench_baseline.py`) report metrics on **three**
populations. These are **not** different coefficient tables — they are one
peptide set split by a quality gate.

The gate: for each peptide observed at all 9 mixing proportions, fit a line
through its 9 points of `m0_ma` (= iso0 / Σ iso0..iso5) vs. nominal heavy
proportion, and take that line's R².

| Population                  | Definition  | Meaning |
|-----------------------------|-------------|---------|
| **All**                     | any R²      | every peptide seen at all 9 proportions — this is the "any R²" set |
| **Curated (R² > 0.95)**     | R² > 0.95   | m0 abundance tracks the mixing proportion cleanly; the NB87a curation gate, and the only peptides used to fit the per-AA coefficient / frozen tables |
| **Uncurated (R² ≤ 0.95)**   | R² ≤ 0.95   | passed "seen at all 9 proportions" but failed the gate — co-eluting / low-SNR / noisy; where peak detection should help most (PROJECT_REVIEW.md, M2 finding 2) |

Two things that are easy to misremember:

- **"Uncurated" is the R² ≤ 0.95 complement, not "any R²".** The any-R²
  population is the **All** row.
- The gate is a **strict `>`** (`bench_aa_coefficients.py`: `r2 > r2_min`,
  `r2_min` default 0.95), so curated is R² > 0.95, not ≥.

`bench_m0_ma_recovery.py` keeps all peptides by running
`load_and_curate(r2_min=-1)` (filter disabled), then tags each row
`is_curated = r2 > 0.95`, so one pass yields all three populations.

## Coverage and weak fractions (`--drop-proportion`)

Curation has a second requirement, separate from the R² gate above: a peptide
must be **observed at every proportion** to enter the curated set. That makes
the curated count hostage to the *weakest acquisition* — one bad LC-MS run
discards every peptide missing from it, regardless of how well-behaved those
peptides are everywhere else.

The `cm` line is the worked example: its `time50` run is weak (~16k vs ~20k
target PSMs), and requiring all 9 proportions bottlenecks cm to 564 curated
peptides. `bench_aa_coefficients.py` and `bench_m0_ma_recovery.py` take
`--drop-proportion PCT` to exclude a named proportion from curation entirely
(coverage then requires only the survivors). For cm, `--drop-proportion 50`
recovers ~1.9k curated peptides and a markedly more stable frozen table — see
`reports/2026-06-23_calibration_benchmark_harness.md` (§ cm-drop50). Pass the **same** `--drop-proportion`
to the freeze inputs and to `bench_m0_ma_recovery.py` so the scored population
matches the table's training population.

## Script catalogue

**Calibration & coefficient training**

| Script | Purpose |
|--------|---------|
| `run_integrate_v0_9_0.py` / `run_integrate_v1_0_0.py` | run `riana integrate` over the calibration set → `*_riana.txt` |
| `run_calibration_benchmark.py` | the standing per-line calibration A/B harness (recovery vs ground truth) |
| `build_ground_truth.py` | build `ground_truth.csv` (proportion ↔ riana filename) |
| `bench_aa_coefficients.py` | NB87a: per-peptide Spep (IsoSpec) + per-AA non-neg regression |
| `build_frozen_tables.py` | bootstrap the frozen per-line `d2o_aa_coefficients_<line>.csv` (OOB R²) |
| `bench_coefficient_tables.py` | compare D₂O labelling-site tables on real in-vivo turnover |
| `bench_o18_coefficients.py` | NB90c: learn the ¹⁸O Spep length-model coefficients |
| `bench_fs_recovery.py` | NB87a step 4: fractional-synthesis recovery vs nominal proportion |
| `bench_fs_method_compare.py` | legacy iso0-only vs new full-envelope FS solver |
| `bench_m0_ma_recovery.py` | observed-vs-predicted m0 / envelope RMSE vs the frozen table |

**Peak detection & integration knobs**

| Script | Purpose |
|--------|---------|
| `bench_peak_boundary.py` | peak-boundary integration comparison (M3 Week 0 stub) |
| `bench_baseline.py` | baseline-subtraction integration comparison (M3 Week 0 stub) |
| `bench_smoothing.py` | Savitzky–Golay window sweep; coefficient / FS-bias shift |
| `bench_n_iso_sweep.py` | sweep N_ISO; test-R² vs isotopomer count |
| `bench_niso_crossover.py` | derive the `--fs` widening crossover (N_ISO vs init-width) |
| `bench_mixing_linearity.py` | model-free mixing-linearity sweep metric |
| `bench_zero_sweep.py` | model-free integration-setting sweep on the 0 % (unlabelled) proportion |
| `run_mixconfirm.py` | integrate several peak-detection configs on the full series |

**Fitting & FS**

| Script | Purpose |
|--------|---------|
| `bench_fit_recovery.py` | fit-module `k_deg₀` recovery via the pseudo-time trick |
| `bench_depth_semantics.py` | `--depth` semantics: rows vs distinct-timepoint qualification |
| `bench_mass_defect_theta.py` | prototype mass-defect θ_ΔS (v1.1.0 item 1b) |
| `bench_fs_rail_drop.py` | FS rail-drop ON/OFF A/B: yield / R² / matched within-protein CV / k (run at `--depth 6`) |
| `bench_fs_rail_threshold.py` | FS rail-drop `fs_rail_hi/lo` threshold sweep (why `1.05/−0.05` is the 1.2.0 default) |
| `bench_fs_rail_singlepoint.py` | FS rail bound on a single-timepoint (TMT) set — confirm the tighter rail is not harmful |

**Match-between-runs (MBR)**

| Script | Purpose |
|--------|---------|
| `bench_mbr_ab.py` | include-MBR vs `--exclude-mbr` A/B — does MBR help the fit? |
| `bench_mbr_quality.py` | do MBR-transferred peaks carry real, sensible signal? |
| `bench_mbr_apex_knobs.py` | sweep apex / extraction knobs on the ac16 calibration MBR |
| `bench_mbr_calibration.py` | ground-truth θ recovery for MBR on the mixing series |
| `bench_missingness.py` | per-precursor missingness across a labelling time series (MBR feasibility) |

**I/O, identity & QC**

| Script | Purpose |
|--------|---------|
| `bench_id_path.py` | ID-path concordance for the M3 I/O layer |
| `bench_rt_alignment.py` | is the quantms mzTab RT aligned across runs? |

**Track D — validation (in-vivo animal)**

| Script | Purpose |
|--------|---------|
| `build_lve_bench_set.py` | freeze the canonical within-protein-θ bench set |
| `bench_within_protein_theta.py` | within-protein fractional-labelling spread (no ground-truth FS) |
| `run_integrate_invivo.py` | integrate one in-vivo acquisition from a quantms mzTab → `_riana.txt` |
