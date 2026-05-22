# Riana calibration benchmarks

Benchmark scripts for the M3 rewrite, scored against the D₂O calibration
mixing series in `../data/calibration_d2o_mixing/`. See `PROJECT_REVIEW.md` §3
for how these gate each M3 milestone.

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
recovers 1817 curated peptides and a markedly more stable frozen table — see
the M2 addendum in `PROJECT_REVIEW.md`. Pass the **same** `--drop-proportion`
to the freeze inputs and to `bench_m0_ma_recovery.py` so the scored population
matches the table's training population.

## Scripts

| Script | Purpose |
|--------|---------|
| `build_ground_truth.py`    | build `ground_truth.csv` (proportion ↔ riana filename) |
| `run_integrate_v0_9_0.py`  | run `riana integrate` → the v0.9.0 `*_riana.txt` baseline |
| `bench_aa_coefficients.py` | NB87a port: per-peptide Spep + per-AA non-neg regression |
| `bench_fs_recovery.py`     | NB87a port: fractional-synthesis recovery vs. nominal proportion |
| `bench_n_iso_sweep.py`     | sweep N_ISO; test-R² vs. isotopomer count |
| `bench_smoothing.py`       | sweep Savitzky–Golay window; coefficient / FS-bias shift |
| `build_frozen_tables.py`   | bootstrap the frozen per-cell-line `d2o_aa_coefficients_<line>.csv` |
| `bench_m0_ma_recovery.py`  | observed-vs-predicted m0 / envelope RMSE vs. the frozen table |
| `bench_peak_boundary.py`   | peak-boundary integration comparison (M3 Week 0 stub) |
| `bench_baseline.py`        | baseline-subtraction integration comparison (M3 Week 0 stub) |
| `bench_fit_recovery.py`    | fit-module `k_deg₀` recovery via the pseudo-time trick |
