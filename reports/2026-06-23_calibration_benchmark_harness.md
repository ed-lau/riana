# Standing calibration benchmark — layout, baselines, and tuning workflow

- **Date:** 2026-06-23
- **Status:** **driver shipped** 2026-06-23 (`tests/benchmark/run_calibration_benchmark.py`);
  reproduces the committed `v1.0.0_fs0123` anchor on all three lines (within±0.05
  ac16 24.8 % / cm 21.5 % / ipsc 28.9 %). Design + baselines below.
- **Purpose:** a persistent, per-cell-type calibration harness so future integration-knob
  and MBR tuning is a one-command A/B against a recorded baseline, not an ad-hoc re-run
- **Lines:** ac16, ipsc, cm (D₂O mixing series, 0→100% with ground-truth proportions)

## Why

We keep re-deriving the same calibration recovery numbers by hand (this session ran the
capture×scoring sweep three times before it stabilised). The mixing series is the one
dataset with **ground-truth θ = mixing proportion**, so it is the right harness for
*every* integration/fit knob (mass tol, peak picker, baseline, N_ISO/`--fs`, MBR). It
deserves a standing home with recorded baselines.

## Current state (fragmented — the thing to unify)

Two integrate paths and two output roots exist, at different configs:

| artifact | path | id path | config | has recovery? |
|---|---|---|---|---|
| `calib_<line>_v1` | `runs/calib_<line>_v1/` | SDRF/mzTab | mass_tol 10 (SDRF) | yes — `bench/` (m0/mA train/test R²) |
| `cal_ac16_mbr` | `runs/cal_ac16_mbr/` | SDRF + `--mbr` | mt10 | no |
| version labels | `tests/data/calibration_d2o_mixing/<line>/integrate_outputs/<label>/` | Percolator | many (v0.9.0, v1.0.0=ms2/mt15, adaptive, apex_tall_w015, ms2_w015, …) | via `bench_*` on demand |

Baseline recovery (`calib_<line>_v1/bench/summary.json`, m0/mA train/test R², n_iso=4,
R²_min 0.95): **ac16 train 0.884 / test 0.894** (1676 peptides). ipsc/cm analogous.

Tools that already do the work (the driver just orchestrates them):
- `run_integrate_v1_0_0.py` — Percolator-path integrate (`--adaptive`/`--ria` added this
  session); `riana integrate --sdrf` — the SDRF-path integrate used by `calib_<line>_v1`.
- `bench_fs_method_compare.py` — |θ−f| recovery vs ground truth, per proportion
  (`--score-channels` added this session).
- `bench_within_protein_theta.py` — within-protein θ-spread (also `--score-channels`).
- `bench_m0_ma_recovery.py` — the m0/mA train/test R² in `calib_<line>_v1/bench/`.

## Baseline — capture×scoring recovery sweep (2026-06-23, this session)

`bench_fs_method_compare.py`, |θ−f| vs ground truth, RIA 0.0598, per-line coefficients
(cm uses drop50 + `--drop-proportion 50`). **within ±0.05 at f=1** (fix arm; capture is
irrelevant so adapt matches):

| line | all | iso0-1 | iso0-3 |
|---|---|---|---|
| ac16 | 25.5% | 27.1% | 27.3% |
| cm   | 21.3% | 25.4% | 23.2% |
| ipsc | 30.4% | 31.9% | 32.7% |

Verdict (see `2026-06-23_adaptive_niso_limited_isotopomer.md`): **iso0-3 (`--fs 0 1 2 3`)
is the robust cross-line choice; adaptive capture is neutral.** Full per-proportion
tables in `runs/_sweep_all.out` (regenerable via `runs/_sweep_all.py`).

## Proposed canonical layout

```
runs/calib_<line>/                      # one dir per cell line
  <config_label>/                       # one per knob setting (baseline, mt15, mbr, fs0123, …)
    <sample>_riana.txt                  # integrate output (9 proportions)
    riana_manifest.tsv
    integrate.log
    recovery/                           # bench_fs_method_compare output
      fs_method_compare_summary.json
    within_protein/                     # bench_within_protein_theta (optional)
    config.json                         # the exact IntegrationConfig used
  BASELINE.md                           # the frozen baseline numbers for this line
```

Canonical baseline config (proposed): the **adopted production defaults** — mass_tol 10,
`peak_rt=apex`, `integration_half_width=0.15`, `--iso 0 1 2 3` capture, fit `--fs 0 1 2 3`.
(The `v1.0.0` ms2/mt15 set stays as the 0.9.0-parity anchor.)

## Tuning / MBR workflow (what the harness buys)

```bash
# 1. baseline (once): integrate + recovery at the canonical config -> runs/calib_<line>/baseline/
# 2. a knob A/B: re-integrate only the changed knob into runs/calib_<line>/<knob>/, recovery,
#    diff the per-proportion within-±0.05 / IQR / bias vs baseline.
# 3. MBR: integrate --mbr into runs/calib_<line>/mbr/ ; compare recovery + curve completeness
#    vs baseline (NB the calibration RT-offset caveat — mzTab RT is OpenMS-aligned ~2 min off
#    the local mzML, so MBR-on-calibration needs the per-run RT correction, not a re-search;
#    see 2026-06-23_adaptive_niso_limited_isotopomer.md §4 + memory mbr_rt_axis_dependency).
```

## Next build — `tests/benchmark/run_calibration_benchmark.py` (spec)

A thin orchestrator (no new science): `--line {ac16,cm,ipsc,all} --label <cfg> [integrate
flags] [--fs N]`. Per line: integrate (or reuse if `<label>/` exists + config matches) →
`bench_fs_method_compare` (recovery) → optional `bench_within_protein_theta` → write the
canonical layout above + a one-line summary appended to `BASELINE.md`. Reuses
`run_integrate_v1_0_0.run_line` and the bench mains; ~80 lines. Deferred to a fresh
context (this session was long); the data + tools already exist, so it is pure plumbing.

## Reproduce (today, without the driver)

```bash
python tests/benchmark/run_integrate_v1_0_0.py --line ac16 --out-label v1.0.0   # fixed
python tests/benchmark/run_integrate_v1_0_0.py --line ac16 --adaptive --out-label adaptive
python runs/_sweep_all.py     # the full {fix,adapt}×{all,iso0-1,iso0-3} recovery sweep
```

---

## Datasets & refreshed current state (merged from PROJECT_REVIEW M2, 2026-06-28)

The calibration datasets and their current-defaults numbers, consolidated here
(was PROJECT_REVIEW §3 M2). These are the per-cell-line **D₂O mixing series** —
lysate from fully-labelled cells (≥10 doublings in 6 % D₂O) mixed with unlabelled
lysate at nine nominal heavy fractions `0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5,
100 %` — plus the AC16 **¹⁸O** mixing series. Per peptide × fraction the expected
isotopomer envelope is computable from sequence + enrichment
(`get_peptide_distribution`), so the headline metric is observed-vs-predicted
recovery, gated by the NB87a curation rule (seen at every proportion, per-peptide
m0/mA-vs-proportion R² > 0.95).

### The lines

| line | cell model | label | JPOST | note |
|---|---|---|---|---|
| **ac16** | AC16 human cardiomyocyte (proliferative) | D₂O 6 % | `JPST002443` | |
| **ipsc** | human iPSC (proliferative) | D₂O 6 % | `JPST003556` | |
| **cm** | iPSC-derived cardiomyocyte (**post-mitotic**) | D₂O 6 % | `JPST003582` | drop50 (below) |
| **ac16 ¹⁸O** | AC16 | H₂¹⁸O ~9 % | — | trains `juber_2026_o18_ac16`; see the ¹⁸O reports |

ac16 and ipsc are proliferative (division dilutes label independently of
turnover); the post-mitotic **cm** line isolates turnover from division — the
biology axis for the cross-line coefficient comparison.

### cm drops its 50 % proportion (`cm_drop50`)

cm's `time50` is a weak acquisition (~16 k vs ~20 k target PSMs; ~5.7 k integrated
peptides). Because curation requires a peptide seen at **every** proportion, that
one bad *run* — not a bad peptide set — bottlenecks 9/9 cm to a third of its
peptides and inflates near-zero-labelling coefficients. The "seen at all" rule is
a *coverage* requirement separate from the R² quality gate, so cm is curated on
the **8 surviving proportions** (`--drop-proportion 50`): this recovers ~3× the
peptides, cuts coefficient bootstrap noise ~36 %, and collapses the unphysical
low-labelling coefficients. `d2o_aa_coefficients_cm_drop50.csv` (→ the shipped
`alamillo_2025_cm`) is the cm reference; the 9/9 table is kept as the A/B record.

### Current numbers (calib_v1 — riana 1.1.0 defaults: ±10 ppm from SDRF, apex / 0.15 min, iso0-5)

Refreshed from `runs/calib_{ac16,cm,ipsc}_v1` + `runs/calib_o18` (replaces the
v0.9.0 baseline table). **OOB R²** = out-of-bag R² of the bootstrap per-AA Spep
fit on the curated set (`build_frozen_tables.py`); **m0 RMSE** = observed-vs-
predicted m0/mA recovery RMSE (`m0_ma_recovery_summary.json`), split by the
R² > 0.95 curation gate.

| line | props (gate) | seen → curated (yield) | OOB R² (95 % CI) | m0 RMSE cur / uncur |
|---|---|---|---|---|
| ac16 D₂O | 9 | 3,165 → 1,676 (53 %) | **0.883** [.853, .911] | 0.017 / 0.065 |
| ipsc D₂O | 9 | 4,840 → 3,129 (65 %) | **0.867** [.842, .891] | 0.023 / 0.080 |
| cm D₂O | 8 (−50 %) | 3,673 → 1,886 (51 %) | **0.856** [.824, .882] | 0.024 / 0.073 |
| ac16 ¹⁸O | 9 | 7,900 → 2,083 (26 %) | **0.892** [.869, .912] | ¹ |

¹ The ¹⁸O line's recovery is validated in `2026-06-24_o18_reverse_model.md` /
`2026-06-26_o18_kinetic_fit.md` (curated tracks the proportions, \|bias\| ~0.03–0.07
mid-range); it has no `m0_ma_recovery` artifact in this harness. Its low curated
yield (26 %) is the 8 h-window AC16 regime, not a label/pipeline defect — see
`2026-06-28_ac16_curation_diagnosis.md`.

**Reading the table:** OOB R² (coefficient-fit generalization) is healthy and
similar across lines (0.86–0.89); the in-vitro **ipsc** line curates best (65 %),
cm next (51 % at 8/9), and the ¹⁸O AC16 line lowest (window-limited). Curated m0
recovery RMSE is tight (≤ 0.024) and ~3× lower than uncurated everywhere — the
R² gate is doing its job. These are the standing baseline any integration-knob /
MBR A/B is scored against (workflow above).
