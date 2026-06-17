# MBR feasibility for the mzTab/DDA path — missingness + RT-alignment measurement

- **Date:** 2026-06-17
- **Branch:** `m3-rewrite`
- **Status:** investigation complete — **GO** for MBR on the mzTab/DDA path (DIA needs none)
- **Benches:** [`tests/benchmark/bench_missingness.py`](../tests/benchmark/bench_missingness.py), [`tests/benchmark/bench_rt_alignment.py`](../tests/benchmark/bench_rt_alignment.py)
- **Roadmap:** `PROJECT_REVIEW.md` → Track A ("Match-between-runs … MEASURED 2026-06-17") + handoff priority #4

## Question

Match-between-runs (MBR) transfers a precursor's identity + retention time from a
run where it was identified into a run where it was *not* picked for MS2, filling
holes in a turnover curve that come from identification gaps rather than real
absence. Before building it for the quantms **mzTab/DDA** path, two questions:

1. **Is it worth it?** How incomplete are DDA turnover curves, and how much of the
   missingness is *recoverable* (a donor run exists) vs genuinely absent? Does DIA
   (DIA-NN, which propagates IDs across runs) already solve this, confirming DDA is
   where MBR pays?
2. **Is it feasible?** MBR needs the acceptor-run RT. Does the quantms mzTab carry
   **RT-aligned** retention times (use the donor RT directly), or **raw** per-run
   RTs (we must align ourselves on shared IDs first)?

## TL;DR

- **DDA curves are badly gappy and the gaps are mostly MBR-recoverable.** Only
  8–15% of precursors span all 12 timepoints; the median precursor is seen in 4–5
  of 12 runs; ~40–45% of the slots belonging to precursors that *have* a donor are
  fillable holes. The **t=0 anchor** (the m0 baseline) is the sharpest loss —
  **70%** of later-seen LVE precursors lack it, because that one acquisition ran
  shallow.
- **DIA does not need MBR** (75% complete curves, 5.6% recoverable gap) — DIA-NN's
  internal propagation already does the job. **DDA is where MBR pays.**
- The DDA↔DIA gap **survives a geometry control** (matched to 3 timepoints, 1
  acquisition per slot: DIA 60% vs DDA 32–44%), so it is the acquisition + ID
  pipeline, not the 12-vs-3 timepoint-count difference.
- **quantms RT *is* aligned** (typical cross-run median |ΔRT| ~4–9 s) **but
  imperfectly** — run-specific residuals reach **15–25 s**, larger than Riana's
  narrow integration window and worst on the same shallow anchor runs MBR most
  needs. The fix is a light **robust per-run RT refinement** on shared IDs (richly
  supported: 3,900–6,100 co-IDs per run pair), not a from-scratch aligner.
- **No architectural block** — see [Architecture](#architecture-no-block).

## Data & method

| set | source (gitignored `runs/`) | acquisition + IDs | geometry |
|-----|------------------------------|-------------------|----------|
| LVE (DDA) | `runs/lve_atr/*_LVE_time*_riana.txt` | Orbitrap DDA, quantms mzTab (q≤0.01 PSMs) | 12 timepoints × 1 rep |
| ATR (DDA) | `runs/lve_atr/*_ATR_time*_riana.txt` | Orbitrap DDA, quantms mzTab | 12 timepoints × 1 rep |
| LV (DIA) | `runs/lve_dia/*_riana.txt` | QEHF DIA, DIA-NN (propagated) | 3 timepoints × 3 bioreps |

- **Identity key** = the `concat` column = `SEQUENCE_charge` (precursor level — what
  MBR transfers; integration needs a specific m/z, so charge states are distinct). A
  sequence-level view is reported for reference but the conclusions are the same.
- `bench_missingness.py` builds a precursor × run **presence matrix** per curve from
  the `_riana.txt` integrate outputs (no fit needed).
- `bench_rt_alignment.py` reads the **raw** quantms mzTab PSM section
  (`data/timeseries_lve_atr/quantms_results/quant_tables/…openms.mzTab`), keeps
  target PSMs (q≤0.01, non-decoy), takes the median RT per `(run, precursor)`, and
  compares co-identified precursors across run pairs.
- **Caveat — not a controlled contrast.** DDA (cardiac LVE/ATR, Orbitrap) and DIA
  (LV *atf6small*, QEHF) are different samples/instruments. The DDA numbers stand on
  their own; DIA is the "what complete looks like" reference, not a paired control.

## Results

### 1. Missingness (precursor = `SEQUENCE_charge`)

| set | tp | precursors (union) | complete curves | median tp/precursor | **t0-anchor loss** | overall gap | **recoverable gap** |
|-----|----|--------------------|-----------------|---------------------|--------------------|-------------|---------------------|
| **LVE** (DDA) | 12 | 15,318 | 2,252 (**14.7%**) | 5/12 | **69.5%** | 50.5% | 39.2% |
| **ATR** (DDA) | 12 | 16,631 | 1,394 (**8.4%**) | 4/12 | 39.9% | 57.5% | 45.3% |
| **LV** (DIA) | 3 | 28,916 | 21,599 (**74.7%**) | 3/3 | 9.8% | 11.9% | 5.6% |

- **complete curves** = precursors present in *every* timepoint.
- **t0-anchor loss** = of precursors seen at any *later* timepoint, the fraction
  missing t=0 (the m0 baseline that pins the kinetic fit).
- **overall gap** = empty (precursor,run) slots / all slots, over every observed precursor.
- **recoverable gap** = empty slots among precursors seen in ≥2 runs (a donor exists)
  ⇒ the MBR-addressable population.

Coverage yield (precursors seen in ≥k of 12 timepoints) — MBR would lift precursors
up these rungs:

| ≥k / 12 | LVE | ATR |
|---------|-----|-----|
| ≥2  | 12,033 | 12,251 |
| ≥3  | 10,431 | 10,251 |
| ≥6  | 7,503  | 6,843  |
| ≥10 | 4,436  | 3,267  |
| =12 | 2,252  | 1,394  |

**The t0 anchor is acquisition-quality driven.** LVE `time00` is the *shallowest*
run in its series (**4,792** precursors vs ~7,500 mid-series), so 70% of
later-seen precursors have no m0 baseline. ATR `time00` is its *deepest* (10,892),
so its t0 loss is only 40%. You can't easily re-shoot a flagship acquisition — but
you can transfer later IDs into it. This is MBR's highest-value case.

### 2. Geometry control (3 timepoints, DIA-matched)

DIA has only 3 timepoints, so completing a 3-slot curve is far easier than a
12-slot one — the raw 75% vs 15% partly reflects that. Re-measuring DDA over 3
timepoints isolates the per-acquisition loss from the curve-length effect:

| set | configuration | complete |
|-----|---------------|----------|
| LVE (DDA) | 3-of-11 non-t0 tp, mean over 165 combos | **43.5%** (range 34–52%) |
| ATR (DDA) | 3-of-11 non-t0 tp, mean over 165 combos | **32.5%** (range 18–51%) |
| LV (DIA) | 3 tp × 3 reps | 74.7% |
| LV (DIA) | 3 tp × **1 rep** (matches DDA's 1 acq/slot) | **60.3%** |

At matched geometry **and** matched shots-per-slot (DIA 1-rep 60% vs DDA 32–44%),
the gap survives. So the missingness is the **acquisition mode + ID pipeline**
(stochastic MS2 + no PSM-level propagation in the mzTab path), **not** curve length.

### 3. RT alignment — does quantms align RT?

Pairwise **median |ΔRT| (seconds)** for co-identified precursors across the LVE
timepoints (identity-line residual, no fit):

```
       t00   t01   t02   t03   t04   t06   t08   t10   t15   t20   t25   t30
 t00     ·  17.4  15.7  16.5   5.4  15.3  16.2  13.5  12.6  15.8  12.8   9.4
 t01  17.4     ·   5.3   4.9  16.3   5.5   4.8   6.6   7.7   5.6   8.2  10.7
 t02  15.7   5.3     ·   4.3  13.7   4.5   5.5   4.9   6.1   5.4   6.1   9.0
 t03  16.5   4.9   4.3     ·  14.6   4.2   4.1   5.2   6.3   4.9   6.5   9.4
 t04   5.4  16.3  13.7  14.6     ·  13.2  14.7  11.7  11.0  14.9  11.6   7.6
 t06  15.3   5.5   4.5   4.2  13.2     ·   4.3   4.1   5.3   4.9   5.5   8.7
 t08  16.2   4.8   5.5   4.1  14.7   4.3     ·   5.1   6.0   4.4   6.6   9.5
 t10  13.5   6.6   4.9   5.2  11.7   4.1   5.1     ·   3.7   5.3   3.8   6.8
 t15  12.6   7.7   6.1   6.3  11.0   5.3   6.0   3.7     ·   5.7   3.7   5.9
 t20  15.8   5.6   5.4   4.9  14.9   4.9   4.4   5.3   5.7     ·   5.6   9.5
 t25  12.8   8.2   6.1   6.5  11.6   5.5   6.6   3.8   3.7   5.6     ·   6.3
 t30   9.4  10.7   9.0   9.4   7.6   8.7   9.5   6.8   5.9   9.5   6.3     ·
```

Representative LVE pairs vs t0 (slope≈1 & intercept≈0 & small RMSE ⇒ aligned):

| pair | n | med\|ΔRT\| | p90 | slope | intercept | r | rmse_id | rmse_fit |
|------|---|-----------|-----|-------|-----------|---|---------|----------|
| t00→t01 | 3,928 | 17.4 | 56.9 | 1.010 | −33.86 | 0.9982 | 57.0 | 55.9 |
| t00→t08 | 3,880 | 16.2 | 51.4 | 1.008 | −28.47 | 0.9984 | 54.3 | 53.5 |
| t00→t30 | 3,892 |  9.4 | 27.5 | 1.003 |  −7.71 | 0.9991 | 39.9 | 39.8 |

**Reading it:**

- **RT is aligned, not raw.** Most run pairs sit within **~4–9 s** median |ΔRT| —
  on the order of a chromatographic peak width, and matching the `PROJECT_REVIEW`
  gotcha note's suspected "~0.1 min" offset. So quantms/ProteomicsLFQ does put runs
  on a common frame; MBR can start from the donor RT.
- **But the alignment is imperfect and run-specific.** In LVE, `t00` and `t04` form
  their own cluster (~5 s from each other, **13–17 s** from everything else); in ATR
  the bad run is `t25` (15–25 s). These residuals **exceed Riana's narrow
  integration window** (`integration_half_width` 0.15 min = ±9 s), and in LVE the
  worst-aligned run (`t00`) is the *same shallow run* that most needs MBR.
- **It's an offset, not a warp.** Slopes are ≈1.00–1.03 and a global linear realign
  barely reduces the residual (`rmse_id` ≈ `rmse_fit`) — the systematic part is a
  near-constant per-run shift, while the large RMSE (40–165 s) is a heavy tail of
  precursor-key outliers (same `SEQUENCE_charge` eluting at ≥2 RTs). So the right
  correction is a **robust per-run offset / LOESS** on shared IDs (or RT-anchor +
  local apex re-find), not a global least-squares line.
- **Self-alignment is well-supported.** Every run pair shares **3,900–6,100**
  co-identified precursors — far more than enough anchors to fit a per-run-pair
  alignment.

## Architecture (no block)

- `records.py` already reserves the provenance value `evidence = "q_value" | "mbr"`.
- An MBR-transferred precursor has **no MS2 scan** in the acceptor run. The DIA
  intake (M6b, shipped 2026-06-11 — one day before this was first scoped) already
  built exactly that extraction path: `io/diann` emits `scan = -1` and
  `core/integration.resolve_rt_anchored_scans` maps an apex RT to the nearest MS1
  scan in the mzML, then the normal scan-based extraction runs. **MBR's extraction
  substrate already exists.**
- The one genuinely new piece is the **donor→acceptor transfer + RT refinement** at
  the intake layer (an ID-assembly concern → Track A). The stale `data/mbr_test`
  fixture (Oct-2022, pre-rewrite) is not reusable as-is.

## Recommended next step — MBR design (to be detailed separately)

Scope it at the **intake layer**, per condition group from the manifest:

1. **Donor pool.** For each `(experiment, condition)` curve group, union the target
   precursors across its runs; a precursor identified in ≥1 run is an MBR donor for
   the runs missing it. (Tunable: require ≥N donor observations to suppress
   one-shot false IDs.)
2. **Per-run RT alignment.** For each acceptor run, fit a **robust** RT map from its
   shared IDs to the donor frame (median offset / LOESS — *not* global linear). The
   measurement above shows this is needed (15–25 s residuals) and feasible
   (thousands of anchors).
3. **Transfer.** For each missing (precursor, acceptor-run), emit a synthetic record
   with the aligned RT, `scan = -1`, and `evidence = "mbr"`; carry the donor's
   charge/sequence/mods.
4. **Extract.** Route through the existing `resolve_rt_anchored_scans` →
   RT-anchored MS1 extraction. Let the apex finder re-center within a tolerance that
   covers the worst-case alignment residual (~25 s).
5. **Control quality.** A transfer-FDR / score gate; flag MBR rows so downstream fit
   and rollup can weight or audit them. Remember **recoverable-gap is an upper bound**
   — MBR transfers an ID+RT, but a genuinely below-LOD precursor still integrates to
   noise; quantify yield-and-verify, don't assume.

## Reproduce

```bash
# inputs live under the gitignored runs/ (regenerable via integrate)
python tests/benchmark/bench_missingness.py            # missingness + geometry control
python tests/benchmark/bench_rt_alignment.py --chamber LVE   # RT alignment matrix
python tests/benchmark/bench_rt_alignment.py --chamber ATR
```
