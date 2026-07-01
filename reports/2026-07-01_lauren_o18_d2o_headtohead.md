# Lauren ¹⁸O/D₂O — second iPSC head-to-head + mesoderm ¹⁸O series

- **Status:** **Integrated + fit + rolled up.** Two new analysis-grade sets harden §5.2
  of the ¹⁸O technical note: a same-cells (SCVI480 iPSC) D₂O↔¹⁸O head-to-head and a
  mesoderm ¹⁸O series. **The trustworthy §5.2 readouts reproduce on a second iPSC line**
  (ranking ρ ≈ 0.67, within-protein geom CV ≈ 0.14); **curation yield is 3–4× lower than
  boomi** — diagnosed (sampling design × R² metric, below). **Data source:** the ¹⁸O run is
  `runs/timeseries_lauren5_7_ipsc_mesoderm_o18` (quantms re-searched without Sage, which hangs
  on the crash-recovered empty-MS1 scans; mzML symlinked from the network drive). Numbers are
  unchanged from the prior `lauren5_7_final_20260630` integration to within rounding.
- **Sibling reports:** [`2026-06-26_o18_kinetic_fit.md`](2026-06-26_o18_kinetic_fit.md)
  (boomi/juber head-to-head), [`2026-06-28_spep_curation_gate.md`](2026-06-28_spep_curation_gate.md).

## TL;DR

- **The method reproduces across iPSC lines.** SCVI480 (lauren) D₂O-vs-¹⁸O ranking
  ρ = **0.67** peptide and within-protein geom CV **≈ 0.14** land on top of the boomi AICS52
  head-to-head (peptide ρ 0.68, CV 0.14) — the readouts §5.2 relies on hold on independent
  cells. Protein-level (linear-φ rollup) ρ is lower for lauren (**0.55** vs boomi 0.66) —
  the sparse-data rollup effect (yield section).
- **D₂O absolute scale is reproducible, ¹⁸O is not** — as §5.2 argued. Cross-line
  per-protein (rollup): D₂O ρ = 0.56, log₂(boomi/lauren) = −0.08 (parity); ¹⁸O ρ = 0.44,
  log₂ = +0.70 (~1.6× shift). The ¹⁸O table is the transferred AC16 one; D₂O uses the
  iPSC-matched `alamillo_2025_ipsc`.
- **Yield — diagnosed.** lauren R²>0.8 yield is **¹⁸O 9.7% / D₂O 14.1%** vs boomi 42.9 / 33.2%,
  at identical CV/ρ — a sampling-design × R²-metric effect (front-loaded 0→7 h grid + flat-curve
  pathology + ~1.5× noise floor), not the estimator and not lost IDs. Full work-up:
  [`2026-07-01_lauren_yield_gap_diagnosis.md`](2026-07-01_lauren_yield_gap_diagnosis.md).
- **This run was only possible after this cycle's integrate fixes** (mass-based
  scan↔precursor guard, per-run failure isolation, git-SHA cache, ~8× MS1 precache,
  empty-scan handling — CHANGELOG 2026-06-30/07-01). The 384-file multi-fraction set
  froze/failed on every one of those before the fixes.

## Data

| dataset | cells | label | RIA (SDRF) | timepoints | window | bioreps × fractions | runs |
|---|---|---|---|---|---|---|---|
| lauren_5_7 (control) | SCVI480 iPSC | ¹⁸O | 0.0794 | 12 | 0→24 h | 2 × 8 | 192 |
| lauren_5_7 (mesoderm) | SCVI480 meso | ¹⁸O | 0.0794 | 12 | 0→24 h | 2 × 8 | 192 |
| lauren_9 | SCVI480 iPSC | D₂O | 0.0598 | 12 | 0→24 h | 2 × 8 | 192 |
| boomi (ref) | AICS52 iPSC | ¹⁸O / D₂O | 0.0897 / 0.0598 | 9 | 0→24 h | 2 × 1 | 18 |

Multi-fraction (8 LC fractions) — a first for Riana; collapsed at fit (`--fraction-collapse
sum`, a no-op on single-fraction boomi/juber). lauren_9's 2 empty-scan fractions were
backfilled via `--resume` after the empty-scan fix (2026-06-30) — **all 192 runs are in this
fit** (verified: 192 `_riana.txt`, 192 integrate manifest rows).

## Curation + within-protein consistency (R²>0.8 + label-aware Spep; o18≥6, D₂O≥8)

| arm | n_fit | yield (R²>0.8) | curated | within-prot geom CV (n_prot) |
|---|---|---|---|---|
| **lauren5_7 iPSC ¹⁸O** | 16,413 | **9.7 %** | 1,400 | **0.159** (142) |
| lauren5_7 mesoderm ¹⁸O | 9,604 | 15.2 % | 1,298 | 0.143 (143) |
| **lauren9 iPSC D₂O** | 31,391 | **14.1 %** | 4,328 | **0.137** (522) |
| boomi iPSC ¹⁸O — ref | 20,021 | 42.9 % | 5,704 | 0.140 (650) |
| boomi iPSC D₂O — ref | 18,448 | 33.2 % | 5,708 | 0.139 (620) |

Within-protein geom CV is **flat at ~0.14–0.16 across every arm and both labels** — the
same value boomi and the §5.2 juber/iPSC sets show. Yield is the outlier (below).

## Head-to-head — same cells (SCVI480 iPSC), D₂O vs ¹⁸O

| metric | lauren SCVI480 | boomi AICS52 (§5.2) |
|---|---|---|
| Spearman(k) — peptide | **0.668** (n=749) | 0.68 (n=3002) |
| Spearman(k) — protein (linear-φ rollup) | **0.55** (n=792) | 0.66 (n=1067) |
| within-prot geom CV: ¹⁸O / D₂O | 0.159 / 0.137 | 0.140 / 0.139 |
| median log₂(k_D₂O / k_¹⁸O) | **+0.17** (≈ parity) | −0.46 (¹⁸O ≈ 1.4× faster) |

Peptide ranking and within-protein CV **reproduce**; protein-rollup ρ softens for lauren
(0.55 vs boomi 0.66) as the rollup pools a broader, noisier protein set on the low-yield data.
The absolute D₂O/¹⁸O scale does not reproduce (+0.17 vs −0.46) — expected, since ¹⁸O uses the
transferred AC16 table and the scale is line-dependent (§5.2). n is smaller here because of the
lower yield, not selection. **Protein-level throughout this report = linear-φ rollup k
(`riana rollup`, R²≥0.8) merged on accession, matching the §5.2 note.**

## Cross-dataset reproducibility (per-protein, same label, SCVI480 vs AICS52)

| pair | peptide ρ | protein ρ (rollup) | median log₂(boomi/lauren), rollup |
|---|---|---|---|
| lauren ¹⁸O vs boomi ¹⁸O | 0.428 (n=346) | 0.44 (n=728) | **+0.70** (~1.6× scale shift) |
| lauren D₂O vs boomi D₂O | 0.584 (n=890) | **0.56** (n=932) | −0.08 (parity) |

The D₂O per-protein rates agree across two independent iPSC lines in both **ranking and
scale**; ¹⁸O agrees on ranking but carries a ~1.6× scale offset from the transferred table.
This is the cleanest statement yet of "D₂O matched-table = quantitative; ¹⁸O transferred =
rank-only" — and it holds between lines, not just within one.

## Yield — DIAGNOSED (sampling design × R² metric, + ~1.5× noise)

lauren yield (9–15 %) sits ~3–4× below boomi (33–43 %) while CV and ρ are unchanged.
Full work-up in [`2026-07-01_lauren_yield_gap_diagnosis.md`](2026-07-01_lauren_yield_gap_diagnosis.md)
(`runs/lauren_yield_diagnostics.py`). Verdict: **it is a sampling-design × metric effect that
drops n without biasing the surviving numbers — not slower turnover and not a failed MS run.**

- **Not turnover** — k distributions comparable (lauren ¹⁸O med 0.037 vs boomi 0.045; lauren
  D₂O 0.039 *faster* than boomi 0.032).
- **Not lost IDs** — lauren has *more* PSMs, proteins, and timepoints (median n_points 18 vs 15).
  Curation is a quantitative-SNR gate, not an identification one — hence 7 k IDs above 3–4× lower
  curation. The "few shared peptides" (749 vs 3002) is the intersection of two ~10–15 % arms
  (arithmetic), not disagreement (ρ still 0.67).
- **Dominant cause** — lauren crams **11 of 12 timepoints into 0→7 h** (FS ≤ 0.23, signal ≈ noise)
  with no sampling on the 7→24 h rising shoulder, so the pinned-exponential R² is depressed even
  for good fits (flat-curve pathology; R²-vs-k arch peaks mid-k, collapses low/high-k in *both*
  sets). Compounded by a **~1.5× higher noise floor** (t0 FS MAD 0.13–0.14 vs boomi 0.09).

### Curation-gate sensitivity — a CI (dk) rescue (Lau 2018)

R² is the wrong *sole* gate for flat / low-dynamic-range curves; the principled fix gates on the
**rate-constant's relative uncertainty** `relunc = (ci_hi − ci_lo)/(2k)` (90 % bootstrap CI
half-width over k — a coefficient of variation of k̂; Lau 2018 / Sadygov d2ome). Second gate =
**R² > 0.8 OR (R² > 0.6 AND relunc < 0.25)** (`runs/gate_comparison.py`):

| arm | yield (strict → rescue) | within-prot geomCV | D₂O↔¹⁸O pep ρ |
|---|---|---|---|
| lauren ¹⁸O | 8.5 → **17.6 %** | 0.159 → 0.213 | — |
| lauren D₂O | 13.8 → **24.7 %** | 0.137 → 0.193 | — |
| lauren head-to-head | — | — | 0.668 → **0.644** (n 749 → 1751) |
| boomi ¹⁸O / D₂O (ref) | 28.5→35.5 / 30.9→42.5 % | 0.140→0.180 / 0.139→0.184 | 0.679 → 0.656 |

The rescue ~**doubles lauren yield**, narrowing the gap to boomi (~2× from ~3.4×), while the
within-protein CV stays clean (< 0.25), the peptide ranking is preserved (ρ drops ≤ 0.03), and the
absolute-scale conclusions are unchanged (log₂ +0.17 → +0.20). It lifts boomi too — a genuine
low-dynamic-range recovery, not a lauren crutch. A *naive* CI gate (relunc < 1, no R² floor)
over-rescues (geomCV → 0.64). Adopting the rescue as a Riana default is a scoped follow-up
(calibrate threshold + within-protein-CV guardrail); the strict R² > 0.8 numbers above remain the
conservative primary.

## Reproduce

```bash
# integrate (both), matched to boomi: iso0–5, ±10 ppm (SDRF), apex/0.15 min, --workers N
riana integrate <mzml_dir> <mzTab> --sdrf <sdrf> --out <run_dir> --workers 8

# fit — depth 6, label-aware Spep default (o18→6, D₂O→8), RIA from manifest
riana fit --manifest runs/timeseries_lauren5_7_ipsc_mesoderm_o18/riana_manifest.tsv \
  --label o18 --coefficients juber_2026_o18_ac16 --depth 6 -W 4 -o runs/timeseries_lauren5_7_ipsc_mesoderm_o18
riana fit --manifest runs/timeseries_lauren9/riana_manifest.tsv \
  --label hw  --coefficients alamillo_2025_ipsc --depth 6 -W 4 -o runs/timeseries_lauren9

# rollup — weighted / unique / linear-simple; lauren_5_7 gets mesoderm-vs-iPSC Δk
riana rollup --manifest runs/timeseries_lauren5_7_ipsc_mesoderm_o18/riana_manifest.tsv \
  --method weighted --model "linear simple" --reference-condition control --min-r2 0.8 \
  -o runs/timeseries_lauren5_7_ipsc_mesoderm_o18
riana rollup --manifest runs/timeseries_lauren9/riana_manifest.tsv \
  --method weighted --model "linear simple" --min-r2 0.8 -o runs/timeseries_lauren9

# evaluation (this report's numbers)
python runs/lauren_headtohead.py
```

Rollups: lauren_5_7 **2,079 proteins** (control + mesoderm, with Δk contrast),
lauren_9 **1,655 proteins**.
