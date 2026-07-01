# Lauren ¹⁸O/D₂O — second iPSC head-to-head + mesoderm ¹⁸O series

- **Status:** **Integrated + fit + rolled up.** Two new analysis-grade sets harden §5.2
  of the ¹⁸O technical note: a same-cells (SCVI480 iPSC) D₂O↔¹⁸O head-to-head and a
  mesoderm ¹⁸O series. **The trustworthy §5.2 readouts reproduce on a second iPSC line**
  (ranking ρ ≈ 0.67, within-protein geom CV ≈ 0.14); **curation yield is 3–4× lower than
  boomi** — the one open item, deferred to a dedicated session.
- **Sibling reports:** [`2026-06-26_o18_kinetic_fit.md`](2026-06-26_o18_kinetic_fit.md)
  (boomi/juber head-to-head), [`2026-06-28_spep_curation_gate.md`](2026-06-28_spep_curation_gate.md).

## TL;DR

- **The method reproduces across iPSC lines.** SCVI480 (lauren) D₂O-vs-¹⁸O ranking
  ρ = **0.67** peptide / **0.63** protein and within-protein geom CV **≈ 0.14** land on
  top of the boomi AICS52 head-to-head (ρ 0.68/0.66, CV 0.14). These are the readouts
  §5.2 relies on — both hold on independent cells.
- **D₂O absolute scale is reproducible, ¹⁸O is not** — as §5.2 argued. Cross-line
  per-protein: D₂O ρ = 0.60, log₂(boomi/lauren) = −0.08 (parity); ¹⁸O ρ = 0.52,
  log₂ = +0.55 (1.5× shift). The ¹⁸O table is the transferred AC16 one; D₂O uses the
  iPSC-matched `alamillo_2025_ipsc`.
- **Open: yield.** lauren R²>0.8 yield is **¹⁸O 9.4% / D₂O 14.1%** vs boomi 42.9 / 33.2%,
  despite identical CV and ρ. Likely data-quality / regime (lower ¹⁸O RIA, crash-recovered
  raw, multi-fraction), not the estimator. Chased next session.
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
| **lauren5_7 iPSC ¹⁸O** | 16,681 | **9.4 %** | 1,394 | **0.159** (140) |
| lauren5_7 mesoderm ¹⁸O | 9,705 | 15.2 % | 1,311 | 0.145 (146) |
| **lauren9 iPSC D₂O** | 31,391 | **14.1 %** | 4,328 | **0.137** (522) |
| boomi iPSC ¹⁸O — ref | 20,021 | 42.9 % | 5,704 | 0.140 (650) |
| boomi iPSC D₂O — ref | 18,448 | 33.2 % | 5,708 | 0.139 (620) |

Within-protein geom CV is **flat at ~0.14–0.16 across every arm and both labels** — the
same value boomi and the §5.2 juber/iPSC sets show. Yield is the outlier (below).

## Head-to-head — same cells (SCVI480 iPSC), D₂O vs ¹⁸O

| metric | lauren SCVI480 | boomi AICS52 (§5.2) |
|---|---|---|
| Spearman(k) — peptide | **0.670** (n=742) | 0.68 (n=3002) |
| Spearman(k) — protein (median-k) | **0.625** (n=661) | 0.66 (n=1267) |
| within-prot geom CV: ¹⁸O / D₂O | 0.159 / 0.137 | 0.140 / 0.139 |
| median log₂(k_D₂O / k_¹⁸O) | **+0.17** (≈ parity) | −0.46 (¹⁸O ≈ 1.4× faster) |

Ranking and within-protein CV **reproduce**; the absolute D₂O/¹⁸O scale does not
(+0.17 vs −0.46) — expected, since ¹⁸O uses the transferred AC16 table and the scale is
line-dependent (§5.2). n is smaller here because of the lower yield, not selection.

## Cross-dataset reproducibility (per-protein, same label, SCVI480 vs AICS52)

| pair | peptide ρ | protein ρ | median log₂(boomi/lauren) |
|---|---|---|---|
| lauren ¹⁸O vs boomi ¹⁸O | 0.416 (n=350) | 0.517 (n=605) | **+0.55** (1.5× scale shift) |
| lauren D₂O vs boomi D₂O | 0.584 (n=890) | **0.604** (n=1081) | −0.08 (parity) |

The D₂O per-protein rates agree across two independent iPSC lines in both **ranking and
scale**; ¹⁸O agrees on ranking but carries a 1.5× scale offset from the transferred table.
This is the cleanest statement yet of "D₂O matched-table = quantitative; ¹⁸O transferred =
rank-only" — and it holds between lines, not just within one.

## The open item — yield

lauren yield (9–15 %) sits ~3–4× below boomi (33–43 %) while CV and ρ are unchanged. Since
low yield leaves the trustworthy readouts intact (it dropped n, not the numbers), it is a
data-quality / regime question, not an estimator one. Leading candidates, next session:

1. **Lower ¹⁸O enrichment** (RIA 0.079 vs boomi's 0.090) — a smaller envelope shift per
   timepoint lowers per-peptide R² (D₂O RIA is identical to boomi, yet D₂O yield is also
   down, so this is at most partial).
2. **Crash-recovered raw quality** — the empty scans were one visible symptom; broader
   per-spectrum noise would depress R² across the board (fits both labels being down).
3. **Multi-fraction collapse** — summing intensities across 8 fractions could add
   cross-fraction variance; boomi/juber never exercised this path.

First diagnostics: per-timepoint curated-FS accumulation and the R²-distribution shape,
lauren vs boomi, per label; and a fraction-collapse-off spot check.

## Reproduce

```bash
# integrate (both), matched to boomi: iso0–5, ±10 ppm (SDRF), apex/0.15 min, --workers N
riana integrate <mzml_dir> <mzTab> --sdrf <sdrf> --out <run_dir> --workers 8

# fit — depth 6, label-aware Spep default (o18→6, D₂O→8), RIA from manifest
riana fit --manifest runs/lauren5_7_final_20260630/riana_manifest.tsv \
  --label o18 --coefficients juber_2026_o18_ac16 --depth 6 -W 4 -o runs/lauren5_7_final_20260630
riana fit --manifest runs/timeseries_lauren9/riana_manifest.tsv \
  --label hw  --coefficients alamillo_2025_ipsc --depth 6 -W 4 -o runs/timeseries_lauren9

# rollup — weighted / unique / linear-simple; lauren_5_7 gets mesoderm-vs-iPSC Δk
riana rollup --manifest runs/lauren5_7_final_20260630/riana_manifest.tsv \
  --method weighted --model "linear simple" --reference-condition control --min-r2 0.8 \
  -o runs/lauren5_7_final_20260630
riana rollup --manifest runs/timeseries_lauren9/riana_manifest.tsv \
  --method weighted --model "linear simple" --min-r2 0.8 -o runs/timeseries_lauren9

# evaluation (this report's numbers)
python runs/lauren_headtohead.py
```

Rollups: lauren_5_7 **2,102 proteins** (control + mesoderm, with Δk contrast),
lauren_9 **1,655 proteins**.
