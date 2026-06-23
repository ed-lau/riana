# Adaptive N_ISO & limited-isotopomer scoring — what earns its keep

- **Date:** 2026-06-23
- **Branch:** `m3-rewrite`
- **Status:** investigation complete — **fit-side narrow scoring (`--fs`) is the keeper; integrate-side adaptive capture (`--iso auto`) is neutral, kept opt-in**
- **Commits:** `92766d3` (B0-B2 adaptive integrate), `7d638e8` (B3 H4′ + B4 scoring), `ad9091f` (bench tooling), `2ecd920` (`--fs` productionized)
- **Benches:** `bench_within_protein_theta.py` (+`--score-channels`), `bench_fs_method_compare.py` (+`--score-channels`, `nansum`), `run_integrate_v1_0_0.py` (+`--adaptive`/`--ria`)
- **Roadmap:** `PROJECT_REVIEW.md` → Track B (adaptive N_ISO); memory `m8_adaptive_niso_robust_envelope`

## Question

The Track-B plan was to make N_ISO **per-peptide and automatic** by running the
IsoSpec forward model *at integrate* (capture each peptide's real envelope width
from the init∪final union ≥1% abundance), rather than a fixed `--iso 0 1 2 3 4 5`.
Two motivations:

1. **Accuracy** — capture wide so the labelled envelope isn't truncated; then fit
   on a clean subset ("integrate wide, fit narrow") to dodge co-eluting
   contaminants in the high isotopomers.
2. **UX / unblock** — the user shouldn't have to pick `--iso` at integrate; the
   isotopomer choice moves to *fit* time (`--fs`), where it belongs.

This session built B0-B4 and asked: **does adaptive capture actually improve the
science, or is the value entirely in the fit-side scoring choice?**

## TL;DR

- **The win is fit-side `--fs` limited-isotopomer scoring, not adaptive capture.**
  A 2×2 control on the ac16 calibration (capture {fixed, adaptive} × scoring {all,
  iso0-3}) shows iso0-3 scoring tightens the recovery core at every mixing
  proportion, while adaptive capture adds **nothing** (`fix·iso0-3 ≈ adapt·iso0-3`,
  fixed marginally better) — because scoring iso0-3 needs only iso0-3, which the
  standard fixed integrate already has.
- **`--fs` is productionized** (`riana fit --fs 0 1 2 3`) with two guards. It runs
  on existing integrate output, no re-integration.
- **`--iso auto` (B1) stays opt-in.** Neutral-to-slightly-negative on θ-spread (LVE),
  geometric k-CV (LVE), and recovery (calibration), and ~3-5× slower to integrate.
  Its remaining justification is UX (parameter-free integrate) + the TMT
  precursor-mass-from-`iso0` path (separate motivation), not N_ISO accuracy.
- **H4′ mix-then-normalize (B3) is correct and a near-no-op at low θ** — load-bearing
  only once you score a narrow subset (the steeper the truncation, the more the
  normalization order matters).
- **The RT↔scan offset is OpenMS alignment, not a `.raw`-vs-`.mzML` artifact** — the
  mzml re-search does **not** fix it (measured identical offsets). Tangential to
  N_ISO but it closes an open MBR question.

## What shipped

- **B0** — `IntegrationConfig.{ria_max, adaptive_iso, iso_abundance_floor, iso_max}`;
  `riana integrate --iso auto` + `--ria`; RIA read from SDRF `characteristics[precursor
  enrichment]`.
- **B1** — `algorithms.isotope_dist.adaptive_channel_masses`: per-peptide N_ISO from
  the init∪final envelope ≥1%, **Commerford upper-bound Spep** for the final-envelope
  width, **init (unlabelled) averaged-isotopolog mass** as the extraction target
  (iso0 == precursor m0; `iso{N}_ppm_error` = drift-from-unlabelled, a future
  orthogonal-θ substrate). NaN-padded ragged output to the run-wide max width.
- **B3** — `solve_fs_d2o` mixes init/final in the **full-cluster basis** → truncates
  to the scoring channels → renormalizes → compares (the H4′ order; the old
  normalize-each-then-mix is only correct when init/final share the in-window mass
  fraction).
- **B4 / `--fs`** — `FitConfig.score_channels`; `riana fit --fs 0 1 2 3` scores the FS
  on the leading subset. **Guard 1** (run-level): required integrate channels are
  conditional on `--fs`, so `integrate --iso 0 1 2 3` + `fit --fs 0 1 2 3` is valid.
  **Guard 2** (per-peptide): a peptidoform whose envelope ends before the subset is
  scored on what it has; a 1-channel peptidoform returns NaN (fixed an old shape-crash).

## Results

### 1. LVE/ATR within-protein A/B (in-vivo, no ground truth → precision proxy)

12,440 (protein, timepoint) cells, RIA 0.046, frozen proteotypic bench set. Metric =
within-protein **θ robust-SD** = median over cells of `1.4826·MAD(θ_peptides)`; θ =
fraction-new ∈ [0,1].

| stratum | fixed | adaptive |
|---|---|---|
| curated (kinetic R²>0.95) | 0.0469 | 0.0461 |
| all | 0.1007 | 0.1072 |
| uncurated | 0.1247 | 0.1330 |

Adaptive ≈ fixed on clean peptides, **worse on the noisy tail** (wide capture imports
high-channel noise on messy peptides). Curated yield also dips (5846→5565 peptides).

**Transferable cross-check — within-protein k-CV** (per-protein spread of fitted
turnover rate; the scale-free metric):

| population | fixed geo CV `1.4826·MAD(ln k)` | adaptive |
|---|---|---|
| curated R²>0.95 | 16.8% | 17.7% |
| R²>0.80 | 25.8% | 25.4% |
| all converged | 46.1% | 49.4% |

Same story as θ-spread, in k-space — so the conclusion is not a metric artifact.
(median k ≈ 0.06/day; LVE curated geo-CV ~17% sits a touch below the ~20-25% typical
for D₂O.) **Note on metrics:** θ robust-SD is in absolute θ units and is confounded by
where θ sits on the curve (k × sampling-time), so it is good for a *within-dataset* A/B
but not cross-dataset; the **geometric k-CV** is unit- and rate-invariant and is the
transferable one.

### 2. Calibration recovery — the decisive 2×2 (ac16, ground-truth mixing proportion)

`bench_fs_method_compare.py`, |θ−f| vs known mixing proportion, RIA 0.0598, ac16
coefficients. **Capture {fixed v1.0.0, adaptive} × scoring {all, iso0-3}.**

within ±0.05 (higher = better):

| f | fix·all | **fix·iso0-3** | adapt·all | adapt·iso0-3 |
|---|---|---|---|---|
| 0   | 19.4% | 22.1% | 18.0% | 21.9% |
| 0.5 | 20.8% | 23.8% | 19.8% | 23.6% |
| 1   | 25.5% | 27.3% | 25.2% | 27.3% |

IQR (lower = tighter core):

| f | fix·all | fix·iso0-3 | adapt·all | adapt·iso0-3 |
|---|---|---|---|---|
| 0 | 0.421 | 0.285 | 0.511 | 0.288 |
| 1 | 0.242 | 0.228 | 0.244 | 0.233 |

- **Scoring iso0-3 helps** (within ±0.05 +2-3 pp, IQR much tighter, bias down: f=1
  +0.025→+0.014). Consistent at every proportion.
- **Adaptive capture adds nothing** — `fix·iso0-3 ≈ adapt·iso0-3` (fixed marginally
  better); `adapt·all` is the worst cell.
- **Caveat — it's a core-tightening, not an RMSE win.** iso0-3 trades slightly fatter
  high-f tails (RMSE +~0.01 at f=1, where some peptides' label genuinely spreads past
  iso3) for a tighter, less-biased core. iso0-1 is too aggressive (loses signal at high
  f); **iso0-3 is the sweet spot.**

**cm replicates ac16** (drop50 config, RIA 0.0598): narrow scoring helps at every f
(within ±0.05 +1.9 pp; IQR f=0 0.300→0.218), and `fix·iso0-3 ≈ adapt·iso0-3` again —
capture adds nothing, `adapt·all` is worst (IQR f=0 0.383). The capture-vs-scoring
conclusion holds on **two** independent cell lines. (ipsc adaptive integrate was
pathologically slow — ~18 min/file on its denser profile mzML — so a third-line
confirmation is deferred; the split is structural and expected to hold.)

### 3. H4′ mix-then-normalize (B3) is a near-no-op at low θ

Re-scoring LVE with the corrected solver moved `fixed·all` 0.1007→0.1026 and left
adaptive unchanged — negligible, as predicted for RIA 4.6%. The correction is
load-bearing only when truncation is severe (narrow `--fs`, high θ), which is exactly
where B4 operates — so **B3 is the prerequisite for B4**, validated together.

### 4. RT↔scan offset is OpenMS alignment, not the search input

Median |mzML-scan-RT − mzTab-RT|, raw-search vs mzml-search mzTab on the *same* local
mzML (ac16):

| file | raw search | mzml search |
|---|---|---|
| 0% (low θ) | 2.16 min | 2.15 min |
| 100% (high θ) | 0.56 min | 0.58 min |

Identical → the offset is ProteomicsLFQ RT alignment, applied by both pipelines.
**Re-searching from mzML does not remove it.** The offset is larger on the *low*-θ file,
so it is not a high-θ effect either (per-proportion ID recovery is flat-to-rising with
proportion). Direct scan-based extraction is immune; only RT-anchored MBR transfers see
it. For the Percolator-path 0.9.0 calibration bench, none of this matters. See memory
`mbr_rt_axis_dependency` (TBD #1 now void).

### 5. Short-peptide stratification — the per-peptide signal is the *wide* tail

Hypothesis: adaptive helps **short** peptides by ending the envelope early vs a fixed
`--iso 0..5/6`. Tested on ac16, |fs−f| stratified by per-peptide N_ISO:

| N_ISO | n_pep | all-6 | iso0-3 | Δ |
|---|---|---|---|---|
| ≤5 (short) | ~1,400 (~5%) | 0.153 | 0.140 | −0.013 |
| 6-9 (bulk) | ~17,000 | 0.10-0.13 | 0.07-0.11 | −0.017 to −0.022 |
| ≥12 (wide) | ~400 | 0.10-0.25 | 0.20-0.33 | **+0.08** |

Three takeaways: (a) short peptides are **rare** at 6% D₂O (median N_ISO 7); (b) narrow
scoring helps **broadly** — long peptides benefit *more*, short peptides stay hardest, so
it is **not** a special short-peptide rescue; (c) the real per-peptide signal is the
opposite end — **global iso0-3 over-truncates the genuinely wide N_ISO≥12 peptides**
(+0.08, cutting real labelled signal). So the per-peptide envelope earns its keep as a
**scoring-width selector** (don't under-score the wide tail), not as a capture-accuracy
lever — see Future work.

## Verdict / roadmap impact

- **Keep & promote `--fs` (B3+B4).** Modest, free, fit-side; works on existing output.
  A default around iso0-3 is defensible for D₂O at 4.6-6% enrichment; document the
  core-vs-tail (IQR vs RMSE) tradeoff.
- **`--iso auto` (B1) stays opt-in / experimental.** No N_ISO-accuracy payoff at these
  enrichments and ~3-5× slower (451 s/file on 182 MB profile mzML vs faster fixed). Its
  live justifications are UX (parameter-free integrate) and the TMT precursor-mass path
  (`track_c_tmt_chemical_mod`), not envelope-width accuracy. Plausible accuracy payoff
  only at genuinely high enrichment where the labelled envelope truly spreads past iso5.

## Future work — per-peptide `--fs` (brainstorm, not built)

Currie et al. used a crude per-peptide rule keyed on D₂O labelling sites (Spep): <15
sites → iso0/iso1 ratio, 15-35 → iso0/iso2, >35 → iso1/iso3 — i.e. shift the scoring
window *up* as labelling capacity grows. Same direction as finding #5 (wide-envelope
peptides need higher/more channels). In the current whole-cluster-RMSE paradigm the
clean analog is a **per-peptide `score_channels = f(N_ISO or Spep)`**: e.g. iso0-3 for
short/medium, iso0-5 for very long peptides — capping where the contamination starts but
not under-scoring the genuinely wide ones. Open questions: (i) extend (keep iso0, add
high channels) vs shift (drop the suppressed iso0 at very high Spep, à la Currie); (ii)
key on N_ISO (already computed by `adaptive_channel_masses`) vs Spep; (iii) whether a
soft per-channel SNR/robust weight beats a hard per-peptide cutoff (the original "robust
matcher" research question). This is the natural next refinement of the flat global `--fs`.

## Reproduce

```bash
# calibration adaptive integrate (Percolator path; fixed v1.0.0 already exists)
python tests/benchmark/run_integrate_v1_0_0.py --line ac16 --adaptive --ria 0.0598 --out-label adaptive

# recovery 2x2 (repeat for v1.0.0 vs adaptive; --score-channels 4 = iso0-3)
python tests/benchmark/bench_fs_method_compare.py \
  --inputs tests/data/calibration_d2o_mixing/ac16/integrate_outputs/v1.0.0 \
  --ground-truth tests/data/calibration_d2o_mixing/ac16/ground_truth.csv \
  --coefficients tests/data/calibration_d2o_mixing/ac16/d2o_aa_coefficients_ac16.csv \
  --ria 0.0598 --score-channels 4 --output-dir runs/_cal_fix_iso03

# LVE within-protein theta + k-CV: integrate fixed vs --iso auto, then
#   bench_within_protein_theta.py --method ... [--score-channels N]
#   riana fit --manifest <run>/riana_manifest.tsv --coefficients commerford -r 0.046
# RT offset: parse PSM spectra_ref scan + retention_time, reconcile via IndexedMzML.rt_for_scans
# (calibration mzML/PSMs are gitignored; regenerate via quantms / snakemake)
```
