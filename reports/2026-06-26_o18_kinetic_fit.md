# ¹⁸O kinetic fit — time-series validation, D₂O head-to-head, R-script reconciliation

- **Date:** 2026-06-26
- **Branch:** `1.1.0`
- **Status:** **Investigation — kinetic `fit --label o18` validated; no reverse-model change required.** Curation + estimator recommendations recorded; one deferred refinement (t0 anchoring).
- **Inputs:** new AC16 + iPSC D₂O-vs-¹⁸O time series (`runs/{juber_ac16,boomi_ipsc}_{o18,d2o}`, gitignored); the legacy R analysis `data/timeseries_juber_ac16_o18/00_Riana_RMD_18` (RIANA v0.8.3 + Percolator); comparison driver `runs/compare_kinetics.py`.
- **Roadmap:** v1.1.0 item #3 — the deferred **kinetic** o18 fit (the reverse model shipped 2026-06-24; this is the time-series turn). See [[o18_reverse_model_design]].

> **Update (2026-06-27) — re-verified after the bootstrap-OOB coefficient refit.** The
> `juber_2026_o18_ac16` table was re-frozen by bootstrap (see the 2026-06-24 report's update)
> and the iPSC/AC16 fits below were re-run. The numbers are **stable**: iPSC FS(24 h) 0.66,
> Spearman(k) 0.66 (peptide & protein), within-protein geomCV 0.149 / 0.146, yield 43 / 33 %;
> AC16 essentially unchanged (ρ 0.65, yield 3.0 / 4.6 %). All curation decisions (R² ≥ 0.8,
> Sₚₑₚ ≥ 5, depth 6) and conclusions stand.

## Question

The reverse model + mixing-calibration fit shipped 2026-06-24, but the **kinetic** `fit --label o18` had never run (no time series). With AC16 + iPSC D₂O-vs-¹⁸O series in hand: does it recover turnover? A first AC16 pass looked alarming (median FS(8h) ≈0.07, 3% R²≥0.8), while the user's legacy R script — on the same AC16 data and a *known-wrong* FS calc — recovered a clean accumulation to ~0.29. Is the ¹⁸O reverse model missing something?

## TL;DR

- **The reverse model is sound.** **iPSC ¹⁸O** (RIA 0.0897, 24 h) recovers turnover cleanly: **36% of peptides R²≥0.8**, monotonic FS 0→**0.67** by 24 h, t½≈14 h. ¹⁸O works wherever there is adequate signal.
- **AC16's weakness is the regime, not the label.** Over 8 h at 5.8% RIA, **D₂O struggles identically** (4.1% R²≥0.8 vs ¹⁸O 3.0%; D₂O curated FS only 0.23 by 8 h). The short window caps *both*.
- **¹⁸O ≈ D₂O on the curated set, and *out-curates* it at iPSC.** Head-to-head (depth 6, curated): iPSC Spearman(k) **0.66 at both peptide and protein level** (the linear-φ rollup with unique parsimony does *not* degrade it); ¹⁸O reaches 43% R²≥0.8 vs D₂O's 33%; within-protein robust geomCV ~0.15, identical between labels (range of the best D₂O tables, [2026-06-25](2026-06-25_d2o_coefficient_tables.md)). The two **agree on ranking** but differ on **absolute scale** (iPSC ~1.5×, ¹⁸O faster) — expected (AC16 ¹⁸O table applied to iPSC; no iPSC ¹⁸O calibration). AC16 stays 8 h-window-limited (protein-level ρ drops to 0.47). Stable across depth 3/4/6.
- **The R-script's 0.29 vs riana's ~0.22 is two parts selection + one part deflation.** The R violin is its R²≥0.8 *kept* subset (selection); on the *same* clean peptides riana gives 0.19–0.22. The residual ~1.3× deflation is a **too-hot reference**: riana's calibration-trained Spep (incl. serine, which the R omits) predicts more iso2 per unit FS than the in-vitro kinetic reality.
- **Clamping FS to [0,1] is neutral** (−1.0 / +0.3 / +1.2 pt across AC16 o18 / AC16 d2o / iPSC o18) — not the lever. The real t0 lever is **anchoring**: the production **simple-exponential** model (FS=1−e^(−k·t)) forces FS(0)=0, but the reference-offset observed FS(0)≈−0.08 is an irreducible residual that tanks the *exponential* R² (3% AC16). A diagnostic `lm(FS~t)` with a **free intercept** absorbs that offset (→12–14%) — that gap is the free baseline, **not** clamping (and a free intercept is the wrong fix: it spends a parameter on noise). The right fix anchors the *reference* (observed-t0 init) with the intercept still pinned — **A/B: +3 pt yield, lower within-protein CV**. All head-to-head numbers below use the **peptide-level simple-exponential** fit (no rollup); the linear figures are diagnostic only.

*(k is per **hour** — the labeling_time axis is hours. AC16 t½≈20 h, iPSC t½≈14 h — both biologically sensible for ~16 h doubling.)*

## Data

| dataset | RIA (SDRF) | timepoints | window | bioreps |
|---|---|---|---|---|
| AC16 ¹⁸O / D₂O | 0.0583 / 0.0598 | 8 | 0→8 h | 1 |
| iPSC ¹⁸O / D₂O | 0.0897 / 0.0598 | 9 | 0→24 h | 2 |

Integrated at the default iso0–5 grid, RIA auto-resolved per-experiment from the manifest (a fit-side fix this session — `fit_project` now reads `precursor_enrichment`; iPSC needs 0.0897, not the 0.06 default).

## iPSC ¹⁸O — the model works

Production fit (`--label o18 --coefficients juber_2026_o18_ac16`, full envelope), median FS per timepoint:

| t (h) | 0 | 2 | 4 | 6 | 8 | 12 | 24 |
|---|---|---|---|---|---|---|---|
| curated FS | 0.00 | 0.08 | 0.17 | 0.27 | 0.34 | 0.45 | **0.67** |

36.1% R²≥0.8 (20% R²≥0.9), curated median k≈0.049/h (t½≈14 h), Spep median 7.1.

## ¹⁸O-vs-D₂O head-to-head (curated R²≥0.8 + Spep≥5, **depth 6**)

Analytical-grade at `--depth 6` (peptidoforms seen at ≥6 distinct timepoints — the cleanest set; depth note below). Peptide-level = simple exponential joined on `concat`; protein-level = linear **φ=log(1−θ)** rollup, **unique parsimony**, `--min-r2 0.8`:

| metric | AC16 (8 h) | iPSC (24 h) |
|---|---|---|
| %R²≥0.8: ¹⁸O / D₂O | 3.3 / 4.6% | **42.9 / 33.0%** |
| Spearman(k) — **peptide**-level | 0.65 (n=320) | **0.66** (n=3479) |
| Spearman(k) — **protein**-level rollup | 0.47 (n=368) | **0.66** (n=1067) |
| within-protein robust geomCV — ¹⁸O / D₂O | 0.158 / 0.169 | **0.150 / 0.146** |
| median log2(k_D₂O / k_¹⁸O) | +0.14 | **−0.54** (¹⁸O ~1.5× faster) |

iPSC is analysis-grade: ¹⁸O *out-curates* D₂O (43% vs 33%), the two **agree on ranking at both peptide and protein level (ρ=0.66/0.66 — the rollup does not degrade it)**, and within-protein consistency is identical (~0.15, the range of the best D₂O tables, [2026-06-25](2026-06-25_d2o_coefficient_tables.md)). AC16 is 8 h-window-limited for *both* labels; its protein-level concordance drops (0.47) as the rollup pools few, noisy curated peptides. The two differ on **absolute scale ~1.5×** (¹⁸O faster, t½ 14 h vs D₂O 20 h at iPSC) — expected: the AC16 ¹⁸O table is applied to iPSC (no iPSC ¹⁸O calibration exists; D₂O uses the iPSC-matched `alamillo_2025_ipsc`), plus genuine label-chemistry differences. **Ranking + within-protein consistency are the trustworthy cross-label readouts; absolute k is not expected to match a priori.**

**Depth robustness.** Stable across `--depth` 3/4/6: peptide Spearman ~0.65 (AC16) / 0.66 (iPSC), bias stable, ¹⁸O ≈ D₂O geomCV at every depth. Deeper coverage *raises* yield (iPSC 36→43%, since ≥6-timepoint peptides are the abundant well-measured ones — sparse peptides *dilute* yield, not inflate it) and tightens within-protein geomCV ~10% (sparse 3-point fits carry unreliable k). So **depth 6 is the analytical-grade default** — it sharpens the metrics without changing a single conclusion.

## Curation-gate sensitivity — a CI (dk) rescue

R² is a *goodness-of-fit* statistic that collapses when a curve is flat (low k or low dynamic
range) even for a well-measured peptide — the issue Lau et al. *Nat Commun* 2018 (and Sadygov's
d2ome) address by gating on the **rate constant's confidence interval** instead. The legacy RMD
optimizer this pipeline descends from (`optim_util.R::fitRiana`) computes exactly that — an
analytical `dk` (SE of k from dA/dk); Riana already emits the bootstrap analogue
(`sd`/`ci_lo`/`ci_hi`, 5–95th pct). Second gate =
**R² > 0.8 OR (R² > 0.6 AND `relunc` < 0.25)**, `relunc = (ci_hi − ci_lo)/(2k)` (90 % CI
half-width over k, a coefficient of variation of k̂). Effect on the head-to-head
(`runs/gate_comparison.py`; strict → rescue):

| dataset | yield ¹⁸O / D₂O | within-prot geomCV ¹⁸O / D₂O | peptide ρ | protein ρ |
|---|---|---|---|---|
| iPSC (boomi) | 28.5→35.5 / 30.9→42.5 % | 0.140→0.180 / 0.139→0.184 | 0.679→0.656 | 0.661→0.607 |
| AC16 (juber) | 5.0→6.4 / 6.0→7.8 % | 0.130→0.175 / 0.136→0.164 | 0.673→0.661 | 0.606→0.620 |

The rescue lifts yield everywhere (~1.2–1.4× iPSC, and ~2× on the sampling-starved lauren SCVI480
set — [2026-07-01](2026-07-01_lauren_yield_gap_diagnosis.md)) while the within-protein CV stays
clean (< 0.25) and the ¹⁸O↔D₂O ranking + scale conclusions are unchanged (peptide ρ drops ≤ 0.03;
AC16 protein ρ even ticks up as more curated peptides stabilize the rollup). A *naive* CI gate
(relunc < 1, no R² floor) over-rescues (geomCV → 0.6). This confirms the ¹⁸O curation is
R²-metric-limited, not estimator-limited; productionizing the CI rescue (calibrated threshold +
CV guardrail) is a scoped curation upgrade that benefits every dataset. All headline numbers
keep the conservative strict R² > 0.8 as primary.

## R-script reconciliation

The R method (`18o_fs_violin.Rmd`): **A_t = iso0/iso2** (the +2 channel single-¹⁸O lands in), mapped onto a linear interpolation between simulated natural and fully-labeled ratios, **clamped [0,1]**, curated to R²≥0.8. The "H4′-wrong" part is interpolating the *ratio* (concavity) rather than mixing intensities then taking the ratio. Reconciliation on riana's current integrate:

- The R violin's **0.29 @ 8 h is its R²≥0.8 *kept* subset** (2066 peptides, selection). The R method's *all-peptides* median is ~0.14; riana-full ~0.09.
- On the **R's exact kept peptides**, riana recovers (clamped): full **0.19**, iso0-3 **0.20**, iso0/iso2-ratio **0.22** — vs R's 0.29. So ~70–75%; the gap is the too-hot reference, not a broken method.

## Three knobs evaluated

**(1) Clamp FS to [0,1] — neutral; keep the permissive diagnostic bounds.** Linear-R² curation yield, unclamped (−0.1,1.2) vs clamped:

| dataset | unclamped | clamped[0,1] | Δ |
|---|---|---|---|
| AC16 ¹⁸O | 12.1% | 11.1% | −1.0 |
| AC16 D₂O | 16.8% | 17.1% | +0.3 |
| iPSC ¹⁸O | 50.9% | 52.1% | +1.2 |

Clamping does not move yield — so it need **not** apply uniformly, and the diagnostic bounds (load-bearing for the calibration NBs' f=0/1 reconstruction check) can stay. The genuine t0 issue is the production simple-exponential forcing FS(0)=0 against a reference-offset observed FS(0)≈−0.08. The fix is **t0 anchoring of the *reference*** — use each peptide's observed t0 envelope as its init (like the empirical-t0 anchoring shipped for `fs_ds`, [[mass_defect_theta_design]]) — **not** a free intercept (a free intercept spends a 2nd parameter on noise and would loosen well-behaved peptides; FS<0 is a non-physical artifact). Anchoring keeps the model 1-parameter / pinned at the origin and fixes the *data*. **A/B (AC16 ¹⁸O, Spep≥5, fixed-intercept exponential): observed-t0 anchor vs theoretical init → R²≥0.8 yield 6.6%→9.6%, within-protein geomCV 0.105→0.037** (CV on small n=7–9 proteins; yield is the robust signal). Worth productionizing for the o18 kinetic path after a full-set / iPSC confirmation.

**(2) Spep cutoff.** Spep distribution (from coefficients): AC16 mean 5.2 / median 4.6 / 90th 9.1; iPSC mean 7.8. AC16 yield + FS(8h) by Spep bin:

| Spep bin | n | %R²≥0.8 | medFS8h |
|---|---|---|---|
| [0,3) | 8757 | 0.7% | 0.00 |
| [3,5) | 11382 | 2.7% | 0.073 |
| [5,7) | 8295 | 3.7% | 0.092 |
| [7,10) | 5617 | 4.5% | 0.113 |
| [10,∞) | 2539 | 6.5% | 0.123 |

[0,3) is dead; monotonic above. **Curate at Spep≥5** (captures 45% of AC16, the median-and-above), per the distribution rather than a raw D/E/N/S count.

**(3) Scoring channels (AC16 ¹⁸O, Spep≥5, clamped) — a bias/variance tradeoff.**

| channels | %R²≥0.8 (yield) | medFS8h (magnitude) |
|---|---|---|
| iso0-2 | 10.3% | **0.170** (least deflated) |
| iso0-3 | 15.3% | 0.150 |
| full (iso0-5) | **18.9%** | 0.134 (most deflated) |

Fewer channels recover more FS magnitude (less pull toward the FS-insensitive iso0/iso1) but are noisier (lower yield); the full envelope is the best *yield* (production default) but most deflated. iso0-3 ≈ middle; **iso0-3 is not a free win over full** — the H4′ normalize-truncate-renormalize already neutralizes the empty iso4/5.

## Recommendations

1. **No reverse-model change** — iPSC validates it; ship the kinetic fit as-is.
2. **Curate metrics** at **R²≥0.8 + Spep≥5**; report within-protein spread as robust geometric CV (`1.4826·MAD(ln k)`, [2026-06-23](2026-06-23_robust_turnover_cv.md)). Accept that low-intensity / ambiguous-proteoform peptides drop out.
3. **Do not blanket-clamp** FS (neutral) and **do not free the intercept** (loosens fits, non-physical). The real low-SNR refinement is **t0 anchoring of the reference** (observed-t0 init, intercept stays pinned) — A/B-positive (+3pt yield); productionize for o18 after a full-set / iPSC confirmation.
4. **Keep the full envelope** as the production scoring default (best yield); the deflation is a known steady-state-calibration-vs-kinetic-regime gap, not a bug. Optionally expose fewer-channel scoring for FS-accuracy-prioritized analyses.
5. **iPSC is the analysis-grade dataset**; AC16's 8 h window limits both labels — use it for method work, not for turnover conclusions.

## Reproduce

```bash
# integrate (per dataset; RIA auto from SDRF, iso0-5)
riana integrate data/timeseries_juber_ac16_o18/mzml <quantms mzTab> \
  --sdrf <sdrf> -W 9 -o runs/juber_ac16_o18
# kinetic fit (RIA auto-resolved per-experiment from the manifest)
riana fit --manifest runs/juber_ac16_o18/riana_manifest.tsv \
  --label o18 --coefficients juber_2026_o18_ac16 -W 4
# head-to-head (curated)
python runs/compare_kinetics.py \
  runs/juber_ac16_o18/riana_fit_peptides.txt runs/juber_ac16_d2o/riana_fit_peptides.txt o18 d2o
```
