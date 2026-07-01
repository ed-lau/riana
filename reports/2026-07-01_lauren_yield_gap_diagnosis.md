# Lauren ¹⁸O/D₂O — diagnosing the curation-yield gap vs boomi

- **Status:** **Diagnosed.** The 3–4× lower R²>0.8 curation yield in the lauren
  (SCVI480) sets vs boomi (AICS52) is **predominantly a sampling-design × metric
  effect, compounded by a ~1.5× higher per-point noise floor — not a failed MS run and
  not slower turnover.** Identification is *higher* in lauren; the loss is entirely at
  the quantitative-fit gate.
- **Parent report:** [`2026-07-01_lauren_o18_d2o_headtohead.md`](2026-07-01_lauren_o18_d2o_headtohead.md)
  (the head-to-head that flagged yield as the open item).
- **Data source:** ¹⁸O run = `runs/timeseries_lauren5_7_ipsc_mesoderm_o18` (Sage-free quantms
  re-search; mzML symlinked from network drive) — numbers unchanged vs the prior integration to
  within rounding.
- **Reproduce:** `python runs/lauren_yield_diagnostics.py`

## TL;DR

The user's puzzle — lauren identifies **more** proteins (~7,000 vs ~5,000) yet curates
**3–4× fewer** for turnover — resolves cleanly: curation is a *quantitative signal-to-noise*
gate, not an identification one. Four candidates tested; verdict for each:

| candidate | verdict | evidence |
|---|---|---|
| **Slower turnover / less proliferation** | **NOT the driver** | k distributions comparable: lauren ¹⁸O med k 0.037 vs boomi 0.045; lauren D₂O 0.039 *faster* than boomi 0.032 |
| **Poor MS run / lost IDs / low depth** | **NOT the driver** | lauren has *more* PSMs, *more* proteins, and *more* timepoints (median n_points 18 vs 15). Yield doesn't rise with depth in lauren (R²>0.8 ~3–5 % at every depth bin) |
| **Sampling design (front-loaded grid) × R² metric** | **DOMINANT** | lauren puts **11 of 12 timepoints in 0→7 h** where FS only reaches 0.23 (signal ≈ noise), then a lone 24 h point. R² is pathologically low when dynamic range ≈ scatter |
| **Higher per-point noise (data quality)** | **SECONDARY (~1.5×)** | t0 FS scatter (pure noise; true FS=0) is **0.13–0.14 in lauren vs 0.09 in boomi** |

## 1. It is not turnover rate

Physically-valid fits (rail-hits dropped), median k_deg [IQR]:

| arm | med k | note |
|---|---|---|
| lauren ¹⁸O | 0.037 [0.021, 0.076] | slightly below boomi ¹⁸O |
| lauren D₂O | 0.039 [0.025, 0.071] | **above** boomi D₂O |
| boomi ¹⁸O | 0.045 [0.031, 0.065] | — |
| boomi D₂O | 0.032 [0.021, 0.047] | — |

Turnover regimes overlap; the "less-proliferation → lower k → fewer curated" hypothesis is
not supported by the data. (SCVI480 is not systematically slower.)

## 2. The R² distribution is wholesale-shifted, not bimodal

Spep-gated, no R²/k filter — median R² [IQR] and pass rates:

| arm | med R² | >0.5 | >0.8 | >0.9 |
|---|---|---|---|---|
| lauren ¹⁸O | **0.26** [−0.09, 0.61] | 33 % | 10 % | 3 % |
| boomi ¹⁸O | **0.77** [0.46, 0.91] | 73 % | 46 % | 27 % |
| lauren D₂O | 0.36 | 40 % | 15 % | 6 % |
| boomi D₂O | 0.69 | 67 % | 36 % | 18 % |

The entire lauren R² distribution sits ~0.4–0.5 lower. It is a systematic reduction in
per-fit goodness, not a subpopulation of failures.

## 3. The mechanism — sampling design meets the flat-curve R² pathology

Median FS per timepoint (control, Spep-gated):

```
lauren ¹⁸O:  0h -.01  0.5h -.01  1h .01  1.5h .02  2h .03  2.5h .04  3h .08  4h .12  5h .15  6h .19  7h .23   ···   24h .54
boomi ¹⁸O:   0h -.04         1h -.01  2h .06        3h .11  4h .15        6h .26  8h .33        12h .43         24h .65
```

- **lauren spends 11 of 12 timepoints in FS 0→0.23** — the low-dynamic-range shoulder — then
  jumps to 0.54 at 24 h with **no sampling between 7 h and 24 h**, exactly where the rising
  exponential is most informative.
- **boomi spreads 9 timepoints across FS 0→0.65**, catching the climb at 6/8/12 h.

R² ≈ 1 − (noise²/signal-variance). When most points sit where signal ≈ noise, R² is low
*even for a perfectly good fit*. This is the classic flat-curve pathology, visible as an
R²-vs-k arch in **both** datasets (R² peaks at mid-k 0.02–0.08 and collapses for slow and
fast peptides) — but lauren's front-loaded grid parks the bulk of its measurements in the
low-SNR region, so far more peptides fall under the gate.

At matched k the effect still shows lauren below boomi (e.g. k∈(0.04,0.08]: R²>0.8 pass rate
14 % lauren vs 69 % boomi), which is the second factor:

## 4. The secondary factor — a ~1.5× higher noise floor

t0 FS scatter is a pure measurement-noise readout (true FS(0)=0 for every peptide):

| arm | MAD(FS) at t0 |
|---|---|
| lauren ¹⁸O / D₂O | **0.134 / 0.140** |
| boomi ¹⁸O / D₂O | **0.088 / 0.095** |

lauren's per-point noise is ~1.5× boomi's — consistent with the crash-recovered raw and/or
added variance from summing 8 LC fractions (a path boomi/juber never exercised). This raises
the noise term in R² across every timepoint and compounds §3. It lowers ¹⁸O's 24 h plateau
too (0.54 vs 0.65), shrinking dynamic range further (partly the lower RIA 0.079 vs 0.090;
D₂O, at identical RIA, reaches 0.58 ≥ boomi's 0.54, so plateau is not the main lever).

## 5. Why the D₂O↔¹⁸O comparison has "few shared peptides"

This is arithmetic, not a separate MS failure: the shared set is the *intersection* of two
independently-curated arms. At lauren's ~10–15 % per-arm yield the intersection is ~1–2 %
(n=742), vs boomi's ~35–45 % giving ~9 % (n=3,002). The low n follows directly from §3–4;
it does not indicate the peptides disagree (ranking ρ is still 0.67).

## 6. The principled fix — a dk / CI rescue (Lau 2018), and its limits

Lau et al. *Nat Commun* 2018 (and later Sadygov's d2ome) recognized that **R² is the wrong
sole gate when the curve is flat**: a slow-but-well-measured peptide has low R² yet a *tight
rate-constant confidence interval*. The RMD R optimizer this pipeline descends from
(`optim_util.R::fitRiana`) computes exactly this — an analytical `dk` (SE of k from
dA/dk_deg). Riana already emits the bootstrap analogue (`sd`, `ci_lo`, `ci_hi`).

Tested a secondary CI gate `rel_unc = (ci_hi−ci_lo)/(2k)` with a within-protein-CV guardrail:

| rule | lauren ¹⁸O yield | within-prot geomCV | boomi ¹⁸O yield |
|---|---|---|---|
| R²>0.8 (current) | 9.3 % | 0.159 | 29.5 % |
| OR(R²>0.6, relunc<0.20) | 14.1 % | 0.195 | 32.3 % |
| OR(R²>0.5, relunc<0.25) | 20.6 % | 0.218 | 37.3 % |
| relunc<1.0 (naive, no floor) | 66 % | **0.64** ✗ | 58 % |

- A **naive** CI gate over-rescues (within-protein CV blows to 0.64 — 4× worse); it pulls in
  genuinely imprecise fits.
- A **calibrated** gate (modest R² floor + tight relunc) roughly **doubles** lauren yield
  while holding within-protein CV near the strict value, and it also lifts boomi — so it
  recovers real low-dynamic-range peptides, not a lauren-only crutch. The D₂O↔¹⁸O peptide
  ranking is preserved (ρ 0.670 → 0.651 at 2.6× more shared peptides).

**Recommendation:** the yield gap is explained and is *not* a data-integrity problem, so no
re-run is needed for the note. The dk/CI rescue is the right long-term curation upgrade
(literature-anchored, benefits every dataset), but it needs a calibrated threshold and a
within-protein-CV guardrail before it becomes a default — a scoped follow-up, not a quick
switch. For the note, the honest framing is: **lauren reproduces the trustworthy readouts
(ranking, within-protein CV) on a second line; its lower curation count is a known
consequence of a front-loaded sampling grid + the R² gate's flat-curve pathology + ~1.5×
noise, all of which drop n without biasing the surviving numbers.**

## Appendix — sampling-design lesson for future ¹⁸O/D₂O runs

For a 0→24 h turnover experiment, spend timepoints on the **rising shoulder** (≈4–16 h),
not the near-flat 0→3 h region. lauren's 0.5 h-spaced early grid buys little (FS < 0.08
through 3 h) at the cost of the 7→24 h gap that most constrains k. boomi's [0,1,2,3,4,6,8,12,24]
is close to log-spaced and curates far better on identical chemistry.
