# Why the Juber AC16 series curates ~10× worse than the Boomi iPSC series

- **Date:** 2026-06-28
- **Branch:** `1.1.0`
- **Status:** **Diagnosed.** AC16's 3–4 % R²>0.8 curation (vs iPSC's 33–43 %) is **not** a
  data-quality or pipeline defect — it is the **8 h labelling window** (primary) compounded by
  AC16's **lower label accumulation** (slower effective turnover; for ¹⁸O, a lower RIA). The
  two factors multiply to the ~7–13× gap. ID depth, intensity, and timepoint ordering are ruled out.
- **Inputs:** `runs/{juber_ac16,boomi_ipsc}_{o18,d2o}` (uniform current-default integration);
  a window-subset experiment fitting the iPSC integrate restricted to its ≤ 8 h timepoints.
- **Curation throughout:** R² > 0.8 + label-aware Spep (¹⁸O ≥ 6 / D₂O ≥ 8) + **depth 6**
  (`reports/2026-06-28_spep_curation_gate.md`).

## Question

The Juber AC16 ¹⁸O/D₂O series curates at **3–4 %** (R² > 0.8) while the Boomi iPSC series
curates at **33–43 %**, despite the ¹⁸O coefficients being *trained* on AC16 and both
datasets having similar (in fact AC16 has greater) protein-ID depth. Why?

## TL;DR

- **AC16 has *more* IDs, not fewer** — 37 k peptidoforms / 3 740 proteins vs iPSC's 18–20 k /
  3 000–3 200. So low intensity / shallow ID (a natural first guess) is **ruled out**.
- **The FS curve is monotonic in every file** — no swapped/mislabelled timepoints.
- **AC16 accumulates far less measurable label by 8 h.** Population median FS at 8 h:
  **AC16 0.14 (D₂O) / 0.07 (¹⁸O) vs iPSC 0.22 / 0.33.** AC16's median fit R² is **negative**
  at every depth — most curves are noise around a tiny rise.
- **Window-subset experiment is decisive.** Truncating the *iPSC* integrate to its ≤ 8 h
  timepoints and re-fitting drops its curation from **33→12 % (D₂O)** and **43→24 % (¹⁸O)** —
  the **window alone costs ~½–⅔ of the yield.** But iPSC-≤8h still beats AC16 **2.7× (D₂O,
  same RIA) to 7× (¹⁸O)** — the **residual is AC16's lower label accumulation** (slower
  effective turnover; ¹⁸O additionally penalised by its lower RIA 0.058 vs 0.090).
- **Net:** AC16's deficit ≈ **window (~2.7×) × label-accumulation (~2.7× D₂O / ~7× ¹⁸O)**.
  Neither factor alone; they compound. The fix is experimental (longer window), not analytical.

## The datasets (uniform current-default integration)

| dataset | RIA | timepoints (h) | bioreps | peptidoforms | proteins | median R² | %R²>0.8 |
|---|---|---|---|---|---|---|---|
| AC16 ¹⁸O | 0.058 | 0,0.5,1,1.5,2,3,6,8 | 1 | 36,590 | 3,740 | **−0.28** | 3.0 |
| AC16 D₂O | 0.060 | 0,0.5,1,1.5,2,3,6,8 | 1 | 37,252 | 3,738 | **−0.09** | 4.2 |
| iPSC ¹⁸O | 0.090 | 0,1,2,3,4,6,8,12,24 | 2 | 20,021 | 3,214 | 0.749 | 42.9 |
| iPSC D₂O | 0.060 | 0,1,2,3,4,6,8,12,24 | 2 | 18,448 | 2,972 | 0.669 | 33.2 |

AC16 IDs **more** peptidoforms and proteins — the deficit is not coverage or intensity.

## FS accumulation — AC16 barely moves in 8 h (population median FS)

| t (h) | 0 | 1 | 2 | 3 | 6 | 8 | 12 | 24 |
|---|---|---|---|---|---|---|---|---|
| AC16 D₂O | −0.06 | −0.04 | −0.03 | 0.00 | 0.13 | **0.14** | — | — |
| iPSC D₂O | −0.04 | −0.02 | 0.01 | 0.06 | 0.16 | **0.22** | 0.36 | 0.54 |
| AC16 ¹⁸O | −0.09 | −0.08 | −0.06 | −0.04 | 0.04 | **0.07** | — | — |
| iPSC ¹⁸O | −0.04 | −0.01 | 0.06 | 0.11 | 0.26 | **0.33** | 0.43 | 0.65 |

Monotonic everywhere (no swaps). By 8 h AC16 reaches only ~0.07–0.14 FS — barely above the
per-point noise — and then the series **stops**, with no 12 h/24 h high-signal points that
carry the iPSC fits.

## Window-subset experiment (the decisive test)

Re-fit the **iPSC** integrate restricted to its ≤ 8 h timepoints (0,1,2,3,4,6,8), at the same
uniform depth-6 + label-Spep curation, and compare to AC16 (8 h) and full iPSC (24 h):

**D₂O** (RIA identical, 0.060 — isolates window + cell turnover from RIA):

| | median R² | %R²>0.8 + Spep≥8 |
|---|---|---|
| AC16 8 h | −0.05 | 3.9 % |
| **iPSC ≤8 h** | **0.31** | **12.2 %** |
| iPSC full 24 h | 0.67 | 30.9 % |

**¹⁸O:**

| | median R² | %R²>0.8 + Spep≥6 |
|---|---|---|
| AC16 8 h | −0.26 | 1.9 % |
| **iPSC ≤8 h** | **0.50** | **22.1 %** |
| iPSC full 24 h | 0.75 | 28.5 % |

Reading the D₂O column (the clean, RIA-matched one): the **window** takes iPSC 30.9 → 12.2 %
(**2.7×**); the **residual cell/turnover** gap is iPSC-≤8h 12.2 % vs AC16 3.9 % (another
**2.7×**). They compound to the ~7× full gap (30.9/3.9). For ¹⁸O the residual is larger (~7×)
because AC16's lower RIA (0.058) shrinks the per-FS envelope shift on top of the slower
accumulation.

## Hypotheses, adjudicated

| hypothesis | verdict |
|---|---|
| ¹⁸O 9 % vs 6 % RIA | **Contributes (¹⁸O only).** The D₂O gap persists at identical RIA, so RIA is not the cause — but the lower AC16 ¹⁸O RIA does enlarge the ¹⁸O-specific gap (residual 7× vs D₂O's 2.7×). |
| Swapped timepoints | **Rejected.** FS is monotonic in every file. |
| Low intensity / poor quant | **Rejected.** AC16 has *more* IDs (37 k pep / 3 740 prot) than iPSC. |
| AC16 not proliferating (no turnover) | **Partly.** AC16's *curated* k ≈ iPSC's (the measurable proteins turn over normally), but the *population* accumulates less FS by 8 h (0.14 vs 0.22, D₂O) — a real but modest (~1.6×) effective-turnover deficit, not "no turnover". |
| 8 h window too short | **Confirmed, primary.** The iPSC→≤8h subset drops curation 2.7× (D₂O) / 1.8× (¹⁸O); AC16's series stops before the high-FS 12–24 h points that carry the iPSC fits. |

> *Caveat:* the iPSC-≤8h vs AC16 residual also carries a **bio-replicate** difference (iPSC 2
> vs AC16 1 → more points per curve, better R²). It was not separated here; the FS-accumulation
> deficit (0.14 vs 0.22) is biorep-independent and is the core cell signal.

## Conclusion & implications

AC16's low curation is **expected and benign** — a short-window × lower-accumulation product,
not a pipeline or quant defect. Consequences:

- **AC16 is a method-development / calibration dataset, not a turnover-conclusions one** (as
  the technical note §5.2 already states). Its head-to-head ρ ≈ 0.66 on the small curated set
  is reassuring but underpowered.
- **The fix is experimental:** extend the AC16 labelling window past 8 h (toward the 24 h that
  makes iPSC analysis-grade), or accept the small curated subset.
- **Nothing to change in the model or curation gates** — the gates behave correctly; AC16
  simply has little signal to gate.

## Reproduce

```bash
# window-subset: fit the iPSC integrate restricted to <=8 h, compare curation
#   (build a manifest with only the integrate rows where labeling_time <= 8, then:)
riana fit --manifest <ipsc_d2o_8h_manifest.tsv> --label hw --coefficients alamillo_2025_ipsc --depth 3
python runs/headtohead.py runs/juber_ac16_o18/riana_fit_peptides.txt runs/juber_ac16_d2o/riana_fit_peptides.txt AC16
```
