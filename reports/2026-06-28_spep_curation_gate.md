# Spep curation gate — fit quality vs labeling-site count (D₂O & ¹⁸O)

- **Date:** 2026-06-28
- **Branch:** `1.1.0`
- **Status:** **Default decided** for the v1.1.0 Spep gate (`fit`/`rollup` `--min-spep`, label-aware,
  sample-tunable): **D₂O ≥ 8, ¹⁸O ≥ 6.** The gate removes physically under-powered peptidoforms — too
  few labelling sites for the envelope to shift measurably as FS goes 0 → 1 — that the R² gate alone
  lets leak through.
- **Inputs (existing fits, each on its appropriate coefficient table):**
  `runs/lve_atr_clean` (cardiac D₂O, in vivo, 0–30 d, RIA 0.046, `deberneh_2025_rss`),
  `runs/boomi_ipsc_d2o` (iPSC D₂O, in vitro, 0–24 h, RIA 0.060, `alamillo_2025_ipsc`),
  `runs/boomi_ipsc_o18` (iPSC ¹⁸O, in vitro, 0–24 h, RIA 0.090, `juber_2026_o18_ac16`).
- **Feeds:** v1.1.0 QoL item #2 (the Spep gate); surfaced during the LVE/ATR run
  ([2026-06-28_lve_atr_d2o](2026-06-28_lve_atr_d2o.md), where Spep ≥ 6 was applied by hand).

## Question

Spep — a peptidoform's number of labelling sites (from the coefficient table) — sets **how far its
isotopomer envelope moves** between unlabelled and fully-labelled. Too few sites and the shift is
below what the MS can resolve, so the fit is unreliable *regardless of its apparent R²*. Where does
that floor sit, and is one default right for both labels? Set the v1.1.0 `--min-spep` defaults from
data, not a guess.

## TL;DR

- **The cutoff is genuinely label-specific.** A `Spep < 8` cut drops only **2 %** of cardiac D₂O's
  R² > 0.8 fits but **56 %** of iPSC ¹⁸O's — and the ¹⁸O ones are mostly *good* (R² > 0.8 holds at
  43–48 % across Spep 4–10). ¹⁸O's **+2 Da per site** makes few sites measurable; D₂O's small
  per-site mass shift does not. → **D₂O ≥ 8, ¹⁸O ≥ 6.**
- **The R² knee is ~Spep 8 for D₂O, ~6 for ¹⁸O.** Below it median R² falls steeply in every dataset;
  above it, in-vivo D₂O plateaus (~0.90) and ¹⁸O actually *peaks then declines* (Spep 6–10 best).
- **The Spep gate is a weak quality lever once R² > 0.8 is applied** — the R² > 0.8 *fraction* of
  retained peptides barely moves with the cutoff (cardiac 68→69 %, iPSC D₂O 35→38 % from Spep 6→12).
  Its real job is removing **spurious-R² leakage** (a handful of noise-driven high-R² low-site fits)
  and **not wasting compute**, not raising bulk quality.
- **It must stay a per-sample knob.** The right floor rises with shorter time series and slower
  turnover (more sites needed to see a small FS change) — iPSC D₂O (24 h window) keeps improving with
  Spep and never plateaus, where cardiac D₂O (30 d) is flat by Spep 12. The defaults are a starting
  point, not a constant.

## Fit quality vs Spep

Median R² and the R² > 0.8 fraction, binned by Spep (each dataset on its own coefficient table):

**D₂O cardiac LVE+ATR** (in vivo, median Spep 22) — climbs, **plateaus ~Spep 12**:

| Spep | n | % set | med R² | %R²>0.8 |
|---|---|---|---|---|
| [4,6) | 78 | 1% | 0.782 | 49% |
| [6,8) | 323 | 2% | 0.816 | 54% |
| [8,10) | 538 | 4% | 0.872 | 63% |
| [10,12) | 768 | 5% | 0.870 | 62% |
| [12,16) | 2,187 | 15% | 0.904 | 67% |
| [16,∞) | 10,857 | 74% | 0.905 | 69% |

**D₂O iPSC** (in vitro, median Spep 14.5) — lower overall (24 h regime); **monotone, no plateau**:

| Spep | n | % set | med R² | %R²>0.8 |
|---|---|---|---|---|
| [4,6) | 707 | 4% | 0.420 | 16% |
| [6,8) | 1,433 | 8% | 0.528 | 20% |
| [8,10) | 1,994 | 11% | 0.615 | 29% |
| [10,12) | 2,190 | 12% | 0.640 | 31% |
| [12,16) | 3,947 | 22% | 0.680 | 35% |
| [16,∞) | 7,793 | 43% | 0.723 | 40% |

**¹⁸O iPSC** (in vitro, median Spep 7.2) — **peaks at Spep 6–10, then declines**:

| Spep | n | % set | med R² | %R²>0.8 |
|---|---|---|---|---|
| [0,4) | 3,097 | 16% | 0.638 | 33% |
| [4,6) | 4,332 | 22% | 0.754 | 43% |
| [6,8) | 4,052 | 20% | 0.785 | 48% |
| [8,10) | 3,254 | 16% | 0.782 | 48% |
| [12,16) | 2,237 | 11% | 0.747 | 43% |
| [16,∞) | 722 | 4% | 0.659 | 34% |

## Recommended defaults

| label | default `--min-spep` | basis |
|---|---|---|
| **D₂O / `hw`** | **8** | R² knee; below it median R² drops ≥ 0.05 in both regimes and the small per-site shift is under-powered. (In-vivo plateau supports up to **12** for strict / clean-only use; raise for short-window or slow-turnover samples.) |
| **¹⁸O / `o18`** | **6** | Quality peaks at Spep 6–10 then declines; the +2 Da/site shift keeps Spep 4–6 measurable, so a D₂O-style 8+ would wrongly cut good fits. Matches the prior ¹⁸O curation (≥ 5). |

**Leakage removed (Spep below the label floor, but R² > 0.8):** cardiac D₂O 215 (2 % of R²>0.8),
iPSC D₂O 420 (7 %) — the noise-driven high-R² low-site fits the gate exists to catch. For ¹⁸O the
same metric at < 8 is 56 %, which is exactly why ¹⁸O's floor is lower, not higher.

## Notes & caveats

- **The AC16 ¹⁸O run is excluded** — it is regime-limited (negative median R² at every Spep over its
  8 h window), so it informs the window, not the Spep floor. Only the good-signal iPSC ¹⁸O and the two
  D₂O sets set these numbers.
- **Spep is coefficient-table-dependent and that is correct.** Each dataset uses its own table
  (`alamillo_2025_ipsc` gives genuinely lower Spep than `deberneh_2025_rss`, matching the lower
  in-vitro incorporation), and the gate acts on the Spep the fit actually used to build the envelope —
  so the floor is a property of the envelope shift, comparable across tables.
- **The Spep gate ≠ the R² gate.** R² > 0.8 already filters bulk quality; Spep adds the *physical*
  prior (can this peptide's envelope move enough to be measured at all), catching fits that pass R²
  on noise. Both are wanted.
- **Sample-tunable by design.** Defaults are a starting point; the report's own numbers show the floor
  should track RIA, time-series length, and turnover rate — hence `--min-spep` is exposed, not fixed.
