# D₂O labelling-site tables — Ilchenko 2019 & Deberneh 2025 vs Commerford 1983

- **Date:** 2026-06-25
- **Branch:** `1.1.0`
- **Status:** experiment complete — **both LC-MS in-vivo tables tighten within-protein k agreement ~18–21% vs Commerford 1983**; shipped as presets (`ilchenko_2019`, `deberneh_2025_rss`); `commerford` renamed `commerford_1983`. They do **not** move the fs↔fs_ds residual → that residual is not the coefficient table.
- **Bench:** `tests/benchmark/bench_coefficient_tables.py`
- **Sources:** Ilchenko/Sadygov 2019 ([PMC8201887](https://pmc.ncbi.nlm.nih.gov/articles/PMC8201887/), Table 2, N_aa); Deberneh 2025 ([PMC12302360](https://pmc.ncbi.nlm.nih.gov/articles/PMC12302360/), Table S1, dataset-1)
- **Roadmap / memory:** the Commerford/LVE residual follow-on of `mass_defect_theta_design`; ties to Track D (within-protein θ).

## Question

The `fs_ds` (mass-defect) work surfaced a residual on in-vivo LVE that an earlier
hypothesis pinned on the **Commerford 1983 tritium** labelling-site table being off for
LC-MS D₂O. Two modern LC-MS-derived per-AA tables exist (Sadygov-lineage nested fits).
Do they (a) tighten the **intensity-FS** itself — measured as within-protein k agreement
— and/or (b) improve the **fs↔fs_ds** cross-check?

## TL;DR

- **Both new tables beat Commerford 1983 on within-protein k consistency** (a protein's
  peptides should agree on one k): median robust geometric k-CV **0.199 → 0.158
  (ilchenko, −21%) / 0.162 (deberneh RSS, −18%)** on the LVE manifest. Strong evidence
  the 1983 tritium values are suboptimal for LC-MS D₂O and the LC-MS tables are better.
- **fs↔fs_ds agreement is unchanged** across tables (bias +0.011 → +0.016/+0.020;
  within-0.1 ~28–30%). So the LVE fs_ds residual (+0.084, report 2026-06-25) is **not**
  the coefficient table — it lies elsewhere (still open).
- **Shipped** `ilchenko_2019` + `deberneh_2025_rss`; renamed `commerford` → `commerford_1983`.
  RSS is the **method the Deberneh authors recommend** (Pearson 0.98 with truth; MPE needs
  2× the spectral accuracy). Default preset unchanged (`commerford_1983`) pending a wider
  cross-tissue check before switching.

## Tables shipped (our `(amino_acid, coefficient)` spec)

- **`ilchenko_2019`** — Table 2 N_aa (incorporated H per AA), M0-M1 isotopomers, WT-mouse
  plasma tryptic peptides, 1 week ²H₂O.
- **`deberneh_2025_rss`** — Table S1 **RSS** column (residual-sum-of-squares nested
  whole-dataset fit, dataset-1), the authors' recommended method.
- `commerford_1983` — the prior `commerford` preset (= Deberneh's "Tritium" reference
  column), renamed for provenance clarity.

The LC-MS tables differ from the 1983 tritium values most at **R, H, K, P, S** (all
lower) and **D** (higher): e.g. Arg 3.34→1.33/1.70, His 2.88→1.37/1.96, Lys 0.54→0.12/0.05,
Pro 2.59→0.78/1.68, Asp 1.89→2.49/2.98. This matches the Deberneh paper's own caveats vs
tritium (Asp computed >55% **higher**, Arg ~40% **lower**, His/Pro/Ser lower; Ala/Gly/Gln/
Glu/Cys agree to ~10%).

## Results

`bench_coefficient_tables.py` on `runs/lve_fixed_ab/riana_manifest.tsv`, RIA 0.045, intensity
FS (`simple`), within-protein robust geometric k-CV = `1.4826·MAD(ln k)` over curated
(R²≥0.9) peptides, median over (protein, condition) groups with ≥3 peptides:

| table | med within-protein k geo-CV | n groups | fs_ds−fs bias | \|Δ\|<0.1 |
|---|---|---|---|---|
| commerford_1983 | 0.199 | 734 | +0.011 | 29.5% |
| **ilchenko_2019** | **0.158** | 661 | +0.020 | 28.3% |
| **deberneh_2025_rss** | **0.162** | 738 | +0.016 | 28.5% |

(18,701 peptides fitted per table.)

- **Within-protein k-CV drops ~18–21%** with either LC-MS table — the headline win. The
  better Spep makes a protein's peptides agree on k.
- **ilchenko** has the lowest CV but **fewer curated groups** (661 vs 738) — a few proteins
  drop below the ≥3-curated-peptide bar (its very low K/T/W coefficients shift some short
  peptides' fits). **deberneh_rss** gives nearly-as-low CV while **retaining the most
  groups**, so it is the more robust practical choice.
- **fs↔fs_ds agreement does not improve** → the mass-defect residual is independent of the
  per-AA table (it is not a Commerford artifact). Open for separate investigation.

## Interpretation

The 1983 Commerford tritium values systematically over-state labile H at R/H/K/P/S for
LC-MS D₂O; the LC-MS in-vivo tables (Ilchenko M0-M1; Deberneh RSS nested fit) reflect the
actual incorporation better, so per-peptide k within a protein converges. Deberneh's own
turnover check (k within ~8% of tritium for abundant proteins despite per-AA differences)
explains why the effect is a **tightening**, not a wholesale k shift. The fs_ds residual
being table-invariant rules out Spep as its cause.

## Recommendation / next

- Keep `deberneh_2025_rss` and `ilchenko_2019` as presets; **deberneh_2025_rss is the
  recommended in-vivo D₂O table** (lowest CV among those that keep full curated coverage,
  and the paper's recommended RSS method).
- Defer switching the **default** preset from `commerford_1983` until a cross-tissue check
  (LVE is one in-vivo dataset; the calibration cell lines have their own trained tables).
- The fs↔fs_ds LVE residual stays open (not the table) — revisit with the RIA / mixture-
  model angle, not the Spep table.

## Reproduce

```bash
python -m tests.benchmark.bench_coefficient_tables \
  runs/lve_fixed_ab/riana_manifest.tsv --ria 0.045 \
  --tables commerford_1983 ilchenko_2019 deberneh_2025_rss
# (calibration mzML/PSMs gitignored; regenerate the manifest via the SDRF integrate path.)
```
