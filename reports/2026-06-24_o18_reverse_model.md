# ¹⁸O (H₂¹⁸O) reverse model + production fit-side

- **Date:** 2026-06-24
- **Branch:** `1.1.0`
- **Status:** **SHIPPED** — `fit --label o18` works end-to-end; coefficient structure decided (DENQ**+S**); two presets bundled (`o18_ac16` in-vitro, `o18_previs` in-vivo).
- **Harness:** [`tests/benchmark/bench_o18_coefficients.py`](../tests/benchmark/bench_o18_coefficients.py) + [`tests/benchmark/_helpers/o18_forward_model.py`](../tests/benchmark/_helpers/o18_forward_model.py) (port of NB90c); guards in [`tests/test_o18.py`](../tests/test_o18.py).
- **Roadmap:** `PROJECT_REVIEW.md` → §3 handoff item #3 (o18 rewrite).

## Question

`riana fit --label o18` errored ("reimplemented post-M4"). Turning it on needed (a) the ¹⁸O **labeling-site model** — how many oxygens per peptide exchange into the H₂¹⁸O pool (the "Spep") — trained from calibration data via the reverse model, and (b) the production envelope + FS solver. Two science questions fell out:

1. **What is the coefficient structure?** ¹⁸O labels the backbone + which side chains? Is NB90c's `D/E/N/Q` model right, or do the hydroxyl residues `S/T/Y` add signal? Is the backbone `L−1`, `L`, or `L+1`?
2. **Does the production path recover ground truth** on the AC16 calibration (0→100% H₂¹⁸O mixing series)?

## TL;DR

- **Structure = `b·(L−1) + c_D·D + c_E·E + c_N·N + c_Q·Q + c_S·S`** (6-param length model). **Serine added** to NB90c's DENQ (strong: +3.4 pt held-out R², 21σ; T/Y empirically null). **Backbone `L−1`** confirmed by the data *and* by Previs (MCP 2009): one ¹⁸O per peptide bond, the terminal amino acid's O back-exchanges in tryptic digest.
- **Trained on AC16** (2083 peptides, **test R² 0.899**): `Spep = 0.159·(L−1) + 1.549·D + 1.217·E + 0.377·N + 0.028·Q + 0.524·S` → bundled `o18_ac16`.
- **`fit --label o18` shipped** — the production envelope matches the NB90c oracle to **5e-5**; **curated recovery on real AC16 tracks the proportions** (|bias| ~0.03–0.07 mid-range, monotonic). The D₂O path is untouched.
- **In-vivo reference bundled** (`o18_previs`): the literature mouse model `1·(L−1) + 2·E + N + Q` (reproduces Previs' own worked example LGEYGFQNAILVR = 16). The in-vivo/in-vitro gap is large (backbone 1.0 vs 0.16).

## The ¹⁸O labeling model — what exchanges

Metabolic H₂¹⁸O labels oxygen incorporated into proteins. From Rachdaoui/Previs (*Mol Cell Proteomics* 2009; 8(12):2653, "Measuring Proteome Dynamics in Vivo"):

- **Backbone: 1 ¹⁸O per peptide bond.** A peptide of L residues has L+1 backbone O = (L−1) internal peptide-bond carbonyls + 2 C-terminal carboxyl O. The C-terminal carboxyl oxygens **back-exchange with the H₂¹⁶O digest**, so the stable count is **L−1** — verbatim, "the number of amino acids in the peptide minus one that is (back-)exchanged from the terminal amino acid during tryptic cleavage."
- **Side chains (in vivo):** glutamate **+2**, asparagine **+1**, glutamine **+1** (via de-/reamination). Aspartate is **omitted** from the explicit count (the discussion expects it ~asparagine).
- **Caveat — non-homogeneous labeling:** backbone peptide-bond O may label to a *greater* degree than side-chain O (side chains lag, labeling only via turnover/reamination). Our solver assumes one RIA for all sites — fine for the steady-state in-vitro calibration (full equilibration), a refinement to watch for non-steady-state in-vivo data.

## The reverse model (harness)

Mirrors the D₂O NB87a trainer (`bench_aa_coefficients.py`). Per peptidoform, the IsoSpec forward model builds the **natural** (θ=0) and **fully-labeled** (at the calibration RIA, with `n` candidate ¹⁸O sites) envelopes — ¹⁸O is a **3-isotope** enriched-O pseudo-element (¹⁶O/¹⁷O/¹⁸O), unlike D₂O's 2-isotope H/D, so it shifts mass **+2 Da/site**. For each peptide a **continuous Spep** is fit by minimizing SSE between the predicted `(1−f)·init + f·final` mixture and the observed iso0–iso5 across all 9 mixing proportions; the per-peptide Spep is then regressed onto sequence features (bounded `lsq_linear`, intercept 0). **N_ISO=6** (the +2 Da label spreads signal into m2/m4/m6, unlike D₂O where iso0–3 suffice).

## Coefficient-structure decision

**Serine helps; threonine/tyrosine don't** (AC16, 5-fold CV; re-fit on the same per-peptide Spep, so cheap):

| model | test R² | 5-fold CV R² |
|---|---|---|
| DENQ (NB90c) | 0.862 | 0.861 ± 0.017 |
| **+ S** | **0.899** | **0.891 ± 0.019** |
| + T | 0.862 | no change |
| + Y | 0.862 | no change |

`c_S = 0.52 ± 0.025` (**21σ**); `c_T, c_Y → 0` exactly. Biochemically coherent — serine sits in one-carbon metabolism (its O exchanges), threonine/tyrosine retain theirs. (Cysteine is excluded — its O is chemically-added carbamidomethyl, not metabolic.)

**Backbone `L−1` beats `L` and `L+1`** (the no-intercept model makes the constant offset a real predictive difference):

| backbone | b | test R² | CV R² |
|---|---|---|---|
| **L−1** | 0.159 | **0.899** | **0.891** |
| L | 0.147 | 0.897 | 0.889 |
| L+1 | 0.137 | 0.894 | 0.886 |

Monotonic; the margin is small, but L−1 wins on every metric — and it is exactly the Previs model.

## Results — AC16 in-vitro coefficients (`o18_ac16`)

`Spep = 0.159·(L−1) + 1.549·D + 1.217·E + 0.377·N + 0.028·Q + 0.524·S`; train/test R² **0.891/0.899**; 2083 curated peptides (≈ NB90c's 2087).

| feature | coef | ±SE |
|---|---|---|
| (L−1) | 0.159 | 0.005 |
| D | 1.549 | 0.026 |
| E | 1.217 | 0.019 |
| N | 0.377 | 0.032 |
| Q | 0.028 | 0.028 |
| S | 0.524 | 0.025 |

**Curated FS recovery** on the real AC16 calibration (R²>0.95, Spep≥3, 1987 peptides), via production `solve_fs_o18`:

| true % | median FS | \|bias\| |
|---|---|---|
| 12.5 | 0.094 | 0.031 |
| 25.0 | 0.213 | 0.037 |
| 37.5 | 0.348 | 0.027 |
| 50.0 | 0.431 | 0.069 |
| 62.5 | 0.598 | 0.027 |
| 75.0 | 0.754 | 0.004 |
| 87.5 | 0.815 | 0.060 |

Monotonic, |bias| ~0.03–0.07 mid-range (the 0/100% boundaries are degenerate and excluded, per NB90c). The *uncurated* median is biased lower — noise, removed by the R²>0.95 gate.

## In-vivo vs in-vitro (`o18_previs`)

In vitro the **side-chain carboxyl/hydroxyl O dominate** (D/E near the 2-O ceiling, S=0.52) while the **backbone barely labels** (b=0.16) — backbone carbonyls label only via slow de-novo synthesis, low in cultured cells. In vivo (Previs) the **backbone dominates** (1 ¹⁸O/bond → b≈1) over weeks. Same sites, opposite emphasis by regime — so o18 (like D₂O) needs **per-regime / per-cell-type tables**. The literature in-vivo model ships as `o18_previs` (`1·(L−1) + 2·E + N + Q`, D=S=0); it reproduces Previs' worked example LGEYGFQNAILVR=16. **Serine is our in-vitro empirical addition** — it is *not* in the Previs enumeration; kept for the fitted tables, S=0 in the in-vivo reference.

## Production (`fit --label o18`)

- **Envelope:** `algorithms/isotope_dist.get_peptide_distribution` `label==3` completed (3-isotope enriched-O); `solve_fs_o18` (H4′ mix-then-normalize, parallel to `solve_fs_d2o`, which is left untouched); `spep_from_length_coefficients`.
- **Dispatch:** `core/fitting` selects the Spep model + FS solver on `config.label`; `load_o18_coefficients` reads the `(feature, coefficient)` table. The CLI requires the matching `--coefficients` format per label.
- **Validated:** production envelope == oracle (5e-5); synthetic FS exact; curated real-data recovery (above); D₂O `test_fitting` green; `tests/test_o18.py` (9 guards, incl. the Previs worked example).
- **Deferred:** the internal `{1,2,3}`→`"D2O"/"O18"` label-taxonomy cleanup (cosmetic); the **kinetic** `fit --label o18` curve awaits an o18 **time series** (the calibration is mixing-proportion, not time) — the planned iPSC d2o-vs-o18 head-to-head.

## Reproduce

```bash
# integrate the 9 o18 calibration runs (mzML + quantms mzTab, via the SDRF)
riana integrate data/calibration_o18/mzml \
  data/calibration_o18/quantms_results/quant_tables/samplesheet_o18.sdrf_openms_design_openms.mzTab \
  --sdrf data/calibration_o18/samplesheet_o18.sdrf.tsv -W 4 \
  -o data/calibration_o18/integrate_outputs
# train the length-model coefficients (N_ISO=6)
python tests/benchmark/bench_o18_coefficients.py \
  --inputs data/calibration_o18/integrate_outputs \
  --ground-truth data/calibration_o18/ground_truth_o18.csv \
  --output-dir data/calibration_o18/o18_coefficients --n-iso 6
# fit with the bundled preset (o18_ac16 in-vitro, or o18_previs in-vivo)
riana fit <riana_or_manifest> --label o18 --coefficients o18_ac16
```
