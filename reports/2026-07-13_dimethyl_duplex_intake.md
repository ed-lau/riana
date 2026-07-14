# Dimethyl duplex — channel→sample intake (1.2.0), validated end-to-end

**Date:** 2026-07-13 · **Branch:** `1.2.0` · **Data:** `data/singlepoint_dimethyld2o`
· **Run:** `runs/dimethyl_duplex` · **Tests:** `tests/test_multiplex.py`

## What shipped

RIANA's first **true (sample-axis) multiplexing** intake. Unlike isobaric TMT (channels
share one MS1 cluster → *merged* into an average), dimethyl shifts the whole precursor
(**+8.0444 Da per site**, sites `S = 1 N-term + #K`), so the channels are separate MS1
envelopes and are now kept as **distinct samples**, each fit at its own precursor
enrichment and compared at rollup as separate conditions.

- **`riana/multiplex.py`** (new) — the multiplexing-label registry: channels ↔ SDRF
  `comment[label]` CV terms ↔ the UNIMOD mod a peptidoform carries, plus the site rule
  and per-residue heavy shift. Dimethyl wired (light 36 / medium 199 / heavy 330);
  **SILAC K+6/R+10 registered as geometry** to prove the heterogeneous-shift model (K and
  R shift differently — the offset is a per-residue sum, not uniform).
- **`riana/constants.py`** — heavy dimethyl **UNIMOD:330** (`mod_atoms` `[0,-2,0,0,0,0]`
  + `mod_fixed_isotopes` `[(2,C13),(6,D)]`) and **UNIMOD:199** (DIMETHYL4, prewired) as
  **pinned-isotope pseudo-elements** — the same machinery TMT introduced. New
  `D_MASS = 2.01410177819` (IsoSpec's built-in ²H mass; NIST agrees to 3e-10). 36/199/330
  whitelisted so their peptidoforms survive intake, and deliberately **not** in
  `CHEMICAL_MODS` (that would merge the channels onto one curve — wrong for a sample-axis mod).
- **`io/sdrf.py`** — a multiplex branch beside the isobaric one: N rows/file allowed, every
  channel row **kept** as its own `RunIdentity` (distinct sample/condition/enrichment), and a
  channel-keyed map `{(stem, (label, channel)): RunIdentity}`.
- **`io/mztab.py`** — each PSM is routed to its channel **by its own dimethyl mod**
  (`multiplex.channel_of`), not by one identity per file.
- **`core/pipeline.py`** — `plan_integration` emits **one RunTask per channel** (shared mzML,
  distinct output stem = the SDRF `source name`).

**No manifest or fit schema change was needed.** `fit_project` already resolves RIA
*per `(experiment, condition)` group* — so once the intake gives each channel its own
condition + enrichment, the existing machinery fits them at their own RIA, and the existing
`linear simple` two-condition Δk works unchanged.

## A quantms 1.8.0 regression found on the way (general, not multiplex)

The dataset was searched with **quantms 1.8.0 (OpenMS 3.6.0)**; every prior set used 1.7.0
(OpenMS 3.5.0). 1.8.0 **stopped writing the optional `opt_global_q-value` column**. RIANA
defaulted every PSM to `q = 1.0`, which can never pass the strict `q < --q_value` gate (capped
at 1.0) → **every PSM silently dropped, whole file unusable.**

The q-value did not disappear, it *moved*: it is `search_engine_score[1]`, whose type the
metadata declares — `MS:1001491 percolator:Q value` in 1.8.0, `MS:1003115 OpenMS target-decoy
q-value` in 1.7.0. (Confirmed against `OpenMS/src/openms/source/FORMAT/MzTab.cpp`, which stamps
`[MS,MS:1001491,percolator:Q value,]`; our file's max q = 9.95e-03 is exactly the 1 % FDR cut.)

`io.mztab._resolve_q_value_column` now prefers the explicit column and falls back to
`search_engine_score[1]` **only when the metadata declares it a q-value** — exact FDR semantics
preserved, and **pre-1.8.0 files are byte-identically unaffected**. This was *not* fixed by
substituting PEP (a different, stricter score) or by loosening `<` to `<=`. Regression tests
cover all three cases (`tests/test_io.py`).

## End-to-end validation (`runs/dimethyl_duplex`)

`integrate -W 6 -q 0.01` → `fit --depth 1 --coefficients deberneh_2025_rss` →
`rollup --model "linear simple" --reference-condition control --test-condition rapamycin`.

**Channel separation is exact.** 3 mzML × 2 channels → **6 per-channel outputs**; each
provenance header carries its own identity (`LIV_ctrl_r1` / control / **0.04614** vs
`LIV_rapa_r1` / rapamycin / **0.05611**). `LIV_ctrl_r1` holds 4209 `[UNIMOD:36]` peptidoforms
and **zero** `[UNIMOD:330]`; `LIV_rapa_r1` holds 7180 heavy and **zero** light.

**Per-channel enrichment is applied** (the load-bearing check):
```
fitting curve condition=control   (8453 rows;  RIA=0.0461 from manifest precursor_enrichment)
fitting curve condition=rapamycin (10908 rows; RIA=0.0561 from manifest precursor_enrichment)
```

**The biology is coherent.** 5779 peptides converged; 339 proteins. Median k: control
**0.178/day**, rapamycin **0.134/day**. Of 77 proteins with a Δk, **62 % are slower under
rapamycin** (median Δk **−0.022/day**); of the 33 significant at `p_adj < 0.05`, **73 % slower**
(most-slowed: P49429, Q9DBF1, P11499, O09173, P35505). mTOR inhibition slowing protein
turnover is the expected direction.

## The heavy channel is degraded — and it is mostly **SNR, not spillover**

Heavy-channel precision is 4–5× worse than light (k_cv median 0.36 vs 0.088), and the loss is
**not confined to S=1**, so the S=1 spillover cannot be the main cause. The direct measurement:

| | light (control) | heavy (rapamycin) |
|---|---|---|
| iso0 intensity (median) | 2.44e4 | **4.78e3** |
| rows with zero envelope | 0.1 % | **8.0 %** |
| paired heavy/light intensity, 1908 shared backbones | — | **0.628** |

The duplex is **not balanced**: on matched peptide backbones the heavy channel carries ~37 %
less signal, with 80× more dropouts. Lower SNR → noisier FS (sd 0.279 vs 0.211) → more FS
rail-drops → fewer surviving replicates (rapamycin `n_points` skews to 0–1; control to 3) →
worse k_cv → fewer proteins clear the rollup. Replicate depth is the dominant k_cv lever
(the single-timepoint finding from the TMT work), and SNR hits **every** S equally — which is
exactly the S-independent pattern observed. **This is a property of the sample/mixing, not a
pipeline defect.**

**The S=1 spillover is nevertheless visible on top of that floor:** within the heavy channel,
k_cv is worst at **S=1 (0.419)** vs **S=2 (0.306)** — the +8 Da-adjacent population, ungated in
this increment, exactly as `reports/2026-07-06_dimethyl_duplex_spillover.md` predicted.

**Does spillover confound the Δk?** No — it works *against* the observed effect. Light
contamination would drag the heavy FS *up* (the contaminating light channel is more labelled),
making rapamycin look **faster**. We observe it slower, so the biological signal survives the
confound rather than being manufactured by it.

## Caveats / next

- **Long S=1 heavy peptides are ungated** — their heavy FS is light-contaminated. The **exact
  forward-model spillover gate** is the next increment (specced in the plan; primitive factored
  from `tests/benchmark/bench_dimethyl_spillover.py`). Do not lean on individual S=1 heavy k yet.
- **The heavy/light imbalance (0.628) should be checked at the bench** — it caps what this
  dataset can deliver regardless of software. Δk *magnitudes* are provisional until the gate
  lands and the imbalance is understood; the *direction* is robust.
- The fit's provenance header records `ria_max 0.06` (the config default) even though each curve
  is correctly fit at its per-group manifest RIA — a cosmetic provenance nit, pre-existing.
