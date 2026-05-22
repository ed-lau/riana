# D₂O mixing calibration dataset

Ground-truth benchmark dataset for RIANA's `integrate` pipeline. See
`PROJECT_REVIEW.md` §M2 for the rationale; the short version is that
this dataset exists so every future change to integration, peak
detection, smoothing, or baseline subtraction can be regression-gated
against a fixed reference.

## What it is

Lysate from cells cultured in 6% D₂O for ≥10 doublings (effectively
complete proteome labelling) was mixed with lysate from unlabelled
cells at 9 nominal heavy fractions:

```
0%, 12.5%, 25%, 37.5%, 50%, 62.5%, 75%, 87.5%, 100%
```

One LC-MS run per fraction per cell line. Three cell lines are treated
**independently** in all downstream analysis (separate `ground_truth.csv`,
separate coefficient tables, separate baseline results):

- `ac16/` — AC16, a transformed human cardiac fibroblast/myocyte hybrid
  line. **Proliferative.** Raw files: JPOST `JPST002443`.
- `ipsc/` — human induced pluripotent stem cells. **Proliferative.**
  Raw files: JPOST `JPST003556`.
- `cm/` — contractile human iPSC-derived cardiomyocytes (iPSC-CM).
  **Post-mitotic.** Raw files: JPOST `JPST003582`.

The third line is deliberate. AC16 and iPSC are both proliferative, and a
proliferating cell dilutes isotopic label through cell division independently
of protein turnover; iPSC-CM is post-mitotic, so it isolates turnover from
division and is the more physiologically relevant test case. Neither AC16 nor
iPSC-CM is primary cardiomyocyte, but the proliferative-vs-post-mitotic
contrast is exactly the axis that lets the cross-line coefficient comparison
separate real biology from integration artifact (`PROJECT_REVIEW.md`, M2
finding 3): per-AA coefficients that track the proliferative/post-mitotic
split are a biology signal; coefficients that scatter without it point at the
integrator.

`comment[file uri]` in the per-line SDRF (`data/calibration_<line>/samplesheet_<line>_alpine.sdrf.tsv`) has the per-fraction download URLs (`https://storage.jpostdb.org/...`).

## What's in each cell-line subdir

| File / dir | Tracked | Notes |
|---|---|---|
| `ground_truth.csv` | yes | `(mzml_filename, nominal_proportion)` from filename + SDRF. |
| `integrate_outputs/snakemake_reference/` | **no** | Per-fraction `_riana.txt` from the original snakemake-era run; port-validation input for `bench_aa_coefficients.py`. |
| `integrate_outputs/v0.9.0/` | **no** | Per-fraction `_riana.txt` from fresh `riana integrate` on master, pinned config. |
| `integrate_outputs/v0.9.0_smoothing_<N>/` | **no** | Per-fraction `_riana.txt` from the `--smoothing` sweep. |
| `benchmark_results/v0.9.0/` | yes | Coefficient tables, FS recovery, N_ISO + smoothing sweeps, summary metrics. Small CSVs/JSON. |
| `benchmark_results/snakemake_reference/` | yes | Same outputs computed on the snakemake-era integrate (AC16 only — the port-validation snapshot). |

`integrate_outputs/` is **gitignored** — collectively ~450 MB of regenerable
intermediate `_riana.txt` files. The committed `benchmark_results/` CSVs are
the v0.9.0 baseline of record. To regenerate the integrate outputs from raw
data: `python tests/benchmark/run_integrate_v0_9_0.py` (needs the mzML files,
see below).

## What is NOT here (heavy inputs)

Raw inputs live at the repo-root `data/calibration_{ac16,ipsc,cm}/`
directory (already gitignored). They are not committed because:

1. mzMLs and mzTabs are large (mzML ~1.5 GB/line, mzTab ~240 MB/line).
2. Raw `.raw` files are already on JPOST under the accessions above and
   are durably citable through ProteomeXchange — no separate Zenodo
   deposit is needed for the raw data.

| Lives at | Approx size | Contents |
|---|---|---|
| `data/calibration_<line>/mzml/` | 1.5 GB / line | 9 mzML.gz files |
| `data/calibration_<line>/snakemake_results/time*/percolator/` | ~22 MB / line | Per-fraction Comet+Percolator PSMs (the ID input RIANA reads). |
| `data/calibration_<line>/quantms_results/quant_tables/*.mzTab` | 242 MB / line | Single 3-engine ConsensusID mzTab covering all 9 runs. Reserved for M3 (mzTab adapter). |
| `data/calibration_<line>/uniprot_human_reviewed.fasta` | ~14 MB / line | Search database. |
| `data/calibration_<line>/samplesheet_<line>_alpine.sdrf.tsv` | ~3 KB / line | SDRF with JPOST URIs and search params. |

`run_integrate_v0_9_0.py` reads PSMs from those paths and writes
outputs into `integrate_outputs/v0.9.0/` here. It runs all three lines by
default, or one with `--line {ac16,ipsc,cm}`.

## Benchmark scripts

In `tests/benchmark/`. See that directory's README for full descriptions.

- `build_ground_truth.py` — emits `ground_truth.csv` from filename + SDRF.
- `bench_aa_coefficients.py` — ports NB87a (`data/notebook/87a_D2O_LearnAALabelingSites_IsoSpec_AC16.ipynb`). Fits per-peptide Spep via IsoSpec forward model, then per-AA non-negative regression. Emits `d2o_amino_acid_coefficients.csv` + `spep_per_peptide.csv`.
- `bench_fs_recovery.py` — using the learned coefficients, solves `fs` per (peptide, fraction) and compares to nominal `proportion/100`. Emits `fs_recovery.csv`.
- `bench_n_iso_sweep.py` — N_ISO ∈ {2..6} sensitivity sweep.
- `bench_smoothing.py` — `riana integrate --smoothing` sweep.

## Why this benchmark is *not* "predicted m0/mA vs observed"

There is no external ground truth for per-peptide isotope envelopes —
we have no animal calibration curve, no externally-measured Spep table.
The only ground truth is the **mixing proportion**.

The benchmark therefore bootstraps: integrate output → per-peptide Spep
via IsoSpec forward model → per-AA coefficients → `fs` solver → compare
to nominal proportion. This is partially circular (a bad integrator can
absorb error into the coefficients and FS-vs-truth still looks fine),
and the escape signal is **stability of the coefficients** across
N_ISO truncation, smoothing settings, and cell line. Instability ⇒
integrate is doing something wrong.

That freeze has since happened (M3 Week 0): `build_frozen_tables.py` bootstraps
a **per-cell-line** frozen coefficient table (`d2o_aa_coefficients_<line>.csv`)
and `bench_m0_ma_recovery.py` scores observed-vs-predicted m0/mA RMSE against
it. A frozen table is a *constant*, so its bias cancels in any method-vs-method
comparison — which is what M3 regression-gating needs. The table is per-line,
not universal: the M2 AC16-vs-iPSC coefficient divergence (up to 0.54) means a
single shared table is not defensible. The `cm/` line adds the post-mitotic
third point to that comparison.
