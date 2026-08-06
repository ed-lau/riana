[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10614233.svg)](https://doi.org/10.5281/zenodo.10614233)

# Riana — Relative Isotope Abundance Analyzer

Riana takes standard mass-spectrometry spectra (mzML) and peptide identifications
(quantms **mzTab** for DDA, **DIA-NN** `report.parquet` for DIA, or Percolator
output) and returns mass isotopomer distributions for protein-turnover analysis.
It then fits kinetic models (one-exponential, Guan, Fornasiero) to the
time-series, rolls peptides up to proteins, and can test **cross-condition
turnover differences** (Δk) with a linearized model. A PySide6 desktop GUI
(`riana gui`) drives the same steps interactively.

Full documentation: <https://ed-lau.github.io/riana/>

## Install

Riana requires Python 3.10 or newer. We recommend a virtual environment.

```bash
pip install riana
```

For a development install from a clone:

```bash
git clone https://github.com/ed-lau/riana
cd riana
pip install -e ".[dev]"
```

## Quickstart

### Project workflow (recommended)

Riana's primary path keys every run off an **SDRF** sample sheet and a quantms
**mzTab** (DDA) or DIA-NN **report.parquet** (DIA, with the `[dia]` extra). It
writes one identity-stamped `<run>_riana.txt` per run plus a `riana_manifest.tsv`
that chains the stages, so `fit` and `rollup` re-group runs from the manifest
rather than from filenames:

```bash
# 1. Integrate every run in the experiment (one mzML in memory per worker)
riana integrate <mzml_dir> report.mzTab --sdrf samplesheet.sdrf.tsv \
    --workers 4 --out ./out

# 2. Fit the kinetic curve per peptidoform (grouped by the manifest)
riana fit --manifest ./out/riana_manifest.tsv \
    --coefficients deberneh_2025_rss --ria 0.06 --depth 3 --out ./out

# 3. Roll peptides up to proteins
riana rollup --manifest ./out/riana_manifest.tsv \
    --parsimony unique --method weighted --out ./out
```

Add `--mbr` to `integrate` to recover time points lost to stochastic MS2
sampling (gated match-between-runs, DDA only). For a two-condition experiment,
`rollup --model "linear simple" --reference-condition <name>` fits turnover in
φ-space and reports a per-protein Δk with Benjamini–Hochberg-adjusted p-values.

### Single-fraction (Percolator) path

Without `--sdrf`, the ID file is read as a Percolator `target.psms.txt` for a
single mzML (a simpler, demoted tier):

```bash
riana integrate <mzml_dir> <percolator_psms.txt> \
    --sample time1 --iso "0 1 2 3 4 5" --q_value 0.01 --mass_tol 25 --out ./out

riana fit ./out/time0_riana.txt ./out/time1_riana.txt ./out/time3_riana.txt \
    --model simple --label hw --coefficients deberneh_2025_rss --ria 0.06 --out ./out
```

Fitting supports both metabolic-water labels — heavy water (D₂O, `--label hw`, the
default) and ¹⁸O water (`--label o18`) — each with its own `--coefficients` table
family. For **D₂O**, pass a per-amino-acid labeling-site table: a bundled preset
(`deberneh_2025_rss` [recommended], `ilchenko_2019`, `commerford_1983`, or the
`alamillo_2025_{ac16,ipsc,cm}` calibration tables) or a path to your own
`(amino_acid, coefficient)` CSV. For **¹⁸O**, pass a length-model table
(`juber_2026_o18_ac16` [in-vitro], `rachdaoui_2009_o18` [in-vivo mouse], or a
`(feature, coefficient)` CSV).

By default integration uses an apex-centred narrow window
(`--integration-half-width 0.15`, dial it to your chromatographic peak width);
power-user dials live under the **Advanced integration** group of
`riana integrate --help`. To reproduce the 0.9.0 fixed-window behaviour, add
`--peak-rt ms2 --integration-half-width 1.0`.

### GUI

`riana gui` (install the `[gui]` extra) opens a PySide6 desktop app with
**Integrate**, **Model** (fit), and **Protein** (rollup) tabs over the same
engine, with interactive chromatogram, fitted-curve, and φ-space plots.

See `riana <command> --help` for the full argument set, or the
[online docs](https://ed-lau.github.io/riana/) for tutorials. (List flags like
`--iso` take a single comma/space-separated token.)

## File formats

Inputs:

- **mzML** (gzipped or plain) — MS1 spectra, streamed one fraction at a time
- **SDRF** `.sdrf.tsv` — the sample sheet that carries run identity (condition,
  biological replicate, labeling time, acquisition, precursor enrichment); the
  primary intake key
- **mzTab** (quantms / OpenMS) — the DDA peptide identifications
- **DIA-NN** `report.parquet` — the DIA identifications (apex RT resolved to the
  nearest MS1 scan); needs the `[dia]` extra
- **Percolator** `target.psms.txt` — single-mzML demoted tier (header-autodetected)

Outputs (each with a provenance header):

- `<run>_riana.txt` — one row per PSM, one column per integrated isotopomer
- `riana_manifest.tsv` — the stage-aware project index (`integrate` / `fit` /
  `rollup` rows) that chains the steps
- `riana_fit_peptides.txt` / `riana_fit_fractions.txt` — per-peptidoform kinetics
  and the per-timepoint fraction-new with prediction intervals
- `riana_rollup_proteins.txt` / `riana_rollup_fractions.txt` — protein-level k
  (and Δk under the `linear simple` model)

## Pipeline

Riana is orchestration-agnostic: search + identification are owned upstream
(e.g. quantms for DDA, DIA-NN for DIA), and Riana is a linear
`integrate → fit → rollup` chain glued by the manifest, which you compose into
whatever workflow already runs them. The bundled Snakemake example was retired in
1.0.0 — drive the subcommands directly, or from your own workflow manager.

## Citation

If you use Riana in published work, please consider citing the following papers:

> Alamillo, L. et al. *Protocol to Measure Protein Half-Life in Cell Culture Using Heavy Water*. STAR Protoc. 2026
> <https://doi.org/10.1016/j.xpro.2026.104426>

> Currie, J. et al. *Improved Method to Determine Protein Turnover Rates with Heavy Water Labeling by Mass Isotopomer Ratio Selection*. J Proteom Res. 2025
> <https://doi.org/10.1021/acs.jproteome.4c01012>

> Alamillo, L. et al. *Deuterium Labeling Enables Proteome Wide Turnover Kinetics in Cell Culture*. Cell Rep Methods. 2025
> <https://doi.org/10.1016/j.crmeth.2025.101104>

> Hammond, D. et al. *Harmonizing Labeling and Analytical Strategies to Obtain Protein Turnover Rates in Intact Adult Animals*. Mol Cell Proteomics. 2022
> <https://doi.org/10.1016/j.mcpro.2022.100252>

## Contributing

Bug reports and pull requests welcome at
<https://github.com/ed-lau/riana/issues>. See `PROJECT_REVIEW.md` for the
current development roadmap.

## License

MIT — see [LICENSE.md](LICENSE.md).
