[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17613314.svg)](https://doi.org/10.5281/zenodo.17613314)

# Riana — Relative Isotope Abundance Analyzer

Riana takes standard mass-spectrometry spectra (mzML) and peptide-spectrum-match
files (Percolator output, or mzTab via quantms) and returns mass isotopomer
distributions, e.g. for protein turnover analysis. It also fits kinetic models
(one-exponential, Guan, Fornasiero) to time-series isotopomer data.

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

### Integrate

Extract isotopomer intensities from one fraction's mzML using a Percolator PSM
table:

```bash
riana integrate <mzml_dir> <percolator_psms.txt> \
    --sample time1 \
    --iso "0 1 2 3 4 5 6" \
    --q_value 0.01 \
    --mass_tol 25 \
    --out ./out
```

By default integration uses an apex-centred narrow window
(`--integration-half-width 0.15`, dial it to your chromatographic peak width).
To reproduce the 0.9.0 fixed-window behaviour, add
`--peak-rt ms2 --integration-half-width 1.0`.

### Fit

Fit a kinetic model across timepoints. Fitting is heavy-water (D₂O) only and
needs a per-amino-acid labeling-site table via `--coefficients` — a bundled
preset (`commerford` literature, or the `ac16` / `ipsc` / `cm` calibration
tables) or a path to your own `(amino_acid, coefficient)` CSV:

```bash
riana fit ./out/time0_riana.txt ./out/time1_riana.txt ./out/time3_riana.txt \
    --model simple \
    --label hw \
    --coefficients commerford \
    --ria 0.06 \
    --depth 3 \
    --out ./out
```

See `riana integrate --help` and `riana fit --help` for the full argument set,
or the [online docs](https://ed-lau.github.io/riana/) for tutorials. (List
flags like `--iso` take a single comma/space-separated token.)

## File formats

- **mzML** (gzipped or plain) — MS1 spectra, parsed with pymzml
- **Percolator** `target.psms.txt` — Crux Percolator or standalone Percolator
  output (auto-detected by header)
- **mzTab** (quantms / OpenMS) — the second ID intake path
- **Output** — tab-delimited `<sample>_riana.txt` with one row per PSM and
  one column per integrated isotopomer, plus a provenance header

## Pipeline

Riana is orchestration-agnostic: search + identification are owned upstream
(e.g. quantms for DDA, DIA-NN for DIA), and Riana is a linear `integrate → fit`
chain you compose into whatever workflow already runs them. The bundled
Snakemake example was retired in 1.0.0 — drive the subcommands directly, or from
your own workflow manager.

## Citation

If you use Riana in published work, please cite:

> Lau, E. *Riana — Relative Isotope Abundance Analyzer*. Zenodo.
> <https://doi.org/10.5281/zenodo.17613314>

## Contributing

Bug reports and pull requests welcome at
<https://github.com/ed-lau/riana/issues>. See `PROJECT_REVIEW.md` for the
current development roadmap.

## License

MIT — see [LICENSE.md](LICENSE.md).
