[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17613314.svg)](https://doi.org/10.5281/zenodo.17613314)

# Riana — Relative Isotope Abundance Analyzer

Riana takes standard mass-spectrometry spectra (mzML) and peptide-spectrum-match
files (Percolator output today; mzTab via quantms in 1.0) and returns mass
isotopomer distributions, e.g. for protein turnover analysis. It also fits
kinetic models (one-exponential, Guan, Fornasiero) to time-series isotopomer
data.

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
    --iso 0 1 2 3 4 5 6 \
    --q_value 0.01 \
    --r_time 0.5 \
    --mass_tol 25 \
    --out ./out
```

### Fit

Fit a kinetic model across timepoints:

```bash
riana fit ./out/time0_riana.txt ./out/time1_riana.txt ./out/time3_riana.txt \
    --model simple \
    --label 1 \
    --ria 0.06 \
    --depth 3 \
    --out ./out
```

See `riana integrate --help` and `riana fit --help` for the full argument set,
or the [online docs](https://ed-lau.github.io/riana/) for tutorials.

## File formats

- **mzML** (gzipped or plain) — MS1 spectra, parsed with pymzml
- **Percolator** `target.psms.txt` — Crux Percolator or standalone Percolator
  output (auto-detected by header)
- **Output** — tab-delimited `<sample>_riana.txt` with one row per PSM and
  one column per integrated isotopomer

mzTab (quantms) intake is planned for the 1.0 release; see `PROJECT_REVIEW.md`.

## Snakemake workflow

A reference Snakemake pipeline (Comet → Percolator → Riana integrate → Riana
fit) is bundled at `workflow/Snakefile`. Edit `config_template.yaml` to point at
your tooling and data, then:

```bash
snakemake -c -s workflow/Snakefile -d out/snakemake_test --configfile config_template.yaml
```

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
