[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10614233.svg)](https://doi.org/10.5281/zenodo.10614233)

# Riana — Relative Isotope Abundance Analyzer

Riana measures protein turnover from mass spectrometry experiments with metabolic
heavy water (D₂O) or ¹⁸O water labeling. It integrates each peptide's MS1 isotopomer
envelope from mzML, using identifications from quantms (mzTab, DDA), DIA-NN
(`report.parquet`, DIA) or Percolator; fits turnover rate constants; rolls peptides
up to proteins; and tests turnover differences between conditions, including across
tissues of the same animals. A desktop app (`riana gui`) runs the same steps.

**Documentation: <https://ed-lau.github.io/riana/>**
— [Quick Start](https://ed-lau.github.io/riana/install.html) ·
[How it works](https://ed-lau.github.io/riana/about.html) ·
[CLI Reference](https://ed-lau.github.io/riana/reference.html) ·
[FAQ](https://ed-lau.github.io/riana/faqs.html) ·
[Change Log](https://ed-lau.github.io/riana/change.html)

## Install

Riana requires Python 3.10 or newer.

```bash
pip install riana            # command line
pip install "riana[gui]"     # + the desktop app
pip install "riana[dia]"     # + DIA-NN report.parquet intake
```

For development, from a clone: `pip install -e ".[dev]"`.

## Quickstart

Describe the experiment in an SDRF sample sheet, then:

```bash
riana integrate <mzml_dir> report.mzTab --sdrf samplesheet.sdrf.tsv --workers 4 --out ./out
riana fit --manifest ./out/riana_manifest.tsv --coefficients deberneh_2025_rss
riana rollup --manifest ./out/riana_manifest.tsv --model "linear simple" \
    --reference-condition control
```

The last step reports each protein's turnover rate per condition and the difference
between conditions (Δk) with adjusted p-values. See the
[Quick Start](https://ed-lau.github.io/riana/install.html) for the SDRF columns,
multi-tissue designs and the outputs, and `riana <command> --help` for every option.

## Citation

If you use Riana in published work, please consider citing:

> Alamillo, L. et al. *Protocol to Measure Protein Half-Life in Cell Culture Using Heavy Water*. STAR Protoc. 2026
> <https://doi.org/10.1016/j.xpro.2026.104426>

> Currie, J. et al. *Improved Method to Determine Protein Turnover Rates with Heavy Water Labeling by Mass Isotopomer Ratio Selection*. J Proteome Res. 2025
> <https://doi.org/10.1021/acs.jproteome.4c01012>

> Alamillo, L. et al. *Deuterium Labeling Enables Proteome Wide Turnover Kinetics in Cell Culture*. Cell Rep Methods. 2025
> <https://doi.org/10.1016/j.crmeth.2025.101104>

> Hammond, D. et al. *Harmonizing Labeling and Analytical Strategies to Obtain Protein Turnover Rates in Intact Adult Animals*. Mol Cell Proteomics. 2022
> <https://doi.org/10.1016/j.mcpro.2022.100252>

## Contributing

Bug reports and pull requests are welcome at <https://github.com/ed-lau/riana/issues>.
`CHANGELOG.md` records what has shipped and `PROJECT_REVIEW.md` the roadmap.
The CLI Reference page is generated from the command line
(`python riana_website/_gen_reference.py`); a test fails if it is out of date.

## License

MIT — see [LICENSE.md](LICENSE.md).
