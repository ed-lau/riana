"""Smoke tests for the M4 Typer CLI (``riana.cli``).

Covers the surfaces that changed in M4: the single (no ``--engine``) integrate
path, the required ``--coefficients`` on fit, bundled-preset resolution, and the
o18 "pending reimplementation" guard. Numerical parity vs 0.9.0 is gated by
``test_integration_port`` (against the committed golden); these just confirm the
CLI wiring.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from riana import __version__
from riana.cli import app

runner = CliRunner()

SAMPLE1 = Path("tests/data/sample1")
PSMS = SAMPLE1 / "percolator.target.psms.txt"
ONE_TIMEPOINT = SAMPLE1 / "sample1_riana.v0_9_0.txt"


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout


def test_no_args_shows_help():
    result = runner.invoke(app, [])
    # no_args_is_help -> usage, exit 0.
    assert "Usage" in result.stdout


def test_integrate_ms2_runs_and_writes_output(tmp_path):
    result = runner.invoke(app, [
        "integrate", str(SAMPLE1), str(PSMS),
        "-s", "sample1", "-o", str(tmp_path),
        "-i", "0 6", "-q", "1.0",
        "--peak-rt", "ms2", "--integration-half-width", "1.0",
        "-m", "50", "-t", "1",
    ])
    assert result.exit_code == 0, result.output
    out = tmp_path / "sample1_riana.txt"
    assert out.exists()
    df = pd.read_csv(out, sep="\t", index_col=0, comment="#")
    for c in ("iso0", "iso6", "concat", "scan", "sample"):
        assert c in df.columns


def test_integrate_sample_must_end_with_digit(tmp_path):
    result = runner.invoke(app, [
        "integrate", str(SAMPLE1), str(PSMS),
        "-s", "sample", "-o", str(tmp_path), "-i", "0 6",
    ])
    assert result.exit_code != 0
    assert "must end with a number" in result.output


# --- M6a SDRF / manifest paths ----------------------------------------------

_BSA_MZTAB = (
    "MTD\tmzTab-version\t1.0.0\n"
    "MTD\tms_run[1]-location\tfile://20180216_BSA.mzML\n\n"
    "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\t"
    "search_engine\tsearch_engine_score[1]\tmodifications\tretention_time\t"
    "charge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\t"
    "start\tend\topt_global_Posterior_Error_Probability_score\t"
    "opt_global_q-value\topt_global_cv_MS:1002217_decoy_peptide\t"
    "opt_global_cv_MS:1000889_peptidoform_sequence\n"
    # retention_time 541.1 s ≈ the BSA mzML's MS1 RT at scan 4408, so the intake
    # scan↔RT guard (DataError on a >2 min median offset) reconciles and passes.
    "PSM\tRHPEYAVSVLLR\t0\tsp|P02769|ALBU_BOVIN\t1\tdb\tnull\t[, , dummy, 1]\t"
    "0.001\tnull\t541.1\t3\t470.6\t470.6\t"
    "ms_run[1]:controllerType=0 controllerNumber=1 scan=4408\tK\tR\t1\t12\t"
    "0.01\t0.001\t0\tRHPEYAVSVLLR\n"
)

_BSA_SDRF = (
    "source name\tcharacteristics[biological replicate]\t"
    "characteristics[precursor enrichment]\tcharacteristics[labeling time]\t"
    "comment[data file]\tcomment[fraction identifier]\t"
    "comment[technical replicate]\t"
    "comment[proteomics data acquisition method]\tfactor value[condition]\n"
    "bsa_t0\t1\t0.06\t0 day\t20180216_BSA.mzML\t1\t1\t"
    "NT=data-dependent acquisition;AC=MS:1003221\tcontrol\n"
)


def test_integrate_sdrf_writes_per_run_output_and_manifest(tmp_path):
    mztab = tmp_path / "bsa.mzTab"
    mztab.write_text(_BSA_MZTAB)
    sdrf = tmp_path / "bsa.sdrf.tsv"
    sdrf.write_text(_BSA_SDRF)
    result = runner.invoke(app, [
        "integrate", str(SAMPLE1), str(mztab),
        "--sdrf", str(sdrf), "-o", str(tmp_path),
        "--peak-rt", "ms2", "--integration-half-width", "1.0", "-m", "50",
    ])
    assert result.exit_code == 0, result.output
    # one <stem>_riana.txt per run + a manifest
    assert (tmp_path / "20180216_BSA_riana.txt").exists()
    assert (tmp_path / "riana_manifest.tsv").exists()
    header = (tmp_path / "20180216_BSA_riana.txt").read_text()
    assert "# sample bsa_t0" in header and "# condition control" in header


def test_fit_requires_exactly_one_input_source(tmp_path):
    # neither positional files nor --manifest
    result = runner.invoke(app, ["fit", "--coefficients", "commerford"])
    assert result.exit_code != 0
    assert "either positional" in result.output


def test_fit_requires_coefficients_for_hw():
    result = runner.invoke(app, ["fit", str(ONE_TIMEPOINT)])
    assert result.exit_code != 0
    assert "--coefficients is required" in result.output


def test_fit_label_o18_is_recognised_but_errors():
    result = runner.invoke(
        app, ["fit", str(ONE_TIMEPOINT), "--label", "o18"]
    )
    assert result.exit_code != 0
    # The engine-side guard message (raised in fit_run), surfaced by Typer.
    msg = (result.output or "") + str(result.exception or "")
    assert "o18" in msg and "hw" in msg


def test_fit_resolves_bundled_preset(tmp_path):
    """--coefficients commerford resolves the bundled table and reaches fit_run.

    The committed single-timepoint golden was integrated with `--iso 0 6`, so
    fit_run's canonical-isotopomer guard fires — which proves the preset loaded
    and the pipeline ran past the `--coefficients` BadParameter (a missing table
    would have errored earlier, before any coefficients were loaded).
    """
    result = runner.invoke(app, [
        "fit", str(ONE_TIMEPOINT),
        "--coefficients", "commerford",
        "-o", str(tmp_path),
    ])
    assert result.exit_code != 0
    msg = (result.output or "") + str(result.exception or "")
    assert "coefficients from commerford" in msg  # preset resolved + loaded
    assert "isotopomers" in msg                    # reached fit_run's guard


def test_fit_and_rollup_via_manifest_chain(tmp_path):
    """fit --manifest writes next to the manifest (ignoring -o) + records fit
    rows; rollup --manifest reads them and records the protein row."""
    from riana.io.manifest import append_manifest, read_manifest
    from tests.test_pipeline import (
        _coeffs,
        _integrate_rows_from_dfs,
        _make_timepoint_dfs,
    )

    coeffs = _coeffs()
    rows = _integrate_rows_from_dfs(tmp_path, _make_timepoint_dfs(coeffs),
                                    condition="control")
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)
    coeff_csv = tmp_path / "coeffs.csv"
    pd.DataFrame({"amino_acid": list(coeffs), "coefficient": list(coeffs.values())}
                 ).to_csv(coeff_csv, index=False)
    elsewhere = tmp_path / "elsewhere"

    r = runner.invoke(app, [
        "fit", "--manifest", str(mf), "--coefficients", str(coeff_csv),
        "-o", str(elsewhere), "-q", "0.05", "-d", "3"])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "riana_fit_peptides.txt").exists()      # next to manifest
    assert (tmp_path / "riana_fit_fractions.txt").exists()
    assert not (elsewhere / "riana_fit_peptides.txt").exists()  # -o ignored
    assert len(read_manifest(mf, stage="fit")) == 2

    r2 = runner.invoke(app, [
        "rollup", "--manifest", str(mf), "--min-peptides", "1"])
    assert r2.exit_code == 0, r2.output
    assert (tmp_path / "riana_protein.txt").exists()
    assert len(read_manifest(mf, stage="protein")) == 1


def test_rollup_requires_exactly_one_input_source(tmp_path):
    result = runner.invoke(app, ["rollup"])  # neither FIT_DIR nor --manifest
    assert result.exit_code != 0
    assert "either FIT_DIR" in result.output


def test_gui_subcommand_is_registered():
    """`riana gui --help` exits 0 — the subcommand is wired (M4 Phase 2).

    The GUI itself is exercised in test_gui.py; here we only confirm the CLI
    surfaces the command without importing PySide6 (the import is lazy).
    """
    result = runner.invoke(app, ["gui", "--help"])
    assert result.exit_code == 0
    assert "GUI" in result.output
