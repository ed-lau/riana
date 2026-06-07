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


def test_gui_subcommand_is_registered():
    """`riana gui --help` exits 0 — the subcommand is wired (M4 Phase 2).

    The GUI itself is exercised in test_gui.py; here we only confirm the CLI
    surfaces the command without importing PySide6 (the import is lazy).
    """
    result = runner.invoke(app, ["gui", "--help"])
    assert result.exit_code == 0
    assert "GUI" in result.output
