# -*- coding: utf-8 -*-

"""Tests for the M4 Phase 2 GUI.

Two layers, matching the GUI's deliberate split:

1. **Qt-free worker tests** — :mod:`riana.gui.tasks` imports no PySide6, so these
   run with no display and assert the GUI's integrate path is numerically the
   *same* as the CLI's (against the committed ``sample1`` golden, the same gate
   ``test_integration_port`` uses).
2. **Headless Qt smoke** — guarded by ``importorskip("PySide6")`` and run under
   the ``offscreen`` platform; build the window and assert the form→config
   marshalling routes through the shared ``IntegrationConfig.__post_init__``
   validator (the "CLI and GUI cannot drift" contract).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Qt (if present) must run headless in CI / on a server.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from riana.config import FitConfig, IntegrationConfig
from riana.core.integration import PeptideTrace
from riana.gui.tasks import (
    extract_trace,
    integrate_fraction,
    plan_sdrf_integration,
    read_psms,
    run_fit,
    run_fit_manifest,
    run_rollup,
)
from riana.io.percolator import file_indices, fraction_psms
from riana.records import Chromatogram, PSMRecord
from tests.test_fitting import (
    _coefficients_for_target_spep,
    _make_synthetic_dfs,
    _spep_by_seq_from_coefficients,
    _TEST_PEPTIDES,
)

SAMPLE1 = Path("tests/data/sample1")
PSMS = SAMPLE1 / "percolator.target.psms.txt"
MZML = SAMPLE1 / "20180216_BSA.mzML.gz"
GOLDEN = SAMPLE1 / "sample1_riana.v0_9_0.txt"

# The 0.9.0-parity config (ms2 + whole-window), matching test_integration_port.
_MS2_CONFIG = IntegrationConfig(
    sample="sample1", isotopomers=(0, 6), q_value=1.0,
    extraction_half_width=1.0, mass_tol_ppm=50, threads=1, forced_mods=(0.0,),
    peak_rt="ms2", baseline_method="none",
)


# --- Qt-free worker tests ---------------------------------------------------- #


def test_read_psms_returns_records():
    psms = read_psms(str(PSMS), "sample1", ())
    assert psms, "no PSMs parsed"
    assert all(isinstance(p, PSMRecord) for p in psms)
    assert all(p.sample == "sample1" for p in psms)


def test_integrate_fraction_matches_cli_golden():
    """The GUI worker's integrate output == the committed CLI/0.9.0 golden.

    Proves the GUI path goes through the identical numeric core, so the two
    surfaces cannot diverge.
    """
    psms = read_psms(str(PSMS), "sample1", ())
    indices = file_indices(psms)
    assert indices == [0], "sample1 should be a single fraction"
    fraction = fraction_psms(psms, 0)

    df, drift = integrate_fraction(_MS2_CONFIG, fraction, str(MZML), "20180216_BSA")
    assert drift is not None and drift.n > 0

    legacy = pd.read_csv(GOLDEN, sep="\t", index_col=0, comment="#")
    merged = df.merge(
        legacy[["concat", "scan", "iso0", "iso6"]].rename(
            columns={"iso0": "iso0_legacy", "iso6": "iso6_legacy"}
        ),
        on=["concat", "scan"],
        how="inner",
    )
    assert len(merged) == len(legacy), "rowcount drift vs golden"
    for new_col, old_col in [("iso0", "iso0_legacy"), ("iso6", "iso6_legacy")]:
        np.testing.assert_allclose(
            merged[new_col].to_numpy(), merged[old_col].to_numpy(),
            rtol=1e-3, atol=1e-6,
            err_msg=f"GUI worker {new_col} drifted vs CLI golden",
        )


def test_extract_trace_returns_chromatograms():
    from riana.core.integration import extract_peptide_trace
    from riana.io.mzml import IndexedMzML

    psms = read_psms(str(PSMS), "sample1", ())
    # Find a peptide that has signal using a single shared reader (cheap), so the
    # path-based worker only has to re-index the mzML once below.
    good_psm = None
    with IndexedMzML(MZML) as mzml:
        for psm in psms[:8]:
            if extract_peptide_trace(_MS2_CONFIG, psm, mzml).chromatograms[0].rt:
                good_psm = psm
                break
    assert good_psm is not None, "no peptide yielded an extractable trace"

    trace = extract_trace(_MS2_CONFIG, good_psm, str(MZML))
    assert isinstance(trace, PeptideTrace)
    assert set(trace.chromatograms) == {0, 6}
    # ms2 path integrates the whole extraction → no narrowed window.
    assert trace.window is None
    for chrom in trace.chromatograms.values():
        assert isinstance(chrom, Chromatogram)
        assert len(chrom.rt) == len(chrom.intensity) > 0


def test_run_fit_worker_fits_synthetic_series(tmp_path):
    """The Model-tab worker reads per-timepoint files + coefficients and fits.

    Uses the same synthetic D₂O series test_fitting validates, written out as
    `_riana.txt` files so the file-reading + coefficient-loading worker path
    (not just fit_run) is exercised.
    """
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)

    paths = []
    for i, df in enumerate(dfs):
        p = tmp_path / f"time{i}_riana.txt"
        df.to_csv(p, sep="\t", index=False)
        paths.append(str(p))
    coeff_csv = tmp_path / "coeffs.csv"
    pd.DataFrame({"amino_acid": list(coeffs), "coefficient": list(coeffs.values())}
                 ).to_csv(coeff_csv, index=False)

    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06, threads=1)
    result = run_fit(config, paths, str(coeff_csv))

    assert len(result) == len(_TEST_PEPTIDES)
    for col in ("k_deg", "R_squared", "t", "fs", "spep"):
        assert col in result.columns
    assert int(result["k_deg"].notna().sum()) >= 1
    # Per-peptide t/fs are retained (these drive the GUI curve plot).
    converged = result[result["k_deg"].notna()].iloc[0]
    assert len(converged["t"]) == len(converged["fs"]) >= config.depth


def test_run_rollup_worker_rolls_fit_outputs_to_proteins(tmp_path):
    """The Protein-tab worker reads a fit dir and produces a protein table.

    Exercises the file-reading worker path (not just rollup_proteins): write the
    `riana fit` outputs, then roll up via the same worker the GUI submits.
    """
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    paths = []
    for i, df in enumerate(dfs):
        p = tmp_path / f"time{i}_riana.txt"
        df.to_csv(p, sep="\t", index=False)
        paths.append(str(p))
    coeff_csv = tmp_path / "coeffs.csv"
    pd.DataFrame({"amino_acid": list(coeffs), "coefficient": list(coeffs.values())}
                 ).to_csv(coeff_csv, index=False)

    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06, threads=1)
    fit_df = run_fit(config, paths, str(coeff_csv))
    fit_dir = tmp_path / "fit"
    fit_dir.mkdir()
    fit_df.to_csv(fit_dir / "riana_fit_peptides.txt", sep="\t", index=True)
    fit_df.attrs["fractions_long"].to_csv(
        fit_dir / "riana_fit_fractions.txt", sep="\t", index=False)

    proteins, points = run_rollup(
        str(fit_dir), "simple", 0.5, 0.05, 10.0, "unique", 1, 3)
    assert {"protein", "k_deg_median", "k_deg_refit"} <= set(proteins.columns)
    # Each synthetic peptide maps to its own protein (proteotypic) -> 5 proteins.
    assert len(proteins) == len(_TEST_PEPTIDES)
    assert int(proteins["k_deg_refit"].notna().sum()) >= 1
    # The collapsed (t, θ) points behind each refit ride alongside (the curve).
    assert len(points) >= 1
    a_t, a_fs = next(iter(points.values()))
    assert len(a_t) == len(a_fs) > 0


@pytest.mark.skipif(not MZML.exists(), reason="sample1 BSA mzML missing")
def test_plan_sdrf_integration_worker_builds_runtasks(tmp_path):
    """The Integrate-tab SDRF planning worker resolves SDRF+mzTab -> RunTasks."""
    from tests.test_pipeline import _BSA_MZTAB, _BSA_SDRF

    mztab = tmp_path / "bsa.mzTab"
    mztab.write_text(_BSA_MZTAB)
    sdrf = tmp_path / "bsa.sdrf.tsv"
    sdrf.write_text(_BSA_SDRF)
    config = IntegrationConfig(
        isotopomers=(0, 1, 2, 3, 4, 5), mass_tol_ppm=25, peak_rt="ms2",
        integration_half_width=1.0, extraction_half_width=1.0,
    )
    tasks = plan_sdrf_integration(config, str(sdrf), str(SAMPLE1), str(mztab))
    assert len(tasks) == 1
    task = tasks[0]
    assert task.stem == "20180216_BSA"
    assert task.identity.sample == "bsa_t0"
    assert task.psms and task.mzml_path.endswith(".mzML.gz")


def test_run_fit_manifest_worker_fits_from_manifest(tmp_path):
    """The Model-tab manifest worker fits curves from a manifest via fit_project."""
    from riana.io.manifest import append_manifest
    from tests.test_pipeline import (
        _coeffs,
        _integrate_rows_from_dfs,
        _make_timepoint_dfs,
    )

    coeffs = _coeffs()
    rows = _integrate_rows_from_dfs(
        tmp_path, _make_timepoint_dfs(coeffs), condition="control")
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)
    coeff_csv = tmp_path / "coeffs.csv"
    pd.DataFrame({"amino_acid": list(coeffs), "coefficient": list(coeffs.values())}
                 ).to_csv(coeff_csv, index=False)

    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06, threads=1)
    result = run_fit_manifest(config, str(mf), str(coeff_csv))
    assert {"k_deg", "condition", "experiment"} <= set(result.columns)
    assert set(result["condition"]) == {"control"}
    # M5 long table rides along for the GUI to write riana_fit_fractions.txt.
    assert "fractions_long" in result.attrs and not result.attrs["fractions_long"].empty


# --- Headless Qt smoke ------------------------------------------------------- #

pytest.importorskip("PySide6")


@pytest.fixture
def main_window(qtbot):
    """A constructed MainWindow with a cheap, unused executor.

    The smoke tests only build widgets and marshal config — they never submit to
    the pool — so a ThreadPoolExecutor stand-in avoids process-spawn overhead.
    """
    from concurrent.futures import ThreadPoolExecutor

    from riana.gui.main_window import MainWindow

    pool = ThreadPoolExecutor(max_workers=1)
    window = MainWindow(pool=pool, default_threads=1)
    qtbot.addWidget(window)
    yield window
    pool.shutdown(wait=False)


def test_window_has_integrate_model_and_protein_tabs(main_window):
    titles = [main_window.tabs.tabText(i) for i in range(main_window.tabs.count())]
    assert titles == ["Integrate", "Model", "Protein"]


def test_protein_tab_build_params_defaults(main_window):
    params = main_window.protein_tab.build_params()
    assert params["parsimony"] == "unique"
    assert params["model"] == "simple"
    assert params["min_peptides"] == 2
    assert params["min_points"] == 3
    assert params["min_r2"] is None  # 0 on the spin -> gate off


def test_protein_tab_plots_refit_curve_on_row_selection(main_window):
    """Selecting a protein row draws the collapsed points + the refit curve."""
    from pyqtgraph import PlotDataItem

    tab = main_window.protein_tab
    result = pd.DataFrame({
        "experiment": [""], "condition": [""], "protein": ["P1"],
        "n_peptides": [3], "k_deg_median": [0.40],
        "k_median_lo": [0.3], "k_median_hi": [0.5],
        "n_points": [4], "k_deg_refit": [0.42],
        "k_refit_lo": [0.3], "k_refit_hi": [0.5], "R_squared_refit": [0.98],
    })
    tab._result_df = result
    tab._points = {("", "", "P1"): ([0.0, 1.0, 2.0, 3.0],
                                    [0.0, 0.3, 0.55, 0.7])}
    tab._last_params = tab.build_params()
    tab.model.set_dataframe(result)

    tab.table.setCurrentIndex(tab.model.index(0, 0))

    curves = [it for it in tab.curve.plot.getPlotItem().items
              if isinstance(it, PlotDataItem)]
    assert len(curves) == 2  # collapsed points scatter + fitted refit line
    assert tab.curve.plot.getPlotItem().titleLabel.text == "P1"


def test_build_config_defaults_round_trip(main_window):
    cfg = main_window.integrate_tab.build_config()
    assert isinstance(cfg, IntegrationConfig)
    assert cfg.peak_rt == "apex"
    assert cfg.isotopomers == (0, 1, 2, 3, 4, 5)
    # apex default: ehw = integration_half_width (0.15) + 0.33 apex offset.
    assert cfg.extraction_half_width == pytest.approx(0.48)


def test_integrate_tab_has_sdrf_and_workers(main_window):
    tab = main_window.integrate_tab
    assert tab.sdrf_edit.text() == ""              # SDRF path (optional) wired
    assert tab.workers_spin.value() == 1           # cross-file workers control


def test_model_tab_has_manifest_field(main_window):
    assert main_window.model_tab.manifest_edit.text() == ""


def test_build_config_surfaces_post_init_validation(main_window):
    """A bad widget value raises the *same* ValueError the CLI surfaces."""
    tab = main_window.integrate_tab
    tab.ihw_edit.setText("-1")  # negative half-width → __post_init__ rejects
    with pytest.raises(ValueError):
        tab.build_config()


def test_build_config_rejects_empty_isotopomers(main_window):
    tab = main_window.integrate_tab
    tab.iso_edit.setText("")
    with pytest.raises(ValueError):
        tab.build_config()


def test_model_tab_build_config_defaults(main_window):
    cfg = main_window.model_tab.build_config()
    assert isinstance(cfg, FitConfig)
    assert cfg.model == "simple"
    assert cfg.label == "hw"
    assert cfg.depth == 3
    assert cfg.ria_max == pytest.approx(0.06)
    # The coefficients combo defaults to the bundled literature preset.
    assert main_window.model_tab.coeff_combo.currentText() == "commerford"


def test_model_tab_plots_fitted_curve_on_row_selection(main_window):
    """Selecting a result row draws the observed points + fitted curve."""
    from pyqtgraph import PlotDataItem

    tab = main_window.model_tab
    rdf = pd.DataFrame(
        {
            "t": [[0.0, 1.0, 2.0, 3.0]],
            "fs": [[0.0, 0.3, 0.55, 0.7]],
            "k_deg": [0.4], "R_squared": [0.98], "sd": [0.02], "spep": [8.0],
            "ci_lo": [0.36], "ci_hi": [0.44], "protein id": ["sp|X|T"],
        },
        index=pd.Index(["PEPTIDEK_2"], name="concat"),
    )
    tab._result_df = rdf
    tab._last_config = tab.build_config()
    tab._populate_results(rdf)

    # Selecting row 0 fires currentRowChanged → the curve handler.
    tab.table.setCurrentIndex(tab.model.index(0, 0))

    curves = [it for it in tab.curve.plot.getPlotItem().items
              if isinstance(it, PlotDataItem)]
    assert len(curves) == 2  # observed scatter + fitted line
    assert tab.curve.plot.getPlotItem().titleLabel.text == "PEPTIDEK_2"
