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
    extraction_half_width=1.0, mass_tol_ppm=50,
    peak_rt="ms2", baseline_method="none",
)


# --- Qt-free worker tests ---------------------------------------------------- #


def test_read_psms_returns_records():
    psms = read_psms(str(PSMS), "sample1")
    assert psms, "no PSMs parsed"
    assert all(isinstance(p, PSMRecord) for p in psms)
    assert all(p.sample == "sample1" for p in psms)


def test_integrate_fraction_matches_cli_golden():
    """The GUI worker's integrate output == the committed CLI/0.9.0 golden.

    Proves the GUI path goes through the identical numeric core, so the two
    surfaces cannot diverge.
    """
    psms = read_psms(str(PSMS), "sample1")
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

    psms = read_psms(str(PSMS), "sample1")
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
                       ria_max=0.06)
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
                       ria_max=0.06)
    fit_df = run_fit(config, paths, str(coeff_csv))
    fit_dir = tmp_path / "fit"
    fit_dir.mkdir()
    fit_df.to_csv(fit_dir / "riana_fit_peptides.txt", sep="\t", index=True)
    fit_df.attrs["fractions_long"].to_csv(
        fit_dir / "riana_fit_fractions.txt", sep="\t", index=False)

    proteins, points = run_rollup(
        str(fit_dir), "simple", 0.5, 0.05, 10.0, "unique", 1, 3)
    assert {"protein", "method", "k_deg", "peptide_median_k"} <= set(proteins.columns)
    # Each synthetic peptide maps to its own protein (proteotypic) -> 5 proteins.
    assert len(proteins) == len(_TEST_PEPTIDES)
    assert int(proteins["k_deg"].notna().sum()) >= 1
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
                       ria_max=0.06)
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
    window = MainWindow(pool=pool)
    qtbot.addWidget(window)
    yield window
    pool.shutdown(wait=False)


def test_window_has_integrate_model_and_protein_tabs(main_window):
    titles = [main_window.tabs.tabText(i) for i in range(main_window.tabs.count())]
    assert titles == ["Integrate", "Model", "Protein"]


def test_hint_bar_mirrors_widget_tooltip(main_window):
    """The fixed hint area shows a control's tooltip on Enter/FocusIn — no hover-hold."""
    from PySide6.QtCore import QEvent

    w = main_window.integrate_tab.iso_edit
    assert w.toolTip()                                   # has a tooltip to mirror
    main_window._hint_filter.eventFilter(w, QEvent(QEvent.Type.Enter))
    assert main_window.hint_label.text() == " ".join(w.toolTip().split())


def test_key_form_fields_have_tooltips(main_window):
    """Coverage guard: the main-path inputs all carry help (feeds the hint bar)."""
    it, mt, pt = (main_window.integrate_tab, main_window.model_tab,
                  main_window.protein_tab)
    for w in (it.mzml_edit, it.id_edit, it.sdrf_edit, it.out_edit, it.qvalue_spin,
              it.peak_rt_combo, it.ihw_edit,
              mt.manifest_edit, mt.coeff_combo, mt.label_combo, mt.depth_spin,
              mt.ria_spin, mt.qvalue_spin,
              pt.manifest_edit, pt.parsimony_combo, pt.min_peptides_spin):
        assert w.toolTip(), f"{type(w).__name__} is missing a tooltip"


def test_protein_tab_build_params_defaults(main_window):
    params = main_window.protein_tab.build_params()
    assert params["parsimony"] == "unique"
    assert params["model"] == "simple"
    assert params["min_peptides"] == 2
    assert params["min_points"] == 3
    assert params["min_r2"] is None  # 0 on the spin -> gate off


def test_protein_tab_is_on_the_manifest_path(main_window):
    """The rollup tab takes a manifest, like Integrate (SDRF) and Model (manifest).

    The GUI runs the SDRF/manifest project path only; the fit-output-directory
    input is retired (CLI-only). `build_params` must surface the manifest, and the
    old ``fit_dir`` key/field must be gone so a stale reference can't creep back.
    """
    tab = main_window.protein_tab
    assert hasattr(tab, "manifest_edit")
    assert not hasattr(tab, "fit_dir_edit")

    tab.manifest_edit.setText("/proj/riana_manifest.tsv")
    params = tab.build_params()
    assert params["manifest"] == "/proj/riana_manifest.tsv"
    assert "fit_dir" not in params


def test_protein_tab_plots_refit_curve_on_row_selection(main_window):
    """Selecting a protein row draws the collapsed points + the refit curve."""
    from pyqtgraph import PlotDataItem

    tab = main_window.protein_tab
    result = pd.DataFrame({
        "experiment": [""], "condition": [""], "protein": ["P1"],
        "method": ["weighted"], "n_peptides": [3], "n_points": [4],
        "k_deg": [0.42], "ci_lo": [0.3], "ci_hi": [0.5], "R_squared": [0.98],
        "peptide_median_k": [0.40],
    })
    tab._result_df = result
    tab._points = {("", "", "P1"): ([0.0, 1.0, 2.0, 3.0],
                                    [0.0, 0.3, 0.55, 0.7])}
    tab._last_params = tab.build_params()
    tab.model.set_dataframe(result)

    tab.table.setCurrentIndex(tab.model.index(0, 0))

    curves = [it for it in tab.curve.plot.items
              if isinstance(it, PlotDataItem)]
    assert len(curves) == 2  # collapsed points scatter + fitted refit line
    assert tab.curve.plot.titleLabel.text == "P1"


def test_protein_tab_linear_model_toggles_and_params(main_window):
    """Selecting 'linear simple' swaps the ODE rate knobs for the φ knobs and
    surfaces them in build_params."""
    tab = main_window.protein_tab
    tab.model_combo.setCurrentText("linear simple")
    assert tab.phi_limit_spin.isVisibleTo(tab)
    assert not tab.kp_spin.isVisibleTo(tab)
    tab.reference_edit.setText("control")
    params = tab.build_params()
    assert params["model"] == "linear simple"
    assert params["phi_limit"] == -4.0
    assert params["reference_condition"] == "control"
    assert "workers" in params


def test_protein_tab_linear_plots_phi_space(main_window):
    """A linear-simple protein row overlays both conditions in φ-space (each a
    points scatter + a through-origin k line) with the Δk in the title."""
    from pyqtgraph import PlotDataItem

    tab = main_window.protein_tab
    tab.model_combo.setCurrentText("linear simple")
    result = pd.DataFrame({
        "experiment": ["e", "e"], "condition": ["control", "atrium"],
        "protein": ["P1", "P1"], "method": ["linear simple"] * 2,
        "n_peptides": [3, 3], "n_points": [5, 5],
        "k_deg": [0.05, 0.10], "ci_lo": [0.04, 0.09], "ci_hi": [0.06, 0.11],
        "R_squared": [0.98, 0.98], "peptide_median_k": [0.05, 0.10],
        "delta_k": [0.05, 0.05], "delta_k_se": [0.005, 0.005],
        "delta_k_p": [1e-9, 1e-9], "delta_k_p_adj": [1e-9, 1e-9],
    })
    tab._result_df = result
    tab._points = {
        ("e", "control", "P1"): ([0.0, 1, 2, 4, 8], [0.0, 0.05, 0.1, 0.18, 0.3]),
        ("e", "atrium", "P1"): ([0.0, 1, 2, 4, 8], [0.0, 0.1, 0.19, 0.33, 0.55]),
    }
    tab._last_params = tab.build_params()
    tab.model.set_dataframe(result)

    tab.table.setCurrentIndex(tab.model.index(0, 0))

    curves = [it for it in tab.curve.plot.items
              if isinstance(it, PlotDataItem)]
    # Two conditions × (points scatter + k line) = 4 data items.
    assert len(curves) >= 4
    assert "Δk" in tab.curve.plot.titleLabel.text


def test_curve_view_draws_ci_ribbon(main_window):
    """plot_fit / plot_linear shade a CI ribbon (FillBetweenItem) when the k CI
    bounds are given, without adding extra data series (PlotDataItem)."""
    from pyqtgraph import FillBetweenItem, PlotDataItem

    cv = main_window.protein_tab.curve
    cv.plot_fit("PEP", [0, 1, 2, 4, 8], [0, 0.3, 0.5, 0.7, 0.85], 0.2,
                "simple", {}, ci_lo=0.15, ci_hi=0.25)
    items = cv.plot.items
    assert sum(isinstance(i, FillBetweenItem) for i in items) == 1
    assert sum(isinstance(i, PlotDataItem) for i in items) == 2  # points + fit line

    cv.plot_linear("P1", {
        "control": ([0, 1, 2, 4, 8], [0, .05, .1, .18, .3], 0.05, 0.04, 0.06),
        "atrium": ([0, 1, 2, 4, 8], [0, .1, .19, .33, .55], 0.10, 0.09, 0.11),
    }, phi_limit=-4.0, delta_k=0.05)
    items = cv.plot.items
    assert sum(isinstance(i, FillBetweenItem) for i in items) == 2  # one per condition


def test_curve_view_splits_mbr_and_folded_points(main_window):
    """plot_fit draws direct / MBR / chemical-fold (Met-Ox, future TMT) points as
    separate series so each provenance is visually distinct."""
    from pyqtgraph import PlotDataItem

    cv = main_window.protein_tab.curve
    t = [0, 1, 2, 4]
    fs = [0.0, 0.3, 0.5, 0.7]
    cv.plot_fit("PEP", t, fs, 0.2, "simple", {},
                evidence=["q_value", "mbr", "q_value", "q_value"],
                metox=[False, False, True, False])
    series = [i for i in cv.plot.items if isinstance(i, PlotDataItem)]
    # direct (2 pts) + MBR (1) + folded (1) + fitted line = 4 PlotDataItems.
    assert len(series) == 4
    names = {i.name() for i in series}
    assert any(n and n.startswith("MBR") for n in names)
    assert any(n and n.startswith("folded") for n in names)


def test_build_config_defaults_round_trip(main_window):
    cfg = main_window.integrate_tab.build_config()
    assert isinstance(cfg, IntegrationConfig)
    assert cfg.peak_rt == "apex"
    assert cfg.isotopomers == (0, 1, 2, 3, 4, 5)
    # apex default: ehw = integration_half_width (0.15) + 0.33 apex offset.
    assert cfg.extraction_half_width == pytest.approx(0.48)
    # The Advanced group is collapsed by default and every dial sits at its
    # IntegrationConfig default — the GUI builds the *same* config the CLI does.
    default = IntegrationConfig()
    for field in ("prominence_k", "width_rel_height", "apex_n_consensus",
                  "smoothing", "smoothing_polyorder", "mass_difference",
                  "ppm_alert", "apex_search_half_width", "write_intensities",
                  "check_scan_id", "scan_precursor_tol_ppm", "mbr",
                  "mbr_min_donor_runs", "mbr_donor_q", "mbr_min_snr",
                  "mbr_min_scans"):
        assert getattr(cfg, field) == getattr(default, field), field


def test_build_config_advanced_widgets_flow_through(main_window):
    """Every Advanced dial reaches the frozen config — the no-drift contract for
    the knobs that used to be CLI-only or hidden entirely."""
    tab = main_window.integrate_tab
    tab.peak_rt_combo.setCurrentText("consensus")
    tab.prominence_spin.setValue(5.0)
    tab.width_rel_spin.setValue(0.5)
    tab.apex_n_spin.setValue(3)
    tab.apex_search_spin.setValue(0.4)
    tab.smoothing_combo.setCurrentText("7")
    tab.smoothing_poly_spin.setValue(3)
    tab.mass_diff_spin.setValue(1.5)
    tab.ppm_alert_spin.setValue(12.0)
    tab.write_intensities_check.setChecked(True)
    tab.id_check.setChecked(False)            # guard off
    tab.precursor_tol_spin.setValue(5.0)
    tab.ext_override_check.setChecked(True)
    tab.ext_spin.setValue(0.9)
    tab.mbr_box.setChecked(True)
    tab.mbr_donor_runs_spin.setValue(3)
    tab.mbr_donor_q_spin.setValue(0.005)
    tab.mbr_snr_spin.setValue(6.0)
    tab.mbr_scans_spin.setValue(5)

    cfg = tab.build_config()
    assert cfg.prominence_k == 5.0
    assert cfg.width_rel_height == 0.5
    assert cfg.apex_n_consensus == 3
    assert cfg.apex_search_half_width == 0.4
    assert cfg.smoothing == 7 and cfg.smoothing_polyorder == 3
    assert cfg.mass_difference == 1.5
    assert cfg.ppm_alert == 12.0
    assert cfg.write_intensities is True
    assert cfg.check_scan_id is False and cfg.scan_precursor_tol_ppm == 5.0
    assert cfg.extraction_half_width == 0.9   # explicit override wins over derive
    assert cfg.mbr is True
    assert cfg.mbr_min_donor_runs == 3 and cfg.mbr_donor_q == 0.005
    assert cfg.mbr_min_snr == 6.0 and cfg.mbr_min_scans == 5


def test_build_config_iso_auto_single_int_and_list(main_window):
    """--iso exposure: 'auto' sets adaptive_iso (+ ria_max from the spin); a single
    index N = the iso0..isoN range; an explicit list stays non-contiguous."""
    tab = main_window.integrate_tab
    tab.iso_edit.setText("auto")
    tab.ria_spin.setValue(0.046)
    cfg = tab.build_config()
    assert cfg.adaptive_iso is True
    assert cfg.ria_max == pytest.approx(0.046)

    tab.iso_edit.setText("3")
    cfg2 = tab.build_config()
    assert cfg2.adaptive_iso is False
    assert cfg2.isotopomers == (0, 1, 2, 3)

    tab.iso_edit.setText("0 6")
    assert tab.build_config().isotopomers == (0, 6)


def test_integrate_tab_has_sdrf_and_workers(main_window):
    tab = main_window.integrate_tab
    assert tab.sdrf_edit.text() == ""              # SDRF path (optional) wired
    assert tab.workers_spin.value() == 1           # cross-file workers control
    assert tab.mass_tol_spin.value() == 10         # default matches the config


def test_integrate_tab_reflects_sdrf_mass_tolerance(main_window, tmp_path):
    """Picking an SDRF reflects its precursor mass tolerance in the spinbox."""
    from tests.test_pipeline import _BSA_SDRF

    header, row = _BSA_SDRF.strip().split("\n")
    sdrf = tmp_path / "s.sdrf.tsv"
    sdrf.write_text(
        header + "\tcomment[precursor mass tolerance]\n" + row + "\t25 ppm\n")
    tab = main_window.integrate_tab
    tab.sdrf_edit.setText(str(sdrf))
    tab._resolve_sdrf_mass_tol()
    assert tab.mass_tol_spin.value() == 25


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
    assert main_window.model_tab.coeff_combo.currentText() == "deberneh_2025_rss"
    # --fs defaults to the full envelope (no limited-isotopomer scoring).
    assert cfg.score_channels is None and cfg.fs_auto is False


def test_model_tab_build_config_fs_scoring(main_window):
    """--fs exposure: 'full envelope' -> None; 'auto' -> fs_auto; 'iso0-N' ->
    score_channels = N+1 (the same count the CLI builds)."""
    tab = main_window.model_tab
    tab.fs_combo.setCurrentText("auto")
    cfg = tab.build_config()
    assert cfg.fs_auto is True and cfg.score_channels is None

    tab.fs_combo.setCurrentText("iso0-3")
    cfg2 = tab.build_config()
    assert cfg2.score_channels == 4 and cfg2.fs_auto is False

    tab.fs_combo.setCurrentText("full envelope")
    cfg3 = tab.build_config()
    assert cfg3.score_channels is None and cfg3.fs_auto is False


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

    curves = [it for it in tab.curve.plot.items
              if isinstance(it, PlotDataItem)]
    assert len(curves) == 2  # observed scatter + fitted line
    # header now carries the protein id alongside the peptide concat (for export).
    title = tab.curve.plot.titleLabel.text
    assert title.startswith("PEPTIDEK_2") and "sp|X|T" in title


def test_curve_view_plot_dmass_draws_per_channel_series(main_window):
    """plot_dmass draws one series per channel; spacing drops iso0 (≡0), absolute
    keeps it; empty arrays fall back to the placeholder."""
    from pyqtgraph import PlotDataItem

    cv = main_window.protein_tab.curve
    t = [0.0, 1.0, 2.0, 4.0]
    # 4 points × 4 channels (mDa). iso0 spacing is 0 by construction.
    dmass = [[0.0, 1.0, 2.0, 3.0], [0.0, 2.0, 4.0, 6.0],
             [0.0, 3.0, 6.0, 9.0], [0.0, 4.0, 8.0, 12.0]]
    dspacing = [[0.0, 0.5, 1.0, 1.5], [0.0, 1.0, 2.0, 3.0],
                [0.0, 1.5, 3.0, 4.5], [0.0, 2.0, 4.0, 6.0]]

    cv.plot_dmass("PEP", t, dmass, dspacing, mode="spacing")
    series = [i for i in cv.plot.items if isinstance(i, PlotDataItem)]
    assert {i.name() for i in series} == {"iso1", "iso2", "iso3"}  # iso0 dropped

    cv.plot_dmass("PEP", t, dmass, dspacing, mode="mass")
    series = [i for i in cv.plot.items if isinstance(i, PlotDataItem)]
    assert {i.name() for i in series} == {"iso0", "iso1", "iso2", "iso3"}

    # No mass-accuracy data -> placeholder, no series.
    cv.plot_dmass("PEP", t, [], [], mode="spacing")
    series = [i for i in cv.plot.items if isinstance(i, PlotDataItem)]
    assert not series

    # Out-of-order input (file/concat order, not ascending x) must be sorted before
    # the connecting line is drawn — else the line zig-zags. x=1.0 arrives early.
    t_unsorted = [0.0, 1.0, 0.25, 0.5]
    ds_unsorted = [[0.0, 0.1], [0.0, 0.4], [0.0, 0.2], [0.0, 0.3]]
    cv.plot_dmass("PEP", t_unsorted, ds_unsorted, ds_unsorted, mode="spacing")
    line = next(i for i in cv.plot.items if isinstance(i, PlotDataItem))
    xs = list(line.getData()[0])
    assert xs == sorted(xs)  # ascending x, no zig-zag


def test_curve_view_fit_overlays_fs_ds(main_window):
    """plot_fit overlays fs_ds as its own series (skipping NaN) when given."""
    from pyqtgraph import PlotDataItem

    cv = main_window.protein_tab.curve
    t = [0.0, 1.0, 2.0, 4.0]
    fs = [0.0, 0.3, 0.5, 0.7]
    fs_ds = [0.0, 0.35, float("nan"), 0.66]  # one NaN point is dropped
    cv.plot_fit("PEP", t, fs, 0.2, "simple", {}, fs_ds=fs_ds)
    series = {i.name(): i for i in cv.plot.items if isinstance(i, PlotDataItem)}
    fsds_series = [n for n in series if n and n.startswith("fs_ds")]
    assert fsds_series, "fs_ds overlay series missing"
    xs, _ = series[fsds_series[0]].getData()
    assert len(xs) == 3  # the NaN point is skipped


def test_model_tab_show_fsds_toggle_overlays(main_window):
    """The Fit-view 'show fs_ds' checkbox overlays the second estimate on selection."""
    from pyqtgraph import PlotDataItem

    tab = main_window.model_tab
    rdf = pd.DataFrame(
        {
            "t": [[0.0, 1.0, 2.0, 3.0]],
            "fs": [[0.0, 0.3, 0.55, 0.7]],
            "fs_ds": [[0.0, 0.33, 0.5, 0.72]],
            "k_deg": [0.4], "R_squared": [0.98], "sd": [0.02], "spep": [8.0],
            "ci_lo": [0.36], "ci_hi": [0.44], "protein id": ["sp|X|T"],
        },
        index=pd.Index(["PEPTIDEK_2"], name="concat"),
    )
    tab._result_df = rdf
    tab._last_config = tab.build_config()
    tab._populate_results(rdf)
    tab.table.setCurrentIndex(tab.model.index(0, 0))

    # Off by default: observed + fit line only.
    assert not tab._show_fsds_check.isChecked()
    n_off = len([i for i in tab.curve.plot.items if isinstance(i, PlotDataItem)])
    tab._show_fsds_check.setChecked(True)
    names = {i.name() for i in tab.curve.plot.items if isinstance(i, PlotDataItem)}
    assert any(n and n.startswith("fs_ds") for n in names)
    assert len([i for i in tab.curve.plot.items
                if isinstance(i, PlotDataItem)]) == n_off + 1


def test_curve_view_plot_dmass_anchor_zeroes_t0(main_window):
    """anchor=True subtracts the unlabelled (t=0) point's per-channel Δ from every
    point, so each channel's series starts at 0."""
    from pyqtgraph import PlotDataItem

    cv = main_window.protein_tab.curve
    t = [0.0, 0.5, 1.0]
    ds = [[0.0, 0.4, 0.8], [0.0, 0.9, 1.6], [0.0, 1.4, 2.4]]  # offset +0.4/+0.8 at t0
    cv.plot_dmass("PEP", t, ds, ds, mode="spacing", anchor=True)
    series = {i.name(): i for i in cv.plot.items if isinstance(i, PlotDataItem)}
    # iso1's first (t=0) value is the anchor → 0 after subtraction.
    xs, ys = series["iso1"].getData()
    assert ys[list(xs).index(0.0)] == pytest.approx(0.0, abs=1e-9)


def test_model_tab_view_toggle_switches_fit_and_dmass(main_window):
    """The Model-tab View toggle re-renders the selected peptide as Fit / Δspacing /
    Δmass without a fresh row selection."""
    from pyqtgraph import PlotDataItem

    tab = main_window.model_tab
    rdf = pd.DataFrame(
        {
            "t": [[0.0, 1.0, 2.0, 3.0]],
            "fs": [[0.0, 0.3, 0.55, 0.7]],
            "dmass": [[[0.0, 1.0, 2.0], [0.0, 2.0, 4.0],
                       [0.0, 3.0, 6.0], [0.0, 4.0, 8.0]]],
            "dspacing": [[[0.0, 0.5, 1.0], [0.0, 1.0, 2.0],
                          [0.0, 1.5, 3.0], [0.0, 2.0, 4.0]]],
            "k_deg": [0.4], "R_squared": [0.98], "sd": [0.02], "spep": [8.0],
            "ci_lo": [0.36], "ci_hi": [0.44], "protein id": ["sp|X|T"],
        },
        index=pd.Index(["PEPTIDEK_2"], name="concat"),
    )
    tab._result_df = rdf
    tab._last_config = tab.build_config()
    tab._populate_results(rdf)
    tab.table.setCurrentIndex(tab.model.index(0, 0))

    # Default view = Fit: observed scatter + fitted line.
    assert tab._fit_radio.isChecked()
    assert len([i for i in tab.curve.plot.items if isinstance(i, PlotDataItem)]) == 2

    # Toggle Δspacing -> per-channel series (iso1, iso2; iso0 dropped).
    tab._spacing_radio.setChecked(True)
    names = {i.name() for i in tab.curve.plot.items if isinstance(i, PlotDataItem)}
    assert names == {"iso1", "iso2"}

    # Toggle Δmass -> all channels including iso0.
    tab._mass_radio.setChecked(True)
    names = {i.name() for i in tab.curve.plot.items if isinstance(i, PlotDataItem)}
    assert names == {"iso0", "iso1", "iso2"}


def test_curve_view_calibration_draws_unit_line(main_window):
    """plot_fit in 'calibration' mode draws the recovery line, a 1:1 ideal
    reference (PlotCurveItem, not a data series), and labels x as mixing proportion."""
    from pyqtgraph import PlotCurveItem, PlotDataItem

    cv = main_window.protein_tab.curve
    f = [0.0, 0.25, 0.5, 0.75]
    fs = [0.0, 0.26, 0.49, 0.74]
    cv.plot_fit("PEP", f, fs, 0.98, "calibration", {})
    series = [i for i in cv.plot.items if isinstance(i, PlotDataItem)]
    # observed scatter + recovery line (no CI here) = 2 data series.
    assert len(series) == 2
    assert any((n := i.name()) and n.startswith("recovery") for i in series)
    # The 1:1 ideal is a PlotCurveItem, deliberately not a data series.
    assert any(isinstance(i, PlotCurveItem) for i in cv.plot.items)
    assert cv.plot.getAxis("bottom").labelText == "Mixing proportion"


def test_model_tab_combo_offers_calibration(main_window):
    items = [main_window.model_tab.model_combo.itemText(i)
             for i in range(main_window.model_tab.model_combo.count())]
    assert "calibration" in items


# --- Track E: sortable tables + graph export --------------------------------- #


def test_dataframe_model_sorts_numeric_by_value():
    """``sort`` orders by the raw numeric column, not the ``:.4g`` display text.

    9 / 10 / 100 sort lexicographically as "10" < "100" < "9"; numerically they
    must come back 9 < 10 < 100. This is why the sort lives in the model over the
    raw frame rather than a proxy comparing display strings.
    """
    from riana.gui.models import DataFrameTableModel
    from PySide6.QtCore import Qt

    model = DataFrameTableModel()
    model.set_dataframe(pd.DataFrame({"k": [10.0, 100.0, 9.0], "tag": ["b", "c", "a"]}))

    model.sort(0, Qt.SortOrder.AscendingOrder)
    assert list(model.dataframe["k"]) == [9.0, 10.0, 100.0]

    model.sort(0, Qt.SortOrder.DescendingOrder)
    assert list(model.dataframe["k"]) == [100.0, 10.0, 9.0]


def test_dataframe_model_data_matches_iat_display():
    """The ndarray-cached ``data()`` renders byte-identically to a ``.iat`` read.

    ``data()`` reads a per-column numpy cache rather than ``DataFrame.iat`` for
    speed on the per-cell-per-repaint hot path; this pins that the cache does not
    change any *rendered* value across the dtypes the result tables carry (float
    with NaN, int, object strings, bool) — before and after a sort reorders the
    backing frame.
    """
    from riana.gui.models import DataFrameTableModel
    from PySide6.QtCore import Qt

    def ref(df, r, c):  # the old iat-based data() body, verbatim
        v = df.iat[r, c]
        return f"{v:.4g}" if isinstance(v, float) else str(v)

    df = pd.DataFrame({
        "k_deg": [0.123456, 1234.5, float("nan"), 9.0],
        "n_points": [3, 12, 7, 100],
        "protein": ["sp|P1|A", "sp|P2|B", "", "sp|P3|C"],
        "converged": [True, False, True, False],
    })
    model = DataFrameTableModel()
    model.set_dataframe(df)

    def assert_all_cells_match():
        for r in range(model.rowCount()):
            for c in range(model.columnCount()):
                got = model.data(model.index(r, c), Qt.ItemDataRole.DisplayRole)
                assert got == ref(model.dataframe, r, c), (r, c, got)

    assert_all_cells_match()
    model.sort(0, Qt.SortOrder.DescendingOrder)  # rebuilds the cache
    assert_all_cells_match()


def test_integrate_scan_spans_match_old_full_scan():
    """The precomputed ``concat -> (min,max)`` map equals the old per-click scan.

    ``_show_chromatogram`` used to derive the scan span with
    ``df[df["concat"] == x]["scan"].astype(int)`` min/max on every selection;
    ``_build_scan_spans`` precomputes the same values once. This pins the
    replacement is behaviour-preserving (incl. duplicate concats and an MBR
    ``scan == -1`` row) and that missing/empty frames degrade to an empty map.
    """
    from riana.gui.integrate_tab import IntegrateTab

    df = pd.DataFrame({
        "concat": ["A_2", "A_2", "B_3", "A_2", "B_3"],
        "scan": [10, -1, 7, 20, 7],  # A_2 carries an MBR -1, as the real frame can
    })
    spans = IntegrateTab._build_scan_spans(df)
    for c in df["concat"].unique():  # reference = the old full-frame scan
        same = df[df["concat"] == c]["scan"].astype(int)
        assert spans[c] == (int(same.min()), int(same.max()))
    assert spans == {"A_2": (-1, 20), "B_3": (7, 7)}

    assert IntegrateTab._build_scan_spans(pd.DataFrame()) == {}
    assert IntegrateTab._build_scan_spans(pd.DataFrame({"concat": ["X"]})) == {}


def test_integrate_tab_filters_and_caps_rows(main_window):
    """The Integrate view is capped to `_MAX_DISPLAY_ROWS` and narrowed by the
    filter box, while the full result is retained off-model.

    A multi-file concat is 10⁵–10⁶ rows; the model must only ever hold a
    filtered + capped view (so sort/selection/memory stay bounded), and the
    filter must reach a peptide the cap left off-screen — matching sequence,
    protein id, or concat, case-insensitively.
    """
    from riana.gui.integrate_tab import IntegrateTab, _MAX_DISPLAY_ROWS

    tab = main_window.integrate_tab
    n = _MAX_DISPLAY_ROWS + 50  # more than the cap for the common sequence
    df = pd.DataFrame({
        "file_idx": [0] * (n + 3),
        "scan": list(range(n)) + [1, 2, 3],
        "sequence": ["PEPTIDEK"] * n + ["RARESEQK"] * 3,
        "protein id": ["sp|P1|COMMON"] * n + ["sp|P2|RARE"] * 3,
        "concat": [f"PEPTIDEK_{i % 4 + 1}" for i in range(n)] + ["RARESEQK_2"] * 3,
    })
    tab._full_df = df
    tab._search_key = tab._build_search_key(df)

    # No filter → capped to the display max, but the label reports the full total.
    tab.filter_edit.setText("")
    tab._apply_filter()
    assert len(tab.model.dataframe) == _MAX_DISPLAY_ROWS
    assert f"{n + 3:,}" in tab.rows_label.text()          # full total, not the cap

    # Filter to the rare protein → 3 matches, all shown, none of the common rows.
    tab.filter_edit.setText("P2")
    tab._apply_filter()
    assert len(tab.model.dataframe) == 3
    assert set(tab.model.dataframe["protein id"]) == {"sp|P2|RARE"}

    # Case-insensitive, and matches the sequence column too.
    tab.filter_edit.setText("rareseq")
    tab._apply_filter()
    assert len(tab.model.dataframe) == 3
    assert set(tab.model.dataframe["sequence"]) == {"RARESEQK"}


def test_all_result_tables_have_sorting_enabled(main_window):
    for tab in (main_window.integrate_tab, main_window.model_tab,
                main_window.protein_tab):
        assert tab.table.isSortingEnabled()


def test_protein_row_selection_maps_through_sorted_order(main_window):
    """After a header sort, view row 0 plots the row now on top — no proxy drift.

    The row handlers index ``dataframe.iloc[row]`` directly, so sorting the
    model's own frame is what keeps the selection→curve mapping correct.
    """
    from PySide6.QtCore import Qt

    tab = main_window.protein_tab
    result = pd.DataFrame({
        "experiment": ["", ""], "condition": ["", ""], "protein": ["P1", "P2"],
        "method": ["weighted", "weighted"], "n_peptides": [3, 3],
        "n_points": [4, 4], "k_deg": [0.2, 0.8], "ci_lo": [0.1, 0.7],
        "ci_hi": [0.3, 0.9], "R_squared": [0.98, 0.97],
        "peptide_median_k": [0.2, 0.8],
    })
    tab._result_df = result
    tab._points = {("", "", "P1"): ([0.0, 1.0], [0.0, 0.2]),
                   ("", "", "P2"): ([0.0, 1.0], [0.0, 0.6])}
    tab._last_params = tab.build_params()
    tab.model.set_dataframe(result)

    k_col = list(result.columns).index("k_deg")
    tab.model.sort(k_col, Qt.SortOrder.DescendingOrder)  # P2 (0.8) to the top
    assert tab.model.dataframe.iloc[0]["protein"] == "P2"

    tab.table.setCurrentIndex(tab.model.index(0, 0))
    assert tab.curve.plot.titleLabel.text == "P2"


def test_curve_view_export_writes_png(main_window, tmp_path):
    from riana.gui.export import export_plot

    tab = main_window.model_tab
    assert not tab.curve.save_button.isEnabled()  # placeholder state
    tab.curve.plot_fit("PEP_2", [0.0, 1.0, 2.0], [0.0, 0.3, 0.5], 0.4,
                       "simple", dict(k_p=0.5, k_r=0.05, r_p=10.0))
    assert tab.curve.save_button.isEnabled()

    png = export_plot(tab.curve.plot, str(tmp_path / "fit.png"))
    assert Path(png).stat().st_size > 0


def test_isotopomer_bar_view_plots_relative_envelope(main_window, tmp_path):
    from riana.gui.export import export_plot

    view = main_window.integrate_tab.isobars
    assert not view.save_button.isEnabled()           # placeholder state
    # raw integrated areas (NaN/zero tolerated) → normalised relative bars
    view.plot_isotopomers("PEPTIDEK_2", [100.0, 50.0, 25.0, float("nan"), 0.0, 10.0])
    assert view.save_button.isEnabled()
    out = export_plot(view.plot, str(tmp_path / "iso.png"))
    assert Path(out).stat().st_size > 0


def test_integrate_row_select_fills_isotopomer_bars(main_window):
    tab = main_window.integrate_tab
    tab.model.set_dataframe(pd.DataFrame({
        "concat": ["PEPK_2"],
        "iso0": [100.0], "iso1": [50.0], "iso2": [25.0],
        "iso3": [12.0], "iso4": [6.0], "iso5": [3.0],
    }))
    tab.isobars.show_placeholder("…")
    tab._show_isotopomers(0)                           # the sync wiring on row-select
    assert tab.isobars.save_button.isEnabled()


def test_chromatogram_view_export_writes_png(main_window, tmp_path):
    from riana.gui.export import export_plot

    trace = PeptideTrace(
        concat="PEPTIDEK_2",
        chromatograms={0: Chromatogram(
            isotopomer=0, target_mz=500.0, mass_tol_ppm=10.0,
            scans=(1, 2, 3), rt=(10.0, 10.1, 10.2), intensity=(5.0, 9.0, 4.0))},
        window=(10.0, 10.2),
    )
    view = main_window.integrate_tab.chromatogram
    assert not view.save_button.isEnabled()  # placeholder state
    view.plot_trace(trace)
    assert view.save_button.isEnabled()

    out = export_plot(view.plot, str(tmp_path / "chrom.png"))
    assert Path(out).stat().st_size > 0


def test_model_tab_manifest_plots_selected_condition_group(main_window):
    """A manifest fit has one result row per (experiment, condition) per peptide.

    Selecting a row must plot *that* group, not always the first — the bug where
    the curve showed a different condition than the row's n_mbr/n_metox counts.
    """
    from pyqtgraph import PlotDataItem

    tab = main_window.model_tab
    rdf = pd.DataFrame(
        {
            "experiment": ["e", "e"], "condition": ["atrium", "control"],
            "k_deg": [0.3, 0.4], "R_squared": [0.97, 0.98],
            "sd": [0.01, 0.01], "spep": [8.0, 8.0],
            "ci_lo": [0.25, 0.35], "ci_hi": [0.35, 0.45], "protein id": ["P", "P"],
            "n_points": [3, 4], "n_mbr": [0, 1], "n_metox": [0, 0], "n_clean": [3, 3],
            "t": [[0.0, 1.0, 2.0], [0.0, 1.0, 2.0, 4.0]],
            "fs": [[0.0, 0.2, 0.4], [0.0, 0.25, 0.45, 0.6]],
            "evidence": [["q_value"] * 3,
                         ["q_value", "mbr", "q_value", "q_value"]],
            "metox": [[False] * 3, [False] * 4],
        },
        index=pd.Index(["PEP_2", "PEP_2"], name="concat"),
    )
    tab._result_df = rdf
    tab._last_config = tab.build_config()
    tab._populate_results(rdf)

    # Both groups are visible and distinguishable by the condition column.
    assert list(tab.model.dataframe["condition"]) == ["atrium", "control"]

    def _series_names():
        return {i.name() for i in tab.curve.plot.items
                if isinstance(i, PlotDataItem)}

    # control (row 1) carries the MBR point …
    tab.table.setCurrentIndex(tab.model.index(1, 0))
    assert any(n and n.startswith("MBR") for n in _series_names())

    # … atrium (row 0) is all-direct — no MBR series (the iloc[0] bug would have
    # shown atrium for both).
    tab.table.setCurrentIndex(tab.model.index(0, 0))
    assert not any(n and n.startswith("MBR") for n in _series_names())

    # The header carries protein + condition so an exported figure self-identifies.
    title = tab.curve.plot.titleLabel.text
    assert title.startswith("PEP_2") and "atrium" in title


def test_curve_view_legend_is_outside_the_plot_viewbox(main_window):
    """The legend lives in its own column, not anchored inside the data ViewBox —
    so its MBR/fold sample glyphs can't be mistaken for plotted points."""
    cv = main_window.model_tab.curve
    cv.plot_fit("PEP", [0, 1, 2], [0.0, 0.3, 0.5], 0.2, "simple", {},
                evidence=["q_value", "mbr", "q_value"])
    # legend's parent is the dedicated legend ViewBox, not the plot's own ViewBox.
    assert cv.plot.legend.parentItem() is cv._legend_vb
    assert cv.plot.legend.parentItem() is not cv.plot.getViewBox()


def test_main_window_sets_app_icon(main_window):
    assert not main_window.windowIcon().isNull()
