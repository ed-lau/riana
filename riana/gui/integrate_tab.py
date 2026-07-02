# -*- coding: utf-8 -*-

"""The Integrate tab: build an :class:`IntegrationConfig` from a form, run the
integration asynchronously, and show progress / log / results / per-fraction
drift, with a click-to-inspect chromatogram.

The GUI runs the **SDRF path** only: the search-ID file is read as the quantms
mzTab (DDA) / DIA-NN report.parquet (DIA) and routed through
:mod:`riana.core.pipeline` — :func:`~riana.core.pipeline.plan_integration` (in a
worker) builds the per-run tasks, this tab dispatches each over its *own* shared
pool, then :func:`~riana.core.pipeline.finalize_run` writes one identity-stamped
``<stem>_riana.txt`` per run and a ``riana_manifest.tsv``. It is **file-parallel**
(the *Workers* control, bounded by one mzML in memory per concurrent run). The
demoted bare-Percolator (no-SDRF) intake stays CLI-only — GUI users are all on
the SDRF path.

The form builds the *same* frozen :class:`~riana.config.IntegrationConfig` the
CLI builds, so its ``__post_init__`` is the single shared validator. CPU work is
awaited on the shared ``ProcessPoolExecutor`` via the Qt-free
:mod:`riana.gui.tasks` workers so the UI stays responsive.
"""

from __future__ import annotations

import asyncio
import os
import re
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtCore import QModelIndex, Qt, QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from qasync import asyncSlot

from riana.config import IntegrationConfig
from riana.core.pipeline import finalize_run
from riana.gui.chromatogram import ChromatogramView, IsotopomerBarView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import (
    extract_trace,
    integrate_fraction,
    load_integrate_results,
    plan_sdrf_integration,
)
from riana.io.manifest import MANIFEST_FILENAME, append_manifest
from riana.records import PSMRecord

# A multi-file integrate concat is 10⁵–10⁶ rows (a 384-file fractionated series is
# millions). Handing all of them to the QTableView makes sort / selection / memory
# the bottleneck even with a fast model, and nobody scrolls millions of rows to
# find a peptide. So the view shows at most this many rows and the user narrows
# with the filter box; the full result is on disk in the per-run ``_riana.txt``.
_MAX_DISPLAY_ROWS = 5000
# Text columns the filter box matches against (case-insensitive substring).
_FILTER_COLUMNS = ("sequence", "protein id", "concat")


class IntegrateTab(QWidget):
    """Form + async runner + results/chromatogram for ``riana integrate``."""

    def __init__(
        self,
        pool,
        status_cb: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__()
        self.pool = pool
        self._status_cb = status_cb or (lambda _msg: None)
        self._cancelled = False
        self._running = False
        # file_idx -> mzML path, for on-demand chromatogram extraction.
        self._fraction_mzml: dict[int, str] = {}
        # concat -> (min_scan, max_scan), precomputed once per result so the
        # chromatogram extraction on a row click is an O(1) lookup, not a
        # full-frame scan of a 10^5-10^6-row concat.
        self._scan_spans: dict[str, tuple[int, int]] = {}
        # The full result stays here; the table model only ever holds a
        # filtered + capped *view* of it (see ``_apply_filter``). ``_search_key``
        # is a lowercased ``sequence\x00protein\x00concat`` column built once so a
        # filter is a single ``str.contains`` rather than three per keystroke.
        self._full_df: pd.DataFrame | None = None
        self._search_key: pd.Series | None = None
        #: whether the Output dir already holds integrate results to display.
        self._results_available = False

        self._build_ui()

    # --- UI construction ---------------------------------------------------- #
    def _build_ui(self) -> None:
        splitter = QSplitter(Qt.Orientation.Horizontal)

        splitter.addWidget(self._build_form())
        splitter.addWidget(self._build_results())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([360, 740])

        outer = QVBoxLayout(self)
        outer.addWidget(splitter)

    def _path_row(self, edit: QLineEdit, on_browse: Callable[[], None]) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        button = QPushButton("Browse…")
        button.clicked.connect(on_browse)
        layout.addWidget(button)
        return row

    def _build_form(self) -> QWidget:
        box = QGroupBox("Integrate")
        form = QFormLayout(box)

        self.mzml_edit = QLineEdit()
        self.mzml_edit.setToolTip(
            "Folder of mzML files to integrate (one per run/timepoint).")
        self.mzml_edit.setPlaceholderText("Folder containing mzML file(s)")
        form.addRow("mzML folder", self._path_row(self.mzml_edit, self._pick_mzml))

        self.id_edit = QLineEdit()
        self.id_edit.setToolTip(
            "Search-ID file: quantms mzTab (DDA) or DIA-NN report.parquet covering every run in the SDRF.")
        self.id_edit.setPlaceholderText("quantms mzTab (DDA) or DIA-NN report.parquet")
        form.addRow("Search ID", self._path_row(self.id_edit, self._pick_id))

        self.sdrf_edit = QLineEdit()
        self.sdrf_edit.setToolTip(
            "SDRF samplesheet — drives the identity/manifest path (timepoints, conditions, RIA, mass tolerance all come from it).")
        self.sdrf_edit.setPlaceholderText("SDRF .tsv — drives the identity / manifest path")
        # Reflect the SDRF's precursor mass tolerance in the spinbox on entry
        # (the spinbox is built below; the slot reads it at call time).
        self.sdrf_edit.editingFinished.connect(self._resolve_sdrf_mass_tol)
        form.addRow("SDRF", self._path_row(self.sdrf_edit, self._pick_sdrf))

        self.iso_edit = QLineEdit("5")
        self.iso_edit.setToolTip(
            "Isotopomers to integrate. A single index N = the m0..mN envelope "
            "(e.g. '5' = m0-m5, the D2O default); an explicit list for a "
            "non-contiguous set (e.g. '0 6'); or 'auto' for adaptive per-peptide "
            "N_ISO from the IsoSpec init+final envelope.")
        form.addRow("Isotopomers", self.iso_edit)

        self.ria_spin = QDoubleSpinBox()
        self.ria_spin.setDecimals(3)
        self.ria_spin.setRange(0.0, 1.0)
        self.ria_spin.setSingleStep(0.005)
        self.ria_spin.setValue(0.06)
        self.ria_spin.setToolTip(
            "Precursor enrichment (RIA max). Used with --iso auto to shape the "
            "fully-labelled envelope; on the SDRF path it is read from "
            "characteristics[precursor enrichment] instead.")
        form.addRow("Precursor enrichment", self.ria_spin)

        self.qvalue_spin = QDoubleSpinBox()
        self.qvalue_spin.setToolTip(
            "Integrate only PSMs with q-value below this (FDR threshold).")
        self.qvalue_spin.setDecimals(4)
        self.qvalue_spin.setRange(0.0, 1.0)
        self.qvalue_spin.setSingleStep(0.001)
        self.qvalue_spin.setValue(0.01)
        form.addRow("Max q-value", self.qvalue_spin)

        self.mass_tol_spin = QSpinBox()
        self.mass_tol_spin.setRange(1, 500)
        self.mass_tol_spin.setValue(
            IntegrationConfig.__dataclass_fields__["mass_tol_ppm"].default)
        self.mass_tol_spin.setSuffix(" ppm")
        self.mass_tol_spin.setToolTip(
            "Auto-filled from the SDRF's comment[precursor mass tolerance] when "
            "you pick an SDRF; override here if you want.")
        form.addRow("Mass tolerance", self.mass_tol_spin)

        self.peak_rt_combo = QComboBox()
        self.peak_rt_combo.setToolTip(
            "Integration-window anchor: apex (spike winner), ms2 (0.9.0 parity), or consensus (median apex over m0–m3; prefer at high D2O).")
        self.peak_rt_combo.addItems(["apex", "ms2", "consensus"])
        form.addRow("Window anchor", self.peak_rt_combo)

        self.ihw_edit = QLineEdit("0.15")
        self.ihw_edit.setToolTip(
            "Integration half-width in RT minutes (dial to your peak width), or 'auto' to detect boundaries.")
        self.ihw_edit.setPlaceholderText("RT min, or 'auto'")
        form.addRow("Integration ½-width", self.ihw_edit)

        self.workers_spin = QSpinBox()
        self.workers_spin.setRange(1, os.cpu_count() or 1)
        self.workers_spin.setValue(1)
        self.workers_spin.setToolTip(
            "Runs/files integrated concurrently (one mzML in memory each; "
            "2-4 suits a many-timepoint series). The parallelism lever.")
        form.addRow("Workers (files)", self.workers_spin)

        self.out_edit = QLineEdit(".")
        self.out_edit.setToolTip(
            "Output directory for the _riana.txt files and the manifest. If it "
            "already holds a project (a riana_manifest.tsv with integrate rows), "
            "'Load results' displays those without re-integrating.")
        # The Output dir doubles as the project locator: point it at an existing
        # project folder and the Load button offers to display its integrate
        # results (the Integrate tab has no manifest field — integrate *creates*
        # the manifest, so a manifest input would invert the data flow).
        self.out_edit.editingFinished.connect(self._check_for_saved_results)
        form.addRow("Output dir", self._path_row(self.out_edit, self._pick_out))

        # The power-user tuning dials — collapsed by default so the everyday
        # form stays simple; mirrors the CLI's "Advanced integration" + "MBR"
        # rich-help panels (all build the same frozen IntegrationConfig).
        form.addRow(self._build_advanced())

        buttons = QHBoxLayout()
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._on_run)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self._on_cancel)
        self.load_button = QPushButton("Load results")
        self.load_button.setEnabled(False)
        self.load_button.setToolTip(
            "Display the integrate results already saved in the Output dir (the "
            "manifest's stage='integrate' rows) without re-integrating. Set the "
            "mzML folder too if you want the chromatogram on row-select.")
        self.load_button.clicked.connect(self._on_load)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        buttons.addWidget(self.load_button)
        form.addRow(buttons)

        self.results_hint = QLabel("")
        self.results_hint.setStyleSheet("color: palette(mid);")
        self.results_hint.setWordWrap(True)
        form.addRow(self.results_hint)

        self.error_label = QLabel("")
        self.error_label.setStyleSheet("color: #b00020;")
        self.error_label.setWordWrap(True)
        form.addRow(self.error_label)

        # The form can get tall once Advanced is expanded — keep it scrollable
        # so Run/Cancel never fall off the bottom of a short window.
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(box)
        return scroll

    def _build_advanced(self) -> QWidget:
        """The collapsible *Advanced* group: every remaining ``IntegrationConfig``
        knob (the tuning dials + the MBR sub-group), so the GUI can build the
        full config the CLI can — no surface drift. Checkable + collapsed by
        default; toggling shows/hides the inner widgets."""
        box = QGroupBox("Advanced…")
        box.setCheckable(True)
        box.setChecked(False)
        outer = QVBoxLayout(box)

        inner = QWidget()
        box.toggled.connect(inner.setVisible)
        inner.setVisible(False)
        outer.addWidget(inner)
        form = QFormLayout(inner)

        # --- peak / baseline dials (moved out of the everyday form) ----------
        self.baseline_combo = QComboBox()
        self.baseline_combo.setToolTip(
            "In-window baseline subtraction: none, noise_floor, snip, or asls.")
        self.baseline_combo.addItems(["none", "noise_floor", "snip", "asls"])
        form.addRow("Baseline", self.baseline_combo)

        self.apex_combo = QComboBox()
        self.apex_combo.addItems(["tallest", "nearest"])
        self.apex_combo.setToolTip("Apex pick rule for peak-rt apex/consensus.")
        form.addRow("Apex selection", self.apex_combo)

        self.apex_search_spin = QDoubleSpinBox()
        self.apex_search_spin.setRange(0.0, 10.0)
        self.apex_search_spin.setDecimals(2)
        self.apex_search_spin.setSingleStep(0.05)
        self.apex_search_spin.setValue(0.25)
        self.apex_search_spin.setSuffix(" min")
        self.apex_search_spin.setToolTip(
            "Half-width bounding the apex search around the PSM/MBR RT prior "
            "(0 = whole extraction).")
        form.addRow("Apex search ½-width", self.apex_search_spin)

        self.prominence_spin = QDoubleSpinBox()
        self.prominence_spin.setRange(0.1, 100.0)
        self.prominence_spin.setDecimals(2)
        self.prominence_spin.setSingleStep(0.5)
        self.prominence_spin.setValue(3.0)
        self.prominence_spin.setToolTip(
            "Apex-finder strictness — clear K·1.4826·MAD(trace); higher = stricter.")
        form.addRow("Prominence k", self.prominence_spin)

        self.width_rel_spin = QDoubleSpinBox()
        self.width_rel_spin.setRange(0.01, 0.99)
        self.width_rel_spin.setDecimals(2)
        self.width_rel_spin.setSingleStep(0.05)
        self.width_rel_spin.setValue(0.05)
        self.width_rel_spin.setToolTip(
            "Apex-height fraction for boundary detection — only used with "
            "Integration ½-width = auto (0.05 = 5% of apex, 0.5 = FWHM).")
        form.addRow("Width rel-height", self.width_rel_spin)

        self.apex_n_spin = QSpinBox()
        self.apex_n_spin.setRange(1, 10)
        self.apex_n_spin.setValue(4)
        self.apex_n_spin.setToolTip(
            "Channels (m0..m{N-1}) the consensus apex pools over — only used "
            "with Window anchor = consensus.")
        form.addRow("Apex N (consensus)", self.apex_n_spin)

        # --- extraction / smoothing / mass ----------------------------------
        ext_row = QWidget()
        ext_layout = QHBoxLayout(ext_row)
        ext_layout.setContentsMargins(0, 0, 0, 0)
        self.ext_override_check = QCheckBox("override")
        self.ext_spin = QDoubleSpinBox()
        self.ext_spin.setToolTip(
            "Extraction half-width (RT min): how much XIC to pull. Only used when 'override' is checked.")
        self.ext_spin.setRange(0.01, 20.0)
        self.ext_spin.setDecimals(2)
        self.ext_spin.setSingleStep(0.05)
        self.ext_spin.setValue(0.5)
        self.ext_spin.setSuffix(" min")
        self.ext_spin.setEnabled(False)
        self.ext_override_check.toggled.connect(self.ext_spin.setEnabled)
        self.ext_override_check.setToolTip(
            "Extraction half-width: how much XIC to pull. Off = derived from "
            "Integration ½-width / Window anchor (the default).")
        ext_layout.addWidget(self.ext_override_check)
        ext_layout.addWidget(self.ext_spin)
        form.addRow("Extraction ½-width", ext_row)

        self.smoothing_combo = QComboBox()
        self.smoothing_combo.addItems(["off", "3", "5", "7", "9", "11"])
        self.smoothing_combo.setToolTip("Savitzky-Golay window (odd ≥ 3); off disables.")
        form.addRow("Smoothing window", self.smoothing_combo)

        self.smoothing_poly_spin = QSpinBox()
        self.smoothing_poly_spin.setRange(2, 5)
        self.smoothing_poly_spin.setValue(2)
        self.smoothing_poly_spin.setToolTip("SG polynomial order; only used when smoothing is on.")
        form.addRow("Smoothing poly-order", self.smoothing_poly_spin)

        self.mass_diff_spin = QDoubleSpinBox()
        self.mass_diff_spin.setRange(0.1, 10.0)
        self.mass_diff_spin.setDecimals(9)
        self.mass_diff_spin.setSingleStep(0.001)
        self.mass_diff_spin.setValue(1.003354835)
        self.mass_diff_spin.setToolTip("Mass step between isotopomers (C13 default).")
        form.addRow("Mass difference", self.mass_diff_spin)

        self.ppm_alert_spin = QDoubleSpinBox()
        self.ppm_alert_spin.setRange(0.1, 1000.0)
        self.ppm_alert_spin.setValue(20.0)
        self.ppm_alert_spin.setSuffix(" ppm")
        self.ppm_alert_spin.setToolTip("Warn when a fraction's median ppm error exceeds this.")
        form.addRow("Drift alert", self.ppm_alert_spin)

        self.write_intensities_check = QCheckBox("write pre-integration trace")
        form.addRow("Intensities", self.write_intensities_check)

        # --- intake scan↔precursor guard ------------------------------------
        self.id_check = QCheckBox("enabled")
        self.id_check.setChecked(True)
        self.id_check.setToolTip(
            "Per-run check that mzTab spectra_ref scans point at the matching "
            "precursor m/z in this mzML (catches a wrong mzML↔mzTab pairing / "
            "quantms filename-prefix scramble; mass-based, immune to RT "
            "alignment). Uncheck only for a run you know is correctly paired.")
        form.addRow("Scan↔precursor guard", self.id_check)

        self.precursor_tol_spin = QDoubleSpinBox()
        self.precursor_tol_spin.setToolTip(
            "Intake guard: per-scan precursor-m/z match tolerance (ppm). Lenient "
            "to drift; a wrong file lands hundreds of ppm off.")
        self.precursor_tol_spin.setRange(1.0, 1000.0)
        self.precursor_tol_spin.setDecimals(1)
        self.precursor_tol_spin.setSingleStep(1.0)
        self.precursor_tol_spin.setValue(10.0)
        self.precursor_tol_spin.setSuffix(" ppm")
        self.id_check.toggled.connect(self.precursor_tol_spin.setEnabled)
        form.addRow("Precursor m/z tolerance", self.precursor_tol_spin)

        # --- MBR sub-group (the checkable box state IS the mbr flag) ----------
        self.mbr_box = QGroupBox("Match-between-runs (MBR)")
        self.mbr_box.setCheckable(True)
        self.mbr_box.setChecked(False)
        self.mbr_box.setToolTip(
            "Transfer a confidently-identified precursor into runs of its "
            "(experiment, condition) curve that missed it. SDRF/mzTab DDA path "
            "only; no-op on DIA.")
        mbr_form = QFormLayout(self.mbr_box)

        self.mbr_donor_runs_spin = QSpinBox()
        self.mbr_donor_runs_spin.setToolTip(
            "MBR: a precursor must be confidently identified in at least this many runs to seed a transfer.")
        self.mbr_donor_runs_spin.setRange(2, 100)
        self.mbr_donor_runs_spin.setValue(2)
        mbr_form.addRow("Min donor runs", self.mbr_donor_runs_spin)

        self.mbr_donor_q_spin = QDoubleSpinBox()
        self.mbr_donor_q_spin.setToolTip(
            "MBR: donor q-value threshold — only IDs at/below this seed transfers.")
        self.mbr_donor_q_spin.setRange(0.0, 1.0)
        self.mbr_donor_q_spin.setDecimals(4)
        self.mbr_donor_q_spin.setSingleStep(0.001)
        self.mbr_donor_q_spin.setValue(0.01)
        mbr_form.addRow("Donor q-value", self.mbr_donor_q_spin)

        self.mbr_snr_spin = QDoubleSpinBox()
        self.mbr_snr_spin.setRange(0.0, 1000.0)
        self.mbr_snr_spin.setDecimals(1)
        self.mbr_snr_spin.setSingleStep(1.0)
        self.mbr_snr_spin.setValue(4.0)
        self.mbr_snr_spin.setToolTip("Apex-SNR floor for transfers (0 = ungated).")
        mbr_form.addRow("Min apex SNR", self.mbr_snr_spin)

        self.mbr_scans_spin = QSpinBox()
        self.mbr_scans_spin.setRange(0, 1000)
        self.mbr_scans_spin.setValue(3)
        self.mbr_scans_spin.setToolTip("Min nonzero scans in the integration window (0 = off).")
        mbr_form.addRow("Min scans", self.mbr_scans_spin)

        form.addRow(self.mbr_box)
        return box

    def _build_results(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        self.progress = QProgressBar()
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.drift_label = QLabel("No run yet.")
        self.drift_label.setWordWrap(True)
        layout.addWidget(self.drift_label)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setFixedHeight(120)
        layout.addWidget(self.log)

        results_split = QSplitter(Qt.Orientation.Vertical)

        # The table, with a filter box + row-count note above it. The model holds
        # only a capped/filtered view (see class note), so the box is how the user
        # reaches a peptide the cap left off-screen.
        table_panel = QWidget()
        table_layout = QVBoxLayout(table_panel)
        table_layout.setContentsMargins(0, 0, 0, 0)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter:"))
        self.filter_edit = QLineEdit()
        self.filter_edit.setClearButtonEnabled(True)
        self.filter_edit.setPlaceholderText("sequence / protein / concat contains…")
        self.filter_edit.setToolTip(
            "Show only rows whose sequence, protein id, or concat contains this "
            "text (case-insensitive). The table shows at most "
            f"{_MAX_DISPLAY_ROWS:,} rows at a time — filter to find a peptide in a "
            "large run. The full result is on disk in each run's _riana.txt.")
        self.filter_edit.textChanged.connect(self._on_filter_text)
        filter_row.addWidget(self.filter_edit, stretch=1)
        self.rows_label = QLabel("")
        self.rows_label.setStyleSheet("color: palette(mid);")
        filter_row.addWidget(self.rows_label)
        table_layout.addLayout(filter_row)

        # Debounce: filtering a millions-row frame is ~1 s, so re-filter only once
        # the user pauses typing, not on every keystroke.
        self._filter_timer = QTimer(self)
        self._filter_timer.setSingleShot(True)
        self._filter_timer.setInterval(250)
        self._filter_timer.timeout.connect(self._apply_filter)

        self.table = QTableView()
        self.model = DataFrameTableModel()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSortingEnabled(True)
        self.table.selectionModel().currentRowChanged.connect(self._on_row_changed)
        table_layout.addWidget(self.table)
        results_split.addWidget(table_panel)

        # Two synced views of the selected peptide: the RT-domain chromatogram and
        # the abundance-domain isotopomer (m0..mN) bar chart, side by side.
        plots = QSplitter(Qt.Orientation.Horizontal)
        self.chromatogram = ChromatogramView()
        self.isobars = IsotopomerBarView()
        plots.addWidget(self.chromatogram)
        plots.addWidget(self.isobars)
        plots.setSizes([460, 280])
        results_split.addWidget(plots)
        results_split.setSizes([320, 300])

        layout.addWidget(results_split, stretch=1)
        return panel

    # --- file pickers ------------------------------------------------------- #
    def _pick_mzml(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select mzML folder")
        if path:
            self.mzml_edit.setText(path)

    def _pick_id(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select search-ID file (quantms mzTab or DIA-NN parquet)",
            filter="Search ID (*.mzTab *.parquet *.txt);;All files (*)"
        )
        if path:
            self.id_edit.setText(path)

    def _pick_sdrf(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select SDRF samplesheet",
            filter="SDRF (*.tsv *.sdrf.tsv);;All files (*)"
        )
        if path:
            self.sdrf_edit.setText(path)
            self._resolve_sdrf_mass_tol()

    def _resolve_sdrf_mass_tol(self) -> None:
        """Reflect the SDRF's ``comment[precursor mass tolerance]`` in the mass
        tolerance spinbox (the user can still override) — the GUI equivalent of
        the CLI's SDRF mass-tolerance resolution. Silent on a missing/unreadable
        SDRF or one without a tolerance; the run path surfaces real errors."""
        path = self.sdrf_edit.text().strip()
        if not path or not Path(path).is_file():
            return
        try:
            from riana.io.sdrf import read_sdrf
            tol = read_sdrf(path).precursor_mass_tol_ppm
        except Exception:
            return
        if tol is not None:
            self.mass_tol_spin.setValue(int(round(tol)))
            self._status_cb(f"mass tolerance set from SDRF: {int(round(tol))} ppm")

    def _pick_out(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.out_edit.setText(path)
            self._check_for_saved_results()

    # --- config marshalling (the shared-validation contract) ---------------- #
    def build_config(self) -> IntegrationConfig:
        """Build the frozen config from the widgets.

        Raises ``ValueError`` on a bad value — from the iso parse or from
        :meth:`IntegrationConfig.__post_init__`, the *same* validator the CLI
        uses. Callers surface the message inline.
        """
        # Mirror riana.cli.integrate: 'auto' = adaptive N_ISO; a single index N =
        # the contiguous iso0..isoN capture; a multi-value list stays explicit
        # (the o18 '0 6' pair).
        iso_text = self.iso_edit.text().strip()
        adaptive_iso = iso_text.lower() == "auto"
        if adaptive_iso:
            isotopomers = (0, 1, 2, 3, 4, 5)
        else:
            parsed = sorted({int(p) for p in re.split(r"[,\s]+", iso_text) if p})
            if not parsed:
                raise ValueError("Isotopomers must be 'auto', a single index N "
                                 "(= iso0..isoN), or an explicit list, e.g. '0 6'.")
            isotopomers = (tuple(range(parsed[0] + 1)) if len(parsed) == 1
                           else tuple(parsed))

        ihw_text = self.ihw_edit.text().strip()
        ihw: float | str = "auto" if ihw_text == "auto" else float(ihw_text)
        peak_rt = self.peak_rt_combo.currentText()
        # Extraction half-width: an explicit Advanced override, else derived
        # exactly as riana.cli.integrate does — ms2 integrates the whole
        # extraction (= ihw); apex/consensus and 'auto' need room for the apex
        # offset (+0.33).
        if self.ext_override_check.isChecked():
            ehw = float(self.ext_spin.value())
        elif peak_rt == "ms2" and ihw != "auto":
            ehw = float(ihw)
        else:
            ehw = (0.33 if ihw == "auto" else float(ihw)) + 0.33

        smoothing_text = self.smoothing_combo.currentText()
        smoothing = None if smoothing_text == "off" else int(smoothing_text)

        return IntegrationConfig(
            isotopomers=isotopomers,
            adaptive_iso=adaptive_iso,
            ria_max=float(self.ria_spin.value()),
            mass_tol_ppm=int(self.mass_tol_spin.value()),
            extraction_half_width=ehw,
            peak_rt=peak_rt,
            integration_half_width=ihw,
            baseline_method=self.baseline_combo.currentText(),
            apex_selection=self.apex_combo.currentText(),
            apex_search_half_width=float(self.apex_search_spin.value()),
            apex_n_consensus=int(self.apex_n_spin.value()),
            prominence_k=float(self.prominence_spin.value()),
            width_rel_height=float(self.width_rel_spin.value()),
            q_value=float(self.qvalue_spin.value()),
            write_intensities=bool(self.write_intensities_check.isChecked()),
            smoothing=smoothing,
            smoothing_polyorder=int(self.smoothing_poly_spin.value()),
            mass_difference=float(self.mass_diff_spin.value()),
            ppm_alert=float(self.ppm_alert_spin.value()),
            check_scan_id=bool(self.id_check.isChecked()),
            scan_precursor_tol_ppm=float(self.precursor_tol_spin.value()),
            mbr=bool(self.mbr_box.isChecked()),
            mbr_min_donor_runs=int(self.mbr_donor_runs_spin.value()),
            mbr_donor_q=float(self.mbr_donor_q_spin.value()),
            mbr_min_snr=float(self.mbr_snr_spin.value()),
            mbr_min_scans=int(self.mbr_scans_spin.value()),
            out_dir=self.out_edit.text().strip() or ".",
        )

    # --- run flow ----------------------------------------------------------- #
    @asyncSlot()
    async def _on_run(self) -> None:
        if self._running:
            return
        self.error_label.setText("")
        self.log.clear()
        self.model.set_dataframe(pd.DataFrame())
        self.chromatogram.show_placeholder("Running…")
        self.isobars.show_placeholder("Running…")
        self._fraction_mzml.clear()
        self._scan_spans = {}
        self._full_df = None
        self._search_key = None
        self.filter_edit.blockSignals(True)  # clear without arming the debounce
        self.filter_edit.clear()
        self.filter_edit.blockSignals(False)
        self.rows_label.setText("")
        self._cancelled = False

        try:
            config = self.build_config()
        except ValueError as exc:
            self._fail(str(exc))
            return

        mzml_dir = self.mzml_edit.text().strip()
        id_path = self.id_edit.text().strip()
        sdrf_path = self.sdrf_edit.text().strip()
        workers = int(self.workers_spin.value())
        if not mzml_dir or not Path(mzml_dir).is_dir():
            self._fail("Select a valid mzML folder.")
            return
        if not id_path or not Path(id_path).is_file():
            self._fail("Select a valid search-ID file (quantms mzTab or DIA-NN "
                       "report.parquet).")
            return
        if not sdrf_path or not Path(sdrf_path).is_file():
            self._fail("Select an SDRF .tsv — the GUI integrates via the SDRF / "
                       "manifest path. (The bare Percolator path is CLI-only.)")
            return

        os.makedirs(config.out_dir, exist_ok=True)
        self._set_running(True)
        try:
            await self._run_sdrf(config, sdrf_path, mzml_dir, id_path, workers)
        except Exception as exc:  # surface worker/IO errors instead of crashing
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            self._set_running(False)

    async def _run_sdrf(self, config, sdrf_path, mzml_dir, mztab_path, workers):
        """SDRF/manifest path — the *same* core/pipeline plan + per-run unit the
        CLI uses, dispatched over this tab's shared pool (no nested pools)."""
        loop = asyncio.get_running_loop()
        self._info(f"planning runs from SDRF {sdrf_path} …")
        tasks = await loop.run_in_executor(
            self.pool, plan_sdrf_integration, config, sdrf_path, mzml_dir,
            mztab_path,
        )
        if not tasks:
            self._fail("No runs to integrate from the SDRF / mzTab.")
            return
        for t in tasks:
            self._fraction_mzml[t.file_idx] = t.mzml_path

        async def run_one(task):
            df, drift = await loop.run_in_executor(
                self.pool, integrate_fraction, config, task.psms,
                task.mzml_path, task.stem,
            )
            return task.file_idx, df, drift

        results = await self._gather_runs(tasks, workers, run_one)
        if self._cancelled:
            self._info("cancelled — outputs not written.")
            return
        # Finalize in file_idx order (deterministic) + append the manifest, the
        # same outputs integrate_project writes.
        out_dir = Path(config.out_dir)
        rows, frames = [], []
        for task in sorted(tasks, key=lambda x: x.file_idx):
            df, drift = results[task.file_idx]
            rows.append(finalize_run(config, task, df, out_dir, mztab_path))
            self._info(f"wrote {rows[-1].output_path}")
            self._update_drift(task.file_idx, drift, config.ppm_alert)
            frames.append(df)
        append_manifest(out_dir / MANIFEST_FILENAME, rows)
        self._info(f"appended {len(rows)} integrate rows to the manifest")
        self._finish_table(frames, config)

    async def _gather_runs(self, items, workers, run_one):
        """Run ``run_one(item)`` over *items*, ≤ *workers* concurrent, updating
        progress as each finishes.

        ``run_one`` returns ``(key, df, drift)``; returns ``{key: (df, drift)}``.
        Parallelism uses this tab's shared pool bounded by a semaphore (one mzML
        per concurrent run). Cancellation is best-effort: in-flight runs finish
        (a pool task can't be killed) but the caller discards the output via
        ``self._cancelled``.
        """
        sem = asyncio.Semaphore(max(1, int(workers)))

        async def _wrapped(item):
            async with sem:
                return await run_one(item)

        self.progress.setRange(0, len(items))
        results: dict = {}
        done = 0
        for coro in asyncio.as_completed([_wrapped(it) for it in items]):
            key, df, drift = await coro      # always await — no orphan coroutines
            results[key] = (df, drift)
            done += 1
            self.progress.setValue(done)
            if not self._cancelled:
                self._info(f"integrated run {key}")
        return results

    @staticmethod
    def _build_scan_spans(df: pd.DataFrame) -> dict[str, tuple[int, int]]:
        """Map each ``concat`` to its ``(min_scan, max_scan)`` in one vectorised
        groupby.

        This replaces the old per-selection ``df[df["concat"] == x]["scan"]``
        full-frame scan (~160 ms per click on a millions-row multi-file concat)
        with an O(1) dict lookup; the (min, max) is identical to what that scan
        produced. Empty / column-less frames give an empty map (the extractor
        then uses its default window). ``scan`` is coerced to numeric to match the
        old ``.astype(int)`` before min/max.
        """
        if df.empty or "concat" not in df.columns or "scan" not in df.columns:
            return {}
        scans = pd.to_numeric(df["scan"], errors="coerce")
        grouped = scans.groupby(df["concat"], sort=False).agg(["min", "max"])
        return {
            str(concat): (int(mn), int(mx))
            for concat, mn, mx in zip(grouped.index, grouped["min"], grouped["max"])
            if pd.notna(mn) and pd.notna(mx)
        }

    def _finish_table(self, frames, config) -> None:
        if not frames:
            return
        df = pd.concat(frames, ignore_index=True)
        self._last_config = config
        self._display_frame(df)
        self._info(f"done — {len(df)} rows across {len(frames)} run(s).")

    def _display_frame(self, df: pd.DataFrame) -> None:
        """Push a result frame into the (filtered + capped) view.

        The shared tail of a fresh run (:meth:`_finish_table`) and a "Load
        results" load: keep the full frame off-model, precompute the per-concat
        scan spans + filter key once, and show the capped view.
        """
        self._full_df = df
        # Precompute each peptide's scan span once, so a row click is an O(1)
        # lookup instead of a `df[df.concat == x]` scan of the whole frame.
        self._scan_spans = self._build_scan_spans(df)
        # Build the filter key once, then push the capped/filtered view into the
        # model (never the full 10⁵–10⁶-row frame).
        self._search_key = self._build_search_key(df)
        self._apply_filter()
        self.chromatogram.show_placeholder(
            "Select a peptide row to view its chromatogram."
        )
        self.isobars.show_placeholder(
            "Select a peptide row to view its isotopomer envelope."
        )

    # --- results filter (view = filtered + capped) -------------------------- #
    @staticmethod
    def _build_search_key(df: pd.DataFrame):
        """A lowercased ``sequence\\x00protein\\x00concat`` column for the filter.

        Concatenating the searchable columns once (NUL-joined so a match can't
        span a boundary) turns each filter into a single ``str.contains`` on this
        key instead of one per column per keystroke — ~3× cheaper on a big frame.
        Returns ``None`` when the frame has none of the filterable columns.
        """
        cols = [c for c in _FILTER_COLUMNS if c in df.columns]
        if df.empty or not cols:
            return None
        key = df[cols[0]].astype(str)
        for c in cols[1:]:
            key = key + "\x00" + df[c].astype(str)
        return key.str.lower()

    def _on_filter_text(self, _text: str) -> None:
        # Debounced: arm the timer; _apply_filter runs when typing pauses.
        self._filter_timer.start()

    def _apply_filter(self) -> None:
        """Push a filtered + row-capped view of the full result into the model."""
        if self._full_df is None:
            return
        full = self._full_df
        query = self.filter_edit.text().strip()
        if query and self._search_key is not None:
            filtered = full[self._search_key.str.contains(
                query.lower(), regex=False, na=False)]
        else:
            filtered = full
        view = filtered.iloc[:_MAX_DISPLAY_ROWS]
        self.model.set_dataframe(view)
        self._update_rows_label(len(full), len(filtered), len(view), bool(query))

    def _update_rows_label(self, total: int, matched: int, shown: int,
                           filtered: bool) -> None:
        if filtered:
            if shown < matched:
                self.rows_label.setText(
                    f"showing {shown:,} of {matched:,} matches "
                    f"({total:,} total) — narrow the filter")
            else:
                self.rows_label.setText(f"{matched:,} of {total:,} rows match")
        elif shown < total:
            self.rows_label.setText(
                f"showing first {shown:,} of {total:,} rows — filter to find a peptide")
        else:
            self.rows_label.setText(f"{total:,} rows")

    def _on_cancel(self) -> None:
        self._cancelled = True
        self.cancel_button.setEnabled(False)
        self._info("cancelling — in-flight runs finish, but outputs are discarded …")

    # --- load prior integrate results (no re-integrate) --------------------- #
    def _check_for_saved_results(self) -> None:
        """Enable Load + hint when the Output dir already holds integrate results.

        The Output dir is the project locator (the Integrate tab has no manifest
        field): a ``riana_manifest.tsv`` with ``stage="integrate"`` rows there
        means this folder is an existing project whose results can be displayed.
        """
        out_dir = self.out_edit.text().strip()
        found = False
        if out_dir and (Path(out_dir) / MANIFEST_FILENAME).is_file():
            try:
                from riana.io.manifest import read_manifest
                found = bool(read_manifest(
                    Path(out_dir) / MANIFEST_FILENAME, stage="integrate"))
            except Exception:
                found = False
        self._results_available = found
        self.load_button.setEnabled(found and not self._running)
        self.results_hint.setText(
            "✓ integrate results found in this folder — Load to view, or Run to "
            "re-integrate" if found else "")

    @asyncSlot()
    async def _on_load(self) -> None:
        """Display the integrate results already saved in the Output dir."""
        if self._running:
            return
        out_dir = self.out_edit.text().strip()
        if not out_dir or not (Path(out_dir) / MANIFEST_FILENAME).is_file():
            self._fail("Set the Output dir to a folder with a riana_manifest.tsv.")
            return
        try:
            config = self.build_config()  # for the chromatogram (a visualization)
        except ValueError as exc:
            self._fail(str(exc))
            return
        self.error_label.setText("")
        self._cancelled = False
        self._set_running(True)
        self.progress.setRange(0, 0)  # busy
        loop = asyncio.get_running_loop()
        try:
            self._info(f"loading saved integrate results from {out_dir} …")
            df = await loop.run_in_executor(
                self.pool, load_integrate_results, out_dir)
            if df.empty:
                self._fail("No integrate rows / readable outputs in the manifest.")
                return
            self._last_config = config
            # Map file_idx → mzML for the chromatogram, if the mzML folder is set
            # (the isotopomer bars + table work without it).
            self._fraction_mzml = self._resolve_fraction_mzml(df)
            self._display_frame(df)
            note = ("" if self._fraction_mzml
                    else "  (set the mzML folder for the chromatogram)")
            self._info(f"loaded {len(df)} rows from saved integrate results.{note}")
        except Exception as exc:  # surface loader/IO errors inline
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self._set_running(False)

    def _resolve_fraction_mzml(self, df: pd.DataFrame) -> dict[int, str]:
        """Map ``file_idx → mzML path`` from the mzML folder, for the chromatogram.

        A loaded frame carries ``file_idx`` + the run ``file`` stem but not the
        mzML path (that lived in the run tasks). If the mzML folder is set and the
        files are there, rebuild the map so the chromatogram works; otherwise it is
        empty and row-select shows the table + isotopomer bars only.
        """
        folder = self.mzml_edit.text().strip()
        if not folder or not {"file_idx", "file"}.issubset(df.columns):
            return {}
        base = Path(folder)
        out: dict[int, str] = {}
        for file_idx, stem in df[["file_idx", "file"]].drop_duplicates().itertuples(
                index=False):
            for ext in (".mzML", ".mzML.gz", ".mzml", ".mzml.gz"):
                cand = base / f"{stem}{ext}"
                if cand.is_file():
                    out[int(file_idx)] = str(cand)
                    break
        return out

    # --- chromatogram on selection ----------------------------------------- #
    def _on_row_changed(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if current.isValid() and not self._running:
            # The isotopomer bars come straight off the row (no mzML), so update
            # them instantly; the chromatogram needs an async extraction.
            self._show_isotopomers(current.row())
            asyncio.ensure_future(self._show_chromatogram(current.row()))

    def _show_isotopomers(self, row: int) -> None:
        df = self.model.dataframe
        if row < 0 or row >= len(df):
            return
        record = df.iloc[row]
        iso_cols = sorted(
            (c for c in df.columns if c.startswith("iso") and c[3:].isdigit()),
            key=lambda c: int(c[3:]))
        if not iso_cols:
            self.isobars.show_placeholder("No isotopomer columns in this result.")
            return
        concat = str(record.get("concat", f"row {row}"))
        self.isobars.plot_isotopomers(concat, [record[c] for c in iso_cols])

    async def _show_chromatogram(self, row: int) -> None:
        df = self.model.dataframe
        if row < 0 or row >= len(df) or not hasattr(self, "_last_config"):
            return
        record = df.iloc[row]
        file_idx = int(record["file_idx"])
        mzml_file = self._fraction_mzml.get(file_idx)
        if mzml_file is None:
            return
        psm = PSMRecord(
            scan=int(record["scan"]),
            charge=int(record["charge"]),
            sequence=str(record["sequence"]),
            peptide_mass=float(record["peptide mass"]),
            sample=str(record["sample"]),
            file_idx=file_idx,
            pep_id=int(record["pep_id"]),
        )
        # The integrator spans all kept scans of this peptide-charge; the kept
        # set IS the displayed frame, so the span is derived from it — but
        # precomputed once per result (``_build_scan_spans``) rather than scanned
        # per click. ``None`` (peptide absent, e.g. a stale selection) lets the
        # extractor fall back to its own default window.
        scan_span = self._scan_spans.get(psm.concat)

        self.chromatogram.show_placeholder(f"Extracting {psm.concat} …")
        loop = asyncio.get_running_loop()
        trace = await loop.run_in_executor(
            self.pool, extract_trace, self._last_config, psm, mzml_file, scan_span
        )
        if trace is None:
            self.chromatogram.show_placeholder(
                f"No extractable signal for {psm.concat}."
            )
        else:
            self.chromatogram.plot_trace(trace)

    # --- small helpers ------------------------------------------------------ #
    def _set_running(self, running: bool) -> None:
        self._running = running
        self.run_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)
        if running:
            self.load_button.setEnabled(False)
        else:  # re-enable Load iff the Output dir has results (post run/load)
            self._check_for_saved_results()

    def _info(self, message: str) -> None:
        self.log.appendPlainText(message)
        self._status_cb(message)

    def _fail(self, message: str) -> None:
        self.error_label.setText(message)
        self._info(f"error: {message}")

    def _update_drift(self, idx: int, drift, ppm_alert: float) -> None:
        if drift is None or drift.n == 0:
            self.drift_label.setText(f"Fraction {idx}: no mass-accuracy data.")
            return
        flag = " ⚠ exceeds alert" if abs(drift.median_ppm) > ppm_alert else " ✓"
        self.drift_label.setText(
            f"Fraction {idx}: median {drift.median_ppm:+.2f} ppm "
            f"(MAD {drift.mad_ppm:.2f}, n={drift.n}); "
            f"suggested shift {drift.suggested_shift_ppm:+.2f} ppm{flag}"
        )
