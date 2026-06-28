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
from PySide6.QtCore import QModelIndex, Qt
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
from riana.gui.chromatogram import ChromatogramView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import (
    extract_trace,
    integrate_fraction,
    plan_sdrf_integration,
)
from riana.io.manifest import MANIFEST_FILENAME, append_manifest
from riana.records import PSMRecord


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
        self.mzml_edit.setPlaceholderText("Folder containing mzML file(s)")
        form.addRow("mzML folder", self._path_row(self.mzml_edit, self._pick_mzml))

        self.id_edit = QLineEdit()
        self.id_edit.setPlaceholderText("quantms mzTab (DDA) or DIA-NN report.parquet")
        form.addRow("Search ID", self._path_row(self.id_edit, self._pick_id))

        self.sdrf_edit = QLineEdit()
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
        self.peak_rt_combo.addItems(["apex", "ms2", "consensus"])
        form.addRow("Window anchor", self.peak_rt_combo)

        self.ihw_edit = QLineEdit("0.15")
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
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        form.addRow(buttons)

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

        # --- intake scan↔RT guard -------------------------------------------
        self.rt_check = QCheckBox("enabled")
        self.rt_check.setChecked(True)
        self.rt_check.setToolTip(
            "Per-run check that mzTab spectra_ref scans reconcile with this "
            "mzML's RTs (catches the quantms filename-prefix scramble). "
            "Uncheck only for a run you know is correctly paired.")
        form.addRow("Scan↔RT guard", self.rt_check)

        self.scan_rt_tol_spin = QDoubleSpinBox()
        self.scan_rt_tol_spin.setRange(0.1, 60.0)
        self.scan_rt_tol_spin.setDecimals(2)
        self.scan_rt_tol_spin.setSingleStep(0.5)
        self.scan_rt_tol_spin.setValue(3.0)
        self.scan_rt_tol_spin.setSuffix(" min")
        self.rt_check.toggled.connect(self.scan_rt_tol_spin.setEnabled)
        form.addRow("Scan↔RT tolerance", self.scan_rt_tol_spin)

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
        self.mbr_donor_runs_spin.setRange(2, 100)
        self.mbr_donor_runs_spin.setValue(2)
        mbr_form.addRow("Min donor runs", self.mbr_donor_runs_spin)

        self.mbr_donor_q_spin = QDoubleSpinBox()
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

        self.table = QTableView()
        self.model = DataFrameTableModel()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSortingEnabled(True)
        self.table.selectionModel().currentRowChanged.connect(self._on_row_changed)
        results_split.addWidget(self.table)

        self.chromatogram = ChromatogramView()
        results_split.addWidget(self.chromatogram)
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
            check_scan_rt=bool(self.rt_check.isChecked()),
            scan_rt_tol_min=float(self.scan_rt_tol_spin.value()),
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
        self._fraction_mzml.clear()
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

    def _finish_table(self, frames, config) -> None:
        if not frames:
            return
        self.model.set_dataframe(pd.concat(frames, ignore_index=True))
        self.chromatogram.show_placeholder(
            "Select a peptide row to view its chromatogram."
        )
        self._last_config = config
        self._info(
            f"done — {len(self.model.dataframe)} rows across {len(frames)} run(s)."
        )

    def _on_cancel(self) -> None:
        self._cancelled = True
        self.cancel_button.setEnabled(False)
        self._info("cancelling — in-flight runs finish, but outputs are discarded …")

    # --- chromatogram on selection ----------------------------------------- #
    def _on_row_changed(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if current.isValid() and not self._running:
            asyncio.ensure_future(self._show_chromatogram(current.row()))

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
        # set IS the displayed frame, so derive the span from it.
        same = df[df["concat"] == psm.concat]["scan"].astype(int)
        scan_span = (int(same.min()), int(same.max()))

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
