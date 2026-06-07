# -*- coding: utf-8 -*-

"""The Model tab: build a :class:`FitConfig` from a form, fit a D₂O time series
asynchronously, and show per-peptide results with a click-to-inspect fitted
kinetic curve.

Mirrors the Integrate tab's pattern: the form builds the *same* frozen
:class:`~riana.config.FitConfig` the CLI builds (shared ``__post_init__``
validation), and the (single, batched) fit runs on the shared
``ProcessPoolExecutor`` via the Qt-free :func:`riana.gui.tasks.run_fit` worker so
the UI stays responsive. The fitted curve is drawn directly from the result
row's ``t`` / ``fs`` / ``k_deg`` (pure math — no worker round-trip).
"""

from __future__ import annotations

import asyncio
import dataclasses
import os
from concurrent.futures import Future
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtCore import QModelIndex, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from qasync import asyncSlot

from riana.config import FitConfig
from riana.core.fitting import available_coefficient_presets
from riana.gui.curve_view import CurveView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import run_fit
from riana.io.writers import make_provenance, write_dataframe_tsv

# Result columns to show in the table (the per-peptide ``t`` / ``fs`` lists are
# kept off-screen and used only for the curve plot).
_DISPLAY_COLS = ["concat", "k_deg", "R_squared", "sd", "ci_lo", "ci_hi",
                 "spep", "protein id"]


class ModelTab(QWidget):
    """Form + async fit runner + results/fitted-curve for ``riana fit``."""

    def __init__(
        self,
        pool,
        default_threads: int = 1,
        status_cb: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__()
        self.pool = pool
        self._status_cb = status_cb or (lambda _msg: None)
        self._running = False
        self._cancelled = False
        self._future: Future | None = None
        self._result_df: pd.DataFrame | None = None
        self._last_config: FitConfig | None = None

        self._build_ui(default_threads)

    # --- UI construction ---------------------------------------------------- #
    def _build_ui(self, default_threads: int) -> None:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_form(default_threads))
        splitter.addWidget(self._build_results())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([360, 740])
        outer = QVBoxLayout(self)
        outer.addWidget(splitter)

    def _build_form(self, default_threads: int) -> QWidget:
        box = QGroupBox("Fit")
        form = QFormLayout(box)

        # Input timepoint files (one _riana.txt per timepoint).
        self.files_list = QListWidget()
        self.files_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.files_list.setFixedHeight(96)
        form.addRow("Timepoint files", self.files_list)
        file_buttons = QHBoxLayout()
        add_btn = QPushButton("Add…")
        add_btn.clicked.connect(self._add_files)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.files_list.clear)
        file_buttons.addWidget(add_btn)
        file_buttons.addWidget(clear_btn)
        form.addRow("", _row(file_buttons))

        # Coefficients: editable combo of bundled presets, or a CSV path.
        self.coeff_combo = QComboBox()
        self.coeff_combo.setEditable(True)
        presets = available_coefficient_presets()
        self.coeff_combo.addItems(presets)
        if "commerford" in presets:
            self.coeff_combo.setCurrentText("commerford")
        self.coeff_combo.lineEdit().setPlaceholderText("preset name or CSV path")
        coeff_row = QHBoxLayout()
        coeff_row.addWidget(self.coeff_combo, stretch=1)
        coeff_browse = QPushButton("Browse…")
        coeff_browse.clicked.connect(self._pick_coefficients)
        coeff_row.addWidget(coeff_browse)
        form.addRow("Coefficients", _row(coeff_row))

        self.model_combo = QComboBox()
        self.model_combo.addItems(["simple", "guan", "fornasiero"])
        form.addRow("Model", self.model_combo)

        self.label_combo = QComboBox()
        self.label_combo.addItems(["hw", "o18"])
        form.addRow("Label", self.label_combo)

        self.ria_spin = QDoubleSpinBox()
        self.ria_spin.setDecimals(3)
        self.ria_spin.setRange(0.0, 1.0)
        self.ria_spin.setSingleStep(0.005)
        self.ria_spin.setValue(0.06)
        form.addRow("RIA max", self.ria_spin)

        self.depth_spin = QSpinBox()
        self.depth_spin.setRange(1, 100)
        self.depth_spin.setValue(3)
        form.addRow("Depth", self.depth_spin)

        self.qvalue_spin = QDoubleSpinBox()
        self.qvalue_spin.setDecimals(4)
        self.qvalue_spin.setRange(0.0, 1.0)
        self.qvalue_spin.setSingleStep(0.001)
        self.qvalue_spin.setValue(0.01)
        form.addRow("Max q-value", self.qvalue_spin)

        self.kp_spin = self._rate_spin(0.5)
        form.addRow("k_p (guan/forn.)", self.kp_spin)
        self.kr_spin = self._rate_spin(0.05)
        form.addRow("k_r (fornasiero)", self.kr_spin)
        self.rp_spin = self._rate_spin(10.0)
        form.addRow("r_p (fornasiero)", self.rp_spin)

        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, os.cpu_count() or 1)
        self.thread_spin.setValue(max(1, min(default_threads, os.cpu_count() or 1)))
        form.addRow("Threads", self.thread_spin)

        self.out_edit = QLineEdit(".")
        out_row = QHBoxLayout()
        out_row.addWidget(self.out_edit, stretch=1)
        out_browse = QPushButton("Browse…")
        out_browse.clicked.connect(self._pick_out)
        out_row.addWidget(out_browse)
        form.addRow("Output dir", _row(out_row))

        buttons = QHBoxLayout()
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._on_run)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self._on_cancel)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        form.addRow(_row(buttons))

        self.error_label = QLabel("")
        self.error_label.setStyleSheet("color: #b00020;")
        self.error_label.setWordWrap(True)
        form.addRow(self.error_label)
        return box

    def _rate_spin(self, value: float) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setDecimals(4)
        spin.setRange(0.0, 1000.0)
        spin.setSingleStep(0.05)
        spin.setValue(value)
        return spin

    def _build_results(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.summary_label = QLabel("No fit yet.")
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setFixedHeight(100)
        layout.addWidget(self.log)

        results_split = QSplitter(Qt.Orientation.Vertical)
        self.table = QTableView()
        self.model = DataFrameTableModel()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.selectionModel().currentRowChanged.connect(self._on_row_changed)
        results_split.addWidget(self.table)

        self.curve = CurveView()
        results_split.addWidget(self.curve)
        results_split.setSizes([320, 300])
        layout.addWidget(results_split, stretch=1)
        return panel

    # --- file / path pickers ----------------------------------------------- #
    def _add_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select integrate output files",
            filter="riana output (*_riana.txt *.txt);;All files (*)",
        )
        existing = {self.files_list.item(i).text()
                    for i in range(self.files_list.count())}
        for p in paths:
            if p not in existing:
                self.files_list.addItem(p)

    def _pick_coefficients(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select coefficient CSV", filter="CSV (*.csv);;All files (*)"
        )
        if path:
            self.coeff_combo.setCurrentText(path)

    def _pick_out(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.out_edit.setText(path)

    def _selected_files(self) -> list[str]:
        return [self.files_list.item(i).text()
                for i in range(self.files_list.count())]

    # --- config marshalling (shared-validation contract) -------------------- #
    def build_config(self) -> FitConfig:
        """Build the frozen :class:`FitConfig`; ``__post_init__`` is the shared
        validator (raises ``ValueError`` on a bad value, surfaced inline)."""
        return FitConfig(
            model=self.model_combo.currentText(),
            label=self.label_combo.currentText(),
            k_p=float(self.kp_spin.value()),
            k_r=float(self.kr_spin.value()),
            r_p=float(self.rp_spin.value()),
            q_value=float(self.qvalue_spin.value()),
            depth=int(self.depth_spin.value()),
            ria_max=float(self.ria_spin.value()),
            threads=int(self.thread_spin.value()),
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
        self.curve.show_placeholder("Fitting…")
        self._result_df = None
        self._cancelled = False

        try:
            config = self.build_config()
        except ValueError as exc:
            self._fail(str(exc))
            return

        files = self._selected_files()
        if not files:
            self._fail("Add at least one integrate output (_riana.txt) file.")
            return

        coefficients = self.coeff_combo.currentText().strip() or None
        # Same CLI-shaped guard riana.cli.fit applies; o18 is caught later in
        # fit_run with its "reimplemented post-M4" message.
        if config.label == "hw" and not coefficients:
            presets = " | ".join(available_coefficient_presets())
            self._fail(
                f"Coefficients are required for label 'hw'. Pick a preset "
                f"({presets}) or a (amino_acid, coefficient) CSV."
            )
            return

        os.makedirs(config.out_dir, exist_ok=True)
        self._set_running(True)
        self.progress.setRange(0, 0)  # busy
        loop = asyncio.get_running_loop()
        try:
            self._info(f"fitting {len(files)} timepoint file(s) …")
            self._future = loop.run_in_executor(
                self.pool, run_fit, config, files, coefficients
            )
            result_df = await self._future

            if self._cancelled:
                self._info("cancelled.")
                return

            self._write_output(config, result_df, files, coefficients)
            self._result_df = result_df
            self._last_config = config
            self._populate_results(result_df)
            self.curve.show_placeholder("Select a peptide row to view its fit.")
            n_fitted = int(result_df["k_deg"].notna().sum())
            n_well = int((result_df["R_squared"] >= 0.9).sum())
            self.summary_label.setText(
                f"{len(result_df)} peptides; {n_fitted} converged; {n_well} R²≥0.9."
            )
            self._info("done.")
        except asyncio.CancelledError:
            self._info("cancelled.")
        except Exception as exc:  # surface worker/IO/fit errors inline
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            self._future = None
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self._set_running(False)

    def _write_output(self, config, result_df, files, coefficients) -> None:
        """Write ``riana_fit_peptides.txt`` exactly as riana.cli.fit does."""
        out_path = Path(config.out_dir) / "riana_fit_peptides.txt"
        provenance = make_provenance(
            dataclasses.asdict(config),
            id_source=",".join(files),
            extra={"model": config.model, "label": config.label,
                   "coefficients": str(coefficients)},
        )
        write_dataframe_tsv(out_path, result_df, provenance, include_index=True)
        self._info(f"wrote {out_path}")

    def _populate_results(self, result_df: pd.DataFrame) -> None:
        display = result_df.reset_index()
        display = display[[c for c in _DISPLAY_COLS if c in display.columns]]
        self.model.set_dataframe(display)

    def _on_cancel(self) -> None:
        self._cancelled = True
        self.cancel_button.setEnabled(False)
        if self._future is not None:
            # Best effort: a ProcessPool future already running can't be killed,
            # but the result is discarded (the _cancelled flag short-circuits).
            self._future.cancel()
        self._info("cancelling …")

    # --- fitted curve on selection (pure, no worker) ------------------------ #
    def _on_row_changed(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if not current.isValid() or self._running or self._result_df is None:
            return
        concat = self.model.dataframe.iloc[current.row()]["concat"]
        if concat not in self._result_df.index:
            return
        row = self._result_df.loc[concat]
        cfg = self._last_config
        kinetic = dict(k_p=cfg.k_p, k_r=cfg.k_r, r_p=cfg.r_p)
        self.curve.plot_fit(
            str(concat), list(row["t"]), list(row["fs"]),
            float(row["k_deg"]), cfg.model, kinetic,
        )

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


def _row(layout) -> QWidget:
    """Wrap a layout in a QWidget so it can be added as a QFormLayout row."""
    w = QWidget()
    layout.setContentsMargins(0, 0, 0, 0)
    w.setLayout(layout)
    return w
