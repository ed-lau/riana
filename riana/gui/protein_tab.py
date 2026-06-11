# -*- coding: utf-8 -*-

"""The Protein tab: roll the `riana fit` outputs up to protein turnover.

Mirrors the Integrate/Model tabs' pattern — a form gathers the rollup
parameters, the (CPU-bound, bootstrapped) rollup runs on the shared
``ProcessPoolExecutor`` via the Qt-free :func:`riana.gui.tasks.run_rollup`
worker (the *same* :func:`riana.core.protein.rollup_proteins` the CLI calls, so
the surfaces cannot diverge), and the per-protein table lands in a view and on
disk as ``riana_rollup_proteins.txt`` (+ ``riana_rollup_fractions.txt``).

Reads a *fit output directory* (the ``riana_fit_peptides.txt`` +
``riana_fit_fractions.txt`` a `riana fit` run wrote). The per-protein refit
curve view is a follow-up — the collapsed ``(t, θ)`` points would have to be
returned explicitly from the worker (``DataFrame.attrs`` does not reliably
survive the pickle back from the pool).
"""

from __future__ import annotations

import asyncio
import os
from concurrent.futures import Future
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtWidgets import (
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
    QSpinBox,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from qasync import asyncSlot
from PySide6.QtCore import QModelIndex, Qt

from riana.gui.curve_view import CurveView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import run_rollup
from riana.io.writers import (
    ESTIMATE_FLOAT_FORMAT,
    make_provenance,
    write_dataframe_tsv,
)

_PEPTIDES_FILE = "riana_fit_peptides.txt"
_FRACTIONS_FILE = "riana_fit_fractions.txt"


class ProteinTab(QWidget):
    """Form + async rollup runner + per-protein results table."""

    def __init__(
        self,
        pool,
        status_cb: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__()
        self.pool = pool
        self._status_cb = status_cb or (lambda _msg: None)
        self._running = False
        self._cancelled = False
        self._future: Future | None = None
        self._result_df: pd.DataFrame | None = None
        self._points: dict = {}
        self._last_params: dict | None = None
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

    def _build_form(self) -> QWidget:
        box = QGroupBox("Protein rollup")
        form = QFormLayout(box)

        self.fit_dir_edit = QLineEdit("")
        self.fit_dir_edit.setPlaceholderText("folder with riana_fit_*.txt")
        dir_row = QHBoxLayout()
        dir_row.addWidget(self.fit_dir_edit, stretch=1)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._pick_fit_dir)
        dir_row.addWidget(browse)
        form.addRow("Fit output dir", _row(dir_row))

        self.method_combo = QComboBox()
        self.method_combo.addItems(["weighted", "pooled"])
        self.method_combo.setToolTip(
            "weighted (default): biorep-aware per-timepoint inverse-variance "
            "collapse, then refit. pooled: all peptide×timepoint points "
            "(pseudoreplication; for comparison).")
        form.addRow("Method", self.method_combo)

        self.parsimony_combo = QComboBox()
        self.parsimony_combo.addItems(["unique", "isoform"])
        form.addRow("Parsimony", self.parsimony_combo)

        self.model_combo = QComboBox()
        self.model_combo.addItems(["simple", "guan", "fornasiero"])
        form.addRow("Model", self.model_combo)

        self.kp_spin = self._rate_spin(0.5)
        form.addRow("k_p (guan/forn.)", self.kp_spin)
        self.kr_spin = self._rate_spin(0.05)
        form.addRow("k_r (fornasiero)", self.kr_spin)
        self.rp_spin = self._rate_spin(10.0)
        form.addRow("r_p (fornasiero)", self.rp_spin)

        self.min_peptides_spin = QSpinBox()
        self.min_peptides_spin.setRange(1, 1000)
        self.min_peptides_spin.setValue(2)
        form.addRow("Min peptides", self.min_peptides_spin)

        self.min_points_spin = QSpinBox()
        self.min_points_spin.setRange(2, 1000)
        self.min_points_spin.setValue(3)
        form.addRow("Min refit points", self.min_points_spin)

        # Optional peptide R² admission gate; 0 = off (inverse-variance only).
        self.min_r2_spin = QDoubleSpinBox()
        self.min_r2_spin.setDecimals(2)
        self.min_r2_spin.setRange(0.0, 1.0)
        self.min_r2_spin.setSingleStep(0.05)
        self.min_r2_spin.setValue(0.0)
        self.min_r2_spin.setToolTip(
            "Peptide R² gate before rollup. 0 = off (the inverse-variance "
            "weighting already down-weights noisy peptides). Slow-turnover "
            "peptides are still admitted via k ≤ 0.025 & SE ≤ 0.05.")
        form.addRow("Min R² (0 = off)", self.min_r2_spin)

        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, os.cpu_count() or 1)
        self.thread_spin.setValue(1)
        self.thread_spin.setToolTip(
            "Worker threads for the per-protein refit (result is identical "
            "regardless of thread count).")
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

        self.summary_label = QLabel("No rollup yet.")
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
    def _pick_fit_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select fit output folder")
        if path:
            self.fit_dir_edit.setText(path)

    def _pick_out(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.out_edit.setText(path)

    # --- params marshalling ------------------------------------------------- #
    def build_params(self) -> dict:
        """Gather the rollup parameters from the form (testable, Qt-free dict)."""
        r2 = float(self.min_r2_spin.value())
        return {
            "fit_dir": self.fit_dir_edit.text().strip(),
            "method": self.method_combo.currentText(),
            "parsimony": self.parsimony_combo.currentText(),
            "model": self.model_combo.currentText(),
            "kp": float(self.kp_spin.value()),
            "kr": float(self.kr_spin.value()),
            "rp": float(self.rp_spin.value()),
            "min_peptides": int(self.min_peptides_spin.value()),
            "min_points": int(self.min_points_spin.value()),
            "min_r2": (r2 if r2 > 0.0 else None),   # 0 = off
            "threads": int(self.thread_spin.value()),
            "out_dir": self.out_edit.text().strip() or ".",
        }

    # --- run flow ----------------------------------------------------------- #
    @asyncSlot()
    async def _on_run(self) -> None:
        if self._running:
            return
        self.error_label.setText("")
        self.log.clear()
        self.model.set_dataframe(pd.DataFrame())
        self.curve.show_placeholder("Rolling up…")
        self._result_df = None
        self._points = {}
        self._cancelled = False

        p = self.build_params()
        fit_dir = Path(p["fit_dir"]) if p["fit_dir"] else None
        if fit_dir is None or not fit_dir.is_dir():
            self._fail("Pick a fit output directory.")
            return
        for name in (_PEPTIDES_FILE, _FRACTIONS_FILE):
            if not (fit_dir / name).exists():
                self._fail(f"{name} not found in {fit_dir}. Run a fit there first.")
                return

        os.makedirs(p["out_dir"], exist_ok=True)
        self._set_running(True)
        self.progress.setRange(0, 0)  # busy
        loop = asyncio.get_running_loop()
        try:
            self._info(f"rolling up {fit_dir} (parsimony={p['parsimony']}) …")
            self._future = loop.run_in_executor(
                self.pool, run_rollup, str(fit_dir), p["model"],
                p["kp"], p["kr"], p["rp"], p["parsimony"],
                p["min_peptides"], p["min_points"], p["min_r2"],
                0.025, 0.05, p["threads"], p["method"],
            )
            result, points = await self._future
            if self._cancelled:
                self._info("cancelled.")
                return
            self._write_output(result, p)
            self._result_df = result
            self._points = points
            self._last_params = p
            self.model.set_dataframe(result)
            self.curve.show_placeholder("Select a protein row to view its refit.")
            n_fit = int(result["k_deg"].notna().sum())
            self.summary_label.setText(
                f"{len(result)} proteins ({p['method']}); "
                f"{n_fit} with a fitted k_deg."
            )
            self._info("done.")
        except asyncio.CancelledError:
            self._info("cancelled.")
        except Exception as exc:  # surface worker/IO/rollup errors inline
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            self._future = None
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self._set_running(False)

    def _write_output(self, result: pd.DataFrame, params: dict) -> None:
        # When the fit dir is a project (carries a manifest), write the rollup
        # outputs there and record stage='rollup' rows (the project chain);
        # otherwise honor the Output dir.
        from riana.core.protein import build_rollup_fractions

        fit_dir = Path(params["fit_dir"])
        manifest = fit_dir / "riana_manifest.tsv"
        out_dir = fit_dir if manifest.exists() else Path(params["out_dir"])
        provenance = make_provenance(
            {k: params[k] for k in (
                "model", "method", "parsimony", "kp", "kr", "rp",
                "min_peptides", "min_points", "min_r2")},
            id_source=params["fit_dir"],
            extra={"method": params["method"], "parsimony": params["parsimony"],
                   "model": params["model"]},
        )
        out_path = out_dir / "riana_rollup_proteins.txt"
        write_dataframe_tsv(out_path, result, provenance, include_index=False,
                            float_format=ESTIMATE_FLOAT_FORMAT)
        self._info(f"wrote {out_path}")
        written = [out_path]

        rollup_fractions = build_rollup_fractions(result)
        if not rollup_fractions.empty:
            frac_path = out_dir / "riana_rollup_fractions.txt"
            write_dataframe_tsv(frac_path, rollup_fractions, provenance,
                                include_index=False,
                                float_format=ESTIMATE_FLOAT_FORMAT)
            self._info(f"wrote {frac_path} ({len(rollup_fractions)} points)")
            written.append(frac_path)

        if manifest.exists():
            from riana.core.pipeline import record_stage_rows
            record_stage_rows(manifest, "rollup", written, result, provenance)
            self._info(f"recorded {len(written)} rollup rows in {manifest}")

    def _on_cancel(self) -> None:
        self._cancelled = True
        self.cancel_button.setEnabled(False)
        if self._future is not None:
            self._future.cancel()
        self._info("cancelling …")

    # --- per-protein refit curve on selection (pure, no worker) ------------- #
    def _on_row_changed(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if not current.isValid() or self._running or self._result_df is None:
            return
        row = self.model.dataframe.iloc[current.row()]
        key = (row["experiment"], row["condition"], row["protein"])
        pts = self._points.get(key)
        k = row.get("k_deg")
        if pts is None or k is None or pd.isna(k):
            self.curve.show_placeholder(f"{row['protein']}: no refit to show.")
            return
        t_list, fs_list = pts
        p = self._last_params or {}
        kinetic = dict(k_p=p.get("kp", 0.5), k_r=p.get("kr", 0.05),
                       r_p=p.get("rp", 10.0))
        self.curve.plot_fit(
            str(row["protein"]), list(t_list), list(fs_list),
            float(k), p.get("model", "simple"), kinetic,
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
