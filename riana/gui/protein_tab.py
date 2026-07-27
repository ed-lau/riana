# -*- coding: utf-8 -*-

"""The Protein tab: roll the `riana fit` outputs up to protein turnover.

Mirrors the Integrate/Model tabs' pattern — a form gathers the rollup
parameters, the (CPU-bound, bootstrapped) rollup runs on the shared
``ProcessPoolExecutor`` via the Qt-free :func:`riana.gui.tasks.run_rollup`
worker (the *same* :func:`riana.core.protein.rollup_proteins` the CLI calls, so
the surfaces cannot diverge), and the per-protein table lands in a view and on
disk as ``riana_rollup_proteins.txt`` (+ ``riana_rollup_fractions.txt``).

Reads a *manifest* (``riana_manifest.tsv`` — the SDRF/project path, the GUI's
only input, matching the Integrate SDRF and Model manifest fields): the fit
outputs are located from its ``stage="fit"`` rows via
:func:`riana.core.pipeline.fit_outputs_from_manifest`, and the rollup is written
next to the manifest with ``stage="rollup"`` rows recorded — so one manifest
drives the whole ``integrate → fit → rollup`` chain, exactly as the CLI's
``rollup --manifest`` does.
"""

from __future__ import annotations

import asyncio
import multiprocessing
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
from riana.gui.progress import ProgressPump
from riana.gui.tasks import load_rollup_results, run_rollup
from riana.io.writers import (
    ESTIMATE_FLOAT_FORMAT,
    make_provenance,
    write_dataframe_tsv,
)

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
        #: whether the selected manifest already has saved rollup results.
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

    def _build_form(self) -> QWidget:
        box = QGroupBox("Protein rollup")
        form = QFormLayout(box)
        self._form = form

        self.manifest_edit = QLineEdit("")
        self.manifest_edit.setToolTip(
            "riana_manifest.tsv from `integrate`/`fit` (the SDRF/project path). "
            "The fit outputs are located from its stage='fit' rows, and the "
            "rollup is written next to the manifest.")
        self.manifest_edit.setPlaceholderText(
            "riana_manifest.tsv (the SDRF path)")
        # Detect already-computed rollup results when a manifest is entered, so the
        # Load button + hint can offer to display them without re-running.
        self.manifest_edit.editingFinished.connect(self._check_for_saved_results)
        # Fill the linear-model Reference/Test condition dropdowns from the
        # manifest's conditions as soon as it is entered (SDRF/project path).
        self.manifest_edit.editingFinished.connect(self._populate_conditions)
        man_row = QHBoxLayout()
        man_row.addWidget(self.manifest_edit, stretch=1)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._pick_manifest)
        man_row.addWidget(browse)
        form.addRow("Manifest", _row(man_row))

        self.method_combo = QComboBox()
        self.method_combo.addItems(["weighted", "pooled"])
        self.method_combo.setToolTip(
            "weighted (default): biorep-aware per-timepoint inverse-variance "
            "collapse, then refit. pooled: all peptide×timepoint points "
            "(pseudoreplication; for comparison).")
        form.addRow("Method", self.method_combo)

        self.parsimony_combo = QComboBox()
        self.parsimony_combo.setToolTip(
            "Protein attribution: unique (single-accession peptides only) or isoform (fold isoform-shared peptides into the canonical entry).")
        self.parsimony_combo.addItems(["unique", "isoform"])
        form.addRow("Parsimony", self.parsimony_combo)

        self.model_combo = QComboBox()
        self.model_combo.addItems(["simple", "guan", "fornasiero", "linear simple"])
        self.model_combo.setToolTip(
            "simple/guan/fornasiero: nonlinear curve_fit. 'linear simple': the "
            "linearized φ=log(1−θ) cross-sample model — per-condition k + a Δk "
            "test (writes delta_k / delta_k_p_adj), plotted in φ-space.")
        self.model_combo.currentTextChanged.connect(self._on_model_changed)
        form.addRow("Model", self.model_combo)

        self.kp_spin = self._rate_spin(0.5)
        form.addRow("k_p (guan/forn.)", self.kp_spin)
        self.kr_spin = self._rate_spin(0.05)
        form.addRow("k_r (fornasiero)", self.kr_spin)
        self.rp_spin = self._rate_spin(10.0)
        form.addRow("r_p (fornasiero)", self.rp_spin)

        # 'linear simple' only: plateau truncation + Δk reference condition.
        self.phi_limit_spin = QDoubleSpinBox()
        self.phi_limit_spin.setDecimals(1)
        self.phi_limit_spin.setRange(-10.0, -1.0)
        self.phi_limit_spin.setSingleStep(0.5)
        self.phi_limit_spin.setValue(-4.0)
        self.phi_limit_spin.setToolTip(
            "['linear simple'] Plateau truncation: drop points with φ=log(1−θ) "
            "at/below this (saturated tail = noise, not slope). −4 ≈ θ 0.98.")
        form.addRow("φ-limit (linear)", self.phi_limit_spin)
        self.reference_combo = QComboBox()
        self.reference_combo.addItem("")  # blank = auto (alphabetically first)
        self.reference_combo.setToolTip(
            "['linear simple'] Δk baseline condition: delta_k = k(test) − "
            "k(reference). Blank = alphabetically first. Auto-filled from the "
            "manifest's conditions.")
        form.addRow("Reference cond. (linear)", self.reference_combo)
        self.test_combo = QComboBox()
        self.test_combo.addItem("")  # blank = auto (only a 2-condition protein gets Δk)
        self.test_combo.setToolTip(
            "['linear simple'] Δk comparison condition (needs a Reference). Pick a "
            "pair to contrast exactly those two even when >2 conditions are present "
            "— an interim before full all-pairwise. NOTE: the joint fit still pools "
            "the residual variance over ALL conditions in the project, so include "
            "only the conditions you mean to compare in the SDRF. Blank = auto (Δk "
            "only for a protein with exactly two conditions).")
        form.addRow("Test cond. (linear)", self.test_combo)

        self.min_peptides_spin = QSpinBox()
        self.min_peptides_spin.setToolTip(
            "Minimum attributed peptides for a protein to be reported.")
        self.min_peptides_spin.setRange(1, 1000)
        self.min_peptides_spin.setValue(2)
        form.addRow("Min peptides", self.min_peptides_spin)

        self.min_points_spin = QSpinBox()
        self.min_points_spin.setToolTip(
            "Minimum collapsed (t, θ) points for the protein refit.")
        self.min_points_spin.setRange(2, 1000)
        self.min_points_spin.setValue(3)
        form.addRow("Min refit points", self.min_points_spin)

        # Peptide-level biological-replicate floor — distinct from "Min refit points"
        # (that is the protein-level (t, θ) floor). The dominant single-timepoint
        # curation lever. 0 = auto (2 for single-timepoint, off for a time series).
        self.min_fit_points_spin = QSpinBox()
        self.min_fit_points_spin.setToolTip(
            "Peptide biological-replicate floor: keep only peptidoforms fit on ≥ N "
            "distinct (biorep, timepoint) points. 0 = auto — 2 for a single-timepoint "
            "experiment (its dominant curation lever, where n_points IS the replicate "
            "count) and off for a time series; set 1 to disable, 3+ to require more "
            "replicates. Distinct from 'Min refit points' above (the protein-level "
            "(t, θ) refit floor).")
        self.min_fit_points_spin.setRange(0, 1000)
        self.min_fit_points_spin.setValue(0)
        form.addRow("Min fit points (0 = auto)", self.min_fit_points_spin)

        # Peptide R² admission gate; default 0.8 (needed for good geom-CV), 0 = off.
        self.min_r2_spin = QDoubleSpinBox()
        self.min_r2_spin.setDecimals(2)
        self.min_r2_spin.setRange(0.0, 1.0)
        self.min_r2_spin.setSingleStep(0.05)
        self.min_r2_spin.setValue(0.8)
        self.min_r2_spin.setToolTip(
            "Peptide R² gate before rollup (default 0.8 — needed for good "
            "within-protein geom-CV; inverse-variance weighting alone "
            "under-curates). Well-measured flat-curve peptides are still rescued "
            "via the Max k_cv gate below. Set 0 to disable the gate.")
        form.addRow("Min R² (0 = off)", self.min_r2_spin)

        # Flat-curve rescue (only with Min R² > 0): admit a low-R² peptide whose
        # rate constant is nonetheless tightly determined — k_cv = relative
        # uncertainty of k̂ = (ci_hi−ci_lo)/(2·|k|), a scale-free CV.
        self.k_cv_spin = QDoubleSpinBox()
        self.k_cv_spin.setDecimals(2)
        self.k_cv_spin.setRange(0.0, 10.0)
        self.k_cv_spin.setSingleStep(0.05)
        self.k_cv_spin.setValue(0.2)
        self.k_cv_spin.setToolTip(
            "Max relative uncertainty of k̂: k_cv = (ci_hi−ci_lo)/(2·|k|), a scale-free "
            "CV — no retuning across time ranges or k units. Its ROLE depends on the "
            "data: for a MULTI-timepoint series it is a SECONDARY flat-curve rescue "
            "(admit a low-R² peptide whose k is nonetheless tight — only when Min R² > 0); "
            "for a SINGLE-timepoint experiment R² is bypassed, so this becomes the "
            "PRIMARY gate (alongside Min fit points). 0 = off.")
        form.addRow("Max k_cv (0 = off)", self.k_cv_spin)

        self.rescue_r2_spin = QDoubleSpinBox()
        self.rescue_r2_spin.setDecimals(2)
        self.rescue_r2_spin.setRange(0.0, 1.0)
        self.rescue_r2_spin.setSingleStep(0.05)
        self.rescue_r2_spin.setValue(0.6)
        self.rescue_r2_spin.setToolTip(
            "R² floor for the Max k_cv rescue — guards against degenerate k≈0 "
            "fits whose CI collapses to a spuriously tight k_cv. 0.6 is "
            "validated; the exact value barely matters above the negative-R² band.")
        form.addRow("Rescue R² floor", self.rescue_r2_spin)

        self.workers_spin = QSpinBox()
        self.workers_spin.setRange(1, os.cpu_count() or 1)
        self.workers_spin.setValue(1)
        self.workers_spin.setToolTip(
            "Worker *processes* for the per-protein refit — the real lever for "
            "the GIL-bound rollup. Result is identical regardless of N. When >1 "
            "the rollup runs off the shared pool (no nested pools).")
        form.addRow("Workers", self.workers_spin)

        self.out_edit = QLineEdit(".")
        self.out_edit.setToolTip(
            "Output directory for riana_rollup_proteins.txt (ignored on the manifest path).")
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
        self.load_button = QPushButton("Load results")
        self.load_button.setEnabled(False)
        self.load_button.setToolTip(
            "Display the rollup results already saved next to this manifest "
            "(the stage='rollup' rows) without re-running. 'Run' recomputes with "
            "the current form settings instead.")
        self.load_button.clicked.connect(self._on_load)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        buttons.addWidget(self.load_button)
        form.addRow(_row(buttons))

        self.results_hint = QLabel("")
        self.results_hint.setStyleSheet("color: palette(mid);")
        self.results_hint.setWordWrap(True)
        form.addRow(self.results_hint)

        self.error_label = QLabel("")
        self.error_label.setStyleSheet("color: #b00020;")
        self.error_label.setWordWrap(True)
        form.addRow(self.error_label)
        self._on_model_changed(self.model_combo.currentText())  # initial visibility
        return box

    def _rate_spin(self, value: float) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setDecimals(4)
        spin.setRange(0.0, 1000.0)
        spin.setSingleStep(0.05)
        spin.setValue(value)
        return spin

    def _on_model_changed(self, model: str) -> None:
        """Show the ODE rate knobs for the nonlinear models, the φ knobs for the
        linear model — they are mutually exclusive."""
        linear = model == "linear simple"
        for w in (self.kp_spin, self.kr_spin, self.rp_spin):
            self._form.setRowVisible(w, not linear)
        for w in (self.phi_limit_spin, self.reference_combo, self.test_combo):
            self._form.setRowVisible(w, linear)

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
        self.table.setSortingEnabled(True)
        self.table.selectionModel().currentRowChanged.connect(self._on_row_changed)
        results_split.addWidget(self.table)

        self.curve = CurveView()
        results_split.addWidget(self.curve)
        results_split.setSizes([320, 300])
        layout.addWidget(results_split, stretch=1)
        return panel

    # --- file / path pickers ----------------------------------------------- #
    def _pick_manifest(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select riana_manifest.tsv",
            filter="Manifest (*.tsv);;All files (*)")
        if path:
            self.manifest_edit.setText(path)
            self._check_for_saved_results()

    def _pick_out(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder")
        if path:
            self.out_edit.setText(path)

    def _populate_conditions(self) -> None:
        """Fill the linear-model Reference/Test dropdowns from the manifest's
        conditions (blank first item = auto). Best-effort and Qt-cheap: a
        missing/malformed manifest just leaves the blank option, and a prior
        still-valid pick is preserved across a re-read."""
        conditions: list[str] = []
        path = self.manifest_edit.text().strip()
        if path and Path(path).is_file():
            try:
                from riana.io.manifest import read_manifest
                conditions = sorted({
                    r.identity.condition for r in read_manifest(path)
                    if r.identity.condition})
            except Exception:
                conditions = []
        for combo in (self.reference_combo, self.test_combo):
            prior = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("")               # blank = auto
            combo.addItems(conditions)
            keep = combo.findText(prior)
            combo.setCurrentIndex(keep if keep >= 0 else 0)
            combo.blockSignals(False)

    # --- params marshalling ------------------------------------------------- #
    def build_params(self) -> dict:
        """Gather the rollup parameters from the form (testable, Qt-free dict)."""
        r2 = float(self.min_r2_spin.value())
        mfp = int(self.min_fit_points_spin.value())   # 0 = auto (None)
        return {
            "manifest": self.manifest_edit.text().strip(),
            "method": self.method_combo.currentText(),
            "parsimony": self.parsimony_combo.currentText(),
            "model": self.model_combo.currentText(),
            "kp": float(self.kp_spin.value()),
            "kr": float(self.kr_spin.value()),
            "rp": float(self.rp_spin.value()),
            "min_peptides": int(self.min_peptides_spin.value()),
            "min_points": int(self.min_points_spin.value()),
            "min_fit_points": (mfp if mfp > 0 else None),   # 0 = auto
            "min_r2": (r2 if r2 > 0.0 else None),   # 0 = off
            "k_cv_max": float(self.k_cv_spin.value()),
            "rescue_r2": float(self.rescue_r2_spin.value()),
            "workers": int(self.workers_spin.value()),
            "phi_limit": float(self.phi_limit_spin.value()),
            "reference_condition": self.reference_combo.currentText().strip() or None,
            "test_condition": self.test_combo.currentText().strip() or None,
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
        manifest = Path(p["manifest"]) if p["manifest"] else None
        if manifest is None or not manifest.is_file():
            self._fail("Pick a riana_manifest.tsv (from `integrate`/`fit` on the "
                       "SDRF path). The legacy fit-directory rollup is CLI-only.")
            return
        if p["test_condition"] and not p["reference_condition"]:
            self._fail("Pick a Reference condition too — the Test condition is the "
                       "Δk comparison measured against a baseline.")
            return
        if (p["test_condition"] and p["reference_condition"]
                and p["test_condition"] == p["reference_condition"]):
            self._fail("Pick two different conditions — Reference and Test are both "
                       f"'{p['reference_condition']}', and a condition compared to "
                       "itself has no Δk.")
            return
        # Locate the fit outputs from the manifest's stage='fit' rows — the same
        # resolver `rollup --manifest` uses. Raises DataError (no fit rows yet /
        # missing outputs) which we surface inline. The manifest read is tiny, so
        # it stays on the UI thread; the heavy rollup goes to the pool below.
        from riana.core.pipeline import fit_outputs_from_manifest
        from riana.exceptions import DataError
        try:
            pep_path, _frac_path = fit_outputs_from_manifest(str(manifest))
        except DataError as exc:
            self._fail(str(exc))
            return
        fit_dir = Path(pep_path).parent

        os.makedirs(manifest.resolve().parent, exist_ok=True)
        self._set_running(True)
        self.progress.setRange(0, 0)  # busy until the first progress update lands
        loop = asyncio.get_running_loop()
        # A Manager queue carries (done, total) back from the worker (pool process
        # or -W thread); a main-thread QTimer drains it into a determinate bar.
        manager = multiprocessing.Manager()
        progress_q = manager.Queue()
        pump = ProgressPump(self.progress, progress_q, parent=self)
        pump.start()
        try:
            self._info(f"rolling up from manifest {manifest} (model={p['model']}, "
                       f"parsimony={p['parsimony']}) …")
            # workers>1 spawns a ProcessPool inside rollup_proteins; run it on a
            # main-process thread (executor=None) so that pool is NOT nested
            # inside a shared-pool worker (which breaks: BrokenProcessPool).
            executor = None if p["workers"] > 1 else self.pool
            self._future = loop.run_in_executor(
                executor, run_rollup, str(fit_dir), p["model"],
                p["kp"], p["kr"], p["rp"], p["parsimony"],
                p["min_peptides"], p["min_points"], p["min_fit_points"], p["min_r2"],
                p["k_cv_max"], p["rescue_r2"], p["method"],
                p["workers"], p["phi_limit"], p["reference_condition"],
                p["test_condition"],
                progress_q,
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
            # Surface the single-timepoint regime — the core logs it, but that INFO is
            # emitted in a worker process and never reaches the GUI. If every rolled
            # protein sits at one labeling timepoint, R² was N/A and curation rode on
            # k_cv + Min fit points.
            single_tp = (
                "n_timepoints" in result.columns and len(result) > 0
                and bool((result["n_timepoints"].dropna() <= 1).all())
            )
            self.results_hint.setText(
                "Single labeling timepoint — R² is not applicable (bypassed); curation "
                "was on Max k_cv + Min fit points (replicate floor). Min R² / Rescue R² "
                "do not apply." if single_tp else "")
            self._info("done.")
        except asyncio.CancelledError:
            self._info("cancelled.")
        except Exception as exc:  # surface worker/IO/rollup errors inline
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            pump.stop()          # final drain while the queue proxy is still live
            manager.shutdown()
            self._future = None
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self._set_running(False)

    def _write_output(self, result: pd.DataFrame, params: dict) -> None:
        # Default / same-folder Output dir updates the project in place next to the
        # manifest (recording stage='rollup' rows — the integrate→fit→rollup
        # chain); a *different* Output dir forks a self-contained derived project
        # there, leaving the input manifest untouched. Mirrors `rollup --manifest`.
        from riana.core.pipeline import record_stage_rows, resolve_manifest_write
        from riana.core.protein import build_rollup_fractions

        manifest = Path(params["manifest"])
        out_dir, target_manifest = resolve_manifest_write(
            manifest, params["out_dir"], "rollup")
        if Path(out_dir).resolve() != manifest.resolve().parent:
            self._info(f"forking a derived project into {out_dir} "
                       f"(input manifest left untouched)")
        provenance = make_provenance(
            {k: params[k] for k in (
                "model", "method", "parsimony", "kp", "kr", "rp",
                "min_peptides", "min_points", "min_r2", "k_cv_max", "rescue_r2",
                "phi_limit", "reference_condition")},
            id_source=params["manifest"],
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

        record_stage_rows(str(target_manifest), "rollup", written, result, provenance)
        self._info(f"recorded {len(written)} rollup rows in {target_manifest}")

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
        p = self._last_params or {}
        if p.get("model") == "linear simple":
            self._plot_linear_row(row, p)
            return
        key = (row["experiment"], row["condition"], row["protein"])
        pts = self._points.get(key)
        k = row.get("k_deg")
        if pts is None or k is None or pd.isna(k):
            self.curve.show_placeholder(f"{row['protein']}: no refit to show.")
            return
        t_list, fs_list = pts[0], pts[1]      # (t, fs[, var, df]) — curve needs t/fs
        kinetic = dict(k_p=p.get("kp", 0.5), k_r=p.get("kr", 0.05),
                       r_p=p.get("rp", 10.0))
        self.curve.plot_fit(
            str(row["protein"]), list(t_list), list(fs_list),
            float(k), p.get("model", "simple"), kinetic,
            ci_lo=_safe_float(row.get("ci_lo")),
            ci_hi=_safe_float(row.get("ci_hi")),
        )

    def _plot_linear_row(self, row, p: dict) -> None:
        """Overlay every condition of the selected protein in φ-space (the Δk
        view): each condition's clearance points + its through-origin k line."""
        exp, prot = row["experiment"], row["protein"]
        sub = self._result_df[
            (self._result_df["experiment"] == exp)
            & (self._result_df["protein"] == prot)
        ]
        per_condition: dict = {}
        for _, r in sub.iterrows():
            pts = self._points.get((exp, r["condition"], prot))
            if pts is None:
                continue
            t_list, theta_list = pts[0], pts[1]      # (t, θ[, var, df])
            per_condition[str(r["condition"])] = (
                list(t_list), list(theta_list), r.get("k_deg"),
                _safe_float(r.get("ci_lo")), _safe_float(r.get("ci_hi")))
        if not per_condition:
            self.curve.show_placeholder(f"{prot}: no points to show.")
            return
        self.curve.plot_linear(
            str(prot), per_condition,
            phi_limit=float(p.get("phi_limit", -4.0)),
            delta_k=row.get("delta_k"), delta_k_p_adj=row.get("delta_k_p_adj"),
        )

    # --- load prior results (no recompute) ---------------------------------- #
    def _check_for_saved_results(self) -> None:
        """Enable the Load button + hint when the manifest already has a rollup.

        A cheap manifest read (``stage="rollup"`` rows) on manifest entry / after
        a run, so a prior rollup can be displayed without recomputing it.
        """
        manifest = self.manifest_edit.text().strip()
        found = False
        if manifest and Path(manifest).is_file():
            try:
                from riana.io.manifest import read_manifest
                found = bool(read_manifest(manifest, stage="rollup"))
            except Exception:
                found = False
        self._results_available = found
        self.load_button.setEnabled(found and not self._running)
        self.results_hint.setText(
            "✓ saved rollup results found — Load to view, or Run to recompute"
            if found else "")

    @asyncSlot()
    async def _on_load(self) -> None:
        """Display the rollup results already saved next to the manifest."""
        if self._running:
            return
        manifest = self.manifest_edit.text().strip()
        if not manifest or not Path(manifest).is_file():
            self._fail("Pick a riana_manifest.tsv first.")
            return
        self.error_label.setText("")
        self._set_running(True)
        self.progress.setRange(0, 0)  # busy
        loop = asyncio.get_running_loop()
        try:
            self._info(f"loading saved rollup results from {manifest} …")
            proteins, points, header = await loop.run_in_executor(
                self.pool, load_rollup_results, manifest)
            # The refit / φ-space curve needs the model that produced the results;
            # take it from the provenance header, the rest from the form.
            params = self.build_params()
            params["model"] = header.get("model", params["model"])
            self._result_df = proteins
            self._points = points
            self._last_params = params
            self.model.set_dataframe(proteins)
            self.curve.show_placeholder("Select a protein row to view its refit.")
            n_fit = int(proteins["k_deg"].notna().sum())
            self.summary_label.setText(
                f"Loaded {len(proteins)} proteins from saved results "
                f"(model={header.get('model', '?')}, "
                f"method={header.get('method', '?')}); {n_fit} with a fitted k_deg.")
            self._info("loaded saved rollup results.")
        except Exception as exc:  # surface loader/IO errors inline
            self._fail(f"{type(exc).__name__}: {exc}")
        finally:
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self._set_running(False)

    # --- small helpers ------------------------------------------------------ #
    def _set_running(self, running: bool) -> None:
        self._running = running
        self.run_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)
        if running:
            self.load_button.setEnabled(False)
        else:  # re-enable Load iff the manifest has saved results (post run/load)
            self._check_for_saved_results()

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


def _safe_float(value) -> float | None:
    """Coerce a (possibly missing / NaN) cell to float, else None."""
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if f != f else f  # drop NaN
