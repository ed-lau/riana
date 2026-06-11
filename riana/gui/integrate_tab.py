# -*- coding: utf-8 -*-

"""The Integrate tab: build an :class:`IntegrationConfig` from a form, run the
integration asynchronously, and show progress / log / results / per-fraction
drift, with a click-to-inspect chromatogram.

Two intake paths, matching the CLI: with an **SDRF** the search-ID file is read
as the quantms mzTab and routed through :mod:`riana.core.pipeline` —
:func:`~riana.core.pipeline.plan_integration` (in a worker) builds the per-run
tasks, this tab dispatches each over its *own* shared pool, then
:func:`~riana.core.pipeline.finalize_run` writes one identity-stamped
``<stem>_riana.txt`` per run and a ``riana_manifest.tsv``. Without an SDRF it is
the demoted single-mzML Percolator path. Both are **file-parallel** (the
*Workers* control, bounded by one mzML in memory per concurrent run) on top of
the per-run *Threads*.

The form builds the *same* frozen :class:`~riana.config.IntegrationConfig` the
CLI builds, so its ``__post_init__`` is the single shared validator. CPU work is
awaited on the shared ``ProcessPoolExecutor`` via the Qt-free
:mod:`riana.gui.tasks` workers so the UI stays responsive.
"""

from __future__ import annotations

import asyncio
import dataclasses
import json
import os
import re
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtCore import QModelIndex, Qt
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

from riana.config import IntegrationConfig
from riana.core.pipeline import finalize_run
from riana.gui.chromatogram import ChromatogramView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import (
    extract_trace,
    integrate_fraction,
    plan_sdrf_integration,
    read_psms,
)
from riana.io.manifest import MANIFEST_FILENAME, append_manifest
from riana.io.mzml import list_mzml_files, mzml_stem
from riana.io.percolator import file_indices, fraction_psms
from riana.io.writers import make_provenance, write_dataframe_tsv
from riana.records import PSMRecord


class IntegrateTab(QWidget):
    """Form + async runner + results/chromatogram for ``riana integrate``."""

    def __init__(
        self,
        pool,
        default_threads: int = 1,
        status_cb: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__()
        self.pool = pool
        self._status_cb = status_cb or (lambda _msg: None)
        self._cancelled = False
        self._running = False
        # file_idx -> mzML path, for on-demand chromatogram extraction.
        self._fraction_mzml: dict[int, str] = {}

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

    def _path_row(self, edit: QLineEdit, on_browse: Callable[[], None]) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        button = QPushButton("Browse…")
        button.clicked.connect(on_browse)
        layout.addWidget(button)
        return row

    def _build_form(self, default_threads: int) -> QWidget:
        box = QGroupBox("Integrate")
        form = QFormLayout(box)

        self.mzml_edit = QLineEdit()
        self.mzml_edit.setPlaceholderText("Folder containing mzML file(s)")
        form.addRow("mzML folder", self._path_row(self.mzml_edit, self._pick_mzml))

        self.id_edit = QLineEdit()
        self.id_edit.setPlaceholderText("percolator psms.txt — or the mzTab when an SDRF is set")
        form.addRow("Search ID", self._path_row(self.id_edit, self._pick_id))

        self.sdrf_edit = QLineEdit()
        self.sdrf_edit.setPlaceholderText("Optional: SDRF .tsv — enables the identity/manifest path")
        form.addRow("SDRF", self._path_row(self.sdrf_edit, self._pick_sdrf))

        self.sample_edit = QLineEdit("time0")
        self.sample_edit.setToolTip(
            "Percolator (no-SDRF) path only — must end in a digit (e.g. time0). "
            "Ignored when an SDRF is set (identity comes from the SDRF).")
        form.addRow("Sample", self.sample_edit)

        self.iso_edit = QLineEdit("0 1 2 3 4 5")
        form.addRow("Isotopomers", self.iso_edit)

        self.qvalue_spin = QDoubleSpinBox()
        self.qvalue_spin.setDecimals(4)
        self.qvalue_spin.setRange(0.0, 1.0)
        self.qvalue_spin.setSingleStep(0.001)
        self.qvalue_spin.setValue(0.01)
        form.addRow("Max q-value", self.qvalue_spin)

        self.mass_tol_spin = QSpinBox()
        self.mass_tol_spin.setRange(1, 500)
        self.mass_tol_spin.setValue(50)
        self.mass_tol_spin.setSuffix(" ppm")
        form.addRow("Mass tolerance", self.mass_tol_spin)

        self.peak_rt_combo = QComboBox()
        self.peak_rt_combo.addItems(["apex", "ms2", "consensus"])
        form.addRow("Window anchor", self.peak_rt_combo)

        self.ihw_edit = QLineEdit("0.15")
        self.ihw_edit.setPlaceholderText("RT min, or 'auto'")
        form.addRow("Integration ½-width", self.ihw_edit)

        self.baseline_combo = QComboBox()
        self.baseline_combo.addItems(["none", "noise_floor", "snip", "asls"])
        form.addRow("Baseline", self.baseline_combo)

        self.apex_combo = QComboBox()
        self.apex_combo.addItems(["tallest", "nearest"])
        form.addRow("Apex selection", self.apex_combo)

        self.ppm_alert_spin = QDoubleSpinBox()
        self.ppm_alert_spin.setRange(0.1, 1000.0)
        self.ppm_alert_spin.setValue(20.0)
        self.ppm_alert_spin.setSuffix(" ppm")
        form.addRow("Drift alert", self.ppm_alert_spin)

        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, os.cpu_count() or 1)
        self.thread_spin.setValue(max(1, min(default_threads, os.cpu_count() or 1)))
        self.thread_spin.setToolTip("Per-run peptide threads (within one mzML).")
        form.addRow("Threads", self.thread_spin)

        self.workers_spin = QSpinBox()
        self.workers_spin.setRange(1, os.cpu_count() or 1)
        self.workers_spin.setValue(1)
        self.workers_spin.setToolTip(
            "Runs/files integrated concurrently (one mzML in memory each; "
            "2-4 suits a many-timepoint series). Distinct from Threads.")
        form.addRow("Workers (files)", self.workers_spin)

        self.out_edit = QLineEdit(".")
        form.addRow("Output dir", self._path_row(self.out_edit, self._pick_out))

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
            self, "Select search-ID file (Percolator psms or mzTab)",
            filter="Search ID (*.txt *.mzTab);;All files (*)"
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
        iso_text = self.iso_edit.text().strip()
        isotopomers = tuple(
            sorted({int(p) for p in re.split(r"[,\s]+", iso_text) if p})
        )
        if not isotopomers:
            raise ValueError("Isotopomers must list at least one value, e.g. '0 6'.")

        ihw_text = self.ihw_edit.text().strip()
        ihw: float | str = "auto" if ihw_text == "auto" else float(ihw_text)
        peak_rt = self.peak_rt_combo.currentText()
        # Derive the extraction half-width exactly as riana.cli.integrate does:
        # ms2 integrates the whole extraction (= ihw); apex/consensus and 'auto'
        # need room for the apex offset (+0.33).
        if peak_rt == "ms2" and ihw != "auto":
            ehw = float(ihw)
        else:
            ehw = (0.33 if ihw == "auto" else float(ihw)) + 0.33

        return IntegrationConfig(
            sample=self.sample_edit.text().strip(),
            isotopomers=isotopomers,
            mass_tol_ppm=int(self.mass_tol_spin.value()),
            extraction_half_width=ehw,
            peak_rt=peak_rt,
            integration_half_width=ihw,
            baseline_method=self.baseline_combo.currentText(),
            apex_selection=self.apex_combo.currentText(),
            q_value=float(self.qvalue_spin.value()),
            ppm_alert=float(self.ppm_alert_spin.value()),
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
            self._fail("Select a valid search-ID file (Percolator psms, or the "
                       "mzTab when an SDRF is set).")
            return
        if sdrf_path and not Path(sdrf_path).is_file():
            self._fail("The SDRF path is set but is not a file.")
            return

        os.makedirs(config.out_dir, exist_ok=True)
        self._set_running(True)
        try:
            if sdrf_path:
                await self._run_sdrf(config, sdrf_path, mzml_dir, id_path, workers)
            else:
                await self._run_percolator(config, mzml_dir, id_path, workers)
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

    async def _run_percolator(self, config, mzml_dir, id_path, workers):
        """Demoted single-mzML Percolator path; now also file-parallel."""
        sample = config.sample
        if not sample or not sample[-1].isdigit():
            self._fail(f"Sample must end with a number (got {sample!r}).")
            return
        mzml_files = list_mzml_files(mzml_dir)
        if not mzml_files:
            self._fail(f"No mzML files in {mzml_dir}.")
            return
        loop = asyncio.get_running_loop()
        self._info(f"reading PSMs from {id_path} …")
        psms = await loop.run_in_executor(self.pool, read_psms, id_path, sample, ())
        indices = file_indices(psms)
        if len(mzml_files) != len(indices):
            self._fail(
                f"mzML count ({len(mzml_files)}) != distinct file_idx count "
                f"({len(indices)}) in the id file."
            )
            return
        runs = []
        for idx in indices:
            mzml_file = os.path.join(mzml_dir, mzml_files[idx])
            self._fraction_mzml[idx] = mzml_file
            runs.append((idx, mzml_stem(mzml_files[idx]), mzml_file,
                         fraction_psms(psms, idx)))

        async def run_one(run):
            idx, label, mzml_file, fraction = run
            df, drift = await loop.run_in_executor(
                self.pool, integrate_fraction, config, fraction, mzml_file, label)
            return idx, df, drift

        results = await self._gather_runs(runs, workers, run_one)
        if self._cancelled:
            self._info("cancelled — outputs not written.")
            return
        frames = []
        for idx, label, _mzml_file, _fraction in runs:
            df, drift = results[idx]
            self._write_outputs(config, df, drift, id_path, label)
            self._update_drift(idx, drift, config.ppm_alert)
            frames.append(df)
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

    def _write_outputs(self, config, df, drift, id_path, label) -> None:
        """Write ``<sample>_riana.txt`` (+ drift sidecar) exactly as the CLI does."""
        out_file = Path(config.out_dir) / f"{config.sample}_riana.txt"
        provenance = make_provenance(
            dataclasses.asdict(config),
            id_source=str(id_path),
            extra={"mzml": label},
        )
        write_dataframe_tsv(out_file, df, provenance, include_index=True)
        if drift is not None:
            drift_path = out_file.with_suffix(".drift.json")
            with drift_path.open("w") as fh:
                json.dump(dataclasses.asdict(drift), fh, indent=2)
        self._info(f"wrote {out_file}")

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
