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
    QButtonGroup,
    QCheckBox,
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
    QRadioButton,
    QSpinBox,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from qasync import asyncSlot

from riana.config import FitConfig
from riana.core.fitting import available_coefficient_presets
from riana.core.fitting import peptide_summary
from riana.gui.curve_view import CurveView
from riana.gui.models import DataFrameTableModel
from riana.gui.tasks import run_fit, run_fit_manifest
from riana.io.writers import (
    ESTIMATE_FLOAT_FORMAT,
    make_provenance,
    write_dataframe_tsv,
)

# Result columns to show in the table (the per-peptide ``t`` / ``fs`` lists are
# kept off-screen and used only for the curve plot). The n_mbr / n_metox / n_clean
# point census surfaces the MBR / Met-Ox composition alongside the kinetics.
# ``experiment``/``condition`` are present only on the manifest (multi-group) path
# — a peptide then has one result row per group, so they must be shown for the
# rows to be distinguishable (and for the curve to plot the selected group, not
# just the first). They are dropped silently on the single-file path.
_DISPLAY_COLS = ["concat", "experiment", "condition",
                 "k_deg", "R_squared", "sd", "ci_lo", "ci_hi",
                 "spep", "n_points", "n_mbr", "n_metox", "n_clean", "protein id"]


class ModelTab(QWidget):
    """Form + async fit runner + results/fitted-curve for ``riana fit``."""

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
        self._last_config: FitConfig | None = None
        #: (concat, result-row) of the table selection, so the Fit/Δ view toggle can
        #: re-render the same peptide without a fresh row-change event.
        self._selected: tuple[str, pd.Series] | None = None

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

        # SDRF/manifest path: when set, fit from the manifest's integrate rows
        # (curves grouped by (experiment, condition) with the timepoint from the
        # SDRF identity) instead of the timepoint file list above.
        self.manifest_edit = QLineEdit()
        self.manifest_edit.setPlaceholderText(
            "Optional: riana_manifest.tsv — overrides the file list above")
        man_row = QHBoxLayout()
        man_row.addWidget(self.manifest_edit, stretch=1)
        man_browse = QPushButton("Browse…")
        man_browse.clicked.connect(self._pick_manifest)
        man_row.addWidget(man_browse)
        form.addRow("Manifest", _row(man_row))

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
        self.model_combo.addItems(["simple", "guan", "fornasiero", "calibration"])
        self.model_combo.setToolTip(
            "Kinetic models (simple/guan/fornasiero) fit k_deg vs labeling time. "
            "'calibration' fits a through-origin FS-vs-mixing-proportion recovery "
            "line (R² = recovery quality) — pick it for a mixing-calibration run, "
            "or use a manifest and it is auto-selected.")
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

        self.fs_combo = QComboBox()
        self.fs_combo.addItems(["full envelope", "auto",
                                "iso0-1", "iso0-2", "iso0-3", "iso0-4", "iso0-5"])
        self.fs_combo.setToolTip(
            "Limited-isotopomer scoring (--fs): fit the FS on a leading channel "
            "subset (iso0-N) to dodge co-eluting contaminants in the high "
            "channels — integrate wide, fit narrow. 'auto' widens per-peptide by "
            "the natural-abundance envelope width (RIA-invariant). Default scores "
            "the full integrated envelope.")
        form.addRow("FS scoring", self.fs_combo)

        self.exclude_mbr_check = QCheckBox("Exclude match-between-runs points")
        self.exclude_mbr_check.setChecked(False)
        self.exclude_mbr_check.setToolTip(
            "Drop MBR-transferred points (evidence='mbr') before fitting. MBR "
            "points are used by default; tick to fit only directly-identified "
            "points — the with/without-MBR A/B.")
        form.addRow("MBR", self.exclude_mbr_check)

        self.kp_spin = self._rate_spin(0.5)
        form.addRow("k_p (guan/forn.)", self.kp_spin)
        self.kr_spin = self._rate_spin(0.05)
        form.addRow("k_r (fornasiero)", self.kr_spin)
        self.rp_spin = self._rate_spin(10.0)
        form.addRow("r_p (fornasiero)", self.rp_spin)

        self.workers_spin = QSpinBox()
        self.workers_spin.setRange(1, os.cpu_count() or 1)
        self.workers_spin.setValue(1)
        self.workers_spin.setToolTip(
            "Worker *processes* for the per-peptide fit — the real lever for the "
            "GIL-bound fit. Result is identical regardless of N. When >1 the fit "
            "runs off the shared pool (no nested pools).")
        form.addRow("Workers", self.workers_spin)

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
        self.table.setSortingEnabled(True)
        self.table.selectionModel().currentRowChanged.connect(self._on_row_changed)
        results_split.addWidget(self.table)

        # Curve panel: a Fit / Δspacing / Δmass view toggle above the plot. Δspacing
        # (the drift-robust M0-internal mass-defect, DeuteRater ΔSₓ) leads the two Δ
        # modes; Fit is the default. The same selected peptide re-renders on toggle.
        curve_panel = QWidget()
        curve_layout = QVBoxLayout(curve_panel)
        curve_layout.setContentsMargins(0, 0, 0, 0)
        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("View:"))
        self.view_group = QButtonGroup(self)
        self._fit_radio = QRadioButton("Fit")
        self._fit_radio.setChecked(True)
        self._spacing_radio = QRadioButton("Δ spacing")
        self._mass_radio = QRadioButton("Δ mass")
        self._spacing_radio.setToolTip(
            "M0-internal mass-defect ΔSₓ over time (drift-robust; DeuteRater). The "
            "labelled neutromers' spacing from M0 widens as deuterium incorporates.")
        self._mass_radio.setToolTip(
            "Absolute accurate-mass shift obs − init reference per channel "
            "(drift-sensitive; instrument-drift / sanity QC).")
        for rb in (self._fit_radio, self._spacing_radio, self._mass_radio):
            self.view_group.addButton(rb)
            mode_row.addWidget(rb)
            rb.toggled.connect(self._on_view_changed)
        self._anchor_check = QCheckBox("anchor t0/f0 (display)")
        self._anchor_check.setToolTip(
            "DISPLAY ONLY — subtract the unlabelled (t/f=0) point's Δ from every "
            "point in the Δ spacing / Δ mass views (needs a t/f=0 point), so they "
            "start at 0 and show the pure labelling signal. Does NOT change the "
            "written fs_ds column, which is always t0/f0-anchored at fit time.")
        # A checkbox needs to re-render on BOTH transitions; the radios' _on_view_changed
        # guards on `checked` (to dodge their paired off/on double-fire), which would
        # swallow the un-check here — so use a dedicated always-render slot.
        self._anchor_check.toggled.connect(self._on_display_toggle)
        mode_row.addWidget(self._anchor_check)
        # Fit view: overlay the orthogonal mass-defect second estimate fs_ds on the
        # main fs curve, so the per-timepoint cross-check is visible at a glance.
        self._show_fsds_check = QCheckBox("show fs_ds (Fit)")
        self._show_fsds_check.setToolTip(
            "Overlay the mass-defect second estimate fs_ds (already t0/f0-anchored) "
            "as hollow green ◇ on the Fit graph, to compare it per-timepoint against "
            "the intensity fs. Applies to the Fit view.")
        self._show_fsds_check.toggled.connect(self._on_display_toggle)
        mode_row.addWidget(self._show_fsds_check)
        mode_row.addStretch(1)
        curve_layout.addLayout(mode_row)
        self.curve = CurveView()
        curve_layout.addWidget(self.curve)
        results_split.addWidget(curve_panel)
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

    def _pick_manifest(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select riana_manifest.tsv",
            filter="Manifest (*.tsv);;All files (*)"
        )
        if path:
            self.manifest_edit.setText(path)

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
        # --fs scoring: 'full envelope' (None), 'auto' (per-peptide widening), or
        # 'iso0-N' -> score_channels = N+1 (the same leading-channel count the CLI
        # builds). Mutually exclusive (auto vs explicit) — FitConfig validates.
        fs_text = self.fs_combo.currentText()
        fs_auto = fs_text == "auto"
        score_channels = (int(fs_text.rsplit("-", 1)[1]) + 1
                          if fs_text.startswith("iso0-") else None)

        return FitConfig(
            model=self.model_combo.currentText(),
            label=self.label_combo.currentText(),
            k_p=float(self.kp_spin.value()),
            k_r=float(self.kr_spin.value()),
            r_p=float(self.rp_spin.value()),
            q_value=float(self.qvalue_spin.value()),
            depth=int(self.depth_spin.value()),
            ria_max=float(self.ria_spin.value()),
            score_channels=score_channels,
            fs_auto=fs_auto,
            workers=int(self.workers_spin.value()),
            exclude_mbr=self.exclude_mbr_check.isChecked(),
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
        self._selected = None
        self._cancelled = False

        try:
            config = self.build_config()
        except ValueError as exc:
            self._fail(str(exc))
            return

        manifest = self.manifest_edit.text().strip()
        files = self._selected_files()
        if manifest:
            if not Path(manifest).is_file():
                self._fail("The manifest path is set but is not a file.")
                return
        elif not files:
            self._fail(
                "Add timepoint _riana.txt files, or set a manifest (SDRF path).")
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
        # workers>1 spawns a ProcessPool inside fit_run; run it on a main-process
        # thread (executor=None) so that pool is NOT nested inside a shared-pool
        # worker (which breaks: BrokenProcessPool).
        executor = None if config.workers > 1 else self.pool
        try:
            if manifest:
                self._info(f"fitting from manifest {manifest} …")
                self._future = loop.run_in_executor(
                    executor, run_fit_manifest, config, manifest, coefficients)
                id_source = manifest
            else:
                self._info(f"fitting {len(files)} timepoint file(s) …")
                self._future = loop.run_in_executor(
                    executor, run_fit, config, files, coefficients)
                id_source = ",".join(files)
            result_df = await self._future

            if self._cancelled:
                self._info("cancelled.")
                return

            self._write_output(config, result_df, id_source, coefficients,
                               manifest or None)
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

    def _write_output(self, config, result_df, id_source, coefficients,
                      manifest=None) -> None:
        """Write ``riana_fit_peptides.txt`` (+ the M5 ``riana_fit_fractions.txt``)
        exactly as riana.cli.fit does. On the manifest path, root the outputs at
        the manifest's folder (ignore the Output dir) and record stage='fit' rows.
        """
        if manifest:
            out_dir = Path(manifest).resolve().parent
            if Path(config.out_dir).resolve() != out_dir:
                self._info(
                    f"manifest path: writing next to the manifest ({out_dir}); "
                    f"ignoring Output dir {config.out_dir}")
        else:
            out_dir = Path(config.out_dir)
        provenance = make_provenance(
            dataclasses.asdict(config),
            id_source=str(id_source),
            extra={"model": config.model, "label": config.label,
                   "coefficients": str(coefficients)},
        )
        out_path = out_dir / "riana_fit_peptides.txt"
        write_dataframe_tsv(out_path, peptide_summary(result_df), provenance,
                            include_index=True,
                            float_format=ESTIMATE_FLOAT_FORMAT)
        self._info(f"wrote {out_path}")
        written = [out_path]

        fractions = result_df.attrs.get("fractions_long")
        if fractions is not None and not fractions.empty:
            frac_path = out_dir / "riana_fit_fractions.txt"
            write_dataframe_tsv(frac_path, fractions, provenance,
                                include_index=False,
                                float_format=ESTIMATE_FLOAT_FORMAT)
            self._info(f"wrote {frac_path} ({len(fractions)} peptide-timepoints)")
            written.append(frac_path)

        if manifest:
            from riana.core.pipeline import record_stage_rows
            record_stage_rows(manifest, "fit", written, result_df, provenance)
            self._info(f"recorded {len(written)} fit rows in {manifest}")

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
            self._selected = None
            return
        sel = self.model.dataframe.iloc[current.row()]
        concat = sel["concat"]
        if concat not in self._result_df.index:
            self._selected = None
            return
        row = self._result_df.loc[concat]
        if isinstance(row, pd.DataFrame):  # manifest path: same peptide, >1 group
            # One result row per (experiment, condition); plot the group the
            # selected table row belongs to, not just the first — else the curve
            # silently shows a different condition than the row's n_* counts.
            row = _match_group(row, sel)
        self._selected = (str(concat), row)
        self._render_curve()

    def _on_view_changed(self, checked: bool) -> None:
        # QButtonGroup toggled fires for both the de-selected and the newly selected
        # button; act only on the on-event (the other carries the same render).
        if checked:
            self._render_curve()

    def _on_display_toggle(self, _checked: bool) -> None:
        # The display checkboxes (anchor, show fs_ds) must re-render on BOTH check and
        # un-check (unlike the radios), so they don't share _on_view_changed's guard.
        self._render_curve()

    def _render_curve(self) -> None:
        """Draw the selected peptide in the view the toggle selects (Fit / Δ)."""
        if self._selected is None:
            return
        concat, row = self._selected
        if self._spacing_radio.isChecked() or self._mass_radio.isChecked():
            mode = "spacing" if self._spacing_radio.isChecked() else "mass"
            dm = row.get("dmass")
            ds = row.get("dspacing")
            x_label = ("Mixing proportion"
                       if self._last_config is not None
                       and self._last_config.model == "calibration"
                       else "Time")
            self.curve.plot_dmass(
                concat, list(row["t"]),
                list(dm) if dm is not None else [],
                list(ds) if ds is not None else [],
                mode=mode,
                protein=row.get("protein id"),
                condition=row.get("condition"),
                x_label=x_label,
                anchor=self._anchor_check.isChecked(),
            )
            return
        cfg = self._last_config
        if cfg is None:
            return
        kinetic = dict(k_p=cfg.k_p, k_r=cfg.k_r, r_p=cfg.r_p)
        ev = row.get("evidence")
        mx = row.get("metox")
        fsds = row.get("fs_ds") if self._show_fsds_check.isChecked() else None
        self.curve.plot_fit(
            concat, list(row["t"]), list(row["fs"]),
            float(row["k_deg"]), cfg.model, kinetic,
            ci_lo=_safe_float(row.get("ci_lo")),
            ci_hi=_safe_float(row.get("ci_hi")),
            evidence=list(ev) if ev is not None else None,
            metox=list(mx) if mx is not None else None,
            protein=row.get("protein id"),
            condition=row.get("condition"),
            fs_ds=list(fsds) if fsds is not None else None,
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


def _match_group(rows: pd.DataFrame, sel: pd.Series) -> pd.Series:
    """Pick the result row whose group matches the selected table row.

    On the manifest path a peptide has one result row per ``(experiment,
    condition)`` group; the curve must plot the group the user selected. Match on
    whichever group keys both frames carry; fall back to the first row if the
    selection is ambiguous (e.g. group columns absent), preserving prior behaviour.
    """
    mask = pd.Series(True, index=rows.index)
    for key in ("experiment", "condition"):
        if key in rows.columns and key in sel.index:
            mask &= rows[key] == sel[key]
    matched = rows[mask]
    return matched.iloc[0] if len(matched) else rows.iloc[0]


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
