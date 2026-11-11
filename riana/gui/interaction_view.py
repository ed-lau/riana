# -*- coding: utf-8 -*-

"""Small-multiple φ-space panels for one protein's experiment × condition interaction.

One panel per experiment (e.g. heart chamber), each overlaying the reference and test
conditions' collapsed rollup points as clearance φ = log(1 − θ) and the joint
random-animal fit's through-origin lines φ = −k·t (``riana_rollup_interaction.txt``).
Shared axes, so a chamber-specific effect reads as one panel's lines diverging while
the others stay together.
"""

from __future__ import annotations

import html
import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from riana.gui.export import save_plot

#: Reference / test colours (the linear view's first two condition colours).
REFERENCE_COLOUR, TEST_COLOUR = "#1f77b4", "#d62728"


class InteractionView(QWidget):
    """A row of φ-space panels, one per experiment, for the selected protein."""

    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.widget = pg.GraphicsLayoutWidget()
        self.widget.setBackground("w")
        layout.addWidget(self.widget)

        controls = QHBoxLayout()
        controls.addStretch(1)
        self.save_button = QPushButton("Save graph…")
        self.save_button.clicked.connect(self._on_save)
        controls.addWidget(self.save_button)
        layout.addLayout(controls)

        self._title = "interaction"
        #: the PlotItems of the current protein, one per experiment (in panel order).
        self.panels: list[pg.PlotItem] = []
        self.show_placeholder("Select a protein row to view it per experiment.")

    def show_placeholder(self, message: str) -> None:
        self.widget.clear()
        self.panels = []
        self.widget.addLabel(html.escape(message), row=0, col=0, color="#666666")
        self.save_button.setEnabled(False)

    def _on_save(self) -> None:
        save_plot(self.widget, self, default_name=self._title)

    def plot_protein(
        self,
        protein: str,
        panels: list[dict],
        *,
        reference: str,
        test: str,
        phi_limit: float = -4.0,
        summary: str = "",
    ) -> None:
        """Draw one panel per experiment.

        ``panels``: ``[{"experiment", "k_ratio", "k_reference", "k_test",
        "points": {condition: (t_list, theta_list)}}]``. Points at or below
        ``phi_limit`` (dropped by the fit) are drawn hollow.
        """
        from riana.core.linear_model import to_phi

        self.widget.clear()
        self.panels = []
        if not panels:
            self.show_placeholder(f"{protein}: no experiments to show.")
            return
        self._title = f"{protein}_interaction"
        self.save_button.setEnabled(True)
        n = len(panels)
        head = f"<b>{html.escape(str(protein))}</b>"
        if summary:
            head += f"&nbsp;&nbsp; {html.escape(summary)}"
        self.widget.addLabel(head, row=0, col=0, colspan=n)

        t_max = max([1.0] + [max(pts[0]) for panel in panels
                             for pts in panel["points"].values() if len(pts[0])])
        grid = np.linspace(0.0, t_max, 60)
        y_seen = [0.0]           # every drawn φ, for one shared y range
        for i, panel in enumerate(panels):
            plot = self.widget.addPlot(row=1, col=i)
            plot.showGrid(x=True, y=True, alpha=0.2)
            plot.setLabel("bottom", "Time")
            # No SI-prefix rescaling: on a labelled axis pyqtgraph would show φ −0.5
            # as "−500 (x0.001)".
            for axis in ("left", "bottom"):
                plot.getAxis(axis).enableAutoSIPrefix(False)
            if i == 0:
                plot.setLabel("left", "φ = log(1 − θ)")
            title = str(panel["experiment"])
            ratio = panel.get("k_ratio")
            if ratio is not None and math.isfinite(ratio):
                title += f"   k ratio {ratio:.2f}"
            plot.setTitle(html.escape(title))
            for cond, colour, k in ((reference, REFERENCE_COLOUR, panel.get("k_reference")),
                                    (test, TEST_COLOUR, panel.get("k_test"))):
                t_list, theta_list = panel["points"].get(cond, ([], []))
                if len(t_list):
                    t = np.asarray(t_list, dtype=float)
                    phi = to_phi(np.asarray(theta_list, dtype=float))
                    y_seen.extend(phi[np.isfinite(phi)].tolist())
                    kept = phi > phi_limit
                    if kept.any():
                        plot.plot(t[kept], phi[kept], pen=None, symbol="o", symbolSize=7,
                                  symbolBrush=colour, symbolPen=colour)
                    if (~kept).any():
                        plot.plot(t[~kept], phi[~kept], pen=None, symbol="o",
                                  symbolSize=7, symbolBrush=None, symbolPen=colour)
                if k is not None and math.isfinite(k):
                    plot.plot(grid, -float(k) * grid, pen=pg.mkPen(colour, width=2))
                    y_seen.append(-float(k) * t_max)
            # The plateau-truncation threshold, kept out of the autoscale so the
            # panels fit the data (it shows once points approach it).
            plot.addItem(pg.InfiniteLine(
                pos=phi_limit, angle=0,
                pen=pg.mkPen("#888888", style=Qt.PenStyle.DashLine)),
                ignoreBounds=True)
            self.panels.append(plot)
        # One explicit, identical range for every panel so they compare directly.
        # (Not pyqtgraph axis links: those sync by screen position, which is not
        # settled while the layout is built, and left the panels clipped.)
        for plot in self.panels:
            plot.setXRange(0.0, t_max, padding=0.03)
            plot.setYRange(min(y_seen), max(y_seen), padding=0.05)
        self.widget.addLabel(
            f"<span style='color:{REFERENCE_COLOUR}'>●</span> {html.escape(reference)}"
            f"&nbsp;&nbsp;<span style='color:{TEST_COLOUR}'>●</span> {html.escape(test)}"
            "&nbsp;&nbsp; lines: joint random-animal fit, φ = −k·t",
            row=2, col=0, colspan=n)
