# -*- coding: utf-8 -*-

"""Embedded pyqtgraph chromatogram view for the Integrate tab.

Plots the per-isotopomer XICs of a selected peptide and shades the integration
window, using the :class:`~riana.core.integration.PeptideTrace` the worker
returns. pyqtgraph (not matplotlib) for fast, interactive pan/zoom; matplotlib
is reserved for static export elsewhere.
"""

from __future__ import annotations

import pyqtgraph as pg
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QVBoxLayout, QWidget

from riana.core.integration import PeptideTrace

# A small qualitative palette cycled across isotopomer traces.
_PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd",
    "#ff7f0e", "#8c564b", "#17becf", "#bcbd22",
]


class ChromatogramView(QWidget):
    """A pyqtgraph plot of one peptide's isotopomer XICs + integration window."""

    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plot = pg.PlotWidget()
        self.plot.setBackground("w")
        self.plot.addLegend()
        self.plot.setLabel("bottom", "Retention time", units="min")
        self.plot.setLabel("left", "Intensity")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        layout.addWidget(self.plot)
        self.show_placeholder("Run an integration, then select a peptide.")

    def show_placeholder(self, message: str) -> None:
        """Clear the plot and show a centred hint."""
        self.plot.clear()
        self.plot.setTitle(message)

    def plot_trace(self, trace: PeptideTrace) -> None:
        """Draw the isotopomer XICs and shade the integrated RT window."""
        self.plot.clear()
        for i, iso in enumerate(sorted(trace.chromatograms)):
            chrom = trace.chromatograms[iso]
            pen = pg.mkPen(color=_PALETTE[i % len(_PALETTE)], width=2)
            self.plot.plot(
                list(chrom.rt), list(chrom.intensity), pen=pen, name=f"iso{iso}"
            )
        if trace.window is not None:
            lo, hi = trace.window
            shade = QColor(120, 120, 120, 50)
            region = pg.LinearRegionItem(values=(lo, hi), movable=False, brush=shade)
            region.setZValue(-10)
            self.plot.addItem(region)
        self.plot.setTitle(trace.concat)
