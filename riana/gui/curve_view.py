# -*- coding: utf-8 -*-

"""Embedded pyqtgraph view of a peptide's kinetic fit (Model tab).

Plots the per-timepoint fraction-new points the fit computed and overlays the
fitted kinetic-model curve. The model functions in :mod:`riana.core.models` are
pure (vectorised over ``t``), so the curve is drawn on the GUI thread directly
from the result row's ``k_deg`` — no worker round-trip needed.
"""

from __future__ import annotations

import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtWidgets import QVBoxLayout, QWidget

from riana.core import models

# Kinetic-model name → function (same mapping fit_run uses).
_MODEL_FNS = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
}


class CurveView(QWidget):
    """A pyqtgraph plot of one peptide's (t, fraction-new) data + fitted curve."""

    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plot = pg.PlotWidget()
        self.plot.setBackground("w")
        self.plot.addLegend()
        self.plot.setLabel("bottom", "Time")
        self.plot.setLabel("left", "Fraction new")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        layout.addWidget(self.plot)
        self.show_placeholder("Run a fit, then select a peptide.")

    def show_placeholder(self, message: str) -> None:
        self.plot.clear()
        self.plot.setTitle(message)

    def plot_fit(
        self,
        concat: str,
        t: list[float],
        fs: list[float],
        k_deg: float,
        model_name: str,
        kinetic_kwargs: dict,
    ) -> None:
        """Scatter the (t, fs) data and overlay the fitted model curve."""
        self.plot.clear()
        if not t:
            self.show_placeholder(f"{concat}: no fitted data points.")
            return

        # Observed fraction-new points.
        self.plot.plot(
            list(t), list(fs), pen=None,
            symbol="o", symbolSize=8, symbolBrush="#1f77b4", name="observed",
        )

        # Fitted model curve on a dense grid, when k_deg converged.
        if k_deg is not None and math.isfinite(k_deg):
            model_fn = _MODEL_FNS.get(model_name, models.one_exponent)
            grid = np.linspace(0.0, max(t), 200)
            curve = model_fn(grid, k_deg=k_deg, a_0=0.0, a_max=1.0, **kinetic_kwargs)
            self.plot.plot(
                grid, np.asarray(curve, dtype=float),
                pen=pg.mkPen("#d62728", width=2), name=f"fit (k_deg={k_deg:.3g})",
            )
        self.plot.setYRange(0.0, 1.0)
        self.plot.setTitle(concat)
