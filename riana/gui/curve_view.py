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
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QVBoxLayout, QWidget

from riana.core import models

# Kinetic-model name → function (same mapping fit_run uses).
_MODEL_FNS = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
}

# Per-condition colours for the linear (φ-space) overlay; cycled if exceeded.
_CONDITION_COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd",
                      "#ff7f0e", "#17becf"]


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
        self.plot.setLabel("left", "Fraction new")
        self.plot.setTitle(concat)

    def plot_linear(
        self,
        protein: str,
        per_condition: dict,
        *,
        phi_limit: float = -4.0,
        delta_k: float | None = None,
        delta_k_p_adj: float | None = None,
    ) -> None:
        """φ-space view for the ``linear simple`` model — overlay each condition.

        ``per_condition`` maps ``condition → (t_list, theta_list, k_deg)``. For
        each condition we plot the clearance φ = log(1 − θ) points (saturated
        points past ``phi_limit`` drawn hollow, as they are dropped from the fit)
        and the fitted through-origin line φ = −k·t. The slope difference between
        the lines *is* the cross-sample Δk, which the title reports.
        """
        from riana.core.linear_model import to_phi

        self.plot.clear()
        if not per_condition:
            self.show_placeholder(f"{protein}: no points to show.")
            return

        t_max = 1.0
        for i, cond in enumerate(sorted(per_condition)):
            t_list, theta_list, k = per_condition[cond]
            if not t_list:
                continue
            colour = _CONDITION_COLOURS[i % len(_CONDITION_COLOURS)]
            t = np.asarray(t_list, dtype=float)
            phi = to_phi(np.asarray(theta_list, dtype=float))
            t_max = max(t_max, float(t.max()))
            kept = phi > phi_limit
            if kept.any():  # points used in the fit — filled
                self.plot.plot(
                    t[kept].tolist(), phi[kept].tolist(), pen=None, symbol="o",
                    symbolSize=8, symbolBrush=colour, symbolPen=colour, name=cond)
            if (~kept).any():  # truncated (saturated) points — hollow
                self.plot.plot(
                    t[~kept].tolist(), phi[~kept].tolist(), pen=None, symbol="o",
                    symbolSize=8, symbolBrush=None, symbolPen=colour)
            if k is not None and math.isfinite(k):
                grid = np.linspace(0.0, t_max, 100)
                self.plot.plot(
                    grid, (-float(k) * grid).tolist(),
                    pen=pg.mkPen(colour, width=2),
                    name=f"{cond} k={k:.3g}")

        # The plateau-truncation threshold, so the dropped points read as "below
        # the line".
        line = pg.InfiniteLine(pos=phi_limit, angle=0,
                               pen=pg.mkPen("#888888", style=Qt.PenStyle.DashLine))
        self.plot.addItem(line)

        self.plot.setLabel("left", "Clearance  φ = log(1 − θ)")
        title = protein
        if delta_k is not None and math.isfinite(delta_k):
            title += f"   Δk={delta_k:+.3g}/day"
            if delta_k_p_adj is not None and math.isfinite(delta_k_p_adj):
                title += f"  (p_adj={delta_k_p_adj:.2g})"
        self.plot.setTitle(title)
        self.plot.enableAutoRange(axis="y")
