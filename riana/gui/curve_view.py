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
from PySide6.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from riana.core import models
from riana.gui.export import save_plot
from riana.gui.plotting import plot_with_external_legend

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
        # Legend in its own column (outside the data area) so its MBR/fold sample
        # glyphs can't be mistaken for plotted points. ``_legend_vb`` must be
        # retained or the anchored legend loses its host.
        self._plot_widget, self.plot, self._legend_vb = plot_with_external_legend()
        self.plot.setLabel("bottom", "Time")
        self.plot.setLabel("left", "Fraction new")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        layout.addWidget(self._plot_widget)

        controls = QHBoxLayout()
        controls.addStretch(1)
        self.save_button = QPushButton("Save graph…")
        self.save_button.clicked.connect(self._on_save)
        controls.addWidget(self.save_button)
        layout.addLayout(controls)

        self._title = "fit"
        self.show_placeholder("Run a fit, then select a peptide.")

    def show_placeholder(self, message: str) -> None:
        self.plot.clear()
        self.plot.setTitle(message)
        self.save_button.setEnabled(False)

    def _on_save(self) -> None:
        save_plot(self.plot, self, default_name=self._title)

    def plot_fit(
        self,
        concat: str,
        t: list[float],
        fs: list[float],
        k_deg: float,
        model_name: str,
        kinetic_kwargs: dict,
        ci_lo: float | None = None,
        ci_hi: float | None = None,
        evidence: list[str] | None = None,
        metox: list[bool] | None = None,
        protein: str | None = None,
        condition: str | None = None,
    ) -> None:
        """Scatter the (t, fs) data and overlay the fitted model curve.

        ``protein`` and ``condition`` (when given) join ``concat`` in the plot
        title so an exported figure identifies which protein and — for a
        multi-condition manifest fit — which condition group it shows.

        When ``ci_lo`` / ``ci_hi`` (the k_deg CI) are given, a shaded ribbon
        between the model curves at those k bounds shows the fit uncertainty.

        Observed points are split by provenance so the user can see how the curve
        was assembled, each as its own series (aligned to ``t``):

        * **direct ID** — blue ○ (the default).
        * **MBR transfer** (``evidence[i] == "mbr"``) — orange △, the transferred
          quant distinct from direct IDs.
        * **chemical-mod fold** (``metox[i]``) — purple ◇, a point consolidated
          into this curve by ``core.fitting._fit_key`` (Met-Ox today; **TMT** once
          it joins ``constants.CHEMICAL_MODS`` — the flag is generic over chemical
          mods, so TMT folds light up here with no further GUI change).

        MBR keeps priority for the symbol (its established orange △); among the
        direct IDs a chemical-fold point takes the purple ◇.
        """
        self.plot.clear()
        if not t:
            self.show_placeholder(f"{concat}: no fitted data points.")
            return
        self._title = "_".join(
            str(p) for p in (concat, condition) if p
        ).replace("|", "_").replace("/", "_")  # export filename
        self.save_button.setEnabled(True)

        ev = list(evidence) if evidence is not None and len(evidence) == len(t) else None
        mx = list(metox) if metox is not None and len(metox) == len(t) else None

        def _point_class(i: int) -> str:
            if ev is not None and ev[i] == "mbr":
                return "mbr"
            if mx is not None and mx[i]:
                return "folded"
            return "direct"

        # (symbol, size, colour, legend label); folded covers Met-Ox + future TMT.
        styles = {
            "direct": ("o", 8, "#1f77b4", "observed"),
            "mbr": ("t", 10, "#ff7f0e", "MBR"),
            "folded": ("d", 10, "#9467bd", "folded (Met-Ox/TMT)"),
        }
        groups: dict[str, list[tuple[float, float]]] = {k: [] for k in styles}
        for i, (ti, fi) in enumerate(zip(t, fs)):
            groups[_point_class(i)].append((ti, fi))
        for cls, pts in groups.items():
            if not pts:
                continue
            symbol, size, colour, label = styles[cls]
            name = label if cls == "direct" else f"{label} ({len(pts)})"
            self.plot.plot(
                [p[0] for p in pts], [p[1] for p in pts], pen=None,
                symbol=symbol, symbolSize=size, symbolBrush=colour,
                symbolPen=colour, name=name,
            )

        # Fitted model curve on a dense grid, when k_deg converged.
        if k_deg is not None and math.isfinite(k_deg):
            model_fn = _MODEL_FNS.get(model_name, models.one_exponent)
            grid = np.linspace(0.0, max(t), 200)

            def _curve(k):
                return np.asarray(
                    model_fn(grid, k_deg=k, a_0=0.0, a_max=1.0, **kinetic_kwargs),
                    dtype=float)

            # CI ribbon: higher k ⇒ higher θ, so ci_hi is the upper edge.
            if (ci_lo is not None and ci_hi is not None
                    and math.isfinite(ci_lo) and math.isfinite(ci_hi)):
                self._ci_band(grid, _curve(ci_lo), _curve(ci_hi), (214, 39, 40))
            self.plot.plot(
                grid, _curve(k_deg),
                pen=pg.mkPen("#d62728", width=2), name=f"fit (k_deg={k_deg:.3g})",
            )
        self.plot.setYRange(0.0, 1.0)
        self.plot.setLabel("left", "Fraction new")
        header = "  •  ".join(str(p) for p in (concat, protein, condition) if p)
        self.plot.setTitle(header)

    def _ci_band(self, x, lower, upper, rgb: tuple) -> None:
        """Shade a confidence ribbon between *lower* and *upper* over *x*.

        Uses ``PlotCurveItem`` boundaries (faint edges) + a translucent
        ``FillBetweenItem`` — deliberately not ``PlotDataItem``, so the band does
        not count as a data series.
        """
        edge = pg.mkPen(rgb + (110,), width=1)
        lo = pg.PlotCurveItem(np.asarray(x), np.asarray(lower, dtype=float), pen=edge)
        hi = pg.PlotCurveItem(np.asarray(x), np.asarray(upper, dtype=float), pen=edge)
        self.plot.addItem(lo)
        self.plot.addItem(hi)
        self.plot.addItem(pg.FillBetweenItem(lo, hi, brush=pg.mkBrush(rgb + (55,))))

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

        ``per_condition`` maps ``condition → (t_list, theta_list, k_deg[, ci_lo,
        ci_hi])``. For each condition we plot the clearance φ = log(1 − θ) points
        (saturated points past ``phi_limit`` drawn hollow, as they are dropped
        from the fit) and the fitted through-origin line φ = −k·t — with a shaded
        ribbon between the lines at the k CI bounds when ``ci_lo`` / ``ci_hi`` are
        given. The slope difference between the lines *is* the cross-sample Δk,
        which the title reports.
        """
        from riana.core.linear_model import to_phi

        self.plot.clear()
        if not per_condition:
            self.show_placeholder(f"{protein}: no points to show.")
            return
        self._title = protein
        self.save_button.setEnabled(True)

        t_max = 1.0
        for i, cond in enumerate(sorted(per_condition)):
            vals = per_condition[cond]
            t_list, theta_list, k = vals[0], vals[1], vals[2]
            ci_lo, ci_hi = (vals[3], vals[4]) if len(vals) >= 5 else (None, None)
            if not t_list:
                continue
            colour = _CONDITION_COLOURS[i % len(_CONDITION_COLOURS)]
            rgb = pg.mkColor(colour).getRgb()[:3]
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
                # CI ribbon: φ = −k·t, so higher k ⇒ lower (more negative) φ;
                # ci_hi is the lower edge, ci_lo the upper.
                if (ci_lo is not None and ci_hi is not None
                        and math.isfinite(ci_lo) and math.isfinite(ci_hi)):
                    self._ci_band(grid, -float(ci_hi) * grid,
                                  -float(ci_lo) * grid, rgb)
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
