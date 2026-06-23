# -*- coding: utf-8 -*-

"""Shared pyqtgraph scaffolding for the embedded plot views.

The one non-obvious piece is :func:`plot_with_external_legend`: pyqtgraph anchors
a legend *inside* the ViewBox by default, where its sample glyphs (an orange ▲ for
MBR, a purple ◇ for a fold) sit on top of the data and read as real plotted
points. We instead reparent the PlotItem's own legend into a dedicated right-hand
column so it never overlaps the data — while staying the PlotItem's legend, so
``plot(..., name=...)`` still auto-populates it and ``plot.clear()`` still empties
it.
"""

from __future__ import annotations

import pyqtgraph as pg


def plot_with_external_legend(background: str = "w"):
    """Return ``(widget, plot_item, legend_vb)`` with the legend outside the plot.

    ``widget`` is a :class:`~pyqtgraph.GraphicsLayoutWidget` (add it to a layout);
    ``plot_item`` is the central :class:`~pyqtgraph.PlotItem` to draw on. The
    caller must keep a reference to ``legend_vb`` (the legend's host ViewBox) so it
    is not garbage-collected out from under the anchored legend.
    """
    widget = pg.GraphicsLayoutWidget()
    widget.setBackground(background)
    plot = widget.addPlot(row=0, col=0)
    legend = plot.addLegend()
    legend_vb = widget.addViewBox(row=0, col=1, enableMouse=False)
    legend_vb.setMaximumWidth(190)
    legend_vb.setMinimumWidth(120)
    legend.setParentItem(legend_vb)
    legend.anchor((0, 0), (0, 0), offset=(5, 5))
    return widget, plot, legend_vb
