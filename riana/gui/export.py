# -*- coding: utf-8 -*-

"""Save an embedded pyqtgraph plot to a PNG (replaces the old ``--plotcurves``).

The old CLI dumped fitted-curve PNGs; the GUI instead lets the user export
whichever plot they are looking at. pyqtgraph's ``ImageExporter`` renders the
live ``PlotItem`` directly — no worker round-trip, since the figure is already
on screen. Raster only: pyqtgraph's ``SVGExporter`` throws on plots carrying
scatter symbols (the observed-point / fold series these views always draw), so
offering an SVG option would crash on the common case.
"""

from __future__ import annotations

from PySide6.QtWidgets import QFileDialog, QWidget


def save_plot(plot_widget, parent: QWidget, *, default_name: str = "graph") -> str | None:
    """Prompt for a PNG path and export *plot_widget*'s ``PlotItem`` to it.

    Returns the written path, or ``None`` if the user cancelled.
    """
    path, _ = QFileDialog.getSaveFileName(
        parent, "Save graph", f"{default_name}.png",
        filter="PNG image (*.png)",
    )
    if not path:
        return None
    if not path.lower().endswith(".png"):
        path += ".png"
    return export_plot(plot_widget, path)


def export_plot(plot_item, path: str) -> str:
    """Export *plot_item*'s whole scene to *path* as a PNG.

    Exports the **scene**, not just the PlotItem, so the legend — which lives in
    its own column outside the data area (see :mod:`riana.gui.plotting`) — is
    captured too. Split out from the dialog so it is testable headless;
    ``ImageExporter`` is imported lazily (it pulls in Qt machinery the core GUI
    does not need until an export is actually requested).
    """
    from pyqtgraph.exporters import ImageExporter

    ImageExporter(plot_item.scene()).export(path)
    return path
