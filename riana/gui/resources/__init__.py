# -*- coding: utf-8 -*-

"""Packaged GUI assets (the app icon) + a loader robust to zip installs."""

from __future__ import annotations


def app_icon():
    """Return the Riana :class:`QIcon`, or an empty icon if the asset is missing.

    Loads the packaged PNG via :mod:`importlib.resources` and through a
    ``QPixmap`` (not a filesystem path) so it works from a zipped/wheel install.
    Never raises — a missing icon must not stop the GUI from launching.
    """
    from PySide6.QtGui import QIcon, QPixmap

    try:
        from importlib.resources import files

        data = files("riana.gui.resources").joinpath("riana.png").read_bytes()
        pixmap = QPixmap()
        pixmap.loadFromData(data, "PNG")
        return QIcon(pixmap)
    except (FileNotFoundError, ModuleNotFoundError, OSError):
        return QIcon()
