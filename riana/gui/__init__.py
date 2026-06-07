# -*- coding: utf-8 -*-

"""Riana PySide6 GUI (M4 Phase 2).

Optional — requires the ``[gui]`` extra (``pip install 'riana[gui]'``: PySide6 +
qasync + pyqtgraph). Reached via the ``riana gui`` subcommand. Importing this
package does **not** import Qt; the entry point :func:`riana.gui.app.run_gui`
imports PySide6 lazily so the core CLI stays lightweight.

The worker functions in :mod:`riana.gui.tasks` are deliberately Qt-free so they
can run in a ``ProcessPoolExecutor`` and be unit-tested without a display.
"""

from __future__ import annotations

__all__ = ["run_gui"]


def run_gui(threads: int = 1) -> int:
    """Launch the GUI. Lazy re-export of :func:`riana.gui.app.run_gui`."""
    from riana.gui.app import run_gui as _run_gui

    return _run_gui(threads=threads)
