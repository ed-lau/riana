# -*- coding: utf-8 -*-

"""Top-level window: a tabbed shell holding the Integrate, Model, and Protein tabs.

Each tab builds the *same* frozen config / calls the *same* core entry point the
CLI does, off the Qt event loop via the shared ``ProcessPoolExecutor``, so the
surfaces cannot diverge: **Integrate** (`integrate_run`), **Model** (`fit_run`),
**Protein** (`rollup_proteins` — Track C). Calibration is folded into the
Integrate tab's results (the per-fraction drift summary), so there is no
separate Calibration tab.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

from PySide6.QtWidgets import QMainWindow, QTabWidget

from riana import __version__
from riana.gui.integrate_tab import IntegrateTab
from riana.gui.model_tab import ModelTab
from riana.gui.protein_tab import ProteinTab


class MainWindow(QMainWindow):
    """The application's main window.

    Args:
        pool: the shared process pool every tab submits CPU work to.
        default_threads: seeds the Integrate form's thread field.
    """

    def __init__(
        self,
        pool: ProcessPoolExecutor,
        default_threads: int = 1,
    ) -> None:
        super().__init__()
        self.pool = pool
        self.setWindowTitle(f"Riana {__version__}")
        self.resize(1100, 720)

        tabs = QTabWidget()
        self.integrate_tab = IntegrateTab(
            pool=pool, default_threads=default_threads, status_cb=self._status
        )
        self.model_tab = ModelTab(
            pool=pool, default_threads=default_threads, status_cb=self._status
        )
        self.protein_tab = ProteinTab(pool=pool, status_cb=self._status)
        tabs.addTab(self.integrate_tab, "Integrate")
        tabs.addTab(self.model_tab, "Model")
        tabs.addTab(self.protein_tab, "Protein")
        self.tabs = tabs
        self.setCentralWidget(tabs)

        self.statusBar().showMessage("Ready.")

    def _status(self, message: str) -> None:
        """Status-bar callback handed to tabs."""
        self.statusBar().showMessage(message)
