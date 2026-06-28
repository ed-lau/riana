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

from PySide6.QtCore import QEvent, QObject
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QTabWidget,
    QWidget,
)

from riana import __version__
from riana.gui.integrate_tab import IntegrateTab
from riana.gui.model_tab import ModelTab
from riana.gui.protein_tab import ProteinTab
from riana.gui.resources import app_icon


class _HintEventFilter(QObject):
    """Mirror the hovered/focused widget's tooltip into a fixed hint label.

    The help text is then visible the moment the mouse (or keyboard focus) lands on
    a control, instead of only after a hover-hold reveals the native tooltip. Reads
    the existing ``toolTip()`` (walking up to the nearest ancestor that has one), so
    there is nothing to keep in sync — every widget that already has a tooltip is
    covered for free. Installed application-wide; it only acts on Enter / FocusIn.
    """

    _TRIGGERS = (QEvent.Type.Enter, QEvent.Type.FocusIn)

    def __init__(self, hint_label: QLabel) -> None:
        super().__init__()
        self._hint = hint_label

    def eventFilter(self, obj, event) -> bool:
        if event.type() in self._TRIGGERS and isinstance(obj, QWidget):
            tip, w = "", obj
            while w is not None and not tip:
                tip = w.toolTip()
                w = w.parentWidget()
            self._hint.setText(" ".join(tip.split()) if tip else "")
        return False  # never consume — purely observational


class MainWindow(QMainWindow):
    """The application's main window.

    Args:
        pool: the shared process pool every tab submits CPU work to.
    """

    def __init__(
        self,
        pool: ProcessPoolExecutor,
    ) -> None:
        super().__init__()
        self.pool = pool
        self.setWindowTitle(f"Riana {__version__}")
        self.setWindowIcon(app_icon())
        self.resize(1100, 720)

        tabs = QTabWidget()
        self.integrate_tab = IntegrateTab(pool=pool, status_cb=self._status)
        self.model_tab = ModelTab(pool=pool, status_cb=self._status)
        self.protein_tab = ProteinTab(pool=pool, status_cb=self._status)
        tabs.addTab(self.integrate_tab, "Integrate")
        tabs.addTab(self.model_tab, "Model")
        tabs.addTab(self.protein_tab, "Protein")
        self.tabs = tabs
        self.setCentralWidget(tabs)

        # A fixed hint area (right of the transient status message): it shows the
        # help of whatever control the mouse/focus is on, so users see hints without
        # hover-holding. Fed by an app-wide event filter that reads each widget's
        # existing tooltip.
        self.hint_label = QLabel("")
        self.hint_label.setStyleSheet("color: palette(mid);")
        self.statusBar().addPermanentWidget(self.hint_label, 1)
        self._hint_filter = _HintEventFilter(self.hint_label)
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self._hint_filter)

        self.statusBar().showMessage("Ready.")

    def _status(self, message: str) -> None:
        """Status-bar callback handed to tabs."""
        self.statusBar().showMessage(message)
