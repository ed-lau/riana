# -*- coding: utf-8 -*-

"""GUI entry point: wire qasync into the Qt event loop and show the window.

``run_gui`` is what the ``riana gui`` subcommand calls. It installs a
:class:`qasync.QEventLoop` as the asyncio event loop so the Integrate tab's
``@asyncSlot`` handlers can ``await`` CPU work submitted to a shared
``ProcessPoolExecutor`` without freezing the UI.
"""

from __future__ import annotations

import asyncio
import sys
from concurrent.futures import ProcessPoolExecutor

import qasync
from PySide6.QtWidgets import QApplication

from riana.gui.main_window import MainWindow


def run_gui(threads: int = 1) -> int:
    """Create the QApplication, run the qasync loop, return the exit code.

    Args:
        threads: seeds the Integrate form's worker-thread field. CPU work runs
            in worker *processes*; this is the per-fraction in-process thread
            count handed to :class:`~riana.config.IntegrationConfig`.
    """
    app = QApplication.instance() or QApplication(sys.argv)

    loop = qasync.QEventLoop(app)
    asyncio.set_event_loop(loop)

    # One pool shared by every tab; spawn-safe because the workers
    # (riana.gui.tasks) are top-level and take only picklable args.
    pool = ProcessPoolExecutor()

    close_event = asyncio.Event()
    app.aboutToQuit.connect(close_event.set)

    window = MainWindow(pool=pool, default_threads=threads)
    window.show()

    try:
        with loop:
            loop.run_until_complete(close_event.wait())
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    return 0
