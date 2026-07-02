# -*- coding: utf-8 -*-

"""Drive a determinate ``QProgressBar`` from a worker's ``(done, total)`` stream.

The fit / rollup run in a worker (a shared-pool *process*, or a main-process
thread when ``-W`` spawns its own pool), so they can't touch Qt. The core loops
already take a ``progress_callback(done, total)``; the GUI passes one that forwards
throttled updates onto a :class:`multiprocessing.Manager` queue (picklable across
the process boundary — a raw ``mp.Queue`` can't be sent as a task arg under the
spawn start method macOS uses). :class:`ProgressPump` is the GUI-side reader: a
``QTimer`` on the **main thread** drains the queue to its latest value and updates
the bar, so nothing but the timer callback touches Qt. The bar stays busy
(indeterminate) until the first ``total`` arrives, then goes determinate.
"""

from __future__ import annotations

import queue as _queue

from PySide6.QtCore import QTimer


class ProgressPump:
    """Poll a ``(done, total)`` queue on the GUI thread and update a progress bar."""

    def __init__(self, bar, progress_queue, *, interval_ms: int = 100, parent=None) -> None:
        self._bar = bar
        self._queue = progress_queue
        self._timer = QTimer(parent)
        self._timer.setInterval(interval_ms)
        self._timer.timeout.connect(self._drain)

    def start(self) -> None:
        self._timer.start()

    def stop(self) -> None:
        """Stop polling and apply any final update still on the queue.

        Call this **before** the queue's Manager is shut down, so the last drain
        still sees a live proxy.
        """
        self._timer.stop()
        self._drain()

    def _drain(self) -> None:
        """Read every pending update, keep the latest, and apply it to the bar.

        Only the most recent ``(done, total)`` matters for a bar, so intermediate
        updates are collapsed. Any queue error (empty, or a torn-down Manager
        proxy) just ends the drain — progress display must never raise.
        """
        latest = None
        while True:
            try:
                latest = self._queue.get_nowait()
            except _queue.Empty:
                break
            except Exception:  # a shut-down Manager proxy, etc. — stop draining
                break
        if latest is not None:
            self.apply(self._bar, latest)

    @staticmethod
    def apply(bar, update) -> None:
        """Set *bar* determinate to ``(done, total)``; no-op on a non-positive total."""
        done, total = update
        if total and total > 0:
            if bar.maximum() != int(total):
                bar.setRange(0, int(total))
            bar.setValue(int(done))
