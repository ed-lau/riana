# -*- coding: utf-8 -*-

"""Minimal, dependency-free progress reporting for the long CLI loops
(``fit`` / ``rollup`` / ``integrate``).

The core functions take an optional ``progress_callback(done, total)`` and call
it as items complete. :class:`ProgressReporter` is the CLI's renderer:

- on a **TTY** it draws a single carriage-return bar on stderr (overwritten in
  place, cleared with a newline at 100 %);
- when **piped / redirected / non-interactive** (or stderr is a file) it instead
  emits a log line every ``log_every`` percent, so logfiles and CI output stay
  clean and readable.

It deliberately avoids a hard dependency on rich/tqdm (neither is installed) and
never fights the logger: during these loops the logger is otherwise quiet, so the
bar owns stderr until the post-loop summary line.
"""

from __future__ import annotations

import logging
import sys
from typing import Callable

ProgressCallback = Callable[[int, int], None]


class ProgressReporter:
    """A ``progress_callback`` that renders to a TTY bar or periodic log lines."""

    def __init__(
        self,
        total: int,
        label: str,
        logger: logging.Logger | None = None,
        *,
        width: int = 30,
        log_every: int = 10,
        stream=None,
    ) -> None:
        self.total = max(0, int(total))
        self.label = label
        self.logger = logger
        self.width = width
        self.log_every = log_every
        self.stream = stream if stream is not None else sys.stderr
        self.tty = bool(getattr(self.stream, "isatty", lambda: False)())
        self._last_pct = -1
        self._done = False
        self._prev_done = 0

    def __call__(self, done: int, total: int | None = None) -> None:
        total = self.total if total is None else total
        if total <= 0:
            return
        done = min(done, total)
        # Multi-phase (e.g. a manifest fit over several condition curves): when the
        # counter restarts, finish the previous bar line and reset for the new phase.
        if done < self._prev_done:
            if self.tty and not self._done:
                self.stream.write("\n")
                self.stream.flush()
            self._last_pct = -1
            self._done = False
        self._prev_done = done
        if self._done:
            return
        pct = int(100 * done / total)
        if self.tty:
            filled = int(self.width * done / total)
            bar = "█" * filled + "·" * (self.width - filled)
            self.stream.write(f"\r{self.label}: |{bar}| {done}/{total} ({pct}%)")
            self.stream.flush()
            if done >= total:
                self.stream.write("\n")
                self.stream.flush()
                self._done = True
        elif pct >= self._last_pct + self.log_every or done >= total:
            self._last_pct = pct
            msg = f"{self.label}: {done}/{total} ({pct}%)"
            if self.logger is not None:
                self.logger.info(msg)
            if done >= total:
                self._done = True

    def close(self) -> None:
        """Finish the bar line if it was left mid-render (e.g. on an early exit)."""
        if self.tty and not self._done:
            self.stream.write("\n")
            self.stream.flush()
            self._done = True


def iter_progress(iterable, total: int, callback: ProgressCallback | None):
    """Yield from *iterable*, calling ``callback(i, total)`` after each item.

    Unifies the serial and process-pool loops: pass a lazy generator (serial) or
    an ``executor.map`` result (pool) — either way progress ticks per completion
    (the pool yields in submission order as chunks finish, so it is real, ordered
    progress). A ``None`` callback makes this a thin pass-through.
    """
    for i, item in enumerate(iterable, 1):
        yield item
        if callback is not None:
            callback(i, total)
