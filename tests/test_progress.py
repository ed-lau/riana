"""Unit tests for the dependency-free progress reporter (riana.progress)."""

from __future__ import annotations

import io
import logging

from riana.progress import ProgressReporter, iter_progress


def _tty_stream() -> io.StringIO:
    s = io.StringIO()
    s.isatty = lambda: True          # force the carriage-return bar path
    return s


def test_iter_progress_passes_through_and_ticks():
    seen = []
    out = list(iter_progress(range(5), 5, lambda d, t: seen.append((d, t))))
    assert out == [0, 1, 2, 3, 4]                 # values untouched
    assert seen == [(1, 5), (2, 5), (3, 5), (4, 5), (5, 5)]


def test_iter_progress_none_callback_is_passthrough():
    assert list(iter_progress(range(3), 3, None)) == [0, 1, 2]


def test_reporter_non_tty_emits_periodic_log_lines():
    buf = io.StringIO()
    log = logging.getLogger("riana.test.progress.nontty")
    log.handlers[:] = [logging.StreamHandler(buf)]
    log.setLevel(logging.INFO)
    rep = ProgressReporter(0, "fit", log, log_every=25)   # non-TTY (StringIO)
    for i in range(1, 9):
        rep(i, 8)
    lines = buf.getvalue().strip().splitlines()
    assert lines == ["fit: 2/8 (25%)", "fit: 4/8 (50%)",
                     "fit: 6/8 (75%)", "fit: 8/8 (100%)"]


def test_reporter_tty_draws_bar_and_finishes_with_newline():
    s = _tty_stream()
    rep = ProgressReporter(0, "fit", None, width=10, stream=s)
    for i in (1, 5, 10):
        rep(i, 10)
    text = s.getvalue()
    assert "\r" in text and "█" in text           # in-place bar
    assert "10/10 (100%)" in text and text.endswith("\n")


def test_reporter_multi_phase_resets_on_counter_restart():
    s = _tty_stream()
    rep = ProgressReporter(0, "fit", None, width=6, stream=s)
    for i in (1, 2):          # phase 1 (e.g. condition A): total 2
        rep(i, 2)
    for i in (1, 2, 3):       # phase 2 (condition B): total 3 — counter restarts
        rep(i, 3)
    # Both phases reach 100% and each finishes its own bar line.
    assert s.getvalue().count("(100%)") == 2
    assert s.getvalue().count("\n") == 2
