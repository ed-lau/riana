# -*- coding: utf-8 -*-

"""Frozen-app entry point for the PyInstaller bundle (see ``riana-gui.spec``).

Behaviour: launched with **no arguments** (a Finder/Explorer double-click of the
bundle) it opens the GUI; launched **with arguments** it delegates to the full
``riana`` Typer CLI, so the same bundle doubles as a portable command-line binary
(``riana-gui integrate ...``, ``riana-gui --version``).

``multiprocessing.freeze_support()`` is mandatory and must run first: the GUI's
``ProcessPoolExecutor`` (and each tab's ``-W`` workers) spawn workers by
re-executing this frozen binary, and ``freeze_support`` is what makes a re-exec run
the worker payload and exit instead of relaunching the app.
"""

from __future__ import annotations

import multiprocessing
import sys


def main() -> None:
    # Finder passes a "-psn_<serial>" process-serial-number arg on some macOS
    # versions; ignore it so a double-launch still counts as "no arguments".
    argv = [a for a in sys.argv[1:] if not a.startswith("-psn_")]
    if not argv:
        from riana.gui.app import run_gui
        raise SystemExit(run_gui())
    from riana.cli import main as cli_main
    cli_main()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
