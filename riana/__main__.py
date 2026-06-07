# -*- coding: utf-8 -*-

"""riana.__main__: executed when the riana package is run as ``python -m riana``."""

from riana.cli import main

# The ``__name__ == "__main__"`` guard is required, not cosmetic: the GUI's
# ProcessPoolExecutor uses the *spawn* start method (macOS default; the only
# Qt-safe one), and each spawned worker re-imports this module. Without the
# guard, the worker would re-run ``main()`` and multiprocessing would raise the
# "attempt to start a new process before bootstrapping has finished" error.
if __name__ == "__main__":
    main()