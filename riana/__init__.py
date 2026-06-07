# -*- coding: utf-8 -*-

""" Riana package.

1.0.0 layout: the science and pipeline live in the ``core/``, ``algorithms/``,
and ``io/`` subpackages with typed records and frozen config (see
PROJECT_REVIEW.md §3). M4 retired the flat 0.9.0 modules and their re-export
shims along with the ``--engine legacy`` path; the typed pipeline is the only
engine, driven by the Typer CLI in :mod:`riana.cli`.
"""

__all__ = [
    # structural layer
    "config",
    "records",
    "algorithms",
    "core",
    "io",
    "cli",
    # shared modules
    "constants",
    "exceptions",
    "logger",
    "params",
    "utils",
    "data",
]

# PEP 440: 1.0.0 is the M1–M4 rewrite release (typed pipeline, streaming I/O,
# peak detection, mzTab intake, Typer CLI, PySide6 GUI). The .dev marker is
# dropped now that M4 is complete; the git tag + Zenodo DOI follow at release.
__version_info__ = ("1", "0", "0")
__version__ = ".".join(__version_info__)
