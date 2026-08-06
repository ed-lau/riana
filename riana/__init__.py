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

# PEP 440. 1.0.0 (released 2026-06-24) = the M1–M4 rewrite + the Track B N_ISO
# line (adaptive N_ISO, H4′ FS solve, --fs limited-isotopomer scoring). 1.1.0
# (released 2026-07-02) added the experimental science: the ¹⁸O rewrite,
# Δmass-over-time GUI + mass-defect→θ fitting. The active dev line is now 1.2.0
# (branch ``1.2.0``): multiplexed (TMT/dimethyl) + single-timepoint labeling, the
# linear-model WLS default + per-point Var(θ) weighting, and per-condition rollup
# curation. Held at a ``.dev`` marker until the 1.2.0 release.
__version_info__ = ("1", "2", "0", "dev0")
__version__ = ".".join(__version_info__)
