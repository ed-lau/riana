# -*- coding: utf-8 -*-

""" Riana package.

Mid-M3 (1.0.0) rewrite. The package is being restructured from the flat 0.9.0
layout into ``core/``, ``algorithms/``, and ``io/`` subpackages with typed
records and config (see PROJECT_REVIEW.md §3). During the rewrite the legacy
flat modules remain as thin re-export shims so the 0.9.0 CLI keeps working as a
live regression gate; the shims are removed with the legacy pipeline in Week 4.
"""

__all__ = [
    # M3 structural layer
    "config",
    "records",
    "algorithms",
    "core",
    "io",
    # stable shared modules
    "constants",
    "exceptions",
    "logger",
    "params",
    "utils",
    # legacy 0.9.0 pipeline (riana_integrate / riana_fit consume the shims;
    # all five are replaced by core/ + io/ + cli.py in Week 4)
    "accmass",
    "fsynthesis",
    "models",
    "main",
    "peptides",
    "project",
    "riana_fit",
    "riana_integrate",
    "spectra",
]

# PEP 440: 1.0.0.dev0 marks the in-progress M3 rewrite. Bumps to 1.0.0 when M3
# (Week 4: CLI + fit) completes.
__version_info__ = ("1", "0", "0", "dev0")
__version__ = ".".join(__version_info__)
