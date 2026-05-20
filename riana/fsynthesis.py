# -*- coding: utf-8 -*-

""" Backwards-compatibility shim.

Fractional-synthesis calculations were lifted to :mod:`riana.core.fsynthesis` in
the M3 (1.0.0) restructure. This re-export keeps the legacy ``riana_fit`` CLI
working while the rewrite is in progress; the shim is removed together with the
legacy pipeline in M3 Week 4. See PROJECT_REVIEW.md §3.
"""

from riana.core.fsynthesis import (  # noqa: F401
    calculate_a0,
    calculate_fs_fine_structure,
    calculate_fs_m0,
    calculate_label_n,
)
