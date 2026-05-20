# -*- coding: utf-8 -*-

""" Backwards-compatibility shim.

Kinetic models were lifted to :mod:`riana.core.models` in the M3 (1.0.0)
restructure. This re-export keeps the legacy ``riana_fit`` CLI working while the
rewrite is in progress; the shim is removed together with the legacy pipeline in
M3 Week 4. See PROJECT_REVIEW.md §3.
"""

from riana.core.models import (  # noqa: F401
    one_exponent,
    plot_model,
    two_compartment_fornasiero,
    two_compartment_guan,
)
