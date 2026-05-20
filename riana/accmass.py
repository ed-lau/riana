# -*- coding: utf-8 -*-

""" Backwards-compatibility shim.

Accurate-mass calculation was lifted to :mod:`riana.algorithms.mass_calc` in the
M3 (1.0.0) restructure. This re-export keeps the legacy ``riana_integrate`` /
``riana_fit`` CLI working while the rewrite is in progress; the shim is removed
together with the legacy pipeline in M3 Week 4. See PROJECT_REVIEW.md §3.
"""

from riana.algorithms.mass_calc import (  # noqa: F401
    _calc_atom_mass,
    _count_residue_atoms,
    calculate_ion_mz,
    count_atoms,
)
