# -*- coding: utf-8 -*-

"""Bundled per-amino-acid D₂O labeling-site coefficient tables.

Each ``<name>.csv`` has columns ``amino_acid, coefficient`` and is consumed by
:func:`riana.core.fitting.load_aa_coefficients` (resolved by stem via
``riana fit --coefficients <name>``):

- ``commerford`` — Commerford, Carsten & Cronkite 1983, labelable hydrogens per
  residue. The literature/mammalian general default.
- ``ac16`` / ``ipsc`` / ``cm`` — calibration-derived tables bootstrapped from the
  M2/M3 D₂O mixing series (``cm`` uses the recommended drop-time50 fit). Use the
  one matching your cell type; re-derive your own for a new cell type.
"""
