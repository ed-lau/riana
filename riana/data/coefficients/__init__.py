# -*- coding: utf-8 -*-

"""Bundled labeling-site coefficient tables (resolved by CSV stem via
``riana fit --coefficients <name>``). Two families, one per ``--label``:

**D₂O (``--label hw``)** — per-AA labile-hydrogen tables, columns
``amino_acid, coefficient``, consumed by
:func:`riana.core.fitting.load_aa_coefficients`:

- ``deberneh_2025_rss`` — in-vivo mouse LC-MS N_EH by the RSS nested-fit. The
  **default/recommended** mammalian table.
- ``ilchenko_2019`` — in-vivo mouse LC-MS M0-M1 N_aa.
- ``commerford_1983`` — Commerford, Carsten & Cronkite tritium literature table.
- ``alamillo_2025_{ac16,ipsc,cm}`` — calibration-derived cell-line tables
  bootstrapped from the D₂O mixing series. Use the one matching your cell type.

**¹⁸O (``--label o18``)** — length-model tables, columns ``feature, coefficient``
(``Spep = b·(L−1) + c_D·D + c_E·E + c_N·N + c_Q·Q + c_S·S``), consumed by
:func:`riana.core.fitting.load_o18_coefficients`:

- ``juber_2026_o18_ac16`` — in-vitro AC16, trained on the ¹⁸O mixing calibration.
- ``rachdaoui_2009_o18`` — in-vivo mouse reference (Rachdaoui/Previs, MCP 2009).

Re-derive your own table for a new cell type / regime.
"""
