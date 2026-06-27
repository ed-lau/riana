# -*- coding: utf-8 -*-

"""Bundled labeling-site coefficient tables (resolved by CSV stem via
``riana fit --coefficients <name>``). Two families, one per ``--label``:

Loaders read only the leading key + ``coefficient`` columns; any further
columns (fit-uncertainty / quality stats) are **ignored at load** but kept in
the file as provenance. Fitted (calibration-derived) tables carry those stats;
literature tables (transcribed from a publication) carry only the two columns.

**D₂O (``--label hw``)** — per-AA labile-hydrogen tables, key columns
``amino_acid, coefficient``, consumed by
:func:`riana.core.fitting.load_aa_coefficients`:

- ``deberneh_2025_rss`` — in-vivo mouse LC-MS N_EH by the RSS nested-fit. The
  **default/recommended** mammalian table. *(literature: 2 columns)*
- ``ilchenko_2019`` — in-vivo mouse LC-MS M0-M1 N_aa. *(literature: 2 columns)*
- ``commerford_1983`` — Commerford, Carsten & Cronkite tritium literature table.
  *(literature: 2 columns)*
- ``alamillo_2025_{ac16,ipsc,cm}`` — calibration-derived cell-line tables
  bootstrapped from the D₂O mixing series (``build_frozen_tables.py``); retrained
  2026-06-27 on the ``runs/calib_{cell}_v1`` integration (current defaults — ±10 ppm,
  apex / 0.15 min, iso0-5; cm drops its weak 50 % fraction). Use the one matching
  your cell type. *(fitted: + ``std_error`` (bootstrap SE), ``oob_r2`` (out-of-bag
  R²), ``ci_lo, ci_hi`` (95% bootstrap CI), ``boot_frac_nonzero``)*

**¹⁸O (``--label o18``)** — length-model tables, key columns ``feature, coefficient``
(``Spep = b·(L−1) + c_D·D + c_E·E + c_N·N + c_Q·Q + c_S·S``), consumed by
:func:`riana.core.fitting.load_o18_coefficients`:

- ``juber_2026_o18_ac16`` — in-vitro AC16, trained on the ¹⁸O mixing calibration.
  *(fitted: + ``std_error, oob_r2, ci_lo, ci_hi, boot_frac_nonzero`` — same
  bootstrap-OOB freeze as the D₂O tables, via ``bench_o18_coefficients.py``, so
  the two labels' uncertainty and R² are directly comparable)*
- ``rachdaoui_2009_o18`` — in-vivo mouse reference (Rachdaoui/Previs, MCP 2009).
  *(literature: ``feature, coefficient, note``)*

Re-derive your own table for a new cell type / regime.
"""
