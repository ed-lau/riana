# -*- coding: utf-8 -*-

""" Utility functions for riana.

``get_peptide_distribution`` was lifted to :mod:`riana.algorithms.isotope_dist`
in the M3 (1.0.0) restructure and is re-exported here for backwards
compatibility. ``strip_concat`` (peptide-identifier cleaning) still lives here;
M3 Week 2 moves it into the typed ``riana.io`` parsers. See PROJECT_REVIEW.md §3.
"""
import re

from riana.algorithms.isotope_dist import get_peptide_distribution  # noqa: F401 — re-export shim


def strip_concat(sequence: str,
                 ) -> str:
    """
    Cleans up concat sequences (peptide_charge) and remove modifications
    to return peptide string for labeling site calculations

    :param sequence:    concat sequence containing charge and modificaitons
    :return:
    """
    # 2021-05-18 strip all N-terminal n from Comet
    sequence = re.sub(r'^n', '', sequence)

    # Strip all modifications
    sequence = re.sub(r'\[.*?\]', '', sequence)

    # Strip the underscore and charge
    sequence = re.sub(r'_[0-9]+', '', sequence)

    return sequence
