# -*- coding: utf-8 -*-

""" Utility functions for riana.

``strip_concat`` (peptide-identifier cleaning) is the live helper here, used by
``core/fitting`` and ``core/fsynthesis`` to normalise concat sequences before
the IsoSpec forward model. (``get_peptide_distribution`` was lifted to
:mod:`riana.algorithms.isotope_dist` in the M3 restructure; its compatibility
re-export was dropped with the legacy modules in M4.)
"""
import re


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
