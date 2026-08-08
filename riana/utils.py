# -*- coding: utf-8 -*-

""" Utility functions for riana.

``strip_concat`` (peptide-identifier cleaning) is the live helper here, used by
``core/fitting`` and ``core/fsynthesis`` to normalise concat sequences before
the IsoSpec forward model. (``get_peptide_distribution`` was lifted to
:mod:`riana.algorithms.isotope_dist` in the M3 restructure; its compatibility
re-export was dropped with the legacy modules in M4.)
"""
import re

from riana import constants


def is_canonical_peptide(sequence: str) -> bool:
    """True iff every residue is one of the 20 standard amino acids.

    Strips ``[UNIMOD:n]`` mod tokens and the ``_charge`` suffix first (via
    :func:`strip_concat`), then checks the bare residues against
    :data:`riana.constants.CANONICAL_AA`. A peptide carrying a non-canonical residue
    (U/O/B/Z/J/X) has no defined atom composition, so it cannot be mass-computed and is
    dropped at intake. Empty / all-token input returns ``False``.
    """
    core = re.sub(r"[^A-Z]", "", strip_concat(sequence).upper())
    return bool(core) and set(core) <= constants.CANONICAL_AA


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
