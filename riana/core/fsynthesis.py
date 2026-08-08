# -*- coding: utf-8 -*-

""" Functions for calculating fractional synthesis rates.

Lifted from ``riana.fsynthesis`` in the M3 (1.0.0) restructure. One science-layer
fix is applied during the lift: ``calculate_a0`` previously tested
``label == 'aa'`` (a string), but callers dispatched with an integer label, so the
amino-acid ``a_0`` branch was unreachable and AA experiments silently used the
natural-abundance baseline (PROJECT_REVIEW.md §2b). With the v1.1.0 label-taxonomy
cleanup the labels are now the strings ``"D2O"`` / ``"O18"`` / ``"AA"``, and the
test is ``label == "AA"``.

The fixed-site-count ``m0`` FS model here (``calculate_a0`` / ``calculate_fs_m0``)
is superseded in production by the IsoSpec forward/solve model
(``riana.core.fitting`` / ``riana.algorithms.isotope_dist``); it is retained only as
the reference baseline that ``tests/benchmark/bench_fs_method_compare`` compares the
production solver against.
"""

from riana.utils import strip_concat
from riana.algorithms.mass_calc import count_atoms
from riana import constants

import numpy as np


def calculate_a0(sequence: str,
                 label: str,
                 ) -> float:
    """
    Calculates the initial isotope enrichment of a peptide prior to heavy water labeling

    :param sequence:    str: concat sequences
    :param label:       str: the labeling chemistry — ``"D2O"`` (heavy water) or
                        ``"O18"`` (¹⁸O); the legacy ``"AA"`` (amino-acid labeling)
                        returns 1, assuming no heavy prior to labeling
    :return:            float: mi at time 0
    """

    # M3 fix (PROJECT_REVIEW.md §2b): was ``label == 'aa'`` — a string test that
    # never matched the integer label callers passed at the time, so AA
    # experiments silently fell through to the natural-abundance branch below.
    if label == "AA":
        return 1

    else:
        sequence = strip_concat(sequence)
        res_atoms = count_atoms(sequence)
        a0 = np.prod([np.power(constants.iso_abundances[i], res_atoms[i]) for i, v in enumerate(res_atoms)])
        # TODO: this should calculate the full isotopic distribution
        return a0


def calculate_fs_m0(a: np.ndarray,
                    seq: str,
                    label: str,
                    ria_max: float,
                    num_labeling_sites: int,
                    ) -> float:
    """
    Calculates fractional synthesis based on a_t, a_0 (initial), and a_max (asymptote)

    :param a:       m_i at a particular time
    :param seq:     the peptide sequence
    :param label:       str: the labeling chemistry — "D2O" (heavy water), "O18" (¹⁸O), or the legacy "AA" (amino-acid labeling)
    :param ria_max: the precursor RIA
    :param num_labeling_sites: the number of labeling sites

    :return:
    """

    # a0 is the initial m_i value before label onset
    a_0 = calculate_a0(seq, label=label)
    # a max is final m_i value at plateau based on labeling site and precursor RIA
    a_max = a_0 * np.power((1 - ria_max), num_labeling_sites)

    # catch errors from no ria or no labeling site
    if a_max - a_0 == 0:
        # repeat an array of 0 if the input is an ndarray, otherwise return 0
        return np.repeat(0, len(a)) if isinstance(a, np.ndarray) else 0
    else:
        return (a-a_0)/(a_max-a_0)
