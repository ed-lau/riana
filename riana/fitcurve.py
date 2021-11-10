# -*- coding: utf-8 -*-

""" Main. """
import logging
import re
import os
import sys
from functools import partial

from riana.project import ReadDirectory
from riana.peptides import ReadPercolator
from riana.spectra import Mzml

from riana import integrate, params, __version__

import tqdm
import pandas as pd

import numpy as np
from riana import accmass, constants


def strip_concat(sequence: str,
                 ) -> str:
    """
    cleans up concat sequences (peptide_charge) and remove modifications
    to return peptide string for labeling site calculations

    :param sequence:    concat sequence containing charge and modificaitons
    :return:
    """
    # 2021-05-18 strip all N-terminal n from Comet
    sequence = re.sub('^n', '', sequence)

    # Strip all modifications
    sequence = re.sub('\\[.*?]', '', sequence)

    # Strip the underscore and charge
    sequence = re.sub('_[0-9]+', '', sequence)

    return sequence


def calculate_a0(sequence: str,
                 ) -> float:
    """
    calculates the initial isotope enrichment of a peptide prior to heavy water labeling

    :param sequence:    str: concat sequences
    :return:            float: mi at time 0
    """

    sequence = strip_concat(sequence)

    res_atoms = accmass.count_atoms(sequence)

    a0 = np.product([np.power(constants.iso_abundances[i], res_atoms[i]) for i, v in enumerate(res_atoms)])

    return a0


def calculate_label_n(sequence:str,
                      ) -> float:

    sequence = strip_concat(sequence)

    return sum([constants.label_hydrogens.get(char) for char in sequence])


def calculate_fs(a, a_0, a_max):
    return (a-a_0)/(a_max-a_0)


def runfit(args):
    print('Fit function not implemented yet.')
    pass
