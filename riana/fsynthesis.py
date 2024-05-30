# -*- coding: utf-8 -*-

""" Functions for calculating fractional synthesis rates. """

from riana.utils import strip_concat, get_peptide_distribution
from riana.accmass import calculate_ion_mz, count_atoms
from riana.constants import H_MASS, label_deuterium_de, label_deuterium_commerford
from riana import accmass, constants

import numpy as np


def calculate_a0(sequence: str,
                 label: str,
                 ) -> float:
    """
    Calculates the initial isotope enrichment of a peptide prior to heavy water labeling

    :param sequence:    str: concat sequences
    :param label:       str: aa, hw, hw_cell, or o18, if aa, return 1 assuming no heavy prior to labeling
    :return:            float: mi at time 0
    """

    if label == 'aa':
        return 1

    else:
        sequence = strip_concat(sequence)
        res_atoms = accmass.count_atoms(sequence)
        a0 = np.product([np.power(constants.iso_abundances[i], res_atoms[i]) for i, v in enumerate(res_atoms)])
        # TODO: this should calculate the full isotopic distribution
        return a0


def calculate_label_n(sequence: str,
                      label: str,
                      aa_res: str = 'K',
                      ) -> int:
    """
    Calculates labeling sites of the peptide sequence in heavy water
    or amino acid labeling

    :param sequence:    the peptide sequence
    :param label:       aa, hw, o18, or hw_cell; if aa, only return the labelable residues
    :param aa_res:      the amino acid being labeled
    :return:
    """

    # strip modification site and charge from concat sequence
    sequence = strip_concat(sequence)

    # if amino acid labeling, return number of labeled residues
    if label == 'aa':
        # return the sum of each residue in the aa_res
        return sum([sequence.count(i) for i in aa_res])

    # if heavy water (in vivo), return the number of labeling site in heavy water labeling in vivo
    elif label == 'hw':
        return int(sum([constants.label_deuterium_commerford.get(char) for char in sequence]))

    # if heavy water cell, return the differential evolution best fit values
    elif label == 'hw_cell':
        return int(sum([constants.label_deuterium_de.get(char) for char in sequence]))

    # else if o18, return the number of labeling sites for o18
    elif label == 'o18':
        return int(sum([constants.label_oxygens.get(char) for char in sequence]) - 1)


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
    :param label:   aa, hw, or o18
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


def calculate_fs_fine_structure(a: np.ndarray,
                                seq: str,
                                label: str,
                                ria_max: float,
                                formula: str = 'm0_m1',
                                ) -> float:
    """
    Calculates fractional synthesis by calculating the complete isotope envelop using a fine structure calculator,
    the number of labeling sites, and an artificial element. The empirical mi is then matched the theoretical mi
    of any isotope ratio in the fine structure calculator given the fractional synthesis rate. This may present a more
    accurate method of calculating fractional synthesis rates than the calculate_fs_m0 method.
    :param a:                   m_i at a particular time
    :param seq:                 the peptide sequence
    :param label:               hw, hw_cell, or o18
    :param ria_max:             the precursor RIA
    :param formula:             the formula to calculate the fractional synthesis rate
    :return:                    the fractional synthesis rate
    """

    def find_nearest_fs(array, value):
        """
        Finds the nearest fractional synthesis rate in the array to the value
        :param array:
        :param value:
        :return:
        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    seq = strip_concat(seq)

    # Calculate the number of labeling sites
    num_labeling_sites = calculate_label_n(sequence=seq,
                                           label=label)

    # print(num_labeling_sites)

    # Prelabeling distribution
    initial = get_peptide_distribution(seq,
                                       label=label,
                                       num_labeling_sites=num_labeling_sites,
                                       )

    # Postlabeling final distribution at 100% FS
    final = get_peptide_distribution(seq,
                                     label=label,
                                     num_labeling_sites=num_labeling_sites,
                                     deuterium_enrichment_level=ria_max)

    # Get the summed probability of each isotopomer
    pep_mass = calculate_ion_mz(seq=seq)

    initial_envelop = [
        sum([p for (m, p) in zip(initial.masses, initial.probs) if np.abs(m - (pep_mass + isotopomer * H_MASS)) <= 0.1])
        for isotopomer in range(0, 8)]

    final_envelop = [
        sum([p for (m, p) in zip(final.masses, final.probs) if np.abs(m - (pep_mass + isotopomer * H_MASS)) <= 0.1]) for
        isotopomer in range(0, 8)]

    # print(f'Initial envelop: {initial_envelop}')
    # print(f'Final envelop: {final_envelop}')

    # For fractional synthesis from 0% to 100%, mix the initial and final envelop
    fs_array = []

    # For 0.01 to 1.00, mix the initial and final envelop
    for fs in np.arange(0, 1.01, 0.01):
        mixed_envelop = [a * (1 - fs) + b * fs for a, b in zip(initial_envelop, final_envelop)]

        if formula == 'm0_m1':
            fs_array.extend([mixed_envelop[0] / mixed_envelop[1]])

        elif formula == 'm0_m2':
            fs_array.extend([mixed_envelop[0] / mixed_envelop[2]])

        elif formula == 'm0_m3':
            fs_array.extend([mixed_envelop[0] / mixed_envelop[3]])

        elif formula == 'm1_m2':
            fs_array.extend([mixed_envelop[1] / mixed_envelop[2]])

        elif formula == 'm1_m3':
            fs_array.extend([mixed_envelop[1] / mixed_envelop[3]])

        elif formula == 'm0_mA':
            fs_array.extend([mixed_envelop[0] / sum(mixed_envelop[0:6])])

        elif formula == 'm1_mA':
            fs_array.extend([mixed_envelop[1] / sum(mixed_envelop[0:6])])

        elif formula == 'Auto':
            if num_labeling_sites < 15:
                fs_array.extend([mixed_envelop[0] / mixed_envelop[1]])
            elif num_labeling_sites > 35:
                fs_array.extend([mixed_envelop[1] / mixed_envelop[3]])
            else:
                fs_array.extend([mixed_envelop[0] / mixed_envelop[2]])

    print(f'FS array: {fs_array}')
    # Then do a reverse lookup of the empirical m_i to get fs.
    predicted_fs = np.array([find_nearest_fs(fs_array, a_i)/100 for a_i in a])
    print(f'Input array: {a}')
    print(f'Predicted FS: {predicted_fs}')

    return predicted_fs



