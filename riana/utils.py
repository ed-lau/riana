# -*- coding: utf-8 -*-

""" Utility functions for riana. """
import re
from riana.accmass import calculate_ion_mz, count_atoms
from riana.constants import H_MASS, label_oxygens, label_deuterium_de, label_deuterium_commerford
import IsoSpecPy
from IsoSpecPy import IsoTotalProb


def strip_concat(sequence: str,
                 ) -> str:
    """
    Cleans up concat sequences (peptide_charge) and remove modifications
    to return peptide string for labeling site calculations

    :param sequence:    concat sequence containing charge and modificaitons
    :return:
    """
    # 2021-05-18 strip all N-terminal n from Comet
    sequence = re.sub('^n', '', sequence)

    # Strip all modifications
    sequence = re.sub('\\[.*?\\]', '', sequence)

    # Strip the underscore and charge
    sequence = re.sub('_[0-9]+', '', sequence)

    return sequence


def get_peptide_distribution(peptide: str,
                             deuterium_enrichment_level: float = None,
                             label: str = 'hw',
                             num_labeling_sites: int = 0,
                             ) -> IsoSpecPy.Iso:

    """
    Calculates the total isotope distribution of a peptide given the peptide sequence and deuterium enrichment level

    :param peptide:                     the peptide sequence
    :param deuterium_enrichment_level:  the deuterium enrichment level of the sample
    :param label:                       the type of labeling, either 'hw'/'hw_cell' for heavy water or 'o18' for oxygen-18
    :param num_labeling_sites:          the number of labeling sites
    :return:                            IsoSpecPy Distribution of atom counts, isotope masses, and isotope probabilities
    """

    # Check that label must be one of hw, hw_cell, or o18
    assert label in ['hw', 'hw_cell', 'o18'], 'Label must be one of hw, hw_cell, or o18'

    if deuterium_enrichment_level is not None:
        assert 0 < deuterium_enrichment_level <= 1, 'Deuterium enrichment level must be greater than 0 and no greater than 1'

    # Get C, H, O, N, S count using the Riana count_atoms function
    peptide_atoms = count_atoms(peptide)
    # print(peptide_atoms)

    # Supply atom counts to IsoSpecPy.IsoParamsFromDict and unpack to get atom counts, isotope masses. and probabilities
    atom_count_list, isotope_mass_list, isotope_probability_list, _ = IsoSpecPy.IsoParamsFromDict(formula={"C": peptide_atoms[0],
                                                                                                           "H": peptide_atoms[1],
                                                                                                           "O": peptide_atoms[2],
                                                                                                           "N": peptide_atoms[3],
                                                                                                           "S": peptide_atoms[4]})

    if label == 'hw' or label == 'hw_cell':
        # Subtract the number of labeling sites from hydrogen, extend the atom count list with accessible deuterium count
        atom_count_list[1] = atom_count_list[1] - num_labeling_sites
        atom_count_list.extend([num_labeling_sites])
        # print(f'Atom count list: {atom_count_list}')
        # Extend the isotope mass list for deuterium, which is the same as hydrogen
        isotope_mass_list.extend([isotope_mass_list[1]])
        # print(f'Isotope mass list: {isotope_mass_list}')

        # Extend the isotope probabilities for labelable hydrogen sites, which is the isotope enrichment level
        # For the unlabeled samples, we should use the background level of 0.0001157
        if deuterium_enrichment_level is None:
            isotope_probability_list.extend([isotope_probability_list[1]])
        else:
            isotope_probability_list.extend([(1-deuterium_enrichment_level, deuterium_enrichment_level)])
            # TODO: include the background deuterium level here too?

    elif label == 'o18':
        atom_count_list[2] = atom_count_list[2] - num_labeling_sites
        atom_count_list.extend([num_labeling_sites])
        # Extend the isotope mass list for O18, which is the same as O16
        isotope_mass_list.extend([isotope_mass_list[2]])

    # print(f'Isotope probability list: {isotope_probability_list}')

    isotope_dist = IsoTotalProb(prob_to_cover=.999,
                       atomCounts=atom_count_list,
                       isotopeMasses= isotope_mass_list,
                       isotopeProbabilities=isotope_probability_list,
                       use_nominal_masses = True)

    return isotope_dist
