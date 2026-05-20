# -*- coding: utf-8 -*-

""" Theoretical isotope distributions for peptides.

Lifted unchanged from ``riana.utils.get_peptide_distribution`` in the M3 (1.0.0)
restructure (PROJECT_REVIEW.md §3). This is the *production* forward model used
by ``riana fit``.

Note for M3 Week 4: the M2 calibration benchmark carries an independent D2O
forward model in ``tests/benchmark/_helpers/forward_model.py`` (envelope binning,
Spep loss, FS solver). Week 4 lifts that solver into this module, but the
benchmark copy stays frozen and independent — the M2 regression gate depends on
its oracle not importing production code. The two also differ in envelope-binning
convention (integer-nominal ±0.5 vs ``H_MASS``-step ±0.1); reconcile in Week 4,
not before.
"""

import IsoSpecPy
from IsoSpecPy import IsoTotalProb

from riana.algorithms.mass_calc import count_atoms


def get_peptide_distribution(peptide: str,
                             deuterium_enrichment_level: float = None,
                             label: int = 1,
                             num_labeling_sites: int = 0,
                             ) -> IsoSpecPy.Iso:

    """
    Calculates the total isotope distribution of a peptide given the peptide sequence and deuterium enrichment level

    :param peptide:                     the peptide sequence
    :param deuterium_enrichment_level:  the deuterium enrichment level of the sample
    :param label:       int: 1=2H_in_vivo, 2=2H_in_vitro, 3=18O, 4=AA, if AA, return 1 assuming no heavy prior to labeling
    :param num_labeling_sites:          the number of labeling sites
    :return:                            IsoSpecPy Distribution of atom counts, isotope masses, and isotope probabilities
    """

    # Check that label must be one of hw, hw_cell, or o18
    assert label in [1, 2, 3], 'Label must be one of 1 (2H_in_vivo), 2 (2H_in_vitro), or 3 (18O)'

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

    if label == 1 or label == 2:
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

    elif label == 3:
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
