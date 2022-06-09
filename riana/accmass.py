# -*- coding: utf-8 -*-

""" functions that calculate accurate mass for peptides """

import re
from riana import constants, params


def _count_residue_atoms(seq: str,
                         iaa: bool = True,
                         ) -> list:
    """
    given an peptide sequence, count the atoms of carbon, hydrogen, oxygen, nitrogen, sulfur
    in the residue
    TODO: add in selenocysteine and allow other modifications

    :param seq:     str: amino acid sequence
    :param iaa:     bool: whether cysteins are modified by iodoacetamide
    :return:        list: atom counts [C, H, O, N, S]
    """

    tot_atoms: list = [0, 0, 0, 0, 0]

    for char in seq:
        try:
            aa_atoms = constants.aa_atoms[char]
            tot_atoms = [tot_atoms[i] + aa_atoms[i] for i in range(len(aa_atoms))]
        except KeyError:
            raise KeyError

    if iaa:
        num_cysteines = seq.count('C')
        mod_atoms = [atom * num_cysteines for atom in constants.mod_atoms['IAA']]
        tot_atoms = [tot_atoms[i] + mod_atoms[i] for i in range(len(tot_atoms))]

    return tot_atoms


def count_atoms(sequence: str,
                iaa: bool = True,
                ) -> list:
    """
    wrapper for _count_residue_atoms that returns the full peptide atom count

    :param sequence:    str: peptide seuence
    :param iaa:         bool: whether to add iaa atoms to cysteines
    :return:            list: atom counts [C, H, O, N, S]
    """

    res_atoms = _count_residue_atoms(sequence, iaa=iaa,  # add iodoacetamide to cysteine
                                     )

    # Add one oxygen and two hydrogen for peptide mass
    terminal_atoms = [0, 2, 1, 0, 0]

    return [res_atoms[i] + terminal_atoms[i] for i, v in enumerate(res_atoms)]


def _calc_atom_mass(atoms: list,
                    ) -> float:
    """
    given a list of atoms [C, H, O, N, S], return accurate mass

    :param atoms:   list [C, H, O, N, S]
    :return:        float accurate monoisotopic mass
    """

    mass_vec = [constants.C_MASS,
                constants.H_MASS,
                constants.O_MASS,
                constants.N_MASS,
                constants.S_MASS]

    # Get dot product between atom list and mass vector
    mass = sum([atoms[i] * mass_vec[i] for i in range(len(atoms))])

    return mass


def calculate_ion_mz(seq: str,
                     ion: str = 'M',
                     charge: int = 0
                     ) -> float:
    """
    given a peptide sequence and ion type, count the number of atoms, accounting for ion
    type and whether cysteines are measured by IAA

    - ion type
    M: full peptide parent ion (with H2O)
    b: b ion (no addition)
    y: y ion (with H2O)

    :param seq: str amino acid sequence with modifications defined by []
    :param ion: str ion type (default: M to return peptide mass)
    :param charge: int numerical charge (default: 0 to return peptide mass)
    :return: float accurate mass
    """

    assert type(charge) == int, "Charge must be integer."

    mass = 0

    # First, strip all mass shifts and add them to the starting mass
    try:
        mods = [float(mod[1:-1]) for mod in re.findall('\\[.*?]', seq)]
    except ValueError:
        raise ValueError('Modification contains string characters.')

    # 2021-11-22 exclude label mass from peptide mass calculation
    mass += sum(m for m in mods if m not in params.label_mass)

    # 2021-05-18 strip all N-terminal n from Comet
    seq = re.sub('^n', '', seq)

    # Strip all modifications
    stripped = re.sub('\\[.*?]', '', seq)

    res_atoms = _count_residue_atoms(stripped,
                                     iaa=params.iaa,  # add iodoacetamide to cysteine
                                     )

    # dictionary for complementary atoms to add to ion types
    comp_atom_dict = {
        'M':  [0, 2, 1, 0, 0],
        'b':  [0, 0, 0, 0, 0],
        'y':  [0, 2, 1, 0, 0],
        'b_': [0, -2, -1, 0, 0],
        'y_': [0, 0, 0, 0, 0],
    }
    comp_atoms = comp_atom_dict[ion]

    ion_atoms = [res_atoms[i] + comp_atoms[i] for i, v in enumerate(res_atoms)]

    mass += _calc_atom_mass(ion_atoms)

    # Return peptide mass if charge is 0
    if charge > 0:
        mz = (mass + constants.PROTON_MASS * charge) / charge
        return mz

    if charge < 0:
        raise ValueError('Negative charges are not supported.')

    return mass
