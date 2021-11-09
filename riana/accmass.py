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

    :param seq: str amino acid sequence
    :param iaa: bool whether cysteins are modified by iodoacetamide
    :return: list [c, h, o, n, s]
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


def _calc_atom_mass(atoms: list,
                    ) -> float:
    """
    given a list of atoms [c, h, o, n, s], return accurate mass
    :param atoms: list [c, h, o, n, s]
    :return: float accurate monoisotopic mass
    """

    mass_vec = [constants.c_mass,
                constants.h_mass,
                constants.o_mass,
                constants.n_mass,
                constants.s_mass]

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

    mass += sum(mods)

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
        mz = (mass + constants.proton_mass * charge)/charge
        return mz

    if charge < 0:
        raise ValueError('Negative charges are not supported.')

    return mass
