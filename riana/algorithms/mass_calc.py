# -*- coding: utf-8 -*-

""" Accurate-mass calculation for peptides.

Lifted unchanged from ``riana.accmass`` in the M3 (1.0.0) restructure — this is
a numerically-sensitive science module, so it is moved verbatim rather than
rewritten (PROJECT_REVIEW.md §3, "lifted unchanged"). The ``params`` global
dependency in :func:`calculate_ion_mz` is carried along deliberately; threading
``iaa`` through ``IntegrationConfig`` is deferred to the Week 3 pipeline rewrite.
"""

import re
from riana import constants, params


def _count_residue_atoms(seq: str,
                         iaa: bool = True,
                         ) -> list:
    """
    given an peptide sequence, count the atoms of carbon, hydrogen, oxygen, nitrogen, sulfur, phosphorus
    in the residue
    TODO: add in selenocysteine and allow other modifications

    :param seq:     str: amino acid sequence
    :param iaa:     bool: whether cysteins are modified by iodoacetamide
    :return:        list: atom counts [C, H, O, N, S, P]
    """

    tot_atoms: list = [0, 0, 0, 0, 0, 0]

    for char in seq:
        try:
            aa_atoms = constants.aa_atoms[char]
            tot_atoms = [tot_atoms[i] + aa_atoms[i] for i in range(len(aa_atoms))]
        except KeyError:
            raise KeyError

    if iaa:
        # Carbamidomethyl (UNIMOD:4) as a fixed mod on every cysteine — the
        # composition now comes from the unified UniMod-keyed ``mod_atoms``
        # table (was a dedicated ``'IAA'`` key). M7 will retire this fixed-mod
        # flag in favour of passing CAM through ``count_atoms(mods=...)`` once
        # the IO layer threads per-cysteine mods.
        num_cysteines = seq.count('C')
        cam_atoms = [atom * num_cysteines for atom in constants.mod_atoms[4]]
        tot_atoms = [tot_atoms[i] + cam_atoms[i] for i in range(len(tot_atoms))]

    return tot_atoms


def count_atoms(sequence: str,
                iaa: bool = True,
                mods: list = (),
                ) -> list:
    """
    wrapper for _count_residue_atoms that returns the full peptide atom count

    :param sequence:    str: peptide seuence
    :param iaa:         bool: whether to add iaa atoms to cysteines
    :param mods:        iterable of UniMod accession ids (ints) for variable
                        modifications carried by this peptidoform; each mod's
                        ``[C, H, O, N, S, P]`` composition (``constants.mod_atoms``)
                        is added to the envelope formula (M7). Empty by default
                        so the bare-sequence path stays byte-identical.
    :return:            list: atom counts [C, H, O, N, S, P]
    """

    res_atoms = _count_residue_atoms(sequence, iaa=iaa,  # add iodoacetamide to cysteine
                                     )

    # Add one oxygen and two hydrogen for peptide mass
    terminal_atoms = [0, 2, 1, 0, 0, 0]

    atoms = [res_atoms[i] + terminal_atoms[i] for i, v in enumerate(res_atoms)]

    # Variable modifications (M7): add each UniMod's atom composition so the
    # IsoSpec envelope reflects the modified peptidoform, not the bare backbone.
    for unimod_id in mods:
        comp = constants.mod_atoms[unimod_id]
        atoms = [atoms[i] + comp[i] for i in range(len(atoms))]

    return atoms


def _calc_atom_mass(atoms: list,
                    ) -> float:
    """
    given a list of atoms [C, H, O, N, S, P], return accurate mass

    :param atoms:   list [C, H, O, N, S, P]
    :return:        float accurate monoisotopic mass
    """

    mass_vec = [constants.C_MASS,
                constants.H_MASS,
                constants.O_MASS,
                constants.N_MASS,
                constants.S_MASS,
                constants.P_MASS]

    # Get dot product between atom list and mass vector
    mass = sum([atoms[i] * mass_vec[i] for i in range(len(atoms))])

    return mass


def calculate_ion_mz(seq: str,
                     ion: str = 'M',
                     charge: int = 0,
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
        mods = [float(mod[1:-1]) for mod in re.findall(r'\[.*?]', seq)]
    except ValueError:
        raise ValueError('Modification contains string characters.')

    # Every bracketed modification mass contributes to the precursor mass.
    mass += sum(mods)

    # 2021-05-18 strip all N-terminal n from Comet
    seq = re.sub(r'^n', '', seq)

    # Strip all modifications
    stripped = re.sub(r'\[.*?]', '', seq)

    res_atoms = _count_residue_atoms(stripped,
                                     iaa=params.iaa,  # add iodoacetamide to cysteine
                                     )

    # dictionary for complementary atoms to add to ion types
    comp_atom_dict = {
        'M':  [0, 2, 1, 0, 0, 0],
        'b':  [0, 0, 0, 0, 0, 0],
        'y':  [0, 2, 1, 0, 0, 0],
        'b_': [0, -2, -1, 0, 0, 0],
        'y_': [0, 0, 0, 0, 0, 0],
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
