# -*- coding: utf-8 -*-

""" hard-coded constants """

"""
Masses

Isotope mass from NIST https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
Electron mass from NIST https://www.physics.nist.gov/cgi-bin/cuu/Value?meu|search_for=electron+mass
Proton mass is different from hydrogen mass or neutron mass: add neutron mass instead of proton mass 
for the isotopes which may make a difference when iso is high enough (e.g., 12 for KK determination)
and account for mass defect

Mass defect of deuterium is 1.007276466621 +  1.00866491595 + 0.000548579907 - 2.01410177812 = 0.00238818435
Mass defect of C13 is 6 * 1.007276466621 + 7 * 1.00866491595 + 6 * 0.000548579907 - 13.00335483507 = 0.1042

"""
PROTON_MASS = 1.007276466621  # proton mass, identical to scipy.constants.physical_constants['proton mass in u']
NEUTRON_MASS = 1.00866491595
ELECTRON_MASS = 0.000548579907

C13_MASSDIFF = 1.003354835  # mass difference of C13 - C12 = 1.003354835
D_MASSDIFF = 1.00627674589  # mass difference of deuterium - protium = 1.00627674589
SILAC_MASSDIFF = 1.001      # Average mass difference in SILAC (mixture of N15 or C13)

# atomic masses.
# TODO: check against scipy.constants and NIST
C_MASS = 12.0000000
H_MASS = 1.00782503223
O_MASS = 15.99491461957
N_MASS = 14.00307400443
S_MASS = 31.9720711744

"""
Amino acids: number of carbon, hydrogen, oxygen, nitrogen, sulfur for amino acids
"""

aa_atoms = {
    'A': [3, 5, 1, 1, 0],
    'C': [3, 5, 1, 1, 1],
    'D': [4, 5, 3, 1, 0],
    'E': [5, 7, 3, 1, 0],
    'F': [9, 9, 1, 1, 0],
    'G': [2, 3, 1, 1, 0],
    'H': [6, 7, 1, 3, 0],
    'I': [6, 11, 1, 1, 0],
    'K': [6, 12, 1, 2, 0],
    'L': [6, 11, 1, 1, 0],
    'M': [5, 9, 1, 1, 1],
    'N': [4, 6, 2, 2, 0],
    'P': [5, 7, 1, 1, 0],
    'Q': [5, 8, 2, 2, 0],
    'R': [6, 12, 1, 4, 0],
    'S': [3, 5, 2, 1, 0],
    'T': [4, 7, 2, 1, 0],
    'V': [5, 9, 1, 1, 0],
    'W': [11, 10, 1, 2, 0],
    'Y': [9, 9, 2, 1, 0],

    # TODO: include non-canonical AA
    'U': [0, 0, 0, 0, 0],  # Selenocysteine
    'X': [0, 0, 0, 0, 0],
    'B': [0, 0, 0, 0, 0],  # Asn or Asp
}

# Carbon, hydrogen, oxygen, nitrogen, sulfur for modifications
mod_atoms = {
    'IAA': [2, 3, 1, 1, 0],
}

# Commerford, Carsten, and Cronkite 1983 Table 1
# Number of labelable hydrogen atoms per amino acids
# Relative specific activity * H/mole
label_hydrogens = {
    'A': 4.00,
    'C': 1.62,
    'D': 1.89,
    'E': 3.95,
    'F': 0.32,
    'G': 2.06,
    'H': 2.88,
    'I': 1.00,
    'K': 0.54,
    'L': 0.69,
    'M': 1.12,
    'N': 1.89,
    'P': 2.59,
    'Q': 3.95,
    'R': 3.34,
    'S': 2.61,
    'T': 0.20,
    'V': 0.56,
    'W': 0.08,
    'Y': 0.42,
}

# Number of oxygen labels should be 1-#aa plus an extra 1 for every S or T
label_oxygens = {
    'A': 1.,
    'C': 1.,
    'D': 1.,
    'E': 1.,
    'F': 1.,
    'G': 1.,
    'H': 1.,
    'I': 1.,
    'K': 1.,
    'L': 1.,
    'M': 1.,
    'N': 1.,
    'P': 1.,
    'Q': 1.,
    'R': 1.,
    'S': 2.,
    'T': 2.,
    'V': 1.,
    'W': 1.,
    'Y': 1.,
}

"""
Calculate natural abundance of isotopes from NIST 
https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
Berglund and Wieser 2019 http://www.iupac.org/publications/pac/83/2/0397/ 
"""
iso_abundances = [0.9893,  # C12
                  0.999885,  # H1
                  0.9975,  # O16
                  0.99636,  # N14
                  0.9499,  # S32
                  ]
