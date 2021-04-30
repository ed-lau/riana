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
proton_mass = 1.007276466621  # proton mass, identical to scipy.constants.physical_constants['proton mass in u']
neutron = 1.00866491595
electron_mass = 0.000548579907

c13_mass_diff = 1.003354835  # mass difference of C13 - C12 = 1.003354835
deuterium_mass_diff = 1.00627674589  # mass difference of deuterium - protium = 1.00627674589

# atomic masses. TODO: check against scipy.constants and NIST
c_mass = 12.0000000
h_mass = 1.00782503223
o_mass = 15.99491461957
n_mass = 14.00307400443
s_mass = 31.9720711744
proton_mass = 1.007276466621

"""
Amino acids

Number of carbon, hydrogen, oxygen, nitrogen, sulfur for amino acids

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

    # Holding off non-canonical AA till later
    'U': [0, 0, 0, 0, 0],  # Selenocysteine
    'X': [0, 0, 0, 0, 0],
    'B': [0, 0, 0, 0, 0],  # Asn or Asp
}

# Carbon, hydrogen, oxygen, nitrogen, sulfur for modifications
mod_atoms = {
    'IAA': [2, 3, 1, 1, 0],
}
