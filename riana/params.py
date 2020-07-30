# -*- coding: utf-8 -*-

""" Parameters """

# integration params
proton_mass = 1.007276466621  # proton mass, identical to scipy.constants.physical_constants['proton mass in u']
# neutron = 1.00866491595
electron_mass = 0.000548579907

# isotope mass from NIST https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
# electron mass from NIST https://www.physics.nist.gov/cgi-bin/cuu/Value?meu|search_for=electron+mass

# proton mass is different from hydrogen mass or neutron mass
# add neutron mass instead of proton mass for the isotopes which may
# make a difference when iso is high enough (e.g., 12 for KK determination)
# and account for mass defect

# mass defect of deuterium is 1.007276466621 +  1.00866491595 + 0.000548579907 - 2.01410177812 = 0.00238818435
# mass defect of C13 is 6 * 1.007276466621 + 7 * 1.00866491595 + 6 * 0.000548579907 - 13.00335483507 = 0.1042
c13_mass_diff = 1.003354835  # mass difference of C13 - C12 = 1.003354835
deuterium_mass_diff = 1.00627674589  # mass difference of deuterium - protium = 1.00627674589

# match between runs
soft_threshold_q = 0.25  # soft threshold for q-value to include peptides if found in other time points
