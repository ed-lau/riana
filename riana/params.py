# -*- coding: utf-8 -*-

""" Parameters """

"""
Peptides

"""
iaa = True  # add iodoacetamide to all cysteine peptide mass automatically
deuterium_mass_defect = False  # Use the mass difference between 2H and H as isotopomer mass difference
labeled_residue = 'K'   # the labeled amino acid in aa labeling
label_mass = 6.0201     # this is the label mass in amino acid
# labeling; will be excluded from peptide mass calculation unlike other PTM masses

"""
Curve-fitting

"""
initial_k_deg = 0.3   # initial value for k_deg in scipy.optimize
max_iter = 1200       # max number of evaluations in scipy.optimize

"""
Match between runs

"""
soft_threshold_q = 0.25  # soft threshold for q-value to include peptides if found in other time points
min_fraction_mbr = 0.25  # minimum portion of samples in which a sample is identified before used for MAR

"""
Testing

"""
dummy = False   # for testing only, do not integrate but simply return 0
