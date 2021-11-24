# -*- coding: utf-8 -*-

""" kinetic models for curve fitting """

import numpy as np


def one_exponent(t,
                 k_deg: float,
                 a_0: float = 0,
                 a_max: float = 1,
                 **_,
                 ) -> float:
    """
    Simple exponential model

    :param t:       time (x)
    :param k_deg:   degradation rate constant
    :param a_0:     initial isotope enrichment
    :param a_max:   asymptotic isotope enrichment
    :return:        float: a_t
    """

    return a_0 + (a_max - a_0) * (1. - np.exp(-k_deg*t))


def two_compartment_guan(t,
                         k_deg: float,
                         a_0: float = 0,
                         a_max: float = 1,
                         k_p: float = 0.5,
                         **_,
                         ) -> float:
    """
    Two-Compartment model with a lumped precursor rate constant k_p to describe precursor lag
    as described in Guan et al.

    :param t:       time (x)
    :param k_deg:   degradation rate constant
    :param a_0:     initial isotope enrichment
    :param a_max:   asymptotic isotope enrichment
    :param k_p:     precursor rate constant
    :return:        float: a_t
    """

    return a_0 + (a_max-a_0) * (1. - (np.exp(-t * k_deg) * k_p - np.exp(-t * k_p) * k_deg) / (k_p - k_deg))


def two_compartment_fornasiero(t,
                               k_deg: float,
                               a_0: float = 0.,
                               a_max: float = 1.,
                               k_p: float = 0.5,
                               k_r: float = 0.1,
                               r_p: float = 10.,
                               **_,
                               ) -> float:
    """
    Two-Compartment model with a precursor rate constant k_p (i.e., b in the Fornasiero et al. paper)
    to describe precursor lag and a global reutilization rate k_r from unlabeled protein degradation (i.e., a)
    as described in Fornasiero et al. The third parameter r describes the ratio of free vs. protein-bound precursor.

    :param t:       time (x)
    :param k_deg:   degradation rate constant
    :param a_0:     initial isotope enrichment
    :param a_max:   asymptotic isotope enrichment
    :param k_p:     precursor accumulation-breakdown constant
    :param k_r:     precursor reutilization rate constant (i.e., a)
    :param r_p:      free/bound precursor ratio
    :return:        float: a_t
    """

    a = k_r
    b = k_p
    r = r_p

    big_c = np.sqrt(-4 * a * b + (a + b + a * r)**2)
    k1 = (a + b + a * r + big_c)/2
    k2 = (a + b + a * r - big_c)/2
    big_a = -1 * (a - b + a * r - big_c)/(big_c * 2)

    tau1 = 1. / k1  # time constants are reciprocals of the rate constants above
    tau2 = 1. / k2

    # note: the below describes the precursor ria over time curve
    # y  = (1. - big_a * np.exp(-k1 * t) - (1. - big_a) * np.exp(-k2 * t))

    x = 1. - (big_a * tau1 * (1. - np.exp(k_deg * t - k1 * t))/(tau1 - (1. / k_deg))) - \
        ((1. - big_a) * tau2 * (1. - np.exp(k_deg * t - k2 * t)) / (tau2 - (1. / k_deg)))

    return a_0 + (a_max-a_0) * (1. - np.exp(-k_deg * t) * x)
