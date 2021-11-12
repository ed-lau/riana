# -*- coding: utf-8 -*-

""" kinetic models for curve fitting """

import numpy as np


def one_exponent(t: float,
                 k_deg: float,
                 a_0: float = 0,
                 a_max: float = 1,
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


def two_compartment_guan(t: float,
                         k_deg: float,
                         a_0: float = 0,
                         a_max: float = 1,
                         k_p: float = 0.5,
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