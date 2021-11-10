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
    :param a_max:  asymptotic isotope enrichment
    :return:
    """

    return a_0 + (a_max - a_0) * (1. - np.exp(-k_deg*t))

