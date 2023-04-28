# -*- coding: utf-8 -*-

""" kinetic models for curve fitting """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

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


def plot_model(protein,
               peptide,
               k_deg,
               r_squared,
               sd,
               t_series,
               fs_series,
               start_time: int = 0,
               end_time: int = 1,
               model_to_use: callable = one_exponent,
               **model_pars,

               ) -> Figure:
    """
    This function plots the fractional synthesis data and the fitted model.
    :param protein:
    :param peptide:
    :param k_deg:
    :param r_squared:
    :param sd:
    :param t_series:
    :param fs_series:
    :param start_time:
    :param end_time:
    :param model_to_use:
    :param model_pars:

    :return:
    """

    # create plot
    fig = Figure(figsize=(5, 3), dpi=100)

    # fig.suptitle('')


    # Produce the fitting trace
    fit_plot = fig.add_subplot(111)
    fit_plot.set_xlabel('Time')
    fit_plot.set_ylabel('Fractional synthesis')
    fit_plot.set_title(f'Protein: {protein} Sequence: {peptide} R2: {np.round(r_squared, 3)}')
    fit_plot.grid()

    fit_plot.plot(t_series, fs_series, '.', label='Fractional synthesis')

    # 2021-11-20 create clipped t,fs series for fs values out of [-0.2, 1.2] and display them as 'x'
    t_clipped = t_series[(fs_series < -0.2) | (fs_series > 1.2)]
    fs_clipped = fs_series[(fs_series < -0.2) | (fs_series > 1.2)]
    fs_clipped = fs_clipped.clip(min=-0.2, max=1.2)
    if len(t_clipped) > 0:
        plt.plot(t_clipped, fs_clipped, 'rx')

    fit_plot.plot(np.array(range(0, int(np.max(t_series)) + 1)),
             model_to_use(t=np.array(range(0, int(np.max(t_series)) + 1)),
                   k_deg=k_deg,
                   a_0=0.,
                   a_max=1.,
                   **model_pars,
                   ),
             'r-', label=f'k_deg={np.round(k_deg, 3)}'
             )

    fit_plot.plot(np.array(range(0, int(np.max(t_series)) + 1)),
             model_to_use(t=np.array(range(0, int(np.max(t_series)) + 1)),
                   k_deg=k_deg + sd,
                   a_0=0.,
                   a_max=1.,
                   **model_pars,
                   ),
             'r--', label=f'Upper={np.round(k_deg + sd, 3)}'
             )

    fit_plot.plot(np.array(range(0, int(np.max(t_series)) + 1)),
             model_to_use(t=np.array(range(0, int(np.max(t_series)) + 1)),
                   k_deg=k_deg ** 2 / (k_deg + sd),
                   a_0=0.,
                   a_max=1.,
                   **model_pars,
                   ),
             'r--', label=f'Lower={np.round(k_deg ** 2 / (k_deg + sd), 3)}'
             )

    fit_plot.legend()
    fit_plot.set_xlim([start_time-1, end_time+1])
    fit_plot.set_ylim([-0.2, 1.2])

    return fig