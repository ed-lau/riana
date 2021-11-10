# -*- coding: utf-8 -*-

""" Main. """
import logging
import re
import os

import tqdm
import pandas as pd
import numpy as np
from scipy import optimize
from functools import partial
import matplotlib.pyplot as plt

from riana import accmass, constants, models, __version__


def strip_concat(sequence: str,
                 ) -> str:
    """
    cleans up concat sequences (peptide_charge) and remove modifications
    to return peptide string for labeling site calculations

    :param sequence:    concat sequence containing charge and modificaitons
    :return:
    """
    # 2021-05-18 strip all N-terminal n from Comet
    sequence = re.sub('^n', '', sequence)

    # Strip all modifications
    sequence = re.sub('\\[.*?]', '', sequence)

    # Strip the underscore and charge
    sequence = re.sub('_[0-9]+', '', sequence)

    return sequence


def calculate_a0(sequence: str,
                 ) -> float:
    """
    calculates the initial isotope enrichment of a peptide prior to heavy water labeling

    :param sequence:    str: concat sequences
    :return:            float: mi at time 0
    """

    sequence = strip_concat(sequence)

    res_atoms = accmass.count_atoms(sequence)

    a0 = np.product([np.power(constants.iso_abundances[i], res_atoms[i]) for i, v in enumerate(res_atoms)])

    return a0


def calculate_label_n(sequence: str,
                      ) -> float:

    sequence = strip_concat(sequence)

    return sum([constants.label_hydrogens.get(char) for char in sequence])


def calculate_fs(a, a_0, a_max):
    return (a-a_0)/(a_max-a_0)


def fit_all(args):
    """
    performs kinetic curve-fitting of integration output using a kinetics model
    returns a csv file containing the peptide concatamer, k_deg, sd, and r2

    :param args:
    :return:
    """

    #
    # parse arguments
    #

    riana_list = args.riana_path        # list of integration output
    model = args.model                  # kinetic model, e.g., 'simple'
    q_threshold = args.q_value          # peptide q value threshold
    t_threshold = args.depth            # minimal number of time points threshold
    ria_max = args.ria                  # final isotope enrichment level (e.g., 0.046)
    outdir = args.out                   # output directory

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #
    # logging
    #
    fit_log = logging.getLogger('riana.fit')
    fit_log.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(outdir, 'riana_fit.log'))
    fh.setLevel(logging.INFO)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add the handlers to the logger
    fit_log.addHandler(fh)
    fit_log.addHandler(ch)
    fit_log.info(args)
    fit_log.info(__version__)

    #
    # read the integration output files in
    #
    rdf = pd.DataFrame()
    for in_file in riana_list:
        rdf = rdf.append(pd.read_table(in_file))

    # filter by percolator q-value
    rdf_filtered = rdf[rdf['percolator q-value'] < q_threshold]

    # filter by number of time points
    concats = rdf_filtered.groupby('concat')['sample'].nunique()
    concats = concats[concats >= t_threshold]
    rdf_filtered = rdf_filtered[rdf_filtered.concat.isin(concats.index)]

    # create output dictionary
    out_dict = {}

    # loop through each qualifying peptide
    for seq in tqdm.tqdm(rdf_filtered.concat.unique()):

        # create subset dataframe with current concat
        y = rdf_filtered.loc[rdf_filtered['concat'] == seq].copy()
        fit_log.info(f'Fitting peptide {seq} with data shape {y.shape}')

        # get t, mi from subset dataframe
        y['mi'] = y['m0'] / (y['m0'] + y['m1'] + y['m2'] + y['m3'] + y['m4'] + y['m5'])
        fit_log.info(y[['sample', 'mi']])
        t = np.array([float(re.sub('[^0-9]', '', time)) for time in y['sample']])
        mi = np.array(y['mi'].tolist())

        # calculate a_0, a_max, and fractional synthesis
        stripped = strip_concat(seq)
        num_labeling_sites = calculate_label_n(seq)
        a_0 = calculate_a0(seq)
        a_max = a_0 * np.power((1 - ria_max), num_labeling_sites)
        fs = calculate_fs(a=mi, a_0=a_0, a_max=a_max)

        fit_log.info(f'concat: {stripped}, n: {num_labeling_sites}, a_0: {a_0}, a_max: {a_max}')
        fit_log.info([t, fs])

        # perform curve-fitting
        popt, pcov = optimize.curve_fit(f=partial(models.one_exponent, a_0=0., a_max=1.),
                                        xdata=t,
                                        ydata=fs,
                                        bounds=([0], [10]),
                                        )
        # calculate standard deviation from the non-linear least square cov
        sd = np.sqrt(np.diag(pcov))[0]

        # calculate residuals and r2
        residuals = fs - models.one_exponent(t, a_max=1., a_0=0., k_deg=popt[0])
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((fs - np.mean(fs)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        fit_log.info(f'Best fit k_deg: {popt[0]}, sd: {sd}, residuals: {residuals}, R2: {r_squared}')

        # create plot
        fig, ax = plt.subplots()
        plt.plot(t, fs, '.', label='Fractional synthesis')
        plt.plot(np.array(range(0, 31)),
                 models.one_exponent(t=np.array(range(0, int(np.max(t)))), k_deg=popt, a_0=0, a_max=1),
                 'r-', label=f'k_deg={np.round(popt[0], 3)}'
                 )

        plt.plot(np.array(range(0, 31)),
                 models.one_exponent(t=np.array(range(0, int(np.max(t)))), k_deg=popt[0] + sd, a_0=0, a_max=1),
                 'r--', label=f'Upper={np.round(popt[0]+sd, 3)}'
                 )

        plt.plot(np.array(range(0, 31)),
                 models.one_exponent(t=np.array(range(0, int(np.max(t)))),
                                     k_deg=popt[0] ** 2 / (popt[0] + sd),
                                     a_0=0,
                                     a_max=1,
                                     ),
                 'r--', label=f'Lower={np.round(popt[0] ** 2 / (popt[0] + sd), 3)}'
                 )
        plt.xlabel('t')
        plt.ylabel('fs')
        plt.title(f'Sequence: {seq} R**2: {np.round(r_squared, 3)}')
        plt.legend()
        plt.xlim([-1, 32])
        plt.ylim([-0.1, 1.1])
        # plt.show()

        plot_dir = os.path.join(outdir, 'curves')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        fig.savefig(os.path.join(plot_dir, f'{seq}.png'))
        plt.close(fig)

        out_dict[seq] = [popt[0], r_squared, sd]

    out_df = pd.DataFrame.from_dict(out_dict, orient='index', columns=['k_deg', 'R_squared', 'sd'])

    out_df.to_csv(os.path.join(outdir, 'riana_fit_peptides.tsv'))

    return True
