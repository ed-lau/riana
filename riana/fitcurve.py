# -*- coding: utf-8 -*-

""" Main. """
import logging
import re
import os
import sys
from functools import partial

from riana.project import ReadDirectory
from riana.peptides import ReadPercolator
from riana.spectra import Mzml

from riana import integrate, params, __version__

import tqdm
import pandas as pd

import numpy as np
from riana import accmass, constants


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


def calculate_label_n(sequence:str,
                      ) -> float:

    sequence = strip_concat(sequence)

    return sum([constants.label_hydrogens.get(char) for char in sequence])


def calculate_fs(a, a_0, a_max):
    return (a-a_0)/(a_max-a_0)


def fit_all(args):
    """
    fits all concatamers

    :param args:
    :return:
    """

    riana_list = args.riana_path

    model = args.model # 'simple'
    q_threshold = args.q_value
    t_threshold = args.depth
    ria_max = 0.047
    outdir = 'out/snakemake'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # %%

    # open the riana files
    rdf = pd.DataFrame()
    for in_file in riana_list:
        rdf = rdf.append(pd.read_table(in_file))


    # filter by percolator q-value
    rdf_filtered = rdf[rdf['percolator q-value'] < q_threshold]

    # filter by number of time points
    concats = rdf_filtered.groupby('concat')['sample'].nunique()
    concats = concats[concats >= t_threshold]
    rdf_filtered = rdf_filtered[rdf_filtered.concat.isin(concats.index)]

    # rdf_filtered = rdf_filtered[rdf_filtered['concat'] in concats.index]
    print(rdf_filtered.sort_values(by='concat').concat)

    # output dictionary
    out_dict = {}

    # Loop through each qualifying peptide
    for seq in tqdm.tqdm(rdf_filtered.concat.unique()):

        y = rdf_filtered.loc[rdf_filtered['concat'] == seq].copy()

        # TODO: this is fit one

        print(seq)
        print(y.shape)

        if y.shape[0] < t_threshold:
            continue

        y['mi'] = y['m0'] / (y['m0'] + y['m1'] + y['m2'] + y['m3'] + y['m4'] + y['m5'])
        print(y[['sample', 'mi']])

        stripped = strip_concat(seq)
        num_labeling_sites = calculate_label_n(seq)
        a_0 = calculate_a0(seq)
        a_max = a_0 * np.power((1 - ria_max), num_labeling_sites)

        print(stripped, num_labeling_sites, a_0, a_max)
        continue
        # test = y[['sample', 'mi']]

        t = np.array([float(re.sub('[^0-9]', '', time)) for time in y['sample']])
        mi = np.array(y['mi'].tolist())

        fs = riana.fitcurve.calculate_fs(a=mi, a_0=a_0, a_max=a_max)

        print(fs)
        print(mi)

        popt, pcov = optimize.curve_fit(f=partial(riana.models.one_exponent, a_0=0., a_max=1.),
                                        xdata=t,
                                        ydata=fs,
                                        bounds=([0], [10]),
                                        )

        print(popt)
        print(popt[0])
        print(np.sqrt(np.diag(pcov)))

        residuals = fs - riana.models.one_exponent(t, a_max=1., a_0=0., k_deg=popt[0])
        ss_res = np.sum(residuals ** 2)

        ss_tot = np.sum((fs - np.mean(fs)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        sd = np.sqrt(np.diag(pcov))[0]
        print(residuals, r_squared)

        fig, ax = plt.subplots()
        plt.plot(t, fs, '.', label='Fractional synthesis')
        plt.plot(np.array(range(0, 31)),
                 riana.models.one_exponent(t=np.array(range(0, 31)), k_deg=popt, a_0=0, a_max=1),
                 'r-', label=f'k_deg={np.round(popt[0], 3)}'
                 )

        plt.plot(np.array(range(0, 31)),
                 riana.models.one_exponent(t=np.array(range(0, 31)), k_deg=popt[0] + sd, a_0=0, a_max=1),
                 'r--', label=f'sd={sd}'
                 )

        plt.plot(np.array(range(0, 31)),
                 riana.models.one_exponent(t=np.array(range(0, 31)), k_deg=popt[0] ** 2 / (popt[0] + sd), a_0=0,
                                           a_max=1),
                 'r--', label=f'sd={sd}'
                 )
        plt.xlabel('t')
        plt.ylabel('fs')
        plt.title(f'Sequence: {seq} R**2: {np.round(r_squared, 3)}')
        plt.legend()
        plt.xlim([-1, 32])
        plt.ylim([-0.1, 1.1])
        # plt.show()
        fig.savefig(f'out/{seq}.png')

        out_dict[seq] = [popt[0], r_squared, sd]

    # %%
    # out_df = pd.DataFrame.from_dict(out_dict, orient='index', columns=['k_deg', 'R_squared', 'sd'])

    # out_df.to_csv("out/out.csv")

    return True