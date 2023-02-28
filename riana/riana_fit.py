# -*- coding: utf-8 -*-

""" Main. """
import logging
import re
import os
import sys
import gc

import tqdm
import pandas as pd
import numpy as np
from scipy import optimize
from functools import partial
import matplotlib.pyplot as plt

from riana import accmass, constants, models, params, __version__
from riana.logger import get_logger



def strip_concat(sequence: str,
                 ) -> str:
    """
    Cleans up concat sequences (peptide_charge) and remove modifications
    to return peptide string for labeling site calculations

    :param sequence:    concat sequence containing charge and modificaitons
    :return:
    """
    # 2021-05-18 strip all N-terminal n from Comet
    sequence = re.sub('^n', '', sequence)

    # Strip all modifications
    sequence = re.sub('\\[.*?\\]', '', sequence)

    # Strip the underscore and charge
    sequence = re.sub('_[0-9]+', '', sequence)

    return sequence


def calculate_a0(sequence: str,
                 label: str,
                 ) -> float:
    """
    Calculates the initial isotope enrichment of a peptide prior to heavy water labeling

    :param sequence:    str: concat sequences
    :param label:       str: aa, hw, or o18, if aa, return 1 assuming no heavy prior to labeling
    :return:            float: mi at time 0
    """

    if label == 'aa':
        return 1

    else:
        sequence = strip_concat(sequence)
        res_atoms = accmass.count_atoms(sequence)
        a0 = np.product([np.power(constants.iso_abundances[i], res_atoms[i]) for i, v in enumerate(res_atoms)])
        # TODO: this should calculate the full isotopic distribution
        return a0


def calculate_label_n(sequence: str,
                      label: str,
                      aa_res: str = 'K',
                      ) -> float:
    """
    Calculates labeling sites of the peptide sequence in heavy water
    or amino acid labeling

    :param sequence:    the peptide sequence
    :param label:       aa, hw, or o18; if aa, only return the labelable residues
    :param aa_res:       the amino acid being labeled
    :return:
    """

    # strip modification site and charge from concat sequence
    sequence = strip_concat(sequence)

    # if amino acid labeling, return number of labeled residues
    if label == 'aa':
        # return the sum of each residue in the aa_res
        return sum([sequence.count(i) for i in aa_res])

    # if d2o, return the number of labeling site in heavy water labeling
    elif label == 'hw':
        return sum([constants.label_hydrogens.get(char) for char in sequence])

    # else if o18, return the number of labeling sites for o18
    else:
        return sum([constants.label_oxygens.get(char) for char in sequence]) - 1.


def calculate_fs(a, a_0, a_max):
    """
    Calculates fractional synthesis based on a_t, a_0 (initial), and a_max (asymptote)

    :param a:       mi at a particular time
    :param a_0:     initial mi value before label onset
    :param a_max:   final mi value at plateau based on labeling site and precursor RIA
    :return:
    """

    # catch errors from no ria or no labeling site
    if a_max - a_0 == 0:
        # repeat an array of 0 if the input is an ndarray, otherwise return 0
        return np.repeat(0, len(a)) if isinstance(a, np.ndarray) else 0
    else:
        return (a-a_0)/(a_max-a_0)


def fit_all(args):
    """
    Performs kinetic curve-fitting of integration output using a kinetics model
    returns a csv file containing the peptide concatamer, k_deg, sd, and r2

    :param args:
    :return:
    """

    #
    # parse arguments
    #

    # get the logger
    logger = get_logger(__name__, args.out)

    riana_list = args.riana_path        # list of integration output
    model_pars = {'k_p': args.kp,
                  'k_r': args.kr,
                  'r_p': args.rp}
    q_threshold = args.q_value          # peptide q value threshold
    t_threshold = args.depth            # minimal number of time points threshold
    ria_max = args.ria                  # final isotope enrichment level (e.g., 0.046)
    outdir = args.out                   # output directory
    label_ = args.label
    aa_res = args.aa

    # select model
    if args.model == 'simple':
        model = models.one_exponent

    elif args.model == 'guan':
        model = models.two_compartment_guan

    elif args.model == 'fornasiero':
        model = models.two_compartment_fornasiero

    else:
        raise Exception('Unknown kinetics model.')

    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    if args.thread:
        try:
            num_threads = min(os.cpu_count() * 4, int(args.thread))

        except ValueError or TypeError:
            num_threads = 1  # os.cpu_count() * 4

    else:
        num_threads = 1  # os.cpu_count() * 4


    # Get logging
    logger.info(args)
    logger.info(__version__)

    #
    # read the integration output files in
    #
    rdf = pd.concat(pd.read_table(in_file) for in_file in riana_list)

    # filter by percolator q-value
    rdf_filtered = rdf[rdf['percolator q-value'] < q_threshold]

    # filter by number of time points
    # concats = rdf_filtered.groupby('concat')['sample'].nunique()
    # concats = concats[concats >= t_threshold]
    # rdf_filtered = rdf_filtered[rdf_filtered.concat.isin(concats.index)]

    # TODO: 2021-12-17 we should filter the concat dataframe here to take only time series from fractions
    #   that contain at least X number of data points.
    rdf_filtered = rdf_filtered.groupby(['concat', 'file_idx']).filter(lambda x: x['sample'].nunique() >= t_threshold)

    # TODO: catch when no peptide reaches the time or q threshold
    if rdf_filtered.shape[0] == 0:
        raise ValueError('No qualifying peptides for fitting. Check q-value or depth threshold.')

    # get the maximum time in this data set
    max_time = np.max([float(re.sub('[^0-9]', '', time)) for time in rdf_filtered['sample']])

    # create loop of concat indices for concurrent.futures
    loop_ = range(len(rdf_filtered.concat.unique()))
    fit_one_partial = partial(fit_one,
                              concat_list=rdf_filtered.concat.unique(),
                              label=label_,
                              filtered_integrated_df=rdf_filtered.copy(),
                              ria_max=ria_max,
                              model_=model,
                              aa_res=aa_res,
                              model_pars=model_pars,
                              )

    # parallel loops using concurrent.futures
    from concurrent import futures
    with futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
        results = list(tqdm.tqdm(ex.map(fit_one_partial, loop_),
                                 total=max(loop_),
                                 desc=f'Fitting peptides to model - {args.model}:'))

    # collect all the dicts, plot out graphs and output final file
    out_dict = {}

    for res in tqdm.tqdm(results, desc=f'Combining results'):
        out_dict = out_dict | res
        # Note this works in python 3.9, but not in 3.8

    if args.plotcurves:
        for res in tqdm.tqdm(results, desc=f'Plotting curves'):

            seq = list(res.keys())[0]
            k_deg, r_squared, sd, t, fs = list(res.values())[0]

            # do not plot and record result if k_deg is nan
            if np.isnan(k_deg):
                # print(k_deg, r_squared, sd, t, fs)
                continue

            # get first protein name
            protein = rdf_filtered[rdf_filtered.concat == seq]['protein id'].iloc[0]
            first_protein = re.sub('(sp|tr)\|.+?\|(.+?)_(MOUSE|HUMAN|RAT).*', '\\2', protein)
            if protein.count(',') > 0:
                first_protein = first_protein + ' (multi)'

            # create plot
            fig, ax = plt.subplots()
            plt.plot(t, fs, '.', label='Fractional synthesis')

            # 2021-11-20 create clipped t,fs series for fs values out of [-0.2, 1.2] and display them as 'x'
            t_clipped = t[(fs < -0.2) | (fs > 1.2)]
            fs_clipped = fs[(fs < -0.2) | (fs > 1.2)]
            fs_clipped = fs_clipped.clip(min=-0.2, max=1.2)
            if len(t_clipped) > 0:
                plt.plot(t_clipped, fs_clipped, 'rx')

            plt.plot(np.array(range(0, int(np.max(t))+1)),
                     model(t=np.array(range(0, int(np.max(t))+1)),
                           k_deg=k_deg,
                           a_0=0.,
                           a_max=1.,
                           **model_pars,
                           ),
                     'r-', label=f'k_deg={np.round(k_deg, 3)}'
                     )

            plt.plot(np.array(range(0, int(np.max(t))+1)),
                     model(t=np.array(range(0, int(np.max(t))+1)),
                           k_deg=k_deg + sd,
                           a_0=0.,
                           a_max=1.,
                           **model_pars,
                           ),
                     'r--', label=f'Upper={np.round(k_deg + sd, 3)}'
                     )

            plt.plot(np.array(range(0, int(np.max(t))+1)),
                     model(t=np.array(range(0, int(np.max(t))+1)),
                           k_deg=k_deg ** 2 / (k_deg + sd),
                           a_0=0.,
                           a_max=1.,
                           **model_pars,
                           ),
                     'r--', label=f'Lower={np.round(k_deg ** 2 / (k_deg + sd), 3)}'
                     )
            plt.xlabel('t')
            plt.ylabel('fs')
            plt.title(f'Protein: {first_protein} Sequence: {seq} R2: {np.round(r_squared, 3)}')
            plt.legend()
            plt.xlim([-1, max_time+1])
            plt.ylim([-0.2, 1.2])

            if k_deg < 0.01:
                speed = "slow"
            elif k_deg > 1:
                speed = "fast"
            else:
                speed = "mid"

            if r_squared >= 0.5:
                quality = "fit"
            else:
                quality = "poor"

            plot_dir = os.path.join(outdir, f'curves_{speed}_{quality}')

            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)

            fig.savefig(os.path.join(plot_dir, f'{first_protein}_{seq}.png'))
            plt.clf()
            plt.close(fig)
            # gc.collect()


    out_df = pd.DataFrame.from_dict(out_dict, orient='index', columns=['k_deg', 'R_squared', 'sd', 't', 'fs'])

    # add back the kinetic parameters
    out_df = out_df.assign(kp=args.kp, kr=args.kr, rp=args.rp, ria_max=args.ria)

    out_df.to_csv(os.path.join(outdir, 'riana_fit_peptides.csv'))

    num_peps_fitted = out_df[out_df['R_squared'] >= 0.9].shape[0]
    logger.info(f'There are {num_peps_fitted} concats with R2 ≥ 0.9')

    return sys.exit(os.EX_OK)


def fit_one(loop_index,
            concat_list: list,
            label: str,
            filtered_integrated_df: pd.DataFrame,
            ria_max: float,
            model_: callable,
            aa_res: str,
            model_pars: dict,
            ):
    """

    :param loop_index:
    :param concat_list:                 list of concatamer
    :param label:                       label type, determines how mi is calculated. currently must be aa, hw, or o18
    :param filtered_integrated_df:
    :param ria_max:
    :param model_:
    :param aa_res:                      which amino acid residue is labeled, e.g., K
    :param model_pars:
    :return:
    """


    # create subset dataframe with current concat
    seq = concat_list[loop_index]
    y = filtered_integrated_df.loc[filtered_integrated_df['concat'] == seq].copy()

    logger = logging.getLogger('riana_fit')
    logger.info(f'Fitting peptide {seq} with data shape {y.shape}')

    '''
    2021-11-22 Get t, mi from subset dataframe
    The current implementation is to avoid setting explicit mi formulae,
    rather it will return m0 over the sum of all m columns that were integrated 
    
    # if label == 'aa':
    #     y['mi'] = y['m0'] / (y['m0'] + y['m6'])
    # else:
    #     y['mi'] = y['m0'] / (y['m0'] + y['m1'] + y['m2'] + y['m3'] + y['m4'] + y['m5'])
    '''
    # TODO: note that this will create a lot of 0.0 fs data points where there was no m0 intensity
    #   which could affect integration results. Ideally we should remove these data points
    #   unless time is 0 and then filter for data series length again according to the data depth
    #   filter.

    y['colsums'] = y.loc[:, y.columns.str.match('^m[0-9]+$')].sum(axis=1)   # sums all m* columns
    y = y.assign(mi=(y.m0 / y.colsums).where(y.m0 != 0, 0))  # avoid division by 0, if m0 is 0, return 0

    logger.info(y[['sample', 'mi']])
    t = np.array([float(re.sub('[^0-9.]', '', time)) for time in y['sample']])
    mi = np.array(y['mi'].tolist())

    # calculate a_0, a_max, and fractional synthesis
    stripped = strip_concat(seq)
    num_labeling_sites = calculate_label_n(seq,
                                           label=label,
                                           aa_res=aa_res)

    a_0 = calculate_a0(seq, label=label)
    a_max = a_0 * np.power((1 - ria_max), num_labeling_sites)
    fs = calculate_fs(a=mi, a_0=a_0, a_max=a_max)

    logger.info(f'concat: {stripped}, n: {num_labeling_sites}, a_0: {a_0}, a_max: {a_max}')
    logger.info([t, fs])

    # TODO: 2021-11-21 remove t/fs data poionts where fs is nan for any reason
    null_result = {seq: [np.nan, np.nan, np.nan, t, fs]}

    # if there is no labeling site, skip
    if num_labeling_sites == 0:
        return null_result

    # perform curve-fitting
    try:
        popt, pcov = optimize.curve_fit(f=partial(model_, a_0=0., a_max=1., **model_pars),
                                        xdata=t,
                                        ydata=fs,
                                        bounds=([1e-4], [10]),
                                        p0=[params.initial_k_deg],
                                        maxfev=params.max_iter,
                                        )
    # catch error when not converging
    except RuntimeError:
        print(RuntimeError)
        print(stripped)
        print(t, fs)
        return null_result

    except ValueError:
        print(ValueError)
        print(stripped)
        print(t, fs)
        return null_result

    # calculate standard deviation from the non-linear least square cov
    sd = np.sqrt(np.diag(pcov))[0]

    # calculate residuals and r2
    residuals = fs - model_(t, a_max=1., a_0=0., k_deg=popt[0], **model_pars)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((fs - np.mean(fs)) ** 2)

    r_squared = np.nan if ss_tot == 0 else 1. - (ss_res / ss_tot)

    logger.info(f'Best fit k_deg: {popt[0]}, sd: {sd}, residuals: {residuals}, R2: {r_squared}')

    return {seq: [popt[0], r_squared, sd, t, fs]}
