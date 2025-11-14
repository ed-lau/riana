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
from riana.utils import strip_concat
from riana.fsynthesis import calculate_label_n, calculate_fs_m0, calculate_fs_fine_structure


def fit_all(args) -> None:
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
    rdf_filtered = rdf[rdf['percolator q-value'] < q_threshold].copy()
    # print(rdf_filtered)

    # filter by number of time points
    # concats = rdf_filtered.groupby('concat')['sample'].nunique()
    # concats = concats[concats >= t_threshold]
    # rdf_filtered = rdf_filtered[rdf_filtered.concat.isin(concats.index)]

    # TODO: 2021-12-17 we should filter the concat dataframe here to take only time series from fractions
    #   that contain at least X number of data points.
    rdf_filtered = rdf_filtered.groupby(['concat', 'file_idx']).filter(lambda x: x['sample'].nunique() >= t_threshold).copy()

    # TODO: catch when no peptide reaches the time or q threshold
    if rdf_filtered.shape[0] == 0:
        raise ValueError('No qualifying peptides for fitting. Check q-value or depth threshold.')

    # get the maximum time in this data set
    max_time = np.max([float(re.sub('[^0-9]', '', time)) for time in rdf_filtered['sample']])

    # create loop of concat indices for concurrent.futures
    loop_ = range(len(rdf_filtered.concat.unique())) # range(200) #
    fit_one_partial = partial(fit_one,
                              concat_list=rdf_filtered.concat.unique(),
                              label=label_,
                              filtered_integrated_df=rdf_filtered.copy(),
                              ria_max=ria_max,
                              model_=model,
                              aa_res=aa_res,
                              fs_fine_structure=args.fs,
                              model_pars=model_pars,
                              )

    # Single threaded loop
    # '''
    if hasattr(args, "gui"):
        logger.info('Running from GUI, using single thread.')
        out_dict = {}
        for i in tqdm.tqdm(loop_, miniters=len(loop_)/100):
            out_dict = out_dict | fit_one_partial(i)

    # '''

    # parallel loops using concurrent.futures
    # '''
    else:
        logger.info(f'Running from command line, using {num_threads} threads.')
        from concurrent import futures

        with futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            results = list(tqdm.tqdm(ex.map(fit_one_partial, loop_),
                                     total=max(loop_),
                                     desc=f'Fitting peptides to model - {args.model}:'))

        # Collect all the dicts, plot out graphs and output final file
        out_dict = {}

        for res in tqdm.tqdm(results, desc=f'Combining results'):
            out_dict = out_dict | res
            #Note this works in python 3.9, but not in 3.8
        # '''

    # Plot the individual curves
    if args.plotcurves:
        for seq, value in tqdm.tqdm(out_dict.items(), desc=f'Plotting curves'):

            #seq = list(res.keys())[0]
            k_deg, r_squared, sd, t, fs = value # list(res.values())[0]

            # do not plot and record result if k_deg is nan
            if np.isnan(k_deg):
                # print(k_deg, r_squared, sd, t, fs)
                continue

            # get first protein name
            protein = rdf_filtered[rdf_filtered.concat == seq]['protein id'].iloc[0]
            first_protein = re.sub('(sp|tr)\|.+?\|(.+?)_(MOUSE|HUMAN|RAT).*', '\\2', protein)
            if protein.count(',') > 0:
                first_protein = first_protein + ' (multi)'

            fig = models.plot_model(protein=first_protein,
                                    peptide=seq,
                                    k_deg=k_deg,
                                    r_squared=r_squared,
                                    sd=sd,
                                    t_series=t,
                                    fs_series=fs,
                                    start_time=0,
                                    end_time=max_time,
                                    model=model,
                                    model_pars=model_pars)

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

    logger.info("Loop finished")

    out_df = pd.DataFrame.from_dict(out_dict, orient='index', columns=['k_deg', 'R_squared', 'sd', 't', 'fs'])

    # Turn the list columns (t, fs) into aggregates separated by commas, first by expanding the lists
    out_df_2 = out_df.drop(columns=['k_deg', 'R_squared', 'sd']).explode(['t', 'fs']).copy()
    out_df_2.index.name = 'concat'

    # Then group by concat and aggregate the lists into a single string
    out_df_2['t'] = pd.to_numeric(out_df_2['t'], downcast='float')
    out_df_2['fs'] = pd.to_numeric(out_df_2['fs'], downcast='float')

    out_df_2 = out_df_2.round(3).copy()
    out_df_2 = out_df_2.groupby('concat').agg(pd.Series.tolist).copy()
    # Merge back to the data frame with k_deg, R_squared, sd
    out_df_2 = out_df_2.merge(out_df[['k_deg', 'R_squared', 'sd']], left_index=True, right_index=True, how='left').copy()

    # Merge with the original dataframe to get the protein id
    right_merge = rdf_filtered[['concat', 'protein id']].copy().drop_duplicates(subset=['concat']).set_index('concat')
    out_df_2 = out_df_2.merge(right_merge, left_index=True, right_index=True, how='left').copy()

    # Add back the kinetic parameters
    # TODO: save the kinetic parameters in a separate file
    out_df_2 = out_df_2.assign(kp=args.kp, kr=args.kr, rp=args.rp, ria_max=args.ria)

    out_df_2.to_csv(os.path.join(outdir, 'riana_fit_peptides.txt'), sep='\t')

    num_peps_fitted = out_df_2[out_df_2['R_squared'] >= 0.9].shape[0]
    logger.info(f'There are {num_peps_fitted} concats with R2 ≥ 0.9')

    # Remove logger
    logger.handlers.clear()

    return None # sys.exit(os.EX_OK)


def fit_one(loop_index,
            concat_list: list,
            label: int,
            filtered_integrated_df: pd.DataFrame,
            ria_max: float,
            model_: callable,
            aa_res: str,
            fs_fine_structure: str,
            model_pars: dict,
            ) -> dict:
    """

    :param loop_index:                  index of concat_list to fit
    :param concat_list:                 list of concatamer
    :param label:       int: 1=2H_in_vivo, 2=2H_in_vitro, 3=18O, 4=AA, if AA, return 1 assuming no heavy prior to labeling
    :param filtered_integrated_df:      dataframe of integrated data
    :param ria_max:                     maximum ria
    :param model_:                      model function
    :param aa_res:                      which amino acid residue is labeled, e.g., K
    :param fs_fine_structure:           whether to use fine structure for fs
    :param model_pars:                  dictionary of model parameters
    :return:
    """


    # create subset dataframe with current concat, taking only unique sample, file_idx pairs
    seq = concat_list[loop_index]
    y = filtered_integrated_df.loc[filtered_integrated_df['concat'] == seq].drop_duplicates(subset=['sample', 'file_idx']).copy()
    # print(seq)
    # print(y)

    logger = logging.getLogger('riana_fit')
    logger.info(f'Fitting peptide {seq} with data shape {y.shape}')
    logger.info('args.fine_structure: ' + str(fs_fine_structure))

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

    t = np.array([float(re.sub('[^0-9.]', '', time)) for time in y['sample']])


    # calculate a_0, a_max, and fractional synthesis
    stripped = strip_concat(seq)
    num_labeling_sites = calculate_label_n(seq,
                                           label=label,
                                           aa_res=aa_res)

    # if label is hw or hw_cell, use the heavy water labeling dictionary
    if (label == 1 or label == 2 or label == 3) and fs_fine_structure is not None:

        logger.info('Using fine structure calculation for FS')

        if fs_fine_structure == "m0_m1":
            y = y.assign(mi=(y.iso0 / y.iso1).where(y.iso0 != 0, 0))

        elif fs_fine_structure == "m0_m2":
            y = y.assign(mi=(y.iso0 / y.iso2).where(y.iso0 != 0, 0))

        elif fs_fine_structure == "m0_m3":
            y = y.assign(mi=(y.iso0 / y.iso3).where(y.iso0 != 0, 0))

        elif fs_fine_structure == "m1_m2":
            y = y.assign(mi=(y.iso1 / y.iso2).where(y.iso1 != 0, 0))

        elif fs_fine_structure == "m1_m3":
            y = y.assign(mi=(y.iso1 / y.iso3).where(y.iso1 != 0, 0))

        elif fs_fine_structure == "m0_mA":
            y = y.assign(mi=(y.iso0 / (y.iso0 + y.iso1 + y.iso2 + y.iso3 + y.iso4 + y.iso5)).where(y.iso0 != 0, 0))

        elif fs_fine_structure == "m1_mA":
            y = y.assign(mi=(y.iso1 / (y.iso0 + y.iso1 + y.iso2 + y.iso3 + y.iso4 + y.iso5)).where(y.iso1 != 0, 0))

        elif fs_fine_structure == "Auto":
            if num_labeling_sites < 15:
                y = y.assign(mi=(y.iso0 / y.iso1).where(y.iso0 != 0, 0))
            elif num_labeling_sites > 35:
                y = y.assign(mi=(y.iso1 / y.iso3).where(y.iso1 != 0, 0))
            else:
                y = y.assign(mi=(y.iso0 / y.iso2).where(y.iso0 != 0, 0))

        mi = np.array(y['mi'].tolist())

        # Use the fine structure calculator with an artificial element to calculate FS in HW labeling
        fs = calculate_fs_fine_structure(a=mi,
                                         seq=seq,
                                         label=label,
                                         ria_max=ria_max,
                                         formula=fs_fine_structure,
                                         )

    else:
        # Use m0/mA to calculate FS if label is aa or if fine structure is not used
        y['colsums'] = y.loc[:, y.columns.str.match('^iso[0-9]+$')].sum(axis=1)  # sums all m* columns
        y = y.assign(mi=(y.iso0 / y.colsums).where(y.iso0 != 0, 0))  # avoid division by 0, if m0 is 0, return 0
        mi = np.array(y['mi'].tolist())
        logger.info("Using fallback m0/mA calculation for FS")
        # logger.info("label:", label)
        # logger.info("fs_fine_structure:", fs_fine_structure)

        fs = calculate_fs_m0(a=mi,
                             seq=seq,
                             label=label,
                             ria_max=ria_max,
                             num_labeling_sites=num_labeling_sites)

    logger.info(f'concat: {stripped}, n: {num_labeling_sites}')
    logger.info([t, fs])

    # TODO: 2021-11-21 remove t/fs data points where fs is nan for any reason
    null_result = {seq: [np.nan, np.nan, np.nan, t, fs]}

    # if there is no labeling site, skip curve fitting because there is no isotopic labeling
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
