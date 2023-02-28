# -*- coding: utf-8 -*-

""" Functions to integrate isotopomer peaks """


import time

import logging
import re
import os
import sys
from functools import partial
import scipy.signal

import numpy as np
import pandas as pd
import tqdm

from .__init__ import __version__
from riana import constants
from riana.peptides import ReadPercolator
from riana.spectra import Mzml
from riana.logger import get_logger

def integrate_all(args) -> None:
    """
    Improved process to integrate for isotope abundance analysis.
    Idea is to loop through the mzML only once - get all the peptides to be integrated first
    Find all the spectrum ID - m/z combinations to be integrated, then integrate them sequentially

    :param args: Arguments from command line
    :return:

    """

    # ---- Get the logger ----

    logger = get_logger(__name__, args.out)
    logger.info(args)
    logger.info(__version__)

    # ---- Read in the Percolator file ----
    mzid = ReadPercolator(path=args.id_path,
                          sample=args.sample,
                          logger=logger,
                          )


    # 2022-10-31 RIANA will now attempt to read the percolator.log.txt file for fraction (file_idx) mzML assignment
    log_path = os.path.join(os.path.dirname(args.id_path.name),
                            'percolator.log.txt')

    # If the log file exists, use it to read the assignment
    if os.path.exists(log_path):
        logger.warning(f'Percolator log file exists at {log_path} '
                       f'and will be used for index assignment.')

        with open(log_path, 'r') as f:
            lines = f.readlines()

        mzml_files = {}
        for line in lines:
            percolator_indices = 'INFO: Assigning index ([0-9]*) to (.*)\.'
            pattern = re.findall(percolator_indices, line)
            if len(pattern) == 1:
                idx, pathname = pattern[0]
                dirname, filename = os.path.split(pathname)
                mzml_files[int(idx)] = re.sub('\.pep\.xml', '', filename)
                # TODO: will probably have to account for .pin or other input to Percolator

    # If the log file does not exist, assign index naively based on sort
    else:
        logger.warning(f'Percolator log file not found at {log_path}; '
                       f'mzml files will be sorted for index assignment.')
        # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
        # Note this will throw an error if not every fraction results in at least some ID, but we will ignore for now.

        # ---- Read in the mzML file ----

        mzml_file_list = [f for f in os.listdir(args.mzml_path) if re.match('^.*.mz[Mm][Ll]', f)]
        # Sort the mzML files by names
        # Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
        # But we will ignore for now
        mzml_file_list.sort()

        # Make dictionary of idx, filename from enumerate
        mzml_files = {}
        for idx, filename in enumerate(mzml_file_list):
            mzml_files[idx] = re.sub('.mz[Mm][Ll](\.gz)?', '', filename)

    # Print mzml files to log
    logger.info(f'mzml file orders: {mzml_files}')

    # Throw an error if there is no mzML file in the mzml directory
    assert len(mzml_files) != 0, '[error] no mzml files in the specified directory'
    # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
    # Note this will break if the last fraction does not contain at least some ID, but we will ignore for now.
    assert len(mzml_files) == max(mzid.indices) + 1, '[error] number of mzml files not matching id list'

    # create overall result files
    overall_integrated_df = pd.DataFrame()   # all integrated isotopomer peak areas
    overall_intensities_df = pd.DataFrame()  # all raw peak intensities

    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in mzid.indices:

        # Progress message
        logger.info('Doing mzml: {0} ({1} of {2})'.format(
            mzml_files[idx],
            str(idx + 1),
            str(len(mzid.indices))))

        # Make a subset dataframe with the current file index (fraction) being considered
        mzid.get_current_fraction_psms(idx)
        mzid.filter_current_fraction_psms(  # lysine_filter=0,
            peptide_q=args.q_value,
            unique_only=args.unique,
            use_soft_threshold=False,   # 20211109 True
            # match_across_runs=False,  # args.mbr
        )

        # --- Read the mzML file of the current fraction into a dictionary and create MS1/MS2 indices ----
        if os.path.exists(os.path.join(args.mzml_path, mzml_files[idx] + '.mzML')):
            mzml_path = os.path.join(args.mzml_path, mzml_files[idx] + '.mzML')
        elif os.path.exists(os.path.join(args.mzml_path, mzml_files[idx] + '.mzML.gz')):
            mzml_path = os.path.join(args.mzml_path, mzml_files[idx] + '.mzML.gz')
        else:
            raise FileNotFoundError(f'Could not find mzML file for index {idx} at {args.mzml_path}')
        try:
            mzml = Mzml(mzml_path)

        except OSError as e:
            sys.exit('[error] failed to load fraction mzml file. ' + str(e.errno))


        #
        # get peak intensity for each isotopomer in each spectrum ID in each peptide
        #
        # TODO: to accommodate multiple PSMs per concat, this should loop through a concat list.
        loop_ = range(len(mzid.curr_frac_filtered_id_df))
        assert len(loop_) > 0, 'No qualified peptide after filtering'

        get_isotopomer_intensity_partial = partial(get_isotopomer_intensity,
                                                   id_=mzid.curr_frac_filtered_id_df.copy(),
                                                   iso_to_do=args.iso,
                                                   mzml=mzml,
                                                   rt_tolerance=args.r_time,
                                                   mass_tolerance=args.mass_tol,
                                                   use_range=True,
                                                   mass_defect=args.mass_defect,
                                                   smoothing=args.smoothing,
                                                   )

        # Single threaded loop
        # '''
        # results = []
        # for i in loop_:
        #     print(i)
        #     results += integrate_one_partial(i)
        # '''

        # For parallelism, use concurrent.futures instead of multiprocessing for higher speed
        # '''
        from concurrent import futures
        with futures.ThreadPoolExecutor(max_workers=args.thread) as ex:
            intensities_results = list(tqdm.tqdm(ex.map(get_isotopomer_intensity_partial, loop_),
                                                 total=max(loop_),
                                                 desc=f'Extracting isotopomer intensities in sample: {args.sample}'
                                                      f' file: {mzml_files[idx]}'))
        # '''


        #
        # append the raw intensities data frame
        #
        intensities_df = pd.DataFrame().append([res for res in intensities_results], ignore_index=True)
        intensities_df['file'] = mzml_files[idx]
        intensities_df['idx'] = idx
        # id_intensities_df = intensities_dfpd.merge(intensities_df, mzid.curr_frac_filtered_id_df, on='pep_id', how='left')

        if len(overall_intensities_df.index) == 0:
            overall_intensities_df = intensities_df
        else:
            overall_intensities_df = overall_intensities_df.append(intensities_df, ignore_index=True)

        # Now perform integration
        integrated_peaks = []
        for int_res in tqdm.tqdm(intensities_results,
                                 desc=f'Integrating peaks in sample: {args.sample} file: {mzml_files[idx]}'
                                 ):
            iso_areas = [int_res['pep_id'][0]]


            for m in ['m' + str(iso) for iso in args.iso]:

                iso_area = np.trapz(int_res[m], x=int_res['rt'])
                iso_areas.append(iso_area)

            integrated_peaks.append(iso_areas)
        #
        # convert the integrated values into a data frame
        #
        df_columns = ['pep_id'] + ['m' + str(iso) for iso in args.iso]

        integrated_df = pd.DataFrame(integrated_peaks, columns=df_columns)
        id_integrated_df = pd.merge(mzid.curr_frac_filtered_id_df, integrated_df, on='pep_id', how='left')
        id_integrated_df['file'] = mzml_files[idx]

        # bind rows of the current result to the sample master
        if len(overall_integrated_df.index) == 0:
            overall_integrated_df = id_integrated_df
        else:
            overall_integrated_df = overall_integrated_df.append(id_integrated_df, ignore_index=True)

    # write the integrated and intensities results
    save_path = os.path.join(args.out, args.sample + '_riana.txt')
    overall_integrated_df.to_csv(path_or_buf=save_path, sep='\t')

    if args.write_intensities:
        save_path = os.path.join(args.out, args.sample + '_riana_intensities.txt')
        overall_intensities_df.to_csv(path_or_buf=save_path, sep='\t')

    tqdm.tqdm.write('Completed.')

    return None


def get_isotopomer_intensity(index: int,
                             id_: pd.DataFrame,
                             iso_to_do: list,
                             rt_tolerance: float,
                             mass_tolerance: int,
                             mzml,
                             use_range: True,
                             mass_defect: str = 'D',
                             smoothing: int = None,
                             ) -> list:
    """
    get all isotopomer mass intensities from ms1 scans within range for one peptide for integration

    :param index: int The row number of the peptide ID table passed from the wrapper.
    :param id_: protein identification dataframe
    :param iso_to_do: list of isotopomers
    :param rt_tolerance: retention time range in minutes
    :param mass_tolerance: relative mass tolerance (to be converted to ppm) e.g., 50
    :param mzml: mzml file object
    :param use_range: boolean whether to use all PSM scans for each peptide for RT determination
    :param mass_defect: whether to use deuterium or c13 for mass defect
    :param smoothing: int number of scans to smooth over
    :return: list of intensity over time [index, pep_id, m0, m1, m2, ...]

    """

    # determine the mass of protons and c13
    proton = constants.PROTON_MASS
    if mass_defect == 'D':
        iso_added_mass = constants.D_MASSDIFF  # see constants for details
    elif mass_defect == 'SILAC':
        iso_added_mass = constants.SILAC_MASSDIFF
    elif mass_defect == "C13":
        iso_added_mass = constants.C13_MASSDIFF
    else:
        raise ValueError(f'Invalid mass defect: {mass_defect}')

    # get peptide mass, scan number, and charge
    peptide_mass = float(id_.loc[index, 'peptide mass'])

    charge = float(id_.loc[index, 'charge'])

    # calculate precursor mass from peptide monoisotopic mass
    peptide_prec = (peptide_mass + (charge * proton)) / charge

    intensity_over_time = []

    # choose nearby scans from the span of scans of all qualifying psms that match to the peptide-z concat
    if use_range:
        concat_id_ = id_[id_['concat'] == id_.loc[index, 'concat']]
        min_scan = min(concat_id_['scan'])
        max_scan = max(concat_id_['scan'])

        peptide_rt_lower = mzml.rt_idx[np.searchsorted(mzml.scan_idx, min_scan, side='left') - 1]
        peptide_rt_upper = mzml.rt_idx[np.searchsorted(mzml.scan_idx, max_scan, side='left') - 1]

        nearby_ms1_scans = mzml.scan_idx[
            np.where((mzml.rt_idx - peptide_rt_lower > -rt_tolerance) & (mzml.rt_idx - peptide_rt_upper < rt_tolerance))]

    # If not use range, simply use the scan number from the psm as the center for +/- rt_tolerance
    else:
        scan_number = int(id_.loc[index, 'scan']) # - 1
        # get retention time from Percolator scan number # 2021-05-19 added -1
        peptide_rt = mzml.rt_idx[np.searchsorted(mzml.scan_idx, scan_number, side='left') - 1]
        assert isinstance(peptide_rt.item(), float), '[error] cannot retrieve retention time from scan number'

        nearby_ms1_scans = mzml.scan_idx[np.abs(mzml.rt_idx - peptide_rt) <= rt_tolerance]

    # loop through each MS1 spectrum within retention time range
    for scan in nearby_ms1_scans:

        try:
            spec = mzml.msdata[np.array(np.where(mzml.scan_idx == scan)).item()]

        except ValueError or IndexError:
            raise Exception('Scan does not correspond to MS1 data.')

        for iso in iso_to_do:

            # set upper and lower mass tolerance bounds
            prec_iso_am = peptide_prec + (iso * iso_added_mass / charge)

            mass_tolerance_ppm = float(mass_tolerance) * 1e-6

            delta_mass = prec_iso_am * mass_tolerance_ppm /2

            # sum intensities within range
            matching_int = np.sum(spec[np.abs(spec[:, 0] - prec_iso_am) <= delta_mass, 1])

            intensity_over_time.append([iso,
                                        mzml.rt_idx[mzml.scan_idx == scan].item(),
                                        matching_int,
            ]
            )


    if not intensity_over_time:
        raise Exception(f'No intensity profile for peptide {peptide_prec}'
                        f' min {min_scan} {peptide_rt_lower} max {max_scan} {peptide_rt_upper}'
                        f' scans {nearby_ms1_scans} iso {iso_to_do}')


    # integrated = [index] + [(id_.loc[index, 'pep_id'])] + integrate_isotope_intensity(np.array(intensity_over_time),
    #                                                                                   iso_to_do=iso_to_do)

    # 2021-11-10 return intensities -- currently this is a pandas df which may be slower...
    intensity_out = pd.DataFrame(intensity_over_time, columns=['iso', 'rt', 'int'])
    intensity_out = intensity_out.pivot(index='rt', columns='iso', values='int')
    intensity_out.columns = ['m' + str(iso) for iso in iso_to_do]
    intensity_out['rt'] = intensity_out.index
    intensity_out = intensity_out.reset_index(drop=True)
    intensity_out['ID'] = index
    intensity_out['pep_id'] = id_.loc[index, 'pep_id']
    intensity_out['concat'] = id_.loc[index, 'concat']


    # Smoothing with Savitzky-Golay filter
    if smoothing is not None:
        for m in ['m' + str(iso) for iso in iso_to_do]:
            if np.sum(intensity_out[m]) > 0:
                intensity_out[m] = scipy.signal.savgol_filter(x=intensity_out[m],
                                                              window_length=smoothing,
                                                              polyorder=1,
                                                              mode='nearest')

    return intensity_out # [integrated, intensity_out]


def integrate_isotope_intensity(intensity_over_time: np.ndarray,
                                iso_to_do: list,
                                ) -> list:
    """
    Given a list of isotopomer intensity over time, give the integrated intensity of each isotopomer

    :param intensity_over_time: Numpy array of isotopomer, time, intensities
    :param iso_to_do: list of isotopomers
    :return: list of itegrated intensity of each isotopomer
    """

    # integrate the individual isotopomers
    iso_intensity = []

    for j in iso_to_do:

        isotopomer_profile = intensity_over_time[intensity_over_time[:, 0] == j]

        if isotopomer_profile.size > 0:

            # use np.trapz rather than scipy.integrate to integrate
            iso_area = np.trapz(isotopomer_profile[:, 2], x=isotopomer_profile[:, 1])

        # if there is no isotopomer profile, set area to 0
        else:
            iso_area = 0

        iso_intensity.append(iso_area)

    if not iso_intensity:
        raise Exception(f'No positive numerical value integrated for isotopmer {intensity_over_time}')

    return iso_intensity
