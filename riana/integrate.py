# -*- coding: utf-8 -*-

""" Functions to integrate isotopomer peaks """


import time

import logging
import re
import os
import sys
from functools import partial

import numpy as np
import pandas as pd
import tqdm

from riana import constants, params, __version__
from riana.project import ReadDirectory
from riana.peptides import ReadPercolator
from riana.spectra import Mzml


def integrate_all(args):
    """
    Improved process to integrate for isotope abundance analysis.
    Idea is to loop through the mzML only once - get all the peptides to be integrated first
    Find all the spectrum ID - m/z combinations to be integrated, then integrate them sequentially

    :param args: Arguments from command line
    :return:

    """

    #
    # Handle command line arguments
    #

    # Use the sample name if supplied, otherwise use the basename of the mzml
    if args.sample is None:
        current_sample = os.path.basename(os.path.normpath(args.mzml_path))
    else:
        current_sample = args.sample

    # File out
    path_to_write = os.path.join(args.out)
    directory_to_write = os.path.dirname(path_to_write)
    os.makedirs(directory_to_write, exist_ok=True)

    # Logging
    main_log = logging.getLogger('riana')
    main_log.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(directory_to_write, f'riana_integrate_{current_sample}.log'))
    fh.setLevel(logging.INFO)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add the handlers to the logger
    main_log.addHandler(fh)
    main_log.addHandler(ch)

    main_log.info(args)
    main_log.info(__version__)

    #
    # Integrates only peptides matched to one protein
    #
    unique_pep = args.unique

    #
    # Convert the to-do isotopomer list option into a list of integers
    #
    iso_to_do = []

    for char in args.iso.split(','):

        try:
            char = int(char)

            # Only doing down to mz + 15.
            if char <= 15:
                iso_to_do.append(int(char))

        except ValueError or TypeError:
            pass

    if not iso_to_do:
        sys.exit('Error: Invalid isotopomer list given.')

    iso_to_do = list(set(iso_to_do))
    iso_to_do.sort()

    #
    # Percolator q value threshold for peptides and proteins
    #
    if args.q_value:

        try:
            q_threshold = float(args.q_value)

        except ValueError or TypeError:
            main_log.warning('Invalid Q value given - using default value.')
            q_threshold = float(1e-2)

    #
    # Retention time cutoff peptides and proteins
    #
    if args.r_time:

        try:
            rt_tolerance = float(args.r_time)

        except ValueError or TypeError:
            main_log.warning('Invalid retention time tolerance given - using default value.')
            rt_tolerance = float(1.0)
    else:
        rt_tolerance = float(1.0)

    #
    # MS1 mass tolerance for integration
    #
    if args.mass_tol:
        try:
            mass_tolerance = float(args.mass_tol) * 1e-6

        except ValueError or TypeError:
            main_log.warning('Invalid mass tolerance given - using default value.')
            mass_tolerance = float(100) * 1e-6
    else:
        mass_tolerance = float(100) * 1e-6

    #
    # Parallelism
    #
    if args.thread:
        try:
            num_threads = max(os.cpu_count() * 4, int(args.thread))

        except ValueError or TypeError:
            main_log.warning('Invalid thread count. Using default.')
            num_threads = 1  # os.cpu_count() * 4

    else:
        num_threads = 1  # os.cpu_count() * 4

    #
    # Read files in
    #

    # This is the directory that holds the entire project
    #20211109 project = ReadDirectory(dir_loc)

    # Get the master peptide ID list
    mzid = ReadPercolator(path=args.id_path,  # project=project,
                          sample=current_sample,
                          directory_to_write=directory_to_write,
                          #percolator_subdirectory=args.percolator
                          )

    mzid.read_psms()

    # TODO: should remove mbr from the main integration function and move to a different script
    # 20211109 mzid.make_master_match_list(  # lysine_filter=0,
    #   peptide_q=q_threshold,
    #    unique_only=unique_pep,
    #    min_fraction=params.min_fraction_mbr)

    # Each subdirectory is a sample
    # 20211109 samples = project.samples
    # Create the grand total out file
    # 20211109 master_df = pd.DataFrame()

    # 20211109 for current_sample in tqdm.tqdm(samples, desc='Processing Sample', total=len(samples)):

    sample_loc = os.path.normpath(args.mzml_path)   # os.path.join(project.path, current_sample, 'mzml')
    assert os.path.isdir(sample_loc), '[error] mzml path not a valid directory'

    # These are not necessary since there is only one sample/percolator per run now
    mzid.get_current_sample_psms(current_sample=current_sample)
    mzid.get_current_sample_mzid_indices()

    #2021-11-04 account for mzml and mzML
    mzml_files = [f for f in os.listdir(sample_loc) if re.match('^.*.mz[Mm][Ll]', f)]

    # Sort the mzML files by names
    # Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
    # TODO: check whether Percolator uses nat sort and implement if so
    mzml_files.sort()

    # Throw an error if there is no mzML file in the mzml directory
    assert len(mzml_files) != 0, '[error] no mzml files in the specified directory'
    # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
    # Note this will break if the last fraction does not contain at least some ID, but we will ignore for now.
    assert len(mzml_files) == max(mzid.indices) + 1, '[error] number of mzml files not matching id list'

    # create overall result files
    overall_integrated_df = pd.DataFrame()   # all integrated isotopomer peak areas
    overall_intensities_df = pd.DataFrame()  # all raw peak intensities

    #

    #
    # Read the mzml files and do integration
    #

    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in mzid.indices:

        # Progress message
        main_log.info('Doing mzml: {0} ({1} of {2})'.format(
            mzml_files[idx],
            str(idx + 1),
            str(len(mzid.indices))))

        # Make a subset dataframe with the current file index (fraction) being considered
        mzid.get_current_fraction_psms(idx)
        mzid.filter_current_fraction_psms(  # lysine_filter=0,
            peptide_q=q_threshold,
            unique_only=unique_pep,
            use_soft_threshold=False,   # 20211109 True
            # match_across_runs=False,  # args.mbr
        )

        try:
            mzml = Mzml(os.path.join(sample_loc, mzml_files[idx]))

        except OSError as e:
            sys.exit('[error] failed to load fraction mzml file. ' + str(e.errno))

        #
        # read the spectra into dictionary and also create MS1/MS2 indices
        #
        mzml.parse_mzml()

        #
        # get peak intensity for each isotopomer in each spectrum ID in each peptide
        #
        loop_ = range(len(mzid.curr_frac_filtered_id_df))

        integrate_one_partial = partial(integrate_one,
                                        id_=mzid.curr_frac_filtered_id_df.copy(),
                                        iso_to_do=iso_to_do,
                                        mzml=mzml,
                                        rt_tolerance=rt_tolerance,
                                        mass_tolerance=mass_tolerance,
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
        with futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            results = list(tqdm.tqdm(ex.map(integrate_one_partial, loop_),
                                    total=max(loop_),
                                    desc=f'Integrating Peaks in Sample {current_sample}'))
        # '''

        #
        # Convert the output_table into a data frame
        #
        df_columns = ['ID', 'pep_id'] + ['m' + str(iso) for iso in iso_to_do]

        integrated_peaks = [res[0] for res in results]

        integrated_df = pd.DataFrame(integrated_peaks, columns=df_columns)
        id_integrated_df = pd.merge(mzid.curr_frac_filtered_id_df, integrated_df, on='pep_id', how='left')
        id_integrated_df['file'] = mzml_files[idx]

        # Bind rows of the current result to the sample master
        if len(overall_integrated_df.index) == 0:
            sample_master_df = id_integrated_df
        else:
            sample_master_df = overall_integrated_df.append(id_integrated_df, ignore_index=True)

    # 2021-05-07 No longer creates subfolder
    # Create subdirectory if not exists
    # os.makedirs(os.path.join(directory_to_write, current_sample), exist_ok=True)
    # save_path = os.path.join(directory_to_write, current_sample, mzml_files[idx] + '_riana.txt')

    save_path = os.path.join(directory_to_write, current_sample + '_riana.txt')
    overall_integrated_df.to_csv(path_or_buf=save_path, sep='\t')

    # Make the soft-threshold data frame. These are the peptides that are ID'ed at 10 times the q-value
    # as the cut-off in this fraction up to q < 0.1, but has q >= q-value cutoff, and furthermore has been
    # consistently identified in the other samples at the same fraction (median fraction) at the q-value cutoff

    """ 20211109 Bind rows of the current sample master to the total (all time point output) master
    if len(master_df.index) == 0:
        master_df = sample_master_df
    else:
        master_df = master_df.append(sample_master_df, ignore_index=True)
    """

    # Write out the total time point output
    # 20211109 master_df.to_csv(path_to_write, sep='\t')

    return sys.exit(os.EX_OK)


def integrate_one(index: int,
                  id_: pd.DataFrame,
                  iso_to_do: list,
                  rt_tolerance: float,
                  mass_tolerance: float,
                  mzml,
                  ) -> list:
    """
    get all isotopomer mass intensities from ms1 scans within range and integrate

    :param index: int The row number of the peptide ID table passed from the wrapper.
    :param id_: protein identification dataframe
    :param iso_to_do: list of isotopomers
    :param rt_tolerance: retention time range in minutes
    :param mass_tolerance: relative mass tolerance (already converted from ppm) e.g., 50e-6
    :param mzml: mzml file object
    :return: list of intensity over time [index, pep_id, m0, m1, m2, ...]

    """

    # bypass all integration for match between runs
    if params.dummy:
        time.sleep(0.025)
        return [index] + [(id_.loc[index, 'pep_id'])] + [0 for _ in iso_to_do]

    # determine the mass of protons and c13
    proton = constants.proton_mass
    if params.deuterium_mass_defect:
        iso_added_mass = constants.deuterium_mass_diff  # see constants for details
    else:
        iso_added_mass = constants.c13_mass_diff

    # get peptide mass, scan number, and charge
    peptide_mass = float(id_.loc[index, 'peptide mass'])

    # 2021-05-01 percolator scan number is likely 1-indexed
    scan_number = int(id_.loc[index, 'scan']) - 1

    charge = float(id_.loc[index, 'charge'])

    # get retention time from Percolator scan number # 2021-05-19 added -1
    peptide_rt = mzml.rt_idx[np.searchsorted(mzml.scan_idx, scan_number, side='left') - 1 ]

    assert isinstance(peptide_rt.item(), float), '[error] cannot retrieve retention time from scan number'

    # calculate precursor mass from peptide monoisotopic mass
    peptide_prec = (peptide_mass + (charge * proton)) / charge

    intensity_over_time = []

    # choose the scan numbers from the index (beware that scan number is 1-indexed)
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
            delta_mass = prec_iso_am*mass_tolerance/2

            # sum intensities within range
            matching_int = np.sum(spec[np.abs(spec[:, 0] - prec_iso_am) <= delta_mass, 1])

            intensity_over_time.append([iso,
                                        mzml.rt_idx[mzml.scan_idx == scan].item(),
                                        matching_int,
            ]
            )

    if not intensity_over_time:
        raise Exception(f'No intensity profile for peptide {prec_iso_am}')

    intensity_array = np.array(intensity_over_time)

    if not intensity_over_time:
        raise Exception(f'Empty intensity over time for peptide {prec_iso_am}')

    integrated = [index] + [(id_.loc[index, 'pep_id'])] + integrate_isotope_intensity(intensity_array,
                                                                                  iso_to_do=iso_to_do)
    # 2021-11-10 return raw too
    return [integrated, intensity_array]


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
