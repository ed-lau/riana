"""
Relative Isotope Abundance Analyzer v.0.4.0. Build Date : : :.
Written by Edward Lau (lau1@stanford.edu) 2016-2018

Example: python riana.py percolator_test/percolator percolator_test/mzml -v 2 -u -i 0,6,12 -q 0.01 -r 0.5 -k 0 -o kk_test


"""

# from pymzid import Mzid
from readpercolator import ReadPercolator
from integrate_mzml import Mzml
# from time import time
import pandas as pd
import scipy
import re
import os
# NEW 2018-09-09 Now using tqdm for progress bar
from tqdm import tqdm

# mzml_loc = 'percolator_test/mzml'
# mzid_loc = 'percolator_test/percolator'
# lysine_filter=0
# qcutoff=0.001
# unique_pep=1
# rt_tolerance=0.5
# iso_to_do = [0,1,2,3,4,5,6]


def integrate(args):
    """
    Improved process to integrate for isotope abundance analysis.
    Idea is to loop through the mzML only once - get all the peptides to be integrated first
    Find all the spectrum ID - m/z combinations to be integrated, then integrate them sequentially

    :param args: Arguments from command line
    :return:

    """

    # Handle command line arguments
    mzid_loc = args.mzid
    mzml_loc = args.mzml

    assert os.path.isdir(mzml_loc), '[error] mzml directory not valid'
    assert os.path.isdir(mzid_loc), '[error] percolator directory not valid'

    # Handle command line options - creating some placeholders for easier testing.
    unique_pep = args.unique

    #
    # Filter for lysine only
    #
    lysine_filter = args.lysine

    try:
        lysine_filter = int(lysine_filter)
        if lysine_filter not in [1, 2, 3]:
            lysine_filter = 0

    except TypeError or ValueError:
        lysine_filter = 0

    print(lysine_filter)


    #
    # Convert the to-do isotopomer list option into a list of integers
    #
    iso_to_do = []
    for char in args.iso.split(','):

        try:
            char = int(char)

            # Only doing down to mz + 12.
            if char <= 12:
                iso_to_do.append(int(char))

        except ValueError or TypeError:
            pass

    if not iso_to_do:
        sys.exit('Error: Invalid isotopomer list given.')

    iso_to_do = list(set(iso_to_do))
    iso_to_do.sort()

    print('Recognized isotopomer list:', str(iso_to_do))

    #
    # Percolator q value cutoff for peptides and proteins
    #
    if args.qvalue:

        try:
            qcutoff = float(args.qvalue)

        except ValueError or TypeError:
            print('Invalid Q value given - using default value.')
            qcutoff = float(1e-2)

    #
    # Retention time cutoff peptides and proteins
    #
    if args.rtime:

        try:
            rt_tolerance = float(args.rtime)

        except ValueError or TypeError:
            print('Invalid retention time tolerance given - using default value.')
            rt_tolerance = float(1.0)
    else:
        rt_tolerance = float(1.0)

    #
    # Open the mzid file
    # Note 2018-09-07 - I am not sure we should use the Mzid module anymore since Percolator doesn't seem to output
    # Mzid files that conform completely to Mzid standards. In particular I wasn't able to locate where in the xml file
    # is the file_idx (mzml file from which the spectrum was found) is encoded. It would be much easier at this point to use
    # the percolator tab-delimited file. Unless I can look at a multi-fraction Percolator run Mzid and modify the Mzid module from there
    #
    """
    try:
        mzid = Mzid(mzid_loc)
        mzid.make_peptide_summary(take_psm_df_cvParams=False)
        mzid.filter_peptide_summary(lysine_filter=lysine_filter,
                                    protein_q=qcutoff,
                                    peptide_q=1,
                                    unique_only=unique_pep,
                                    require_protein_id=True)

    except OSError as e:
        sys.exit('Failed to load mzid file. ' + str(e.errno))
    """

    try:
        mzid = ReadPercolator(mzid_loc)
        mzid.get_mzid_indices()
        mzid.filter_id_df(lysine_filter=lysine_filter,
                          protein_q=qcutoff,
                          peptide_q=qcutoff,
                          unique_only=unique_pep,
                          require_protein_id=True)

    except OSError as e:
        sys.exit('Failed to load mzid file. ' + str(e.errno))

    # Check that the number of mzMLs in the mzML folder is the same as the maximum of the ID file's file_idx column.
    # Note this will throw an error if not every fraction results in at least some ID, but we will ignore for now.

    mzml_files = [f for f in os.listdir(mzml_loc) if re.match('^.*.mzML', f)]

    # Sort the mzML files by names
    # Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
    # But we will ignore for now
    mzml_files.sort()

    # Throw an error if there is no mzML file in the mzml directory
    assert len(mzml_files) != 0, '[error] no mzml files in the specified directory'
    assert len(mzml_files) == max(mzid.indices) + 1, '[error] number of mzml files not matching id list'

    #
    # Read the mzml files
    #

    # For each file index (fraction), open the mzML file, and create a subset Percolator ID dataframe
    for idx in mzid.indices:

        # Verbosity 0 progress message
        print('[verbosity 0] doing mzml:', mzml_files[idx],
              '(' + str(idx + 1), 'of', str(len(mzid.indices)) + ')', sep=' ')

        # Make a subset dataframe with the current file index (fraction) being considered
        fraction_id_df = mzid.subset_id_df(idx)

        # Arrange the PSM rows by scan number
        fraction_id_df = fraction_id_df.sort_values(by='scan').reset_index(drop=True)

        # Add a dummy peptide ID for this fraction only
        fraction_id_df = fraction_id_df.assign(pep_id=fraction_id_df.index)

        try:
            mzml = Mzml(os.path.join(mzml_loc, mzml_files[idx]))
        except OSError as e:
            sys.exit('[error] failed to load fraction mzml file. ' + str(e.errno))

        #
        # Determine how many peptides to integrate. If test_run is on, only do 50 peptides to save time.
        #
        """
        end = len(fraction_id_df)

        if args.test:
            end = min(50, len(fraction_id_df))

        if end == 0:
            return sys.exit(os.EX_CONFIG)

        if args.verbosity == 2:
            print(fraction_id_df)

        print(end)
        """

    #
    # Step 1: Get the list of peptides and their mass, isotopomers, and scan number to integrate
    #

        input_table = []

        for i in fraction_id_df.pep_id:

            if args.verbosity >= 1 and i % 10 == 0:
                print("verbosity 1: now getting input list " + str(i))

            # Get ALL the MS1 scans to integrate for this particular peptide
            to_do = mzml.get_scans_to_do(int(fraction_id_df.loc[i, 'scan']), rt_tolerance)

            # Append the input table with ID, peptide_id, acc, peptide_seq, all scans to integrate, and z
            for scan_id, rt in to_do:

                input = [i,
                         fraction_id_df.loc[i, 'protein id'],
                         fraction_id_df.loc[i, 'sequence'],
                         fraction_id_df.loc[i, 'charge'],
                         scan_id,
                         rt,
                         fraction_id_df.loc[i, 'spectrum precursor m/z'],
                                    ]

                input_table.append(input)

        #
        # Just for tidiness, convert the input_table to a dataframe
        #
        df_columns = ['pep_id', 'acc', 'seq', 'z', 'scan_id', 'rt', 'calc_mz']

        in_df = pd.DataFrame(input_table, columns=df_columns)

        # Arrange by scan number so the mzml iterator doesn't have to do read the files over and over
        in_df = in_df.sort_values(by='scan_id').reset_index()

        #
        # Step 2: Get peak intensity for each isotopomer in each spectrum ID in each peptide
        #

        counter = len(in_df)
        # Prepare an empty list to encompass the output from each row of the mzidentml summary
        output_table = []

        #t1 = time()

        print('Integrating peak intensities from spectra.')

        for i in tqdm(range(counter)):

            out = []

            # If verbose, print out message and time for every 10 peptide-scan combinations
            """
            if args.verbosity >= 1 and i % 10 == 0:
                print('verbosity 1: integrating peptide-scan combination ' + str(i) + ' of ' + str(counter))

                t2 = time()
                print('verbosity 1: elapsed time of current fraction: ' + str(round((t2 - t1)/60, 2)) + ' minutes.')
                avg_time = (t2-t1)/(i+1)
                remaining_time = ((counter - i) * avg_time) / 60
                print('verbosity 1: estimated remaining time of current fraction: ' + str(round(remaining_time, 2)) + ' minutes.')
            """

            # Print out the accurate mass and scan number currently being integrated
            if args.verbosity == 2:
                print('verbosity 2: integrating mz: ', str(round(float(in_df.loc[i, 'calc_mz']), 4)),
                      ' scan id: ', str(round(in_df.loc[i, 'scan_id'], 2)))

            # Get the intensities of all isotopomers from the spectrum ID
            # NB: the calc_mz from the mzML file is the monoisotopic m/z
            iso = mzml.get_isotope_from_scan_id(peptide_am=float(in_df.loc[i, 'calc_mz']),
                                                z=float(in_df.loc[i, 'z']),
                                                spectrum_id=in_df.loc[i, 'scan_id'],
                                                iso_to_do=iso_to_do)

            if iso:
                out.append(i)
                out.append(in_df.loc[i, 'pep_id'])
                out.append(in_df.loc[i, 'acc'])
                out.append(in_df.loc[i, 'seq'])
                out.append(in_df.loc[i, 'z'])
                out.append(in_df.loc[i, 'scan_id'])
                out.append(in_df.loc[i, 'rt'])
                out.append(in_df.loc[i, 'calc_mz'])

                for j in range(0, len(iso_to_do)):
                    out.append(iso[j][2])

                #print(out)

                output_table.append(out)

        df_columns = ['ID', 'pep_id', 'acc', 'seq', 'z', 'scan_id', 'rt', 'calc_mz']

        for each_iso in iso_to_do:
            df_columns.append('m' + str(each_iso))

        out_df = pd.DataFrame(output_table, columns=df_columns)

        # Arrange by scan number so the mzml iterator doesn't have to do read the files over and over
        out_df = out_df.sort_values(by=['pep_id', 'ID']).reset_index()

        #
        # Step 3: Integrate each peptide (intensity over retention time)
        #
        cols_to_integrate = ['m' + str(iso) for iso in iso_to_do]
        peptides_to_integrate = list(set(list(out_df.pep_id)))

        integrated_output = []

        for peptide in peptides_to_integrate:

            # Subsetting the dataframe to consider only the current pep_id
            subset_df = out_df[out_df['pep_id'] == peptide]
            rt = subset_df['rt']

            # Integrate each isotopomer over the retention time!
            integrated = subset_df.groupby('pep_id')[cols_to_integrate].agg(lambda x: scipy.integrate.trapz(x, rt)).values.tolist()[0]
            integrated = [peptide] + integrated
            integrated_output.append(integrated)

        integrated_df = pd.DataFrame(integrated_output, columns=['pep_id']+cols_to_integrate)

        integrated_out_df = pd.merge(fraction_id_df, integrated_df, on='pep_id', how='left')


        #
        # Convert the output_table into a data frame
        #
        # Create directory if not exists
        os.makedirs(args.out, exist_ok=True)
        save_path = os.path.join(args.out, mzml_files[idx] + '_integrated_out.txt')

        integrated_out_df.to_csv(save_path, sep='\t')

    return sys.exit(os.EX_OK)



#
# Code for running main with parsed arguments from command line
#

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='RIAna v.0.4.0 integrates the relative abundance of'
                                                 'isotopomers')


    parser.add_argument('mzid', help='path to folder containing the search result file (percolator tab delimited)')

    parser.add_argument('mzml', help='path to folder containing mzml files')

    parser.add_argument('-i', '--iso', help='isotopes to do, separated by commas, e.g., 0,1,2,3,4,5 [default: 0,6]',
                              default='0,6')

    parser.add_argument('-u', '--unique', action='store_true', help='integrate unique peptides only')

    parser.add_argument('-k', '--lysine',
                        help='lysine mode, 0=No filter, 1=1 K, 2=1 or more K, 3=KK only [default = 0]',
                        type=int,
                        choices=[0, 1, 2, 3],
                        default=0)

    """
    parser.add_argument('-t', '--test', action='store_true',
                        help='test mode: integrates only first 50 qualifying peptides')
    """

    parser.add_argument('-q', '--qvalue',
                        help='integrate only peptides with q value below this threshold[default: 1e-2]',
                        type=float,
                        default=1e-2)

    parser.add_argument('-r', '--rtime', help='retention time (in minutes, both directions) tolerance for integration',
                        type=float,
                        default=1.0)

    parser.add_argument('-o', '--out', help='name of the output directory [default: riana_out]',
                        default='riana_out')

    parser.add_argument('-v', '--verbosity',
                        help='verbosity of error messages level. 0=quiet, 1=default, 2=verbose',
                        type=int,
                        choices=[0, 1, 2],
                        default=1)

    parser.set_defaults(func=integrate)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)

