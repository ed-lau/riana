"""
Relative Isotope Abundance Analyzer v.0.6.0. Build Date : : :.
Written by Edward Lau (lau1@stanford.edu) 2016-2019

Example: python riana.py data/percolator_test/percolator data/percolator_test/mzml -v 2 -u -i 0,6,12 -q 0.01 -r 0.5 -k 3 -t 2 -m 100 -o out_test

"""

# from pymzid import Mzid
from readpercolator import ReadPercolator
from parse_mzml import Mzml
from integrate_peaks import Peaks
import pandas as pd
import re
import os
from multiprocessing import cpu_count


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
    vb = args.verbosity

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

    if vb > 0:
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
    # MS1 mass tolerance for integration
    #
    if args.masstolerance:
        try:
            mass_tolerance = float(args.masstolerance) * 1e-6

        except ValueError or TypeError:
            print('Invalid mass tolerance given - using default value.')
            mass_tolerance = float(100) * 1e-6
    else:
        mass_tolerance = float(100) * 1e-6

    #
    # Multi-threading
    #
    if args.thread:
        try:
            num_thread = max(cpu_count()-1, int(args.thread))

        except ValueError or TypeError:
            print('Invalid thread count, using default: 1.')
            num_thread = 4
    else:
        num_thread = 4

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
        print('[verbosity 0] doing mzml: {0} ({1} of {2})'.format(
            mzml_files[idx],
            str(idx + 1),
            str(len(mzid.indices))))

        # Make a subset dataframe with the current file index (fraction) being considered
        mzid.subset_id_df(idx)

        try:
            mzml = Mzml(os.path.join(mzml_loc, mzml_files[idx]))
        except OSError as e:
            sys.exit('[error] failed to load fraction mzml file. ' + str(e.errno))

        #
        # Read the spectra into dictionary and also create MS1/MS2 indices as a Peaks object
        #
        mzml.parse_mzml()
        peaks = Peaks(mzml.msdata, mzml.rt_idx, mzml.mslvl_idx)

        #
        # Link ID file, iso_to_do, and rt_tolerance to mzML
        #
        peaks.associate_id(mzid.fraction_id_df)
        peaks.set_iso_to_do(iso_to_do)
        peaks.set_rt_tolerance(rt_tolerance)
        peaks.set_mass_tolerance(mass_tolerance)

        #
        # Get peak intensity for each isotopomer in each spectrum ID in each peptide
        #
        output_table = peaks.get_isotopes_from_amrt_multiwrapper(num_thread=num_thread)

        df_columns = ['ID', 'pep_id'] + ['m' + str(iso) for iso in iso_to_do]
        out_df = pd.DataFrame(output_table, columns=df_columns)

        integrated_out_df = pd.merge(peaks.id, out_df, on='pep_id', how='left')

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

    parser = argparse.ArgumentParser(description='RIANA.py v.0.5.0 integrates the relative abundance of'
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

    parser.add_argument('-m', '--masstolerance', help='mass tolerance in ppm for integration [default 100 ppm]',
                        type=float,
                        default=100)

    parser.add_argument('-t', '--thread', help='thread (default = 4)',
                        type=float,
                        default=4.0)

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

    #gc.enable()
    #gc.set_debug(gc.DEBUG_LEAK)

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
