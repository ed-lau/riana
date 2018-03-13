"""
Relative Isotope Abundance Analyzer v.0.3.0. Build Date : : :.
Written by Edward Lau (edward.lau@me.com) 2016-2018

Example
"""

from pymzid import Mzid
from integrate_mzml import Mzml
from time import time
import pandas as pd
import scipy
import sys
import os


# mzml_loc = '/Users/edwardlau/Desktop/Heart_FASP_1_small/Heart_FASP_1_small.mzML'
# mzid_loc = '/Users/edwardlau/Desktop/Heart_FASP_1_small/percolator.target.mzid'
# lysine_filter=1
# qcutoff=0.01
# unique_pep=1
# rt_tolerance=1
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

    # Handle command line options
    out_loc = args.out
    test_run = args.test
    unique_pep = args.unique

    #
    # Filter for lysine only
    #
    lysine_filter = args.lysine

    try:
        lysine_filter = int(lysine_filter)
        if lysine_filter not in [1, 2]:
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
    #
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

    #
    # Filter only for peptides with percolator posterior error probability below cutoff
    #

    #
    # Filter only for peptides in a protein
    #
    #a.pep_summary_df = a.pep_summary_df.loc[lambda x: x.uniprot.astype(str) != "nan", :]


    #
    # Read the mzml files
    #
    try:
        mzml = Mzml(mzml_loc)
    except OSError as e:
        sys.exit('Failed to load mzml file. ' + str(e.errno))

    #
    # Determine how many peptides to integrate. If test_run is on, only do 5 peptides to save time.
    #
    end = len(mzid.filtered_pep_summary_df)

    if test_run:
        end = min(2, len(mzid.filtered_pep_summary_df))

    if end == 0:
        return sys.exit(os.EX_CONFIG)

    print(mzid.filtered_pep_summary_df)

    print(end)

    #
    # Step 1: Get the list of peptides and their mass, isotopomers, and scan number to integrate
    #

    input_table = []

    for i in range(end):
        print(i)

        # Get ALL the MS1 scans to integrate for this particular peptide
        to_do = mzml.get_scans_to_do(int(mzid.filtered_pep_summary_df.loc[i, 'spectrum_id']), rt_tolerance)

        # Append the input table with ID, peptide_id, acc, peptide_seq, all scans to integrate, and z
        for scan_id, rt in to_do:

            input = [i,
                     mzid.filtered_pep_summary_df.loc[i, 'pep_id'],
                     mzid.filtered_pep_summary_df.loc[i, 'acc'],
                     mzid.filtered_pep_summary_df.loc[i, 'seq'],
                     mzid.filtered_pep_summary_df.loc[i, 'z'],
                     scan_id,
                     rt,
                     mzid.filtered_pep_summary_df.loc[i, 'calc_mz'],
                                ]

            input_table.append(input)

    #
    # Just for tidiness, convert the input_table to a dataframe
    #
    df_columns = ['ID', 'pep_id', 'acc', 'seq', 'z', 'scan_id', 'rt', 'calc_mz']

    in_df = pd.DataFrame(input_table, columns=df_columns)

    # Arrange by scan number so the mzml iterator doesn't have to do read the files over and over
    in_df = in_df.sort_values(by='scan_id').reset_index()

    #
    # Step 2: Get peak intensity for each isotopomer in each spectrum ID in each peptide
    #

    counter = len(in_df)
    # Prepare an empty list to encompass the output from each row of the mzidentml summary
    output_table = [[] for i in range(0, counter)]

    t1 = time()
    print('Extracting intensities from spectra...')



    for i in range(counter):

        out = []

        print('Integrating peptide-scan combination ' + str(i) + ' of ' + str(counter))

        if i % 10 == 0:
            t2 = time()
            print('Done. Extracting time: ' + str(round(t2 - t1, 2)) + ' seconds.')
            avg_time = (t2-t1)/(i+1)
            remaining_time = ((counter - i)/avg_time ) / 60
            print('Remaining time ' + str(remaining_time) + ' minutes.')

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

            output_table[i] = out

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

    integrated_out_df = pd.merge(mzid.filtered_pep_summary_df, integrated_df, on='pep_id', how='left')


    #
    # Convert the output_table into a data frame
    #


    integrated_out_df.to_csv(out_loc, sep='\t')

    return sys.exit(os.EX_OK)


#
# Code for running main with parsed arguments from command line
#

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='RIAna v.0.3.0 integrates the relative abundance of'
                                                 'isotopomers')


    parser.add_argument('mzid', help='path to mzid file')
    parser.add_argument('mzml', help='path to mzml file')

    parser.add_argument('-i', '--iso', help='isotopes to do, separated by commas, e.g., 0,1,2,3,4,5 [default: 0,6]',
                              default='0,6')
    parser.add_argument('-u', '--unique', action='store_true', help='integrate unique peptides only')
    parser.add_argument('-k', '--lysine',
                        help='lysine mode, 0=do not filter, 1=only one lysine, 2=any lysine [default = 0',
                        choices=[0,1,2],
                        default=0)
    parser.add_argument('-t', '--test', action='store_true',
                        help='test mode: integrates only first 5 qualifying peptides')
    parser.add_argument('-q', '--qvalue', help='integrate only peptides with q value below this threshold[default: 1e-2]',
                        default=1e-2)
    parser.add_argument('-r' '--rtime', help='retention time (in minutes, both directions) tolerance for integration',
                        default=1)
    parser.add_argument('-o', '--out', help='name of the output files [default: ria.txt]',
                        default='ria.txt')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose error messages')


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

