"""
Relative Isotope Abundance Analyzer v.0.2.0. Build Date : : :.
Written by Edward Lau (edward.lau@me.com) 2016-2017


Usage:
    riana.py --help
    riana.py integrate <mzid> <mzml>  [--iso=0,6 --lys=0 --unique --test --qvalue=1e-2 --rt=1.0 --out=ria.txt ]
    riana.py integrate_fast <mzid> <mzml>  [--iso=0,6 --lys=0 --unique --test --qvalue=1e-2 --rt=1.0 --out=ria.txt ]
    riana.py merge <index>

Options:
    -h --help            Show this screen.
    -v --version         Show version.
    --iso=0,6            Isotopes to do separated by commas e.g., 0,1,2,3,4,5 [default: 0,6].
    -u --unique          Integrate unique peptides only
    -K --lys=0           Filter peptides with number of lysines (0=do not filter, 1=only one lysine,
                         2=any number of lysines [default: 0]
    -t --test            Integrate only the first 5 qualifying peptides for testing purpose
    -q --qvalue=1e-2     Integrate only peptides/proteins with percolator q value below this threshold [default: 1e-2]
    --rt=1.0             Retention time (in minutes, both directions), tolerance for integration
    --out=ria.txt        Name of output file [default: ria.txt]

Example:
    riana.py integrate_fast percolator.target.mzid Heart_FASP_1.mzML --iso 0,6,12 -t -K 0 --qvalue 5e-3
    riana.py merge ria_index.txt

"""

from docopt import docopt
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
    DEPRECATED
    Main loop for isotope abundance analysis

    :param args: Arguments from command line
    :return:

    """

    # Handle command line arguments
    mzid_loc = args['<mzid>']
    mzml_loc = args['<mzml>']

    # Handle command line options
    out_loc = args['--out']
    test_run = args['--test']
    unique_pep = args['--unique']

    #
    # Filter for lysine only
    #
    lysine_filter = args['--lys']

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
    for char in args['--iso'].split(','):

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
    if args['--qvalue']:

        try:
            qcutoff = float(args['--qvalue'])

        except ValueError or TypeError:
            print('Invalid Q value given - using default value.')
            qcutoff = float(1e-2)

    #
    # Retention time cutoff peptides and proteins
    #
    if args['--rt']:

        try:
            rt_tolerance = float(args['--rt'])

        except ValueError or TypeError:
            print('Invalid retention time tolerance given - using default value.')
            rt_tolerance = float(1.0)
    else:
        rt_tolerance = float(1.0)

    #
    # Open the mzid file. Note that take_psm_df_cvParams should be set to False to avoid conflict
    # between the identical "Percolator:Q value" cvParam names inside both PSM and Peptide fields.
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
    #a.pep_summary_df = a.pep_summary_df.loc[lambda x: x.acc.astype(str) != "nan", :]


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
        end = min(5, len(mzid.filtered_pep_summary_df))

    if end == 0:
        return sys.exit(os.EX_CONFIG)

    print(mzid.filtered_pep_summary_df)

    print(end)
    #
    # Main loop through the mzidentml summary to pass tasks to the integration worker
    #


    #
    # Prepare an empty list to encompass the output from each row of the mzidentml summary
    #
    output_table = [[] for i in range(0, end)]

    for i in range(end):

        out = []


        print('Integrating peptide ' + str(i) + ' of ' + str(end))

        rt = mzml.get_rt_from_scan(float(mzid.filtered_pep_summary_df.loc[i, 'spectrum_id']))
        z = float(mzid.filtered_pep_summary_df.loc[i, 'z'])

        iso = mzml.get_isotopes_from_amrt(peptide_am=float(mzid.filtered_pep_summary_df.loc[i, 'calc_mz']),
                                          peptide_rt=rt,
                                          z=z,
                                          rt_tolerance=rt_tolerance,
                                          iso_to_do=iso_to_do)

        if iso:
            out.append(i)
            out.append(mzid.filtered_pep_summary_df.loc[i, 'acc'])
            out.append(mzid.filtered_pep_summary_df.loc[i, 'seq'])
            out.append(mzid.filtered_pep_summary_df.loc[i, 'z'])
            out.append(rt)
            out.append(mzid.filtered_pep_summary_df.loc[i, 'calc_mz'])

            for j in range(0, len(iso_to_do)):
                out.append(iso[j][1])

            print(out)

            output_table[i] = out



    #
    # Convert the output_table into a data frame
    #
    df_columns = ['ID', 'acc', 'peptide', 'z', 'rt', 'calc_mz',]

    for each_iso in iso_to_do:
        df_columns.append('m' + str(each_iso))

    out_df = pd.DataFrame(output_table, columns=df_columns)

    out_df.to_csv(out_loc, sep='\t')

    return sys.exit(os.EX_OK)



def integrate_fast(args):
    """
    Second try: improve the process to integrate for isotope abundance analysis.
    Idea is to loop through the mzML only once - get all the peptides to be integrated first
    Find all the spectrum ID - m/z combinations to be integrated, then integrate them sequentialy

    :param args: Arguments from command line
    :return:

    """

    # Handle command line arguments
    mzid_loc = args['<mzid>']
    mzml_loc = args['<mzml>']

    # Handle command line options
    out_loc = args['--out']
    test_run = args['--test']
    unique_pep = args['--unique']

    #
    # Filter for lysine only
    #
    lysine_filter = args['--lys']

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
    for char in args['--iso'].split(','):

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
    if args['--qvalue']:

        try:
            qcutoff = float(args['--qvalue'])

        except ValueError or TypeError:
            print('Invalid Q value given - using default value.')
            qcutoff = float(1e-2)

    #
    # Retention time cutoff peptides and proteins
    #
    if args['--rt']:

        try:
            rt_tolerance = float(args['--rt'])

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



def merge(args):
    print('Merging feature to be implemented in coming updates.')


#
# Docopt
#
if __name__ == "__main__":
    args = docopt(__doc__, version='Relative Isotope Analyser v.0.2.0')
    print(args)
    if args['integrate']:
        integrate(args)
    elif args['integrate_fast']:
        integrate_fast(args)
    elif args['merge']:
        merge(args)


