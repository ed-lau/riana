"""
Relative Isotope Abundance Analyzer v.0.2.0
Edward Lau

Usage:
    riana.py --help
    riana.py integrate <mzid> <mzml>  [--thread=2 --iso=0,6 --lys=0 --test --qvalue=1e-2 --out=ria.txt]
    riana.py merge <index>

Options:
    -h --help            Show this screen.
    -v --version         Show version.
    -h --thread=2        Number of threads (max: 8) [default: 2]
    --iso=0,6            Isotopes to do separated by commas e.g., 0,1,2,3,4,5 [default: 0,6].
    -K --lys=0           Filter peptides with number of lysines (0=do not filter, 1=only one lysine,
                         2=any number of lysine [default: 0]
    -t --test            Integrate only the first 5 qualifying peptides for testing purpose
    -q --qvalue=1e-2              Integrate only peptides with percolator PEP below this threshold [default: 1e-2]
    --out=ria.txt        Name of output file [default: ria.txt]

Example:
    riana.py integrate percolator.target.mzid Heart_FASP_1.mzML --iso 0,6,12 -t -K 0 --qvalue 5e-3
    riana.py merge ria_index.txt

"""

from docopt import docopt
from parse_mzid import Mzid
from integrate_mzml import Mzml
import pandas as pd
import sys
import os
from celery import Celery

# Configure Celery
app = Celery('riana',
             backend='rpf://',
             broker='amqp://guest@localhost//')

def integrate(args):
    """
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

    #
    # Number of threads
    #

    num_threads = args['--thread']

    try:
        num_threads = int(num_threads)
        if num_threads < 1:
            num_threads = 1
        num_threads = max(num_threads, 8)
    except TypeError or ValueError:
        num_threads = 2

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
            print('Invalid PEP given - using default value.')
            qcutoff = float(1e-2)




    #
    # Open the mzid file
    #
    try:
        a = Mzid(mzid_loc)
        a.make_peptide_summary()
        a.filter_peptide_summary(lysine_filter=lysine_filter,
                                 protein_q=qcutoff,
                                 peptide_q=qcutoff)

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
        m = Mzml(mzml_loc)
    except OSError as e:
        sys.exit('Failed to load mzml file. ' + str(e.errno))

    #
    # Determine how many peptides to integrate. If test_run is on, only do 5 peptides to save time.
    #
    end = len(a.pep_summary_df)
    if test_run:
        end = min(5, len(a.pep_summary_df))

    if end == 0:
        return sys.exit(os.EX_CONFIG)

    #
    # Main loop through the mzidentml summary to pass tasks to the integration worker
    #

    pool = Pool(processes=num_threads)

    #
    # Prepare an empty list to encompass the output from each row of the mzidentml summary
    #
    output_table = [[] for i in range(0, end)]


    for i in range(0, end):

        print('Doing peptide ' + str(i) + ' of ' + str(len(a.filtered_pep_summary_df)))

        output_table[i] = integrate_worker.delay(a, m, i, iso_to_do)



    #
    # Convert the output_table into a data frame
    #
    df_columns = ['ID', 'uniprot', 'peptide', 'z', 'rt', 'calc_mz',]

    for each_iso in iso_to_do:
        df_columns.append('m' + str(each_iso))

    out_df = pd.DataFrame(output_table, columns=df_columns)

    out_df.to_csv(out_loc, sep='\t')

    return sys.exit(os.EX_OK)



@app.task
def integrate_worker(mzid, mzml, i, iso_to_do):
    """
    Worker for core integration: receive the index on the Mzid peptide summary files,
    and then pass it to the Mzml class to pull the integrated isotopomer abundances
    and return as a single list. The returned list will be bound in the main integrate function
    to a nested list that will then be converted to a data frame.

    :param mzid:
    :param mzml:
    :param i:
    :param iso_to_do:

    :return: out - a list containing isotopomer abundances.
    """
    out = []



    rt = mzml.get_rt_from_scan(float(mzid.filtered_pep_summary_df.loc[i, 'spectrum_id']))
    z = float(mzid.pep_summary_df.loc[i, 'z'])

    iso = mzml.get_isotopes_from_amrt(peptide_am=float(mzid.filtered_pep_summary_df.loc[i, 'calc_mz']),
                                      peptide_rt=rt,
                                      z=z,
                                      rt_tolerance=1.0,
                                      iso_to_do=iso_to_do)

    out.append(i)
    out.append(mzid.filtered_pep_summary_df.loc[i, 'uniprot'])
    out.append(mzid.filtered_pep_summary_df.loc[i, 'seq'])
    out.append(mzid.filtered_pep_summary_df.loc[i, 'z'])
    out.append(rt)
    out.append(mzid.filtered_pep_summary_df.loc[i, 'calc_mz'])

    for j in range(0, len(iso_to_do)):
        out.append(iso[j][1])

    print(out)

    return out


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
    if args['merge']:
        merge(args)


