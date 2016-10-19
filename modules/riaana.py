"""
Relative Isotope Abundance Analyzer

Usage:
    riaana.py <mzid_file> <mzml_file> <out_file> --isotope_to_do=<iso> [--lys_filter=F --test_run=F]

Options:
    -h --help                   Show this screen.
    -v --version                Show version.
    --isotope_to_do=<iso>       Isotopes to do separated by commas e.g., 0,1,2,3,4,5 [default: 0,6].
    --lys_filter=F              Integrate only peptides with one lysine [default: F]
    --test_run=F                Integrate only 1 peptide for testing

Example:
    riaana.py /Users/edwardlau/Desktop/crux-3.0.Darwin.i386/mouse_aa_heart_fasp_1/percolator.target.mzid /Volumes/NewVolumn/Heart_FASP_1.mzML out --isotope_to_do 0,6,12



"""

from docopt import docopt
from parse_mzid import Mzid
from integrate_mzml import Mzml
import pandas as pd


mzid_loc = "/Users/edwardlau/Desktop/crux-3.0.Darwin.i386/mouse_aa_heart_fasp_1/percolator.target.mzid"
#mzml_loc = "/Volumes/Data/proteomics/Liverpool_mouse_turnover/amino_acid_labeling/heart/Heart_FASP_1.mzML"
mzml_loc = "/Volumes/New Volume/Heart_FASP_1.mzML"
iso_to_do = [0,6]

def main(args):

    out_loc = args['<out_file>']

    mzid_loc = args['<mzid_file>']

    mzml_loc = args['<mzml_file>']

    lysine_only = (args['--lys_filter'] == 'T')
    print(lysine_only)

    test_run = (args['--test_run'] == 'T')

    # Convert the arguments into a list of integers
    iso_to_do = []
    for char in args['--isotope_to_do'].split(','):
        try:
            iso_to_do.append(int(char))
        except:
            pass

    print(iso_to_do)
    print(type(iso_to_do))

    if not len(iso_to_do) > 0:
        return 'Isotopomers must be supplied as list of integer'

    print(iso_to_do)

    a = Mzid(mzid_loc)
    a.parse_file()
    a.read_dbs()
    a.read_peptide()
    a.read_pe()
    a.read_psm()
    a.read_protein()
    a.make_summary()


    # Filter only for peptides with percolator PEP < 0.1
    a.pep_summary_df = a.pep_summary_df.loc[lambda x: x.percolator_PEP.astype(float) < 0.1, :]

    # Filter only for peptides in a protein
    #a.pep_summary_df = a.pep_summary_df.loc[lambda x: x.uniprot.astype(str) != "nan", :]


    # Reset index
    a.pep_summary_df = a.pep_summary_df.reset_index()

    # Write up CSV
    # import sys
    # b.to_csv(sys.stdout)
    # b.to_csv("~/Desktop/b.csv")


    #
    # Read the mzml files
    #

    m = Mzml(mzml_loc)

    # Prepare an empty list to encompass the output from each row of the mzidentml summary
    output_table = []
    #for j in range(0,len(iso_to_do)):
    #    output_table.append([])


    # Main loop through the mzidentml summary


    end = len(a.pep_summary_df)
    if test_run:
        end = 1

    for i in range(0, end):

        print('Doing peptide ' + str(i) + ' of ' + str(len(a.pep_summary_df)))


        # If the lysine filter is on and the peptide has no or more than one lysine, skup
        if lysine_only and a.pep_summary_df.loc[i, 'seq'].count('K') != 1:
            print('Skipping peptide with K != 1')
            continue

        iso = m.get_isotopes_from_amrt(peptide_am=float(a.pep_summary_df.loc[i, 'calc_mz']),
                                       peptide_rt=m.get_rt_from_scan(float(a.pep_summary_df.loc[i, 'spectrum_id'])),
                                       z=float(a.pep_summary_df.loc[i, 'z']),
                                       rt_tolerance=1.0,
                                       iso_to_do=iso_to_do)

        print(iso)

        out  = []

        out.append(i)
        out.append(a.pep_summary_df.loc[i, 'seq'])
        out.append(a.pep_summary_df.loc[i, 'z'])
        out.append(m.get_rt_from_scan(float(a.pep_summary_df.loc[i, 'spectrum_id'])))
        out.append(a.pep_summary_df.loc[i, 'calc_mz'])

        for j in range(0, len(iso_to_do)):
             out.append(iso[j][1])

        output_table.append(out)

        print(out)

    print(output_table)
    out_df = pd.DataFrame(output_table)
    print(out_df)

    out_df.to_csv(out_loc, sep='\t')



if __name__ == "__main__":
   args = docopt(__doc__, version='RIA Analyzer v.0.2.0')
   print(args)
   main(args)

