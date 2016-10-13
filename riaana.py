"""
Relative Isotope Abundance Analyzer

Usage:
    riaana.py <mzid_file> <mzml_file> <out_file> --isotope_to_do=<iso>

Options:
    -h --help           Show this screen.
    -v --version        Show version.
    --isotope_to_do=<iso>     Isotopes to do separated by commas e.g., 0,1,2,3,4,5 [default: 0,6].

Example:
    riaana.py /Users/edwardlau/Desktop/crux-3.0.Darwin.i386/mouse_aa_heart_fasp_1/percolator.target.mzid /Users/edwardlau/Desktop/Heart_FASP_1.mzML out --isotope_to_do 0,6,12



"""

from docopt import docopt
from parse_mzid import Mzid
from integrate_mzml import Mzml


mzid_loc = "/Users/edwardlau/Desktop/crux-3.0.Darwin.i386/mouse_aa_heart_fasp_1/percolator.target.mzid"
mzml_loc = "/Users/edwardlau/Desktop/Heart_FASP_1.mzML"
iso_to_do = [0,6]

def main(args):

    output = args['<out_file>']

    mzid_loc = args['<mzid_file>']

    mzml_loc = args['<mzml_file>']

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



    output_table = []

    for j in range(0,len(iso_to_do)):
        output_table.append([])


    for i in range(1, len(a.pep_summary_df)):

        print('Doing peptide ' + str(i) + ' of ' + str(len(a.pep_summary_df)))

        if a.pep_summary_df.loc[i, 'seq'].count('K') != 1:
            print('Skipping peptide with K != 1')

            for j in range(0, len(iso_to_do)):
                output_table[j].append(-1)

            continue

        iso = m.get_isotopes_from_amrt(peptide_am=float(a.pep_summary_df.loc[i, 'calc_mz']),
                                       peptide_rt=m.get_rt_from_scan(float(a.pep_summary_df.loc[i, 'spectrum_id'])),
                                       z=float(a.pep_summary_df.loc[i, 'z']),
                                       rt_tolerance=1.0,
                                       iso_to_do=iso_to_do)

        for j in range(0, len(iso_to_do)):
            output_table[j].append(iso[j][1])

        print(iso)



    # from parse_mzid import SearchResult








if __name__ == "__main__":
   args = docopt(__doc__, version='RIA Analyzer v.0.1.0')
   print(args)
   main(args)

# python3 main.py rmats out_file
