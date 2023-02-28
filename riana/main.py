# -*- coding: utf-8 -*-

""" Main. """

import os
import argparse
from riana import riana_integrate, riana_fit, riana_preprocess, __version__
from typing import List


# Check thread count is below cpu count
class CheckThreadCount(argparse.Action):
    def __call__(self, parser, namespace, values: int, option_string=None):
        if values > os.cpu_count():
            raise argparse.ArgumentTypeError("Thread count must be lower than CPU count")
        setattr(namespace, self.dest, values)


class StoreUniqueSortedIsotopomers(argparse.Action):
    """Checks that the list of arguments contains no duplicates, then stores"""
    def __call__(self, parser, namespace, values: List[int], option_string=None):
        if len(values) > len(set(values)):
            raise argparse.ArgumentError(
                self,
                "You cannot specify the same value multiple times. "
                + f"You provided {values}",
            )
        values.sort()
        setattr(namespace, self.dest, values)


class CheckSampleNameEndsWithNumber(argparse.Action):
    """ Check sample name contains a number then stores"""
    def __call__(self, parser, namespace, values, option_string=None):
        if not values[-1].isdigit():
            raise argparse.ArgumentError(
                self,
                "Sample name must end with a number. "
                + f"You provided {values}",
            )
        setattr(namespace, self.dest, values)



class CheckReadableDir(argparse.Action):
    """ Class to check if directory is readable. """
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


class CheckRTime(argparse.Action):
    """ Class to check r_time is a float between 0 and 10. """
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            values = float(values)
        except ValueError:
            raise argparse.ArgumentTypeError("%r for r_time not a floating-point literal" % (values,))

        if values < 0.0 or values > 10.0:
            raise argparse.ArgumentTypeError("%r for r_time not in range [0.0, 10.0]" % (values,))
        setattr(namespace, self.dest, values)


class CheckQValue(argparse.Action):
    """ Class to check that q values are between 0 and 1. """
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            values = float(values)
        except ValueError:
            raise argparse.ArgumentTypeError("%r for q_value not a floating-point literal" % (values,))

        if values < 0.0 or values > 1.0:
            raise argparse.ArgumentTypeError("%r for q_value not in range [0.0, 1.0]" % (values,))
        setattr(namespace, self.dest, values)

# Check that each character in the -aa argument is a valid amino acid
class CheckAminoAcids(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        for aa in values:
            if aa not in aa_list:
                raise argparse.ArgumentTypeError("Amino acid %r not recognised" % (aa,))
        setattr(namespace, self.dest, values)


# ---- Code for running main with parsed arguments from command line ----
def main():
    """ Main entry point for the riana script """
    parser = argparse.ArgumentParser(description='Riana integrates the relative abundance of'
                                                 ' isotopomers in mass spectrometry data and performs'
                                                 'kinetics modeling',
                                     epilog='For more information, see GitHub repository at '
                                               'https://github.com/ed-lau/riana',
                                     )

    parser.add_argument('-v', '--version',
                        action='version',
                        version='riana {version}'.format(version=__version__))

    # Sub-commands
    subparsers = parser.add_subparsers(help='Type riana function -h for individual help messages',
                                       title='Functions',
                                       description='Riana has the following sub-commands:',
                                       )

    parser_preprocess = subparsers.add_parser('preprocess',
                                          help='Preprocesses Percolator input files',
                                          description='Preprocesses Percolator input files',
                                          epilog='For more information, see GitHub repository at '
                                                 'https://github.com/ed-lau/riana',
                                            )

    parser_integrate = subparsers.add_parser('integrate',
                                             help='Integrates isotopomer abundance over retention time',
                                             description='Integrates isotopomer abundance over retention time',
                                             epilog='For more information, see GitHub repository at '
                                                    'https://github.com/ed-lau/riana',
                                            )

    parser_fit = subparsers.add_parser('fit',
                                       help='Fit to kinetic models *under development*')

    #
    # Arguments for preprocess sub-command
    #

    parser_preprocess.add_argument('id_files',
                                   type=argparse.FileType('r'),
                                   nargs='+',
                                   help='<required> path to Percolator input file(s)')

    parser_preprocess.add_argument('-o', '--out',
                                   help= 'path to output directory',
                                   default = '.')

    parser_preprocess.set_defaults(func=riana_preprocess.main)
    #
    # Arguments for integrate subcommand
    #
    parser_integrate.add_argument('mzml_path',
                                  type=str,
                                  help='<required> path to folder containing the mzml files',
                                  action=CheckReadableDir,
                                  )

    parser_integrate.add_argument('id_path',
                                  type=argparse.FileType('r'),
                                  help='<required> path to the percolator output psms.txt file',
                                  )

    parser_integrate.add_argument('-s', '--sample',
                                  help='sample name to override mzml folder name, must end with a number, e.g., time1',
                                  type=str,
                                  default='time0',
                                  action=CheckSampleNameEndsWithNumber,
                                  )

    parser_integrate.add_argument('-i', '--iso',
                                  help='isotopomer(s) to integrate',
                                  default='0 6',
                                  nargs='+',
                                  type=int,
                                  choices=range(0, 21),
                                  metavar='[0-20]',
                                  action=StoreUniqueSortedIsotopomers,
                                  )

    parser_integrate.add_argument('-u', '--unique',
                                  action='store_true',
                                  help='integrate unique peptides only')

    parser_integrate.add_argument('-t', '--thread',
                                  help='number of threads for concurrency [default: 1]',
                                  type=int,
                                  default=1,
                                  action=CheckThreadCount,
                                  )

    parser_integrate.add_argument('-o', '--out',
                                  help='path to the output directory [default: riana]',
                                  action=CheckReadableDir,
                                  default='.',
                                  )

    parser_integrate.add_argument('-q', '--q_value',
                                  help='integrate only peptides with q value below this threshold [default: 1e-2]',
                                  metavar="FDR[0,1]",
                                  type=float,
                                  action=CheckQValue,
                                  default=1e-2)

    parser_integrate.add_argument('-r', '--r_time',
                                  help='retention time (in minutes, both directions) tolerance for integration',
                                  type=float,
                                  action=CheckRTime,
                                  default=1.0)

    parser_integrate.add_argument('-w', '--write_intensities',
                                  action='store_true',
                                  help='also write pre-integration intensities into a result file')

    parser_integrate.add_argument('-m', '--mass_tol',
                                  help='<integer> mass tolerance in ppm for integration [default 50 ppm]',
                                  type=int,
                                  choices=range(1, 101),
                                  metavar='[1-100]',
                                  default=50)

    parser_integrate.add_argument('-S', '--smoothing',
                                  help='smoothing window size for integration',
                                  type=int,
                                  choices=range(3, 18, 2),
                                  )

    parser_integrate.add_argument('-D', '--mass_defect',
                                  type=str,
                                  choices=['D', 'C13', 'SILAC'],
                                  default='D',
                                  help='mass defect type [default: D]')

    parser_integrate.set_defaults(func=riana_integrate.integrate_all)

    #
    # Arguments for fit subcommand
    #
    parser_fit.add_argument('riana_path',
                            nargs='+',
                            type=str,
                            help='<required> paths to one or more integrate out text '
                                 'files (note: the sample field must include numericals '
                                 'corresponding to time units (e.g., time0, time 6)',
                            )

    parser_fit.add_argument('-m', '--model',
                            type=str,
                            choices=['simple', 'guan', 'fornasiero'],
                            default='simple',
                            help='kinetic models for fitting, currently only the simple '
                                 'exponential model is implemented [default: simple]',
                            )

    parser_fit.add_argument('-l', '--label',
                            type=str,
                            choices=['aa', 'hw', 'o18'],
                            default='hw',
                            help='labeling type [default: hw]')

    parser_fit.add_argument('-a', '--aa',
                            type=str,
                            default='K',
                            help='which amino acid residue is label carrying [default: K]',
                            action=CheckAminoAcids,
                            )

    parser_fit.add_argument('--kp',
                            help='for two-compartment models, the precursor rate constant [default: 0.5]',
                            type=float,
                            default=0.5)

    parser_fit.add_argument('--kr',
                            help='for the fornasiero model, the reutilization rate constant [default: 0.05]',
                            type=float,
                            default=0.05)

    parser_fit.add_argument('--rp',
                            help='for the fornasiero model, '
                                 'the ratio of protein bound to free precursors [default: 10]',
                            type=float,
                            default=10)

    parser_fit.add_argument('-q', '--q_value',
                            help='fits only peptide data points with q value below this threshold [default: 1e-2]',
                            metavar="FDR[0,1]",
                            type=float,
                            action=CheckQValue,
                            default=1e-2)

    parser_fit.add_argument('-d', '--depth',
                            help='fits only peptides identified in at least this many samples [default: 6]',
                            type=int,
                            default=3)

    parser_fit.add_argument('-r', '--ria',
                            help='final isotope enrichment levels, if known [default: 0.5]',
                            type=float,
                            default=0.5)

    parser_fit.add_argument('-o', '--out', help='path to the output directory [default: riana]',
                            default='riana')

    parser_fit.add_argument('-p', '--plotcurves',
                            action='store_true',
                            help='plot fitted curves')

    parser_fit.add_argument('-t', '--thread',
                            help='number of threads for concurrency [default: 1]',
                            type=int,
                            default=1)

    parser_fit.set_defaults(func=riana_fit.fit_all)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # gc.enable()
    # gc.set_debug(gc.DEBUG_LEAK)

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
