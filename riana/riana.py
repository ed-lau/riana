# -*- coding: utf-8 -*-

""" Main. """

from riana import integrate,fitcurve, __version__


#
# Code for running main with parsed arguments from command line
#
def main():
    import argparse

    parser = argparse.ArgumentParser(description='Riana integrates the relative abundance of'
                                                 'isotopomers')

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    # Sub-commands
    subparsers = parser.add_subparsers(help='Type riana function -h for individual help messages',
                                       title='Functions',
                                       description='Riana has the following sub-commands:',
                                       )
    parser_integrate = subparsers.add_parser('integrate',
                                             aliases=['int'],
                                             help='Integrates isotopomer abundance over retention time')
    parser_fit = subparsers.add_parser('fit',
                                       help='Fit to kinetic models. Note implemented yet.')

    # Arguments for integrate subcommand
    parser_integrate.add_argument('dir',
                                  type=str,
                                  help='path to folder containing the mzml files')

    parser_integrate.add_argument('-p', '--percolator',
                                  help='path to the percolator output psms.txt file',
                                  type=str,
                                  default='percolator')

    parser_integrate.add_argument('-s', '--sample',
                                  help='sample name to override mzml folder name, must include numbers, e.g., time1',
                                  type=str,
                                  default=None)

    parser_integrate.add_argument('-i', '--iso',
                                  help='isotopes to do, separated by commas, e.g., 0,1,2,3,4,5 [default: 0,6]',
                                  default='0,6')

    parser_integrate.add_argument('-d', '--deuterium',
                                  action='store_true',
                                  help='experimental feature: use mass defect for deuterium.')

    parser_integrate.add_argument('-u', '--unique',
                                  action='store_true',
                                  help='integrate unique peptides only')

    parser_integrate.add_argument('-t', '--thread',
                                  help='number of threads for concurrency; leave as 0 for auto (default = 0)',
                                  type=int,
                                  default=0)

    parser_integrate.add_argument('-o', '--out', help='path to the output directory [default: riana]',
                                  default='riana')

    parser_integrate.add_argument('-q', '--qvalue',
                                  help='integrate only peptides with q value below this threshold[default: 1e-2]',
                                  type=float,
                                  default=1e-2)

    parser_integrate.add_argument('-r', '--rtime',
                                  help='retention time (in minutes, both directions) tolerance for integration',
                                  type=float,
                                  default=1.0)

    parser_integrate.add_argument('-m', '--mass_tol',
                                  help='mass tolerance in ppm for integration [default 50 ppm]',
                                  type=float,
                                  default=50)

    parser_integrate.set_defaults(func=integrate.integrate_all)

    # Arguments for fit subcommand
    parser_fit.set_defaults(func=fitcurve.runfit)

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
