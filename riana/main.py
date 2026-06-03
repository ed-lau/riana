# -*- coding: utf-8 -*-

""" Main. """

import os
import argparse
from riana import riana_integrate, riana_fit, __version__
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

class StoreUniqueIgnoredMods(argparse.Action):
    """Checks that the list of arguments contains no duplicates, then stores"""
    def __call__(self, parser, namespace, values: List[float], option_string=None):
        if len(values) > len(set(values)):
            raise argparse.ArgumentError(
                self,
                "You cannot specify the same value multiple times. "
                + f"You provided {values}",
            )
        setattr(namespace, self.dest, values)

class StoreUniqueForcedMods(argparse.Action):
    """Checks that the list of arguments contains no duplicates
    If forced_mods is empty, adds 0 to the list, then stores"""

    def __call__(self, parser, namespace, values: List[float], option_string=None):
        if len(values) > len(set(values)):
            raise argparse.ArgumentError(
                self,
                "You cannot specify the same value multiple times. "
                + f"You provided {values}",
            )
        # Sort the values
        values.sort()
        # If values are empty, add 0. Otherwise, add 0 to the front
        if not values:
            values = [0]
        else:
            values.insert(0, 0) # Add 0 to the front of the list

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

    parser_integrate = subparsers.add_parser('integrate',
                                             help='Integrates isotopomer abundance over retention time',
                                             description='Integrates isotopomer abundance over retention time',
                                             epilog='For more information, see GitHub repository at '
                                                    'https://github.com/ed-lau/riana',
                                            )

    parser_fit = subparsers.add_parser('fit',
                                       help='Fit to kinetic models *under development*')

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
                                  default=[0, 6],
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
                                  help='path to the output directory [default: .]',
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
                                  help='<integer> mass tolerance half-width in ppm for integration '
                                       '(integrate ±N ppm around the theoretical m/z) [default 50 ppm]',
                                  type=int,
                                  choices=range(1, 501),
                                  metavar='[1-500]',
                                  default=50)

    parser_integrate.add_argument('-S', '--smoothing',
                                  help='smoothing window size for integration',
                                  type=int,
                                  choices=range(3, 18, 2),
                                  )

    parser_integrate.add_argument('-D', '--mass_difference',
                                  type=float,
                                  default=1.003354835,
                                  help='mass difference between isotopomers [default: 1.003354835]')

    parser_integrate.add_argument('-X', '--ignored_mods',
                                  help='modification(s) to ignore in the search result for'
                                       'calculating true peptide mass. This must match'
                                       'the exact string in the search engine output'
                                       'e.g., 6.02 for SILAC (default: [])',
                                  default=[],
                                  nargs='+',
                                  type=float,
                                  action=StoreUniqueIgnoredMods,
                                  )

    parser_integrate.add_argument('-F', '--forced_mods',
                                  help='modification(s) to always add to each peptide during integration'
                                       'to create separate clusters of isotopomers'
                                       'this is useful for SILAC experiments, e.g., 6.0201 for SILAC',
                                  default=[0],
                                  nargs='+',
                                  type=float,
                                  action=StoreUniqueForcedMods,
                                  )




    parser_integrate.add_argument('--engine',
                                  choices=['legacy', 'new'],
                                  default='legacy',
                                  help='integration engine. "legacy" (default) is the '
                                       'preserved 0.9.0 path; "new" routes through '
                                       'riana.core.integration.integrate_run (M3 Week 3). '
                                       '--engine new emits the Phase D mass-accuracy '
                                       'columns and a per-fraction drift JSON sidecar; '
                                       '--engine legacy preserves the 0.9.0 schema bit-for-bit.')

    parser_integrate.set_defaults(func=_integrate_dispatch)

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
                            type=int,
                            choices=[1, 2, 3, 4],
                            default=1,
                            help='labeling types (1: Deuterium in vivo, 2: Deuterium in vitro, '
                                 '3: Oxygen-18, 4: Amino acid labeling) '
                                 ' [default: 1]')

    parser_fit.add_argument('-a', '--aa',
                            type=str,
                            default='K',
                            help='which amino acid residue(s) are label carrying, e.g., KR [default: K]',
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

    parser_fit.add_argument('-o', '--out', help='path to the output directory [default: .]',
                            default='.')

    parser_fit.add_argument('-p', '--plotcurves',
                            action='store_true',
                            help='plot fitted curves')

    parser_fit.add_argument('-f', '--fs',
                            type=str,
                            choices=['m0_m1', 'm0_m2', 'm0_m3', 'm0_mA', 'm1_m3', 'm1_m2', 'm1_mA', 'Auto'],
                            default=None,
                            help='calculate fractional synthesis using fine structure isotopomers')

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


def _integrate_dispatch(args: argparse.Namespace) -> None:
    """``riana integrate`` entry point — chooses legacy vs new engine."""
    engine = getattr(args, 'engine', 'legacy')
    if engine == 'new':
        _integrate_new(args)
    else:
        riana_integrate.integrate_all(args)


def _integrate_new(args: argparse.Namespace) -> None:
    """``--engine new`` adapter: build typed inputs, fan out per fraction.

    This is the M3 Week 3 CLI wiring (Phase E). It mirrors
    :func:`riana_integrate.integrate_all`'s per-fraction loop but reads PSMs
    via :func:`io.percolator.read_percolator`, opens mzMLs via
    :class:`io.mzml.IndexedMzML`, and integrates through
    :func:`core.integration.integrate_run`. Output is the same
    ``<sample>_riana.txt`` schema (with the Phase D mass-accuracy columns
    appended) plus a per-fraction ``<sample>.drift.json`` sidecar.

    Week 4 replaces argparse with typer/click and re-homes this adapter to
    ``cli.py``; the function lives here for now because the legacy ``main.py``
    is still the live regression gate.
    """
    import dataclasses
    import json
    import os as _os
    import re
    from pathlib import Path

    from riana.config import IntegrationConfig
    from riana.core.integration import integrate_run
    from riana.io.mzml import IndexedMzML
    from riana.io.percolator import file_indices, fraction_psms, read_percolator
    from riana.logger import get_logger

    logger = get_logger(__name__, args.out)
    logger.info('engine=new (M3 Week 3)')
    logger.info(__version__)

    # Build the IntegrationConfig from argparse defaults. The CLI surface
    # doesn't expose peak_method / baseline_method / smoothing_polyorder /
    # ppm_alert yet — they ride at IntegrationConfig defaults
    # (peak_method='fixed_window'; see commit 3d8c715 for why). When Phase C v2
    # lands a stable detected pipeline, the CLI flags follow.
    config = IntegrationConfig(
        sample=args.sample,
        isotopomers=tuple(args.iso),
        mass_tol_ppm=int(args.mass_tol),
        r_time=float(args.r_time),
        q_value=float(args.q_value),
        unique_only=bool(args.unique),
        write_intensities=bool(args.write_intensities),
        smoothing=args.smoothing,
        mass_difference=float(args.mass_difference),
        ignored_mods=tuple(args.ignored_mods),
        forced_mods=tuple(args.forced_mods),
        threads=int(args.thread),
        out_dir=args.out,
    )

    # Read PSMs once; the per-fraction loop filters by file_idx.
    psms_path = args.id_path
    if hasattr(psms_path, 'name'):
        psms_path = psms_path.name
    if hasattr(args.id_path, 'close'):
        args.id_path.close()
    all_psms = read_percolator(psms_path, sample=args.sample,
                               ignored_mods=tuple(args.ignored_mods))

    # mzML directory layout: mirror the legacy resolution (sort-by-name, accept
    # .mzML or .mzML.gz). The legacy `percolator.log.txt` path is unimplemented
    # under the new engine — Week 3 bench inputs come from a tempdir with a
    # single symlinked mzML per fraction (see run_integrate_v0_9_0), so the
    # sort order is unambiguous.
    mzml_files = sorted(
        f for f in _os.listdir(args.mzml_path)
        if re.match(r'^.*\.mz[Mm][Ll](\.gz)?$', f)
    )
    if not mzml_files:
        raise FileNotFoundError(
            f'No mzML files in {args.mzml_path}'
        )
    indices = file_indices(all_psms)
    if len(mzml_files) != len(indices):
        raise ValueError(
            f'mzML count ({len(mzml_files)}) != distinct file_idx count '
            f'({len(indices)}) in {psms_path}'
        )

    for idx in indices:
        mzml_basename = re.sub(r'\.mz[Mm][Ll](\.gz)?$', '', mzml_files[idx])
        mzml_path = _os.path.join(args.mzml_path, mzml_files[idx])
        logger.info(f'integrating fraction {idx}: {mzml_basename}')

        fraction = fraction_psms(all_psms, idx)
        with IndexedMzML(mzml_path) as mzml:
            df = integrate_run(
                config, fraction, mzml, file_label=mzml_basename
            )

        out_file = Path(args.out) / f'{args.sample}_riana.txt'
        df.to_csv(out_file, sep='\t')

        drift = df.attrs.get('drift_summary')
        if drift is not None:
            drift_path = out_file.with_suffix('.drift.json')
            with drift_path.open('w') as f:
                json.dump(dataclasses.asdict(drift), f, indent=2)
        logger.info(f'wrote {out_file} (+ drift sidecar)')

    logger.info('engine=new: done')
    logger.handlers.clear()
