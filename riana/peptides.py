# -*- coding: utf-8 -*-

""" Methods to parse individual peptides from percolator. """

import logging
import os
import re

import pandas as pd

from riana import accmass
from riana.exceptions import DataError


class ReadPercolator(object):
    """
    This class reads the Percolator tab-delimited files for riana. Supports multiple fractions

    """

    def __init__(self,
                 path: str,
                 sample,
                 _ignored_mods,
                 logger: logging.Logger,

                 ):
        """
        :param path: path of the input
        :param sample: sample name

        """

        self.path = path  # 20211109 project.path
        self.logger = logger
        # 20211109 self.samples = project.samples
        # 20211109 self.percolator_subdirectory = percolator_subdirectory
        self.sample = sample
        self.ignored_mods = _ignored_mods

        self.id_df = pd.DataFrame()
        self.curr_sample_id_df = pd.DataFrame()
        self.curr_frac_id_df = pd.DataFrame()
        self.curr_frac_filtered_id_df = pd.DataFrame()

        # for match between run function
        self.match_across_runs_master = pd.DataFrame()

        self.fraction_id_df = pd.DataFrame()

        self.read_psms() # Read the Percolator File
        self.indices = []
        self.get_current_sample_mzid_indices() # Get the indices of the fractions





    def read_psms(self):
        """
        Reads in the Percolator tab delimited file and return a pandas data frame

        :return:
        """

        # 20211109 all_psms = pd.DataFrame()

        # 20211109 for sample in self.samples:

        # 20211109 sample_loc = os.path.join(self.path, sample, self.percolator_subdirectory)
        # 20211109 assert os.path.isdir(sample_loc), '[error] project sample subdirectory not valid'

        self.logger.info('Reading Percolator file at {0}'.format(self.path.name))

        # Detect Crux vs standalone Percolator by header inspection — do not rely
        # on an exception falling through pandas parsing.
        try:
            with open(self.path.name, 'r') as f:
                header_line = f.readline()
        except OSError as e:
            raise DataError(f'Failed to load percolator file {self.path}: {e}') from e

        header_cols = header_line.rstrip('\n').split('\t')
        is_crux = 'spectrum precursor m/z' in header_cols and 'sequence' in header_cols

        if is_crux:
            try:
                self.id_df = pd.read_csv(filepath_or_buffer=self.path, sep='\t')
            except OSError as e:
                raise DataError(f'Failed to load percolator file {self.path}: {e}') from e

            # 2021-12-21 calculate peptide mass because the crux peptide mass column does not include cysteine IAA mass
            self.id_df['peptide mass'] = [
                accmass.calculate_ion_mz(seq, ignored_mods=self.ignored_mods)
                for seq in self.id_df['sequence']
            ]
            self.logger.info(
                f'Crux Percolator file detected via header. '
                f'Loaded {len(self.id_df)} rows.'
            )
        else:
            self.logger.info('Standalone Percolator file detected via header.')
            self._read_standalone_percolator()

        self.id_df.loc[:, 'sample'] = self.sample
        self.id_df['concat'] = self.id_df['sequence'].map(str) + '_' + self.id_df['charge'].map(str)

        self.logger.info('Percolator file for {0} has size {1}'.format(self.sample,
                                                                       self.id_df.shape))

        return True

    def _read_standalone_percolator(self):
        """Parse a standalone Percolator psms.txt (variable-column protein-id rows).

        Expects MSFragger PSMId format ``filename.scan.scan.charge_index``.
        """

        with open(self.path.name, 'r') as f:
            f_ln = f.readlines()

        # The percolator output has different number of columns per row because proteins are separated by tabs.
        # Read the first half of the table without the protein IDs.
        self.id_df = pd.DataFrame([ln.split('\t')[0:5] for ln in f_ln[1:]])
        self.id_df.columns = ['PSMId', 'score', 'percolator q-value', 'percolator PEP', 'peptide']
        self.id_df['percolator q-value'] = self.id_df['percolator q-value'].astype(float)
        self.id_df['percolator PEP'] = self.id_df['percolator PEP'].astype(float)

        # Strip flanking residues from the peptide column: "X.SEQUENCE.Y" -> "SEQUENCE".
        self.id_df['sequence'] = [pep[2:-2] for pep in self.id_df['peptide']]

        # Some MSFragger pep.xml output ships a charge digit in the peptide column
        # before the closing flanker; only strip it when the trailing token is a
        # bare digit (i.e. "1"-"9"), to avoid mangling sequences ending in a
        # modification annotation that happens to be numeric.
        self.id_df['sequence'] = [
            seq[:-1] if (len(seq) >= 2 and seq[-1].isdigit() and not seq[-2].isdigit())
            else seq
            for seq in self.id_df['sequence']
        ]

        self.id_df['flanking aa'] = [pep[0] + pep[-1] for pep in self.id_df['peptide']]

        self.id_df['spectrum precursor m/z'] = 0
        self.id_df['percolator score'] = 0
        self.id_df['spectrum neutral mass'] = 0
        self.id_df['distinct matches/spectrum'] = 0
        self.id_df['peptide mass'] = [
            accmass.calculate_ion_mz(seq, ignored_mods=self.ignored_mods)
            for seq in self.id_df['sequence']
        ]

        self.id_df['protein id'] = [','.join(ln.rstrip().split('\t')[5:]) for ln in f_ln[1:]]

        # PSMId format expected: "filename.scan.scan.charge_index"
        # Parse via single regex against the documented format rather than
        # repeated split/index gymnastics.
        psmid_re = re.compile(r'^(?P<file>.+)\.(?P<scan>\d+)\.\d+\.(?P<charge>\d+)_\d+$')
        parsed = self.id_df['PSMId'].astype(str).str.extract(psmid_re)
        if parsed.isnull().any().any():
            bad = self.id_df.loc[parsed.isnull().any(axis=1), 'PSMId'].head(3).tolist()
            raise DataError(
                'Standalone Percolator PSMId did not match expected MSFragger '
                f"format 'filename.scan.scan.charge_index'. Examples: {bad}"
            )

        self.id_df['charge'] = parsed['charge'].astype(int)
        self.id_df['scan'] = parsed['scan'].astype(int)
        self.id_df['file_name'] = parsed['file'].apply(os.path.basename)

        # TODO: read the Percolator log file to get the actual index assignment
        # and use file names to open the mzml instead.
        sorted_index = sorted(set(self.id_df['file_name']))
        self.id_df['file_idx'] = self.id_df['file_name'].apply(sorted_index.index)

        self.id_df = self.id_df[[
            'file_idx',
            'scan',
            'charge',
            'spectrum precursor m/z',
            'spectrum neutral mass',
            'peptide mass',
            'percolator score',
            'percolator q-value',
            'percolator PEP',
            'distinct matches/spectrum',
            'sequence',
            'protein id',
            'flanking aa',
        ]]



    def get_current_sample_mzid_indices(self):
        """
        Get all the file indices in the current sample's Percolator results file.

        :return:
        """

        self.indices = list(set(self.id_df['file_idx']))

        return True

    def get_current_fraction_psms(self,
                                  idx):
        """
        Return the PSMS only within a fraction

        :param idx:
        :return:
        """

        if idx not in set(self.id_df['file_idx']):
            raise DataError(
                f'Fraction index {idx} not present in Percolator output '
                f'(known indices: {sorted(set(self.id_df["file_idx"]))}).'
            )

        self.curr_frac_id_df = self.id_df.query('file_idx == @idx')

        # TODO: 2021-12-17 we should change the behavior here to get all qualifying scans
        # Remove duplicate sequence/z, keeping the one with the lowest q-value only.
        # self.curr_frac_id_df = self.curr_frac_id_df.sort_values('percolator q-value').drop_duplicates(
        #    subset=['concat'])

        # Arrange the PSM rows by scan number
        self.curr_frac_id_df = self.curr_frac_id_df.sort_values(by='scan').reset_index(drop=True)

        # Add a dummy peptide ID for this fraction only
        self.curr_frac_id_df = self.curr_frac_id_df.assign(pep_id=self.curr_frac_id_df.index)

        return True

    def filter_current_fraction_psms(self,
                                     # lysine_filter=0,
                                     peptide_q: float = 1e-2,
                                     unique_only=False,
                                     # require_protein_id=False,
                                     use_soft_threshold: bool = True,
                                     match_between_runs: bool = False,
                                     ):

        """
       Filter the percolator id data frame by:
        - peptides that belong to any protein identified at a protein Q value
        - peptides below a certain peptide Q value
        - peptides containing certain number of lysines

        # 2017-04-06 Note the require_protein_id flag doesn't work for MSGF+ test files at the moment
        # because it seems the MSGF+ mzID files have no ProteinDetectioNList fields but instead
        # store the protein accessions inside <DBSequence>. Turn the flag off when doing MSGF+.

        :param peptide_q: Peptide-level Q value from command line argument
        :param unique_only: Only doing unique peptides
        :param require_protein_id: Require protein IDs (to filter out some mascot rows with no protein fields)
        :param use_soft_threshold:
        :param match_across_runs:
        :return: True
        """

        # filter by peptide q, unique, and lysine
        self.curr_frac_filtered_id_df = self.filter_df_by_args(self.curr_frac_id_df,
                                                               peptide_q=peptide_q,
                                                               # lysine_filter=lysine_filter,
                                                               unique_only=unique_only)

        self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.reset_index(drop=True)
        self.curr_frac_filtered_id_df['evidence'] = 'q_value'


        return True

    @staticmethod
    def filter_df_by_args(df,
                          peptide_q,
                          # lysine_filter,
                          unique_only):
        """
        # Filter the msater peptide ID list by peptide Q value and peptide uniqueness
        :param df:
        :param peptide_q:
        :param unique_only:
        :return:
        """

        try:
            df = df.loc[lambda x: x['percolator q-value'].astype(float) < peptide_q, :].reset_index(drop=True)

        except KeyError:
            pass

        # if the lysine filter is on and the peptide has no or more than one lysine, skip
        # if lysine_filter == 1:
        #     df = df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) == 1, :].reset_index(drop=True)

        # Get only the peptides associated with one and only one proteins
        if unique_only:
            df = df[(df['protein id'].str.count(',') == 0)].reset_index(drop=True)

        return df
