# -*- coding: utf-8 -*-

""" Methods to parse individual peptides from percolator. """

import pandas as pd
import os
import sys
import logging
import numpy as np

from riana import accmass, params

from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.metrics import r2_score


class ReadPercolator(object):
    """
    This class reads the Percolator tab-delimited files for riana. Supports multiple fractions

    """

    def __init__(self,
                 path: str,
                 sample,
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

        self.logger.info('Reading Percolator file at {0}'.format(self.path.name)) # 20211109 sample_loc))

        # Opening the Percolator tab delimited target.psms.file
        # List all files in the percolator directory ending with target.psms.txt.
        # 20211109 id_files = [f for f in os.listdir(sample_loc) if f.endswith('target.psms.txt')]
        # 20211109 assert len(id_files) <= 1, '[error] check percolator output directory has 1 *.target.psms.txt'

        # If there is no percolator.target.psms.txt, read psms.txt
        # 20211109 if len(id_files) == 0:
        # 20211109     id_files = [f for f in os.listdir(sample_loc) if f.endswith('psms.txt') and
        # 20211109                 not f.startswith('decoy')]
        # 20211109    assert len(id_files) == 1, '[error] check percolator output directory has 1 psms.txt'

        try:
            # Read the Percolator psms file.
            self.id_df = pd.read_csv(filepath_or_buffer=self.path,  # 20211109 os.path.join(sample_loc, id_files[0]),
                                sep='\t')

            # Test for Crux Percolator file
            test_col = len(self.id_df['spectrum precursor m/z'])

            # 2021-12-21 calculate peptide mass because the crux peptide mass column does not include cysteine IAA mass
            self.id_df['peptide mass'] = [accmass.calculate_ion_mz(seq) for seq in self.id_df['sequence']]
            self.logger.info(f'Crux Percolator file detected through the presence of spectrum precursor '
                        f'm/z column with length {test_col}')

        except OSError as e:
            sys.exit('Failed to load mzid file. ' + str(e.errno))

        # Try reading the Standalone Percolator psms file
        except (pd.errors.ParserError, KeyError) as e:
            self.logger.info("Standalone Percolator file detected")

            with open(os.path.join(self.path), 'r') as f:  # 20211109 sample_loc, id_files[0]
                f_ln = f.readlines()

            # The percolator output has different number of columns per row because proteins are separated by tabs
            # Read the first half of the table without the protein IDs
            self.id_df = pd.DataFrame([ln.split('\t')[0:5] for ln in f_ln[1:]])
            self.id_df.columns = ['PSMId', 'score', 'percolator q-value', 'percolator PEP', 'peptide']
            self.id_df['percolator q-value'] = self.id_df['percolator q-value'].astype(float)
            self.id_df['percolator PEP'] = self.id_df['percolator PEP'].astype(float)

            # Create a sequence column for compatibility with Crux percolator
            self.id_df['sequence'] = [pep[2:-2] for pep in self.id_df['peptide']]
            # Create a flanking aa column for compatibility with Crux percolator
            self.id_df['flanking aa'] = [pep[0] + pep[-1] for pep in self.id_df['peptide']][1]

            # Create dummy columns for compatibility
            self.id_df['spectrum precursor m/z'] = 0
            self.id_df['percolator score'] = 0
            self.id_df['spectrum neutral mass'] = 0
            self.id_df['distinct matches/spectrum'] = 0
            self.id_df['peptide mass'] = [accmass.calculate_ion_mz(seq) for seq in self.id_df['sequence']]

            # Then read in the protein names and join them by comma instead of tab
            self.id_df['protein id'] = [','.join(ln.rstrip().split('\t')[5:]) for ln in f_ln[1:]]

            # Split the PSMId column to create file_idx, scan, and charge.
            self.id_df['charge'] = [psm.split('_')[-2] for psm in self.id_df['PSMId']]
            self.id_df['charge'] = self.id_df['charge'].astype(int)
            self.id_df['scan'] = [psm.split('_')[-3] for psm in self.id_df['PSMId']]
            self.id_df['scan'] = self.id_df['scan'].astype(int)
            # The file name is the underscore('_') split until the last 3 parts, then rejoined by underscore
            # in case there are underscores in the filename. We then remove everything
            # We then remove all directories to get base name
            self.id_df['file_name'] = [os.path.basename('_'.join(psm.split('_')[:-3])) for psm in self.id_df['PSMId']]

            # Get the sorted file names, hopefully this is the same index as the Crux Percolator output
            # TODO: Read the Percolator log file to get actual index and use file names to open the mzml instead
            sorted_index = sorted(set(self.id_df['file_name']))
            self.id_df['file_idx'] = self.id_df['file_name'].apply(sorted_index.index)

            self.d_df = self.id_df[['file_idx',
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

        self.id_df.loc[:, 'sample'] = self.sample
        self.id_df['concat'] = self.id_df['sequence'].map(str) + '_' + self.id_df['charge'].map(str)

        self.logger.info('Percolator file for {0} has size {1}'.format(self.sample,
                                                                       self.id_df.shape))

        # 20211109 all_psms = all_psms.append(id_df, sort=False)
        # 20211109 self.logger.info('Master Percolator file has size {0}'.format(all_psms.shape))

        return True



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

        assert idx in list(set(self.id_df['file_idx']))

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
