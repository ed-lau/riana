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

        """20211109 all soft-thresholds and match between runs will be moved 
        #
        # soft_threshold - data frame has peptide/z with q above the q-cutoff but within the soft-threshold
        # and are additionally identified at the q-cutoff in at least 50% samples.
        #
        soft_threshold_q = params.soft_threshold_q

        
        # subset the match across run master df into current fraction index only
        idx = self.curr_frac_id_df['file_idx'][0]
        assert idx in list(set(self.curr_sample_id_df['file_idx']))

        
        match_candidates = self.match_across_runs_master[(self.match_across_runs_master.file_idx == idx)]. \
            reset_index(drop=True).drop(columns=['num_samples', 'file_idx'])

        ## Note we also need to filter the master list by the same selected criteria (unique, lysine)
        ## Might be best to make a new function to do that.

        # if using soft-threshold, get additional peptides
        if use_soft_threshold:

            # get peptides whose q-values are above first q-value threshold but below the soft threshold
            curr_frac_soft_id_df = self.curr_frac_id_df.loc[
                                   lambda x: x['percolator q-value'].astype(float) >= peptide_q, :] \
                .reset_index(drop=True)

            curr_frac_soft_id_df = curr_frac_soft_id_df.loc[
                                   lambda x: x['percolator q-value'].astype(float) < soft_threshold_q, :] \
                .reset_index(drop=True)

            # inner join with total match candidates (identified in certain number of samples)
            curr_frac_soft_id_df = curr_frac_soft_id_df.merge(match_candidates[['concat']],
                                                              on='concat',
                                                              how='inner')

            # append the soft threshold peptides to the filtered id list and label their evidence level
            curr_frac_soft_id_df['evidence'] = 'soft'
            self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.append(curr_frac_soft_id_df, sort=False). \
                reset_index(drop=True)

        if match_between_runs:
            # This should be moved to a separate method or script to first generate a master list based on the
            # full set of Percolator input. The integration function should then integrate everything that is filtered.

            # The Match Across Runs candidates are the peptide/z concats that are identified in multiple
            # samples but are not in the current fraction being considered
            mar_candidates = [p for p in match_candidates['concat'].tolist() if p
                              not in self.curr_frac_filtered_id_df['concat'].tolist()]

            # this currently only works with consistent 2d chromatographic fractions
            # (fraction 0 in one sample corresponds
            # to the same 1st dimension separation to the second sample, etc.)
            pred_df = pd.DataFrame(index=mar_candidates)

            sample_names = list(set(self.master_id_df['sample']))
            current_sample = self.curr_frac_filtered_id_df['sample'][0]

            # loop through each sample
            for i, each_sample in enumerate(sample_names):

                # skip current sample
                if each_sample == current_sample:
                    continue

                # otherwise, get a dataframe of the good peptides in the other sample
                other_sample_df = self.master_id_df.query('sample == @each_sample & file_idx == @idx')
                other_sample_df = other_sample_df.sort_values(by=['percolator q-value']). \
                    drop_duplicates('concat').reset_index(drop=True)

                other_sample_df = self.filter_df_by_args(other_sample_df, peptide_q=peptide_q,
                                                    # lysine_filter=lysine_filter,
                                                    unique_only=unique_only)
                other_sample_df = other_sample_df[['concat', 'scan']]

                compare_df = other_sample_df.merge(self.curr_frac_filtered_id_df[['concat', 'scan']],
                                                   on='concat',
                                                   how='inner')

                compare_df = compare_df.sort_values(by=['scan_x'])

                # Fit regression model
                X = np.array(compare_df.scan_x).reshape(-1, 1)
                y = np.array(compare_df.scan_y)

                regr = BaggingRegressor(DecisionTreeRegressor(max_depth=5),
                                        n_estimators=25, random_state=np.random.RandomState(1))

                regr.fit(X, y)

                # Predict
                y_hat = regr.predict(X)

                # plt.figure()
                # plt.scatter(X, y, c="k", label="training samples")
                # plt.plot(X, y_hat, c="r", label="n_estimators=25", linewidth=1)
                # plt.xlabel("data")
                # plt.ylabel("target")
                # plt.title("Boosted Decision Tree Regression")
                # plt.show()

                # Log the X, y, and y_hat
                self.match_logger.debug('{0} {1} {2} X {3}'.format(current_sample,
                                                                   each_sample,
                                                                   idx,
                                                                   np.array(compare_df.scan_x).tolist()))
                self.match_logger.debug('{0} {1} {2} y {3}'.format(current_sample,
                                                                   each_sample,
                                                                   idx,
                                                                   np.array(compare_df.scan_y).tolist()))
                self.match_logger.debug('{0} {1} {2} y_hat {3}'.format(current_sample,
                                                                       each_sample,
                                                                       idx,
                                                                       np.round(y_hat).tolist()))

                # Log the R2 of the predicted scans
                self.logger.info('Matched chromatograph in same '
                                 'frac in {0}. Bagging regr R2: {1}'.format(each_sample,
                                                                            np.round(r2_score(y, y_hat), 3)))

                # Left join the list of match across runs candidates to their scan to the other sample being considered
                match_df = pd.DataFrame({'concat': mar_candidates})
                match_df = match_df.merge(other_sample_df, how='left', on='concat')
                scan_list = np.array(match_df.scan).reshape(-1, 1)

                # replace with predicted values if there are non nan values
                if len(scan_list[~np.isnan(scan_list)]) > 0:
                    scan_list[~np.isnan(scan_list)] = regr.predict(scan_list[~np.isnan(scan_list)].reshape(-1, 1))
                    pred_df[each_sample] = scan_list

            # remove from the predicted value rows that are completely empty
            # then create a output dataframe mimicking the main percolator output
            pred_df = pred_df.dropna(axis=0, how='all')
            mar_out_df = pd.DataFrame(
                {'file_idx': self.curr_frac_filtered_id_df.file_idx[0],
                 'scan': pred_df.median(axis=1),#.astype(int),  # Median value of all predicted scans
                 'charge': [int(b.split('_')[1]) for b in pred_df.index],
                 'spectrum precursor m/z': np.array(
                     match_candidates.set_index('concat').loc[pred_df.index,]['spectrum precursor m/z']),
                 'spectrum neutral mass': 0,
                 'peptide mass': np.array(match_candidates.set_index('concat').loc[pred_df.index,]['peptide mass']),
                 'percolator score': 1,
                 'percolator q-value': 0,
                 'percolator PEP': 0,
                 'distinct matches/spectrum': 0,
                 'sequence': [str(b.split('_')[0]) for b in pred_df.index],
                 'protein id': 'NA',
                 'flanking aa': 'NA',
                 'sample': current_sample,
                 'concat': pred_df.index,
                 'pep_id': ['matched_' + str(a) for a in range(0, len(pred_df.index))],
                 'evidence': 'match_between_runs'
                 })

            # 2020-07-31 convert scan to nearest integer
            mar_out_df.scan = mar_out_df.scan.astype(np.int)

            self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.append(mar_out_df, sort=False). \
                reset_index(drop=True)
        """
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
