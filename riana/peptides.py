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
                 path,
                 sample,
                 # project,
                 directory_to_write,
                 # percolator_subdirectory,
                 ):
        """
        :param project: path of the input
        :param directory_to_write: path of the output directory


        """

        self.path = path  # 20211109 project.path
        # 20211109 self.samples = project.samples
        # 20211109 self.percolator_subdirectory = percolator_subdirectory
        self.sample = sample

        # logging
        self.logger = logging.getLogger('riana.read_id')

        # 20211109 match between run logs
        # self.match_logger = logging.getLogger('riana.match_across_run')
        # fh = logging.FileHandler(os.path.join(directory_to_write, f'riana_peptides_{self.sample}.log'))
        # fh.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # fh.setFormatter(formatter)
        # self.match_logger.addHandler(fh)

        self.master_id_df = pd.DataFrame()
        self.curr_sample_id_df = pd.DataFrame()
        self.curr_frac_id_df = pd.DataFrame()
        self.curr_frac_filtered_id_df = pd.DataFrame()

        # for match between run function
        self.match_across_runs_master = pd.DataFrame()

        self.indices = []
        self.fraction_id_df = pd.DataFrame()

    def read_psms(self):
        """
        Reads in the Percolator tab delimited file and return a pandas data frame

        :return:
        """

        # 20211109 all_psms = pd.DataFrame()

        # 20211109 for sample in self.samples:

        # 20211109 sample_loc = os.path.join(self.path, sample, self.percolator_subdirectory)
        # 20211109 assert os.path.isdir(sample_loc), '[error] project sample subdirectory not valid'

        self.logger.info('Reading Percolator file at {0}'.format(self.path)) # 20211109 sample_loc))

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
            id_df = pd.read_csv(filepath_or_buffer=self.path,  # 20211109 os.path.join(sample_loc, id_files[0]),
                                sep='\t')
            # Test for Crux Percolator file
            test_col = id_df['spectrum precursor m/z']

            # 2021-12-21 calculate peptide mass because the crux peptide mass column does not include cysteine IAA mass
            id_df['peptide mass'] = [accmass.calculate_ion_mz(seq) for seq in id_df['sequence']]

        except OSError as e:
            sys.exit('Failed to load mzid file. ' + str(e.errno))

        # Try reading the standalone percolator psms file
        except (pd.errors.ParserError, KeyError) as e:
            self.logger.info("Standalone Percolator file detected")

            with open(os.path.join(self.path), 'r') as f:  # 20211109 sample_loc, id_files[0]
                f_ln = f.readlines()

            # The percolator output has different number of columns per row because proteins are separated by tabs
            # Read the first half of the table without the protein IDs
            id_df = pd.DataFrame([ln.split('\t')[0:5] for ln in f_ln[1:]])
            id_df.columns = ['PSMId', 'score', 'percolator q-value', 'percolator PEP', 'peptide']
            id_df['percolator q-value'] = id_df['percolator q-value'].astype(float)
            id_df['percolator PEP'] = id_df['percolator PEP'].astype(float)

            # Create a sequence column for compatibility with Crux percolator
            id_df['sequence'] = [pep[2:-2] for pep in id_df['peptide']]
            # Create a flanking aa column for compatibility with Crux percolator
            id_df['flanking aa'] = [pep[0] + pep[-1] for pep in id_df['peptide']][1]

            # Create dummy columns for compatibility
            id_df['spectrum precursor m/z'] = 0
            id_df['percolator score'] = 0
            id_df['spectrum neutral mass'] = 0
            id_df['distinct matches/spectrum'] = 0
            id_df['peptide mass'] = [accmass.calculate_ion_mz(seq) for seq in id_df['sequence']]

            # Then read in the protein names and join them by comma instead of tab
            id_df['protein id'] = [','.join(ln.rstrip().split('\t')[5:]) for ln in f_ln[1:]]

            # Split the PSMId column to create file_idx, scan, and charge.
            id_df['charge'] = [psm.split('_')[-2] for psm in id_df['PSMId']]
            id_df['charge'] = id_df['charge'].astype(int)
            id_df['scan'] = [psm.split('_')[-3] for psm in id_df['PSMId']]
            id_df['scan'] = id_df['scan'].astype(int)
            # The file name is the underscore('_') split until the last 3 parts, then rejoined by underscore
            # in case there are underscores in the filename. We then remove everything
            # We then remove all directories to get base name
            id_df['file_name'] = [os.path.basename('_'.join(psm.split('_')[:-3])) for psm in id_df['PSMId']]

            # Get the sorted file names, hopefully this is the same index as the Crux Percolator output
            # TODO: Read the Percolator log file to get actual index and use file names to open the mzml instead
            sorted_index = sorted(set(id_df['file_name']))
            id_df['file_idx'] = id_df['file_name'].apply(sorted_index.index)

            id_df = id_df[['file_idx',
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

        id_df.loc[:, 'sample'] = self.sample
        id_df['concat'] = id_df['sequence'].map(str) + '_' + id_df['charge'].map(str)

        self.logger.info('Percolator file for {0} has size {1}'.format(self.sample,
                                                                       id_df.shape))

        # 20211109 all_psms = all_psms.append(id_df, sort=False)
        # 20211109 self.logger.info('Master Percolator file has size {0}'.format(all_psms.shape))

        # 20211109 self.master_id_df = all_psms.copy()
        self.master_id_df = id_df.copy()

        return True

    def make_master_match_list(self,
                               peptide_q=1e-2,
                               # lysine_filter=0,
                               unique_only: bool = True,
                               min_fraction=0.25):
        """
        get the master list of the peptides that are consistently found at a threshold across all samples
        in the project

        :param peptide_q:
        :param min_fraction:
        :return:
        """

        # check if the master peptide ID list has been read first
        if self.master_id_df.empty:
            self.logger.error('Peptide ID list has not been read yet')
            raise Exception

        else:
            self.match_across_runs_master = self.filter_df_by_args(self.master_id_df,
                                                                   peptide_q=peptide_q,
                                                                   # lysine_filter=lysine_filter,
                                                                   unique_only=unique_only)

        # count the number of samples in which the peptide/z is found, also take the most common file_idx
        self.match_across_runs_master = self.match_across_runs_master.groupby(['concat']) \
            .agg({'sample': 'nunique',
                  'file_idx': 'median',
                  'spectrum precursor m/z': 'median',
                  'peptide mass': 'median',
                  }).reset_index(drop=False)

        # rename the n_unique summary column into num_samples
        self.match_across_runs_master.rename(columns={'sample': 'num_samples'}, inplace=True)

        # force the median fraction number into an integer
        self.match_across_runs_master['file_idx'] = pd.to_numeric(self.match_across_runs_master['file_idx'],
                                                                  downcast='integer')

        # take only those peptide/z that are found in at least min_fraction * total samples
        self.match_across_runs_master = self.match_across_runs_master.loc[
                                        lambda x: x['num_samples'].astype(float) >= len(self.samples) * min_fraction, :]

        return True

    def get_current_sample_psms(self,
                                current_sample):
        """
        Get all the PSMs that are associated with the current sample being analyzed
        This is not necessary as of 2021-11-09 since there is only one sample/percolator per run
        :return:
        """

        assert current_sample in list(set(self.master_id_df['sample']))

        self.curr_sample_id_df = self.master_id_df.query('sample == @current_sample')  # .reset_index(drop=True)

        return True

    def get_current_sample_mzid_indices(self):
        """
        Get all the file indices in the current sample's Percolator results file.

        :return:
        """
        self.indices = list(set(self.curr_sample_id_df['file_idx']))

        return True

    def get_current_fraction_psms(self,
                                  idx):
        """
        Return the PSMS only within a fraction

        :param idx:
        :return:
        """

        assert idx in list(set(self.curr_sample_id_df['file_idx']))

        self.curr_frac_id_df = self.curr_sample_id_df.query('file_idx == @idx')

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
