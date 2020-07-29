# -*- coding: utf-8 -*-

""" Methods to parse individual peptides from percolator. """

import pandas as pd
import os
import sys
import logging
import numpy as np

from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.metrics import r2_score


class ReadPercolator(object):
    """
    This class reads the Percolator tab-delimited files for riana. Supports multiple fractions

    """

    def __init__(self, project, directory_to_write,
                 ):
        """
        :param path: path of the folder of the file to be loaded, e.g., "~/Desktop/example/"


        """

        self.path = project.path
        self.samples = project.samples

        # Logging

        self.logger = logging.getLogger('riana.read_id')

        self.match_logger = logging.getLogger('riana.match_across_run')
        fh = logging.FileHandler(os.path.join(directory_to_write, 'riana_match_across_run.log'))
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self.match_logger.addHandler(fh)

        self.master_id_df = pd.DataFrame()
        self.curr_sample_id_df = pd.DataFrame()
        self.curr_frac_id_df = pd.DataFrame()
        self.curr_frac_filtered_id_df = pd.DataFrame()

        # For match between run function
        self.match_across_runs_master = pd.DataFrame()

        self.indices = []
        self.fraction_id_df = pd.DataFrame()

    def read_all_project_psms(self):
        """
        Reads in the all the Percolator tab delimited files in a project and return a pandas data frame

        :return:
        """

        all_psms = pd.DataFrame()

        for sample in self.samples:

            sample_loc = os.path.join(self.path, sample)
            assert os.path.isdir(sample_loc), '[error] project sample subdirectory not valid'

            self.logger.info('Reading Percolator file at {0}'.format(sample_loc))

            # Opening the Percolator tab delimited target.psms.file
            # List all files in the percolator directory ending with target.psms.txt.
            id_files = [f for f in os.listdir(sample_loc) if f.endswith('target.psms.txt')]
            assert len(id_files) == 1, '[error] check percolator output directory has 1 *.target.psms.txt'

            try:
                # Read the Percolator psms file.
                id_df = pd.read_csv(filepath_or_buffer=os.path.join(sample_loc, id_files[0]),
                                    sep='\t')

            except OSError as e:
                sys.exit('Failed to load mzid file. ' + str(e.errno))

            id_df.loc[:, 'sample'] = sample
            id_df['concat'] = id_df['sequence'].map(str) + '_' + id_df['charge'].map(str)

            self.logger.info('Percolator file for {0} has size {1}'.format(sample,
                                                                            id_df.shape))

            all_psms = all_psms.append(id_df, sort=False)

            self.logger.info('Master Percolator file has size {0}'.format(all_psms.shape))

        self.master_id_df = all_psms.copy()

        return True

    def make_master_match_list(self,
                               peptide_q=1e-2,
                               lysine_filter=0,
                               unique_only=True,
                               min_fraction=0.25):
        """
        Get a list of the peptides that are consistently found at a threshold across all samples
        in the project

        :param peptide_q:
        :param min_fraction:
        :return:
        """

        # Check if the master peptide ID list has been read first
        if self.master_id_df.empty:
            self.logger.error('Peptide ID list has not been read yet')
            sys.exit()

        elif not self.master_id_df.empty:
            self.match_across_runs_master = filter_df_by_args(self.master_id_df,
                                                              peptide_q=peptide_q,
                                                              lysine_filter=lysine_filter,
                                                              unique_only=unique_only)

        # Count the number of samples in which the peptide/z is found, also take the most common file_idx
        self.match_across_runs_master = self.match_across_runs_master.groupby(['concat']) \
            .agg({'sample': 'nunique',
                  'file_idx': 'median',
                  'spectrum precursor m/z': 'median',
                  'peptide mass': 'median',
                  }).reset_index(drop=False)

        # Rename the n_unique summary column into num_samples
        self.match_across_runs_master.rename(columns={'sample': 'num_samples'}, inplace=True)

        # Force the median fraction number into an integer
        self.match_across_runs_master['file_idx'] = pd.to_numeric(self.match_across_runs_master['file_idx'],
                                                                  downcast='integer')

        # Take only those peptide/z that are found in at least min_fraction * total samples
        self.match_across_runs_master = self.match_across_runs_master.loc[
                                        lambda x: x['num_samples'].astype(float) >= len(self.samples) * min_fraction, :]

        return True

    def get_current_sample_psms(self,
                                current_sample):
        """
        Get all the PSMs that are associated with the current sample being analyzed
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

        # Remove duplicate sequence/z, keeping the one with the lowest q-value only.
        self.curr_frac_id_df = self.curr_frac_id_df.sort_values('percolator q-value').drop_duplicates(
            subset=['concat'])

        # Arrange the PSM rows by scan number
        self.curr_frac_id_df = self.curr_frac_id_df.sort_values(by='scan').reset_index(drop=True)

        # Add a dummy peptide ID for this fraction only
        self.curr_frac_id_df = self.curr_frac_id_df.assign(pep_id=self.curr_frac_id_df.index)

        return True

    def filter_current_fraction_psms(self,
                                     lysine_filter=0,
                                     protein_q=1e-2,
                                     peptide_q=1e-2,
                                     unique_only=False,
                                     # require_protein_id=False,
                                     use_soft_threshold=True,
                                     match_across_runs=False,
                                     ):

        """
       Filter the percolator id data frame by:
        - peptides that belong to any protein identified at a protein Q value
        - peptides below a certain peptide Q value
        - peptides containing certain number of lysines

        # 2017-04-06 Note the require_protein_id flag doesn't work for MSGF+ test files at the moment
        # because it seems the MSGF+ mzID files have no ProteinDetectioNList fields but instead
        # store the protein accessions inside <DBSequence>. Turn the flag off when doing MSGF+.

        :param lysine_filter: Lysine filter from command line argument
        :param protein_q: Protein-level Q value from command line argument
        :param peptide_q: Peptide-level Q value from command line argument
        :param unique_only: Only doing unique peptides
        :param require_protein_id: Require protein IDs (to filter out some mascot rows with no protein fields)
        :param use_soft_threshold:
        :param match_across_runs:
        :return: True
        """

        #
        # Filter peptides by Protein Q value. To do. Have to read the target.protein.txt
        #
        """
        try:
            self.filtered_protein_df = self.protein_df.loc[lambda x: x.percolator_Q_value.astype(float) < protein_q, :]
            #self.filtered_protein_df = self.filtered_protein_df.reset_index()
        except:
            print('No filtering by protein Q value done.')
            self.filtered_protein_df = self.protein_df
        """

        #
        # Filter peptides by peptide-level Q-value filter
        #

        #
        # Filter peptides by optional lysine filter
        #

        # 2013-04-06 Again I forgot why we only chose the first five columns here. Resetting to all columns for now.
        """
        if require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='inner') # NB: 2017-04-05 might need 'inner' here for riana to work
        elif not require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='left')  # NB: 2017-04-05 might need 'inner' here for riana to work
        """

        #
        # Get the protein Uniprot accession via PE and then via DBS
        #
        """
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.pe_df, how='left')
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.dbs_df, how='left')
        self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)
        """

        self.curr_frac_filtered_id_df = filter_df_by_args(self.curr_frac_id_df,
                                                          peptide_q=peptide_q,
                                                          lysine_filter=lysine_filter,
                                                          unique_only=unique_only)

        self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.reset_index(drop=True)
        self.curr_frac_filtered_id_df['evidence'] = 'q_value'

        #
        # soft_threshold - data frame has peptide/z with q above the q-cutoff but within the soft-threshold
        # and are additionally identified at the q-cutoff in at least 50% samples.
        #

        soft_threshold_q = 0.25

        # Subset the match across run master df into current fraction only
        idx = self.curr_frac_id_df['file_idx'][0]
        assert idx in list(set(self.curr_sample_id_df['file_idx']))

        match_candidates = self.match_across_runs_master[(self.match_across_runs_master.file_idx == idx)]. \
            reset_index(drop=True).drop(columns=['num_samples', 'file_idx'])

        ## Note we also need to filter the master list by the same selected criteria (unique, lysine)
        ## Might be best to make a new function to do that.

        # If using soft-threshold, get additional peptides
        if use_soft_threshold:
            #
            curr_frac_soft_id_df = self.curr_frac_id_df.loc[
                                   lambda x: x['percolator q-value'].astype(float) >= peptide_q, :] \
                .reset_index(drop=True)

            curr_frac_soft_id_df = curr_frac_soft_id_df.loc[
                                   lambda x: x['percolator q-value'].astype(float) < soft_threshold_q, :] \
                .reset_index(drop=True)

            # Do an inner join
            curr_frac_soft_id_df = curr_frac_soft_id_df.merge(match_candidates[['concat']],
                                                              on='concat',
                                                              how='inner')

            curr_frac_soft_id_df['evidence'] = 'soft'

            self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.append(curr_frac_soft_id_df, sort=False). \
                reset_index(drop=True)

        if match_across_runs:

            # The Match Across Runs candidates are the peptide/z concats that are identified in multiple
            # samples but are not in the current fraction being considered
            mar_candidates = [p for p in match_candidates['concat'].tolist() if p
                              not in self.curr_frac_filtered_id_df['concat'].tolist()]

            # This currently only works if you have consistent chromatographic fractions
            # (fraction 0 in one sample corresponds
            # to the same 1st dimension separation to the second sample, etc.)

            pred_df = pd.DataFrame(index=mar_candidates)

            sample_names = list(set(self.master_id_df['sample']))
            current_sample = self.curr_frac_filtered_id_df['sample'][0]

            # Loop through each sample
            for i, each_sample in enumerate(sample_names):

                # Skip current sample
                if each_sample == current_sample:
                    continue

                # Otherwise, get a data frame of good peptides in the other sample
                other_sample_df = self.master_id_df.query('sample == @each_sample & file_idx == @idx')
                other_sample_df = other_sample_df.sort_values(by=['percolator q-value']). \
                    drop_duplicates('concat').reset_index(drop=True)

                other_sample_df = filter_df_by_args(other_sample_df, peptide_q=peptide_q,
                                                    lysine_filter=lysine_filter,
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
                                                                       idx,np.round(y_hat).tolist()))

                # Log the R2 of the predicted scans
                self.logger.info('Matched chromatogr. in same '
                                 'frac in {0}. Bagging regr R2: {1}'.format(each_sample,
                                                                            np.round(r2_score(y, y_hat), 3)))

                # Left join the list of match across runs candidates to their scan to the other sample being considered
                match_df = pd.DataFrame({'concat': mar_candidates})
                match_df = match_df.merge(other_sample_df, how='left', on='concat')
                scan_list = np.array(match_df.scan).reshape(-1, 1)

                # Replace with predicted values
                scan_list[~np.isnan(scan_list)] = regr.predict(scan_list[~np.isnan(scan_list)].reshape(-1, 1))
                pred_df[each_sample] = scan_list

            # TODO: This gives an astype error (NaN for all scans?)
            mar_out_df = pd.DataFrame(
                {'file_idx': self.curr_frac_filtered_id_df.file_idx[0],
                 'scan': pred_df.median(axis=1).astype(int),  # Median value of all predicted scans
                 'charge': [int(b.split('_')[1]) for b in mar_candidates],
                 'spectrum precursor m/z': np.array(
                     match_candidates.set_index('concat').loc[mar_candidates,]['spectrum precursor m/z']),
                 'spectrum neutral mass': 0,
                 'peptide mass': np.array(match_candidates.set_index('concat').loc[mar_candidates,]['peptide mass']),
                 'percolator score': 1,
                 'percolator q-value': 0,
                 'percolator PEP': 0,
                 'distinct matches/spectrum': 0,
                 'sequence': [str(b.split('_')[0]) for b in mar_candidates],
                 'protein id': 'NA',
                 'flanking aa': 'NA',
                 'sample': current_sample,
                 'concat': mar_candidates,
                 'pep_id': ['matched_' + str(a) for a in range(0, len(mar_candidates))],
                 'evidence': 'match_across_runs'
                 })

            self.curr_frac_filtered_id_df = self.curr_frac_filtered_id_df.append(mar_out_df, sort=False). \
                reset_index(drop=True)

        return True


def filter_df_by_args(df,
                      peptide_q,
                      lysine_filter,
                      unique_only):
    """
    # Filter the msater peptide ID list by peptide Q value
    :param df:
    :param peptide_q:
    :param lysine_filter:
    :param unique_only:
    :return:
    """

    try:
        df = df.loc[lambda x: x['percolator q-value'].astype(float) < peptide_q, :].reset_index(drop=True)

    except KeyError:
        pass

        # If the lysine filter is on and the peptide has no or more than one lysine, skip
    if lysine_filter == 1:
        df = df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) == 1, :].reset_index(drop=True)

    elif lysine_filter == 2:
        df = df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) > 0, :].reset_index(drop=True)

    elif lysine_filter == 3:
        df = df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) == 2, :].reset_index(drop=True)

    # Get only the peptides associated with one and only one proteins
    if unique_only:
        df = df[(df['protein id'].str.count(',') == 0)].reset_index(drop=True)

    return df
