#
"""
The goal of this notebook is to re-implement match between runs for RIANA.
Ideally we want to work with Percolator output files only.
So this function will take in a list of Percolator files (from all time points)
It will compare the identification list of different samples against each other
And then it will append each Percolator output files with new rows.
The new rows will be the peptides that are found in other samples but not in the current sample.
The new rows will be appended to the bottom of the Percolator output file

python -m riana preprocess tests/data/mbr_test/ctrl/percolator.target.psms.txt tests/data/mbr_test/treatment1/percolator.target.psms.txt tests/data/mbr_test/treatment2/percolator.target.psms.txt

"""

# 20211109 mzid.make_master_match_list(  # lysine_filter=0,
#   peptide_q=q_threshold,
#    unique_only=unique_pep,
#    min_fraction=params.min_fraction_mbr)

# This file will do an MBR by:
# 1. Reading in the mzid files
# 2. Creating a master match list
# 3. Creating a master peptide list
# 4. Match retention time of peptides using a spline
# 5. Calculate MBR score

# TODO: should remove mbr from the main integration function and move to a different script

# Each subdirectory is a sample
# 20211109 samples = project.samples
# Create the grand total out file
# 20211109 master_df = pd.DataFrame()

# 20211109 for current_sample in tqdm.tqdm(samples, desc='Processing Sample', total=len(samples))

import logging
import argparse
import pandas as pd
from riana.peptides import ReadPercolator
from riana.logger import get_logger

class MatchBetweenRun(object):
    def __init__(self,
                 logger: logging.Logger,
                 list_of_percolator_files,
                 q_values: float = 0.01,
                 min_fraction: float = 0.35,
                 ):
        """
        This class will take in a list of Percolator output files
        :param list_of_percolator_files:  list of Percolator output files
        :param q_values:  q-value threshold
        :param min_fraction:  minimum fraction of samples that the peptide must be found in
        """

        self.list_of_percolator_files = list_of_percolator_files
        self.q_value_threshold = q_values
        self.min_fraction = min_fraction
        self.overall_df = pd.DataFrame()

        self.logger = logger
        self.make_overall_match_list()

    def make_overall_match_list(self):
            """
            get the overall list of the peptides that are consistently found at a threshold across all samples
            in the project

            :return:
            """
            # self.overall_df = pd.concat([pd.read_csv(f, sep='\t').assign(File=f) for f in self.list_of_percolator_files])
            self.overall_df = pd.concat(
                [ReadPercolator(f, sample=f.name, logger=self.logger) for f in self.list_of_percolator_files])
            # Filter by q-value
            # self.overall_df = self.overall_df[self.overall_df['q-value'] <= self.q_value_threshold]

            print(self.overall_df)
            return None

            pass

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
                     'scan': pred_df.median(axis=1),  # .astype(int),  # Median value of all predicted scans
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



    #
    # # count the number of samples in which the peptide/z is found, also take the most common file_idx
    # self.match_across_runs_master = self.match_across_runs_master.groupby(['concat']) \
    #     .agg({'sample': 'nunique',
    #           'file_idx': 'median',
    #           'spectrum precursor m/z': 'median',
    #           'peptide mass': 'median',
    #           }).reset_index(drop=False)
    #
    # # rename the n_unique summary column into num_samples
    # self.match_across_runs_master.rename(columns={'sample': 'num_samples'}, inplace=True)
    #
    # # force the median fraction number into an integer
    # self.match_across_runs_master['file_idx'] = pd.to_numeric(self.match_across_runs_master['file_idx'],
    #                                                           downcast='integer')
    #
    # # take only those peptide/z that are found in at least min_fraction * total samples
    # self.match_across_runs_master = self.match_across_runs_master.loc[
    #                                 lambda x: x['num_samples'].astype(float) >= len(
    #                                     self.samples) * min_fraction, :]

    # return True

def main(args: argparse.Namespace) -> None:

    logger = get_logger(__name__, out_path=args.out)
    mbr = MatchBetweenRun(logger = logger,
                          list_of_percolator_files=args.id_files,)

    return None

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

"""


# # count the number of samples in which the peptide/z is found, also take the most common file_idx
# self.match_across_runs_master = self.match_across_runs_master.groupby(['concat']) \
#     .agg({'sample': 'nunique',
#           'file_idx': 'median',
#           'spectrum precursor m/z': 'median',
#           'peptide mass': 'median',
#           }).reset_index(drop=False)
#
# # rename the n_unique summary column into num_samples
# self.match_across_runs_master.rename(columns={'sample': 'num_samples'}, inplace=True)
#
# # force the median fraction number into an integer
# self.match_across_runs_master['file_idx'] = pd.to_numeric(self.match_across_runs_master['file_idx'],
#                                                           downcast='integer')
#
# # take only those peptide/z that are found in at least min_fraction * total samples
# self.match_across_runs_master = self.match_across_runs_master.loc[
#                                 lambda x: x['num_samples'].astype(float) >= len(self.samples) * min_fraction, :]
#
# return True
# 20211109 match between run logs
# self.match_logger = logging.getLogger('riana.match_across_run')
# fh = logging.FileHandler(os.path.join(directory_to_write, f'riana_peptides_{self.sample}.log'))
# fh.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# fh.setFormatter(formatter)
# self.match_logger.addHandler(fh)

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