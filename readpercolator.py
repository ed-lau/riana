import pandas as pd
import os

class ReadPercolator(object):
    """
    This class reads the Percolator tab-delimited files for riana. Supports multiple fractions

    """


    def __init__(self, path):
        """
        :param path: path of the folder of the file to be loaded, e.g., "~/Desktop/example/"


        """
        self.path = os.path.join(path)
        self.id_df = self.read_target_psms()
        self.indices = self.get_mzid_indices()

    def read_target_psms(self):
        """
        Reads in the Percolator tab delimited file to return a pandas data frame

        :return:
        """

        # Opening the Percolator tab delimited target.psms.file
        # List all files in the percolator directory ending with target.psms.txt.
        id_files = [f for f in os.listdir(self.path) if f.endswith('target.psms.txt')]

        assert len(id_files) == 1, '[error] check percolator output directory has 1 *.target.psms.txt'

        # Read the Percolator psms file.
        id_df = pd.read_csv(filepath_or_buffer=os.path.join(self.path, id_files[0]),
                            sep='\t')

        return id_df

    def get_mzid_indices(self):
        # Get all the file indices in the Percolator results file.
        return list(set(self.id_df['file_idx']))


    def filter_id_df(self, lysine_filter=0, protein_q=1e-2, peptide_q=1e-2, unique_only=True, require_protein_id=False):
        """
       Filter the percolator id dataframe by:
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
        try:
            self.id_df = self.id_df.loc[lambda x: x['percolator q-value'].astype(float) < peptide_q, :]
            self.id_df = self.id_df.reset_index(drop=True)

        except:
            pass

        #
        # Filter peptides by optional lysine filter
        #

        # If the lysine filter is on and the peptide has no or more than one lysine, skup
        if lysine_filter == 1:
            try:
                self.id_df = self.id_df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) == 1, :]
                self.id_df = self.id_df.reset_index(drop=True)

            except:
                pass

        if lysine_filter == 2:
            try:
                self.id_df = self.id_df.loc[lambda x: x.sequence.apply(lambda y: y.count('K')) > 0, :]
                self.id_df = self.id_df.reset_index(drop=True)

            except:
                pass

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

        # Get only the peptides associated with one and only one proteins
        if unique_only:
            self.id_df = self.id_df[(self.id_df['protein id'].str.count(',') == 0)]
            self.id_df = self.id_df.reset_index(drop=True)
            print('filtered by unique peptide')

        self.id_df = self.id_df.reset_index(drop=True)

        return True

    def subset_id_df(self, idx):
        """
        Return only a fraction

        :param idx:
        :return:
        """

        return self.id_df[self.id_df['file_idx'] == idx]
