import pandas as pd
import os


class ReadLipid(object):
    """
    This class reads tab-delimited inclusion lists riana. Supports multiple fractions

    """

    def __init__(self, path):
        """
        :param path: path of the folder of the file to be loaded, e.g., "~/Desktop/example/"


        """
        self.path = os.path.join(path)
        self.id_df = self.read_inclusion()
        self.indices = self.get_mzid_indices()
        self.fraction_id_df = pd.DataFrame()

    def read_inclusion(self):
        """
        Reads in the inclusion tab delimited file to return a pandas data frame

        :return:
        """

        # Opening the Percolator tab delimited target.psms.file
        # List all files in the percolator directory ending with target.psms.txt.
        id_files = [f for f in os.listdir(self.path) if f.endswith('inclusion.txt')]

        assert len(id_files) == 1, '[error] check percolator output directory has 1 *.inclusion.txt'

        # Read the Percolator psms file.
        id_df = pd.read_csv(filepath_or_buffer=os.path.join(self.path, id_files[0]),
                            sep='\t')

        return id_df

    def get_mzid_indices(self):
        # Get all the file indices in the Percolator results file.
        return list(set(self.id_df['file_idx']))


    def subset_id_df(self, idx):
        """
        Return only a fraction, also get only the peptide with the lowest Q

        :param idx:
        :return:
        """

        self.fraction_id_df = self.id_df[self.id_df['file_idx'] == idx].drop_duplicates(
            subset=['id', 'charge'])

        # Arrange the PSM rows by scan number
        self.fraction_id_df = self.fraction_id_df.sort_values(by='rtime').reset_index(drop=True)

        # Add a dummy peptide ID for this fraction only
        self.fraction_id_df = self.fraction_id_df.assign(pep_id=self.fraction_id_df.index)

        return True
