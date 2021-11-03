import unittest
import os
from tempfile import TemporaryDirectory

from riana.project import ReadDirectory
from riana.peptides import ReadPercolator


class MzmlTest(unittest.TestCase):
    """
    Test cases involving reading mzml files
    """

    def setUp(self):
        """
        Read directory and percolator file and store them in class.
        """

        self.dir_path = 'tests/data'
        self.sample_dir = ReadDirectory(self.dir_path)

        with TemporaryDirectory() as temp_dir:
            self.perc = ReadPercolator(self.sample_dir, temp_dir, percolator_subdirectory=".")

        pass

    def tearDown(self):
        """

        :return:
        """

        pass


    def test_files_exist(self):
        """ Tests that mzML exists. """

        self.assertTrue(os.path.exists(os.path.join(self.dir_path, 'sample1', '20180216_BSA.mzML.gz')))
        self.assertTrue(os.path.exists(os.path.join(self.dir_path, 'sample1', 'percolator.target.psms.txt')))


    def test_that_directory_is_read(self):
        """ Reads in sample directory and assert there is one sample. """

        self.assertEqual(len(self.sample_dir.samples), 1)

    def test_that_percolator_can_be_read(self):
        """ Load the percolator file inside the directory"""

        self.perc.read_all_project_psms()
        self.assertFalse(self.perc.master_id_df['file_idx'].empty)

    def test_percolator_first_sample(self):
        """ Get the first sample from the percolator file and open the mzML"""

        self.perc.read_all_project_psms()

        self.perc.get_current_sample_psms(self.perc.samples[0])
        self.perc.get_current_sample_mzid_indices()

        self.assertEqual(len(self.perc.indices), 1)

