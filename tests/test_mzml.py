import unittest
import os
from tempfile import TemporaryDirectory

from riana.read_directory import ReadDirectory
from riana.read_peptide import ReadPercolator

class MzmlTest(unittest.TestCase):
    """
    Test cases involving reading mzml files
    """

    def setUp(self):
        """



        :return:
        """
        global dir_path

        dir_path = 'tests/data'
        pass

    def tearDown(self):
        """

        :return:
        """

        pass


    def test_that_mzML_exists(self):
        """ Tests that mzML exists. """

        self.assertTrue(os.path.exists(os.path.join(dir_path, 'sample1', '20180216_BSA.mzML.gz')))

    def test_that_perc_exists(self):
        """ Tests that percolator exists. """

        self.assertTrue(os.path.exists(os.path.join(dir_path, 'sample1', 'percolator.target.psms.txt')))


    def test_that_directory_is_read(self):
        """ Reads in sample directory and assert there is one sample. """

        sample_dir = ReadDirectory(dir_path)

        self.assertEqual(len(sample_dir.samples), 1)

    def test_that_percolator_can_be_read(self):
        """ Load the percolator file inside the directory"""

        sample_dir = ReadDirectory(dir_path)

        with TemporaryDirectory() as temp_dir:
            perc = ReadPercolator(sample_dir, temp_dir)

        perc.read_all_project_psms()

        self.assertFalse(perc.master_id_df['file_idx'].empty)
