import logging
import os
import unittest

from riana.peptides import ReadPercolator


class MzmlTest(unittest.TestCase):
    """
    Test cases involving reading mzML files and Percolator output.
    """

    @classmethod
    def setUpClass(cls):
        cls.dir_path = os.path.join('tests', 'data', 'sample1')
        cls.psms_path = os.path.join(cls.dir_path, 'percolator.target.psms.txt')

    def setUp(self):
        self._psms_handle = open(self.psms_path, 'r')
        self.addCleanup(self._psms_handle.close)

    def _make_percolator(self):
        return ReadPercolator(
            path=self._psms_handle,
            sample='sample1',
            _ignored_mods=[],
            logger=logging.getLogger('test_mzml'),
        )

    def test_files_exist(self):
        self.assertTrue(os.path.exists(os.path.join(self.dir_path, '20180216_BSA.mzML.gz')))
        self.assertTrue(os.path.exists(self.psms_path))

    def test_percolator_loads_with_expected_columns(self):
        perc = self._make_percolator()
        expected_cols = {
            'file_idx', 'scan', 'charge', 'sequence', 'percolator q-value',
            'peptide mass', 'protein id', 'concat', 'sample',
        }
        self.assertTrue(expected_cols.issubset(perc.id_df.columns))
        self.assertFalse(perc.id_df.empty)

    def test_percolator_indices_populated(self):
        perc = self._make_percolator()
        self.assertGreater(len(perc.indices), 0)
        self.assertIn(perc.indices[0], set(perc.id_df['file_idx']))

    def test_get_current_fraction(self):
        perc = self._make_percolator()
        perc.get_current_fraction_psms(perc.indices[0])
        self.assertFalse(perc.curr_frac_id_df.empty)
        # The pep_id helper column should be set up for downstream integration.
        self.assertIn('pep_id', perc.curr_frac_id_df.columns)
