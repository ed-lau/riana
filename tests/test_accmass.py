import unittest

from riana import accmass


class PeptideMassTest(unittest.TestCase):
    """
    Test cases involving calculating peptide masses
    """

    def setUp(self):
        """

        :return:
        """

        pass

    def tearDown(self):
        """

        :return:
        """

        pass


    def test_peptide_mass(self):
        """ Tests several peptides against Expasy/Percolator """

        self.assertAlmostEqual(accmass.calculate_ion_mz('RHPEYAVSVLLR'), 1438.8045, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('TVMENFVAFVDK'), 1398.6853, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('QNQEYQVLLDVRARLEGEI'), 2272.1811, places=3)

    def test_modified_peptides(self):
        """ Tests peptide modification bracket """

        self.assertAlmostEqual(accmass.calculate_ion_mz('RHPEYA[0]VSVLLR'), 1438.8045, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('RHPEYA[0.00]VSVLLR'), 1438.8045, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('RHPEYA[0]VS[0]VLLR'), 1438.8045, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('RHPEYA[0]VSVL[0]LR'), 1438.8045, places=3)
        self.assertAlmostEqual(accmass.calculate_ion_mz('TDEM[15.99]AHFDRERIPER'), 1916.8750, places=3)

        # There should be no character inside the modification
        with self.assertRaises(ValueError):
            accmass.calculate_ion_mz('RHPE[C57]YALR', charge=-1)

    def test_charge(self):
        """ Tests that charge is positive integer """
        with self.assertRaises(ValueError):
            accmass.calculate_ion_mz('RHPEYALR', charge=-1)

        with self.assertRaises(AssertionError):
            accmass.calculate_ion_mz('RHPEYALR', charge=True)


