"""

Parse MZML v.0.5.0. Build Date : : :.

Note on test file created by:

    $ msconvert
    hw_h_6_small.raw
    --filter "peakPicking vendor"
    --filter "scanTime [2100,2400]"
    --filter "threshold count 150 most-intense"

    $ crux tide-index
    20171120_Sp_Mm_Canonical_16935e.fasta
    SpMmPartial1miss
    --digestion partial-digest
    --missed-cleavages 1

    $ crux tide-search
    --compute-sp T
    ./hw_h_6_small.mzML.gz
    SpMmPartial1miss
    --precursor-window 20
    --precursor-window-type ppm
    --mz-bin-width 0.02
    --mz-bin-offset 0.0
    --isotope-error 1,2

    $crux percolator
    --protein T
    --fido-empirical-protein-q T
    tide-search.target.txt
    --mzid-output T



"""

import unittest
import pymzml as mz

class Mzml(object):

    def __init__(
            self,
            path
    ):
        """
        This class reads mzml files using pymzml and integrates based on the parsed mzid
        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        """

        self.path = path
        self.msdata = {}
        self.rt_idx = {}
        self.mslvl_idx = {}

    def parse_mzml(self):
        """
        Read the mzml file and create data dictionary
        :return:
        """

        #
        run = mz.run.Reader(self.path,
                            MS_precision={
                                1: 20e-6,
                                2: 20e-6
                            })

        for n, spec in enumerate(run):

            if n % 1000 == 0:
                print(
                    'Loading spectrum {0} at retention time {scan_time:1.2f}'.format(
                        spec.ID,
                        scan_time=spec.scan_time_in_minutes()
                    )
                )

            self.mslvl_idx[n + 1] = spec.ms_level
            self.rt_idx[n + 1] = spec.scan_time_in_minutes()

            if spec.ms_level == 1:
                self.msdata[n + 1] = spec.peaks("centroided")

        print(
            'Parsed {0} spectra from file {1}'.format(
                n + 1,
                self.path)
            )

        return True


class TestParseMzml(unittest.TestCase):
    """Tests for `primes.py`."""

    def test_read_small(self):
        """
        Test whether the small file can be read
        """

        Mzml("_testdata/hw_h_6_small.mzML.gz")
        return


if __name__ == '__main__':
    unittest.main()

