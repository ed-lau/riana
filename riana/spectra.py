# -*- coding: utf-8 -*-

""" Methods to parse individual mzml. """

import pymzml as mz
import numpy as np
import logging


class Mzml(object):

    def __init__(
            self,
            path,
    ):
        """
        This class reads mzml files using pymzml and integrates based on the parsed mzid
        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        """

        self.path = path
        self.msdata = None
        self.rt_idx = None
        self.scan_idx = None

        self.logger = logging.getLogger('riana.mzml')
        self.logger.info('Reading mzML at {0}'.format(self.path))

        self.parse_mzml()

    def __repr__(self):
        """ repr """
        return str(self.path)

    def __str__(self):
        """ repr """
        return str(self.path)

    def parse_mzml(self):
        """
        Read the mzml file and create data dictionary
        :return:
        """

        run = mz.run.Reader(self.path,
                            MS_precision={
                                1: 1e-6,
                                2: 10e-6
                            })

        rt_idx = []
        msdata = []
        scan_numbers = []

        for n, spec in enumerate(run):

            # 2020-07-28: take only MS1 level data
            if spec.ms_level == 1:
                scan_numbers.append(n + 1)
                rt_idx.append(spec.scan_time_in_minutes())
                msdata.append(spec.peaks("centroided"))  # a dict of np.ndarray objects

        self.msdata = msdata
        self.rt_idx = np.array(rt_idx, dtype=np.float64)
        self.scan_idx = np.array(scan_numbers, dtype=np.int64)

        self.logger.info(
            'Parsed {0} spectra from file {1}'.format(
                n + 1,
                self.path,
            )
            )

        return True
