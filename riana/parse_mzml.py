# -*- coding: utf-8 -*-

""" Methods to parse individual mzml. """


import pymzml as mz
import numpy as np
import logging


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
        self.msdata = None
        self.rt_idx = None
        self.mslvl_idx = None
        self.scan_idx = None

        self.logger = logging.getLogger('riana.mzml')
        self.logger.info('Reading mzML at {0}'.format(self.path))

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

        #
        run = mz.run.Reader(self.path,
                            MS_precision={
                                1: 20e-6,
                                2: 20e-6
                            })

        mslvl_idx = np.zeros(shape=run.get_spectrum_count(), dtype=np.int)
        rt_idx = np.zeros(shape=run.get_spectrum_count())
        msdata = []
        scan_numbers = np.zeros(shape=run.get_spectrum_count(), dtype=np.int)

        for n, spec in enumerate(run):

            # Check for retention time
            #if n % 1000 == 0:
            #    print(
            #        'Loading spectrum {0} at retention time {scan_time:1.2f}'.format(
            #            spec.ID,
            #            scan_time=spec.scan_time_in_minutes()
            #        )
            #    )

            mslvl_idx[n] = spec.ms_level
            scan_numbers[n] = n+1
            rt_idx[n] = spec.scan_time_in_minutes()

            #if spec.ms_level == 1:
            msdata.append(spec.peaks("centroided"))  # a dict of np.ndarray objects


        self.msdata = msdata
        self.rt_idx = rt_idx
        self.mslvl_idx = mslvl_idx
        self.scan_idx = scan_numbers

        self.logger.info(
            'Parsed {0} spectra from file {1}'.format(
                n + 1,
                self.path)
            )

        return True

