# -*- coding: utf-8 -*-

""" Methods to parse individual mzml. """


import pymzml as mz
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
        self.msdata = {}
        self.rt_idx = {}
        self.mslvl_idx = {}

        self.logger = logging.getLogger('riana.mzml')
        self.logger.info('Reading mzML at {0}'.format(self.path))

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

            # Check for retention time
            #if n % 1000 == 0:
            #    print(
            #        'Loading spectrum {0} at retention time {scan_time:1.2f}'.format(
            #            spec.ID,
            #            scan_time=spec.scan_time_in_minutes()
            #        )
            #    )

            self.mslvl_idx[n + 1] = spec.ms_level
            self.rt_idx[n + 1] = spec.scan_time_in_minutes()

            if spec.ms_level == 1:
                self.msdata[n + 1] = spec.peaks("centroided")

        self.logger.info(
            'Parsed {0} spectra from file {1}'.format(
                n + 1,
                self.path)
            )

        return True

