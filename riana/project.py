# -*- coding: utf-8 -*-

""" Methods to read in project directory. """

import os


class ReadDirectory(object):
    """
    This class reads the Project Directory and find the Percolator and mzML files inside.
    This is needed to create a structure in order to do match between run functionality.

    """

    def __init__(self, path):
        """
        :param path: path of the project folder, e.g., "~/Desktop/example/"


        """
        self.path = os.path.join(path)
        self.samples = self._get_sample_list()

    def _get_sample_list(self):
        """
        Get the list of samples and check that each is a directory
        :return: list names of subfolders
        """
        # Excluded hidden folders
        sample_list = [s.name for s in os.scandir(self.path) if s.is_dir()] # [s for s in sorted(os.listdir(self.path)) if not s.startswith('.')]

        for sample in sample_list:
            sample_loc = os.path.join(self.path, sample)
            assert os.path.isdir(sample_loc), '[error] project sample subdirectory not valid'

        return sample_list
