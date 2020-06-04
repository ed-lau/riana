import unittest
import os
import ftplib

from tqdm import tqdm

class MzmlTest(unittest.TestCase):
    """
    Test cases involving reading mzml files
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
        # target = 'TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML'
        #
        # os.remove(os.path.join('riana', 'tests', '_testdata', target))


    def test_that_directory_exists(self):
        """

        :return:
        """

        pass

    # def test_that_mzml_downloads(self):
    #     """
    #
    #     Write a test for downloading mzml from PXD and ship test percolator file..
    #
    #
    #     :return:
    #     """
    #
    #     ftp = ftplib.FTP("ftp.pride.ebi.ac.uk")
    #     ftp.login(user='', passwd='')
    #     ftp.cwd('/pride/data/archive/2012/03/PXD000001/')
    #
    #     target = 'TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML'
    #
    #     size = ftp.size(target)
    #
    #     host_file = os.path.join('riana', 'tests', '_testdata', target)
    #
    #     try:
    #         with open(host_file, 'wb') as local_file:
    #
    #             with tqdm(total=size,
    #                       unit='B', unit_scale=True, unit_divisor=1024,
    #                       ) as pbar:
    #
    #                 pbar.set_description("Downloading test data from ProteomeXchange")
    #
    #                 def callback(data):
    #                     pbar.update(len(data))
    #                     local_file.write(data)
    #
    #                 ftp.retrbinary('RETR {}'.format(target), callback)
    #
    #     except ftplib.error_perm:
    #         print('error FTP')
    #         pass
    #
    #     ftp.quit()
    #
    #     # return True
    #
    #     pass
