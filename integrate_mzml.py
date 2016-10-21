"""

Integrate MZML v.0.1.0. Build Date : : :.
Written by Edward Lau (edward.lau@me.com)

Classes that concern parsing mzIdentML files and creating summary tables

"""

import pymzml as mz
import os.path
import pandas as pd
import scipy.integrate
from time import time

class Mzml(object):
    def __init__(self, path):
        """

        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        """


        self.path = os.path.join(path)
        self.msdata = mz.run.Reader(self.path, MS1_Precision=20e-6, MSn_Precision=20e-6)
        self.index = self.make_index()



    def make_index(self):
        """

        :return: index: a dictionary ms1 scan number vs. rt
        """

        # Index retention time; turn this into a dictionary please.
        index = {}

        for spectrum in self.msdata:

            # # Print progress every 1000 spectra
            try:
                if spectrum['id'] % 1000 == 0:
                    print('Indexing ' + str(spectrum['id']) + ' of ' +
                          str(self.msdata.getSpectrumCount()) + ' spectra.')

            except:
                pass

            # Only indexing MS1 scans
            if spectrum['ms level'] == 1:
                index[spectrum['id']] = spectrum['MS:1000016']

        return index


    def get_rt_from_scan(self, peptide_scan):

        # Some spectral properties: 'id', 'ms level', 'total ion current'
        # NB: According to the mzml OBO, MS:1000016 is the retention time
        return self.msdata[peptide_scan]['MS:1000016']


    def get_isotopes_from_amrt(self, peptide_am, peptide_rt, z, rt_tolerance, iso_to_do):

        if self.index == {}:
            print('No index found: creating new index.')
            self.make_index()

        print(peptide_am, peptide_rt,)

        timeDependentIntensities = []

        # Choose the scan numbers from the index
        nearbyScans = []
        for scan_id, scan_rt in self.index.items():
            if abs(scan_rt - peptide_rt) <= rt_tolerance:
                nearbyScans.append([scan_id, scan_rt])

        t1 = time()
        print('Extracting intensities from spectra...')

        # Loop through each spectrum, check if it is an MS1 spectrum, check if it is within 1 minute of retention time
        for nearbyScan_id, nearbyScan_rt in nearbyScans:
            spectrum = self.msdata[nearbyScan_id]

            #Loop through every isotope in the to-do list
            for i in iso_to_do:

                matchList = spectrum.hasPeak(peptide_am + (i*1.003/z))

                if matchList:
                    for mz, I in matchList:
                        timeDependentIntensities.append([nearbyScan_rt, i, I, mz])

        t2 = time()
        print('Done. Extracting time: ' + str(round(t2 - t1, 2)) + ' seconds.')

        #print(timeDependentIntensities)

        print('Integrating...')
        # Integrate the individual isotopomers
        allIso = []

        for j in iso_to_do:

            isotopomer_profile = []

            for scan, i, I, mz in timeDependentIntensities:
                if i == j:
                    isotopomer_profile.append([scan, i, I, mz])

            # If there is no isotopomer profile, set area to 0
            if isotopomer_profile:
                iso_df = pd.DataFrame(isotopomer_profile)
                iso_area = scipy.integrate.trapz(iso_df[2], iso_df[0])
                # Remove all negative areas
                iso_area = max(iso_area, 0)

            else:
                iso_area = 0

            allIso.append([j, iso_area])

        return allIso

#
#   For doctest
#
if __name__ == '__main__':
    import doctest
    doctest.testmod()