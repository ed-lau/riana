"""

Integrate MZML v.0.1.0. Build Date : : :.
Written by Edward Lau (edward.lau@me.com)


"""

import pymzml as mz
import os.path
import pandas as pd
import scipy.integrate
from time import time

class Mzml(object):
    def __init__(self, path):
        """
        This class reads mzml files using pymzml and integrates based on the parsed mzid

        #######
        Examples for mzml
        #######

        import pymzml as mz

        data = "mzml/mhcc_d01_A10.mzML.gz"

        # Reading the data into a Reader object
        spectra = mz.run.Reader(data)

        # Get the number of spectra
        spectra.getSpectrumCount()

        # Get info of spectra, including ID of each individual spectrum
        spectra.info
        spectra.info['filename']

        # Loop through all the spectra in the Reader object, get their ID
        ms1_list = []
        ms2_list = []

        for spectrum in spectra:

        if spectrum['ms level'] == 1:
            ms1_list.append(spectrum['id'])

        if spectrum['ms level'] == 2:
            ms2_list.append(spectrum['id'])


        # Get all peaks from a spectrum

        spectra[9093].peaks

        # Sort the spectra by intensity
        def getKey(item):
            return item[1]

        sorted_spectrum = sorted(spectra[9093].peaks, reverse=True, key = getKey)

        # To get only the top 150 peaks:
        sorted_spectrum[:150]

        # Write all ms1 file
        run2 = mz.run.Writer(filename = 'write_test.mzML', run=spectra , overwrite = True)

        for spectrum in spectra:

        if spectrum['ms level'] == 1:
            print(spectrum['id'])
            run2.addSpec(spectrum)

        run2.save()




        #######
        #######
        #######


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
        i = 0
        for spectrum in self.msdata:
            i += 1

            # # Print progress every 1000 spectra
            try:
                if i % 1000 == 0:
                    print('Indexing ' + str(i) + ' of ' +
                          str(self.msdata.getSpectrumCount()) + ' spectra (ID: ' + str(spectrum['id']) + ').' )

            except:
                pass

            # Only indexing MS1 scans
            if spectrum['ms level'] == 1:
                index[spectrum['id']] = spectrum['MS:1000016']

        return index


    def get_rt_from_scan(self, peptide_scan):
        """
        Given the scan number, return the retention time

        :param peptide_scan: the peptide scan number
        :return: the retention time
        """

        # Some spectral properties: 'id', 'ms level', 'total ion current'
        # NB: According to the mzml OBO, MS:1000016 is the retention time
        return self.msdata[peptide_scan]['MS:1000016']


    def get_isotopes_from_amrt(self, peptide_am, peptide_rt, z, rt_tolerance, iso_to_do):
        """

        :param peptide_am:
        :param peptide_rt:
        :param z:
        :param rt_tolerance:
        :param iso_to_do:
        :return:
        """

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
            try:
                spectrum = self.msdata[nearbyScan_id]

            except KeyError:
                print("Key not found")
                return []


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