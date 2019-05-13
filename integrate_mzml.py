"""

Integrate MZML v.0.1.0. Build Date : : :.
Written by Edward Lau (edward.lau@me.com) 2016-2017


"""

import pymzml as mz
import pandas as pd
import scipy.integrate
import gc
from pathos.multiprocessing import ProcessingPool as Pool

class Mzml(object):
    def __init__(self, path):
        """
        This class reads mzml files using pymzml and integrates based on the parsed mzid
        :param path: path of the mzml file to be loaded, e.g., "~/Desktop/example.mzml"
        """

        self.path = path
        self.msdata = {}
        self.rt_idx = {}
        self.mslvl_idx = {}
        self.id = pd.DataFrame()
        self.iso_to_do = []
        self.rt_tolerance = 30
        self.njobs = 10
        self.result = []


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

            if n %% 100 == 0:
                print(
                    'Spectrum {0}, MS level {ms_level} @ RT {scan_time:1.2f}'.format(
                        spec.ID,
                        ms_level=spec.ms_level,
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



    def set_iso_to_do(self, iso_to_do):
        """
        Setter
        :param iso_to_do:
        :return:
        """

        self.iso_to_do = iso_to_do


    def set_rt_tolerance(self, rt_tolerance):
        """
        Setter
        :param rt_tolerance:
        :return:
        """

        self.rt_tolerance = rt_tolerance

    def associate_id(self, id_df):
        """
        Link the id file to this mzml

        :param id_df:
        :return:
        """

        self.id = id_df

    # def make_index(self):
    #     """
    #     Generate two indices:
    #     MS1 index: a dictionary of ms1 scan number vs. rt
    #     MS2 index: a dictionary of ms2 scan number vs. rt
    #
    #
    #     :return: True
    #     """
    #
    #     # Index retention time; turn this into a dictionary please.
    #     i = 0
    #     for spectrum in self.msdata:
    #         i += 1
    #
    #         # # Print progress every 1000 spectra
    #         try:
    #             if i % 1000 == 0:
    #                 print('Indexing ' + str(i) + ' of ' +
    #                       str(self.msdata.getSpectrumCount()) + ' spectra (ID: ' + str(spectrum['id']) + ').' )
    #
    #         except:
    #             pass
    #
    #         # Only indexing MS1 and MS2 scans
    #         if spectrum['ms level'] == 1:
    #             self.ms1_index[spectrum['id']] = spectrum['MS:1000016']
    #         if spectrum['ms level'] == 2:
    #             self.ms2_index[spectrum['id']] = spectrum['MS:1000016']
    #
    #     return True


    # def get_rt_from_scan(self, peptide_scan):
    #     """
    #     For the deprecated integrate function
    #     Given the scan number, return the retention time
    #
    #     :param peptide_scan: the peptide scan number
    #     :return: the retention time
    #     """
    #
    #     # Some spectral properties: 'id', 'ms level', 'total ion current'
    #     # NB: According to the mzml OBO, MS:1000016 is the retention time
    #     return self.msdata[peptide_scan]['MS:1000016']

    # def get_scans_to_do(self, peptide_scan, rt_tolerance):
    #     """
    #     For the new integrate_fast function
    #     Given the scan number, return all the scan IDs to integrate
    #
    #     :param peptide_scan:    MS2 scan number
    #     :param rt_tolerance:    Retention time tolerance in min
    #     :return: the ID of the scans to be integrated
    #     """
    #
    #     peptide_rt = self.ms2_index[peptide_scan]
    #
    #     # if self.ms2_index == {}:
    #     #     print('No index found: creating new index.')
    #     #     self.make_index()
    #
    #     # Choose the scan numbers from the index
    #     nearbyScans = []
    #     for scan_id, scan_rt in self.ms1_index.items():
    #         if abs(scan_rt - peptide_rt) <= rt_tolerance:
    #             nearbyScans.append([scan_id, scan_rt])
    #
    #     return nearbyScans

    #N = 100
    #pbar = tqdm(total=N)
    #res = [None] * N  # result list of correct size

    #def update_bar(self):
        # note: input comes from async `wrapMyFunc`
     #   self.res[i] = ans  # put answer into correct index of result list
      #  self.pbar.update()

    def get_isotopes_from_scan_id_multiwrapper(self):

        self.map = Pool().map
        print("multiwrapper")
        result = self.map(self.multithread_test, range(10))
        return result

    def multithread_test(self, index):
        return index*2

    def get_isotope_from_scan_id_wrapper(self, index):
        """
        Wrapper for the get_isotope_from_scan_id() function below

        :param in_df:
        :return:
        """
        #print("wrapper", index)
        #x = self.get_isotope_from_scan_id(peptide_am=float(self.id.loc[index, 'calc_mz']),
        #                                     z=float(self.id.loc[index, 'z']),
        #                                     spectrum_id=self.id.loc[index, 'scan_id'],
        #                                     iso_to_do=self.iso_to_do)

        self.result = self.get_isotopes_from_amrt_2(peptide_am=float(self.id.loc[index, 'peptide mass']), # spectrum precursor m/z'
                                        peptide_scan=int(self.id.loc[index, 'scan']),
                                        z=float(self.id.loc[index, 'charge']),
                                        rt_tol=self.rt_tolerance,
                                        iso_to_do=self.iso_to_do)

        #gc.collect()

        return self.result

    # def get_isotope_from_scan_id(self, peptide_am, z, spectrum_id, iso_to_do):
    #     """
    #     For the new integrate_fast function, get isotope intensities of a scan
    #     given a peptide m/z and RT combination
    #     :param peptide_am:  Peptide accurate mass
    #     :param z:           Peptide charge
    #     :param spectrum_id: Scan number?
    #     :param iso_to_do:   List of isotopomers to integrate
    #     :return:
    #     """
    #
    #     timeDependentIntensities = []
    #
    #     # Get the spectrum based on the spectrum number
    #     try:
    #         spectrum = self.msdata[spectrum_id]
    #
    #     except KeyError:
    #
    #         print('[error] spectrum index out of bound')
    #         return []
    #
    #     # 2018-09-07 Need to catch a number of errors of XML tree not
    #     # Being able to read the spectrum object returned by pymzml
    #
    #     except xml.etree.ElementTree.ParseError:
    #
    #         print('[warning] XML eTree does not appear to be able to read this spectrum',
    #                   '(scan number:', str(spectrum_id) + ')', sep=' ')
    #
    #         return []
    #
    #     # If it still fails, use the non string version and hope it doesn't freeze
    #     except BaseException as err:
    #         print(err)
    #         spectrum = self.msdata[spectrum_id]
    #
    #     # 2019-05-08 the new pymzml has a weird bug where it would crash when loading some spectra if I specify
    #     # the spectrum number as an int; this problem goes away if the identifier is a str but sometimes it retrieves
    #     # the wrong spectrum, e.g., for the liverpool test if I retrieve '4698' it gets the spectrum 14698 instead
    #     # Since this seems to happen less frequently than the freezing scenario I will use the str for now
    #
    #     if not spectrum.ID == spectrum_id:
    #         print('[warning] pyMZML failed to retrieve correct spectrum',
    #               '(scan number:', str(spectrum_id) + ')', sep=' ')
    #
    #         return []
    #
    #     assert spectrum.ms_level == 1, '[error] specified spectrum is not a parent ion scan'
    #
    #     #Loop through every isotope in the to-do list
    #     for i in iso_to_do:
    #
    #         iso_mz = peptide_am + ((i * 1.003) / z)
    #
    #         matchList = spectrum.has_peak(iso_mz)
    #
    #         if matchList:
    #             for mz, I in matchList:
    #                 timeDependentIntensities.append([spectrum_id, i, I, mz])
    #         else:
    #             timeDependentIntensities.append([spectrum_id, i, 0, iso_mz])
    #
    #
    #     return timeDependentIntensities


    def get_isotopes_from_amrt_2(self, peptide_am, peptide_scan, z, rt_tol, iso_to_do):

        peptide_rt = self.rt_idx.get(peptide_scan)
        # Calculate precursor mass from
        peptide_prec = (peptide_am + (z * 1.007825)) / z
        #peptide_prec = peptide_am

        timeDependentIntensities = []

        # Choose the scan numbers from the index; one-line list comprehension?
        nearby_scans = [[i, rt] for i, rt in self.rt_idx.items()
                        if abs(rt - peptide_rt) <= rt_tol and self.mslvl_idx[i] == 1]

        # Loop through each spectrum, check if it is an MS1 spectrum, check if it is within 1 minute of retention time
        for nearbyScan_id, nearbyScan_rt in nearby_scans:
            # Get the spectrum based on the spectrum number

            for iso in iso_to_do:

                peptide_prec_isotopomer_am = peptide_prec + (iso * 1.007825 / z)
                upper = peptide_prec_isotopomer_am + (peptide_prec_isotopomer_am * 50e-6)
                lower = peptide_prec_isotopomer_am - (peptide_prec_isotopomer_am * 50e-6)

                matching_int = sum([I for mz_value, I in self.msdata.get(nearbyScan_id) if upper > mz_value > lower])

                timeDependentIntensities.append([nearbyScan_rt, iso, matching_int, peptide_prec_isotopomer_am])


            # gc.collect()

        # Integrate the individual isotopomers
        allIso = []

        for j in iso_to_do:

            isotopomer_profile = [[scan, i, I, mz_value] for scan, i, I, mz_value in timeDependentIntensities if i == j]

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

    def get_isotopes_from_amrt(self, peptide_am, peptide_scan, z, rt_tol, iso_to_do):
        """
         For the deprecated integrate function, get isotope intensities of all relevant scans
         given a peptide m/z and RT combination, then integrate over time

        :param peptide_am:
        :param peptide_rt:
        :param z:
        :param rt_tolerance:
        :param iso_to_do:
        :return:
        """

        # if self.ms2_index == {}:
        #     print('No index found: creating new index.')
        #     self.make_index()

        #print(peptide_am, peptide_scan,)

        peptide_rt = self.rt_idx.get(peptide_scan)
        # Calculate precursor mass from
        peptide_prec = (peptide_am + z*1.007825)/z

        timeDependentIntensities = []

        # Choose the scan numbers from the index; one-line list comprehension?
        nearby_scans = [[i, rt] for i, rt in self.rt_idx.items()
                        if abs(rt - peptide_rt) <= rt_tol and self.mslvl_idx[i] == 1]


        # Loop through each spectrum, check if it is an MS1 spectrum, check if it is within 1 minute of retention time
        for nearbyScan_id, nearbyScan_rt in nearby_scans:
            # Get the spectrum based on the spectrum number
            try:
                spectrum_peaks = self.msdata.get(nearbyScan_id)
                spectrum = mz.spec.Spectrum(measured_precision=20e-6)
                spectrum.measured_precision = 20e-6
                spectrum.set_peaks(spectrum_peaks, "raw")
                spectrum.set_peaks(spectrum_peaks, "centroided")


            except KeyError:

                print('[error] spectrum index out of bound')

                continue

            # # 2018-09-07 Need to catch a number of errors of XML tree not
            # # Being able to read the spectrum object returned by pymzml
            #
            # except xml.etree.ElementTree.ParseError:
            #
            #     print('[warning] XML eTree does not appear to be able to read this spectrum',
            #           '(scan number:', str(nearbyScan_id) + ')', sep=' ')
            #
            #     continue
            #
            # # If it still fails, use the non string version and hope it doesn't freeze
            # except BaseException as err:
            #     print(err)
            #     spectrum = self.msdata[nearbyScan_id]
            #
            # # 2019-05-08 the new pymzml has a weird bug where it would crash when loading some spectra if I specify
            # # the spectrum number as an int; this problem goes away if the identifier is a str but sometimes it retrieves
            # # the wrong spectrum, e.g., for the liverpool test if I retrieve '4698' it gets the spectrum 14698 instead
            # # Since this seems to happen less frequently than the freezing scenario I will use the str for now
            #
            # if not spectrum.ID == nearbyScan_id:
            #     print('[warning] pyMZML failed to retrieve correct spectrum',
            #            '(scan number:', str(nearbyScan_id) + ')', sep=' ')

            #   continue

            # assert spectrum.ms_level == 1, '[error] specified spectrum is not a parent ion scan'

            #Loop through every isotope in the to-do list
            for i in iso_to_do:

                matching_peaks = spectrum.has_peak(peptide_prec + (i*1.007825/z))

                if matching_peaks:
                    for mz_value, I in matching_peaks:
                        timeDependentIntensities.append([nearbyScan_rt, i, I, mz_value])

            del spectrum
            #gc.collect()

        # Integrate the individual isotopomers
        allIso = []

        for j in iso_to_do:

            isotopomer_profile = [[scan, i, I, mz_value] for scan, i, I, mz_value in timeDependentIntensities if i == j]

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
