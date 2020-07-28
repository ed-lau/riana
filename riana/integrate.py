import os
import pandas as pd
import numpy as np
import scipy.integrate

def integrate_one(index,
                  id,
                  iso_to_do,
                  rt_tolerance,
                  mass_tolerance,
                  mzml,
                  ):
    """
    Wrapper for the get_isotope_from_amrt() function below

    :param index: int The row number of the peptide ID table passed from the wrapper.
    :return: list [index, pep_id, m0, m1, m2, ...]

    Given peptide accurate mass and retention time and charge, find all the isotopic peaks intensity at each
    scan within the retention time window

    :param peptide_am: float Accurate peptide mass
    :param peptide_scan: int Scan number
    :param am_is_mz: bool Is the accurate mass actually prec m/z (for lipid)
    :param scan_is_rt: bool Is the scan number actually rtime (for lipid)
    :param z: int Peptide charge
    :return: List of intensity over time
    """


    peptide_mass = float(id.loc[index, 'peptide mass'])
    scan_number = int(id.loc[index, 'scan'])
    charge = float(id.loc[index, 'charge'])

    # print('Peptide mass', peptide_mass)
    # print('scan number', scan_number)
    # print('charge', charge)


    # Isotope mass from NIST https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
    # Electron mass from NIST https://www.physics.nist.gov/cgi-bin/cuu/Value?meu|search_for=electron+mass

    # 2019-05-29 proton mass is different from hydrogen mass is different fron neutron mass
    # Perhaps we should add neutron mass instead of proton mass for the isotopes which may
    # make a difference when iso is high enough (e.g., 12 for KK determination)
    # in the future we may have to account for mass defects
    proton = 1.007276466621  # 1.007825
    # The above can also be accessed through scipy.constants.physical_constants['proton mass in u'] but
    # we will hard code for now

    # neutron = 1.00866491595
    # electron = 0.000548579907

    # Mass defect of deuterium is 1.007276466621 +  1.00866491595 + 0.000548579907 - 2.01410177812 = 0.00238818435
    # Mass defect of C13 is 6 * 1.007276466621 + 7 * 1.00866491595 + 6 * 0.000548579907 - 13.00335483507 = 0.1042
    # Mass difference of C13 - C12 = 1.003354835
    # Mass difference of deuterium - protium = 1.00627674589
    iso_added_mass = 1.003354835

    # Get retention time from Percolator scan number
    peptide_rt = mzml.rt_idx[mzml.scan_idx == scan_number]

    assert isinstance(peptide_rt.item(), float), '[error] cannot retrieve retention time from scan number'


    # Calculate precursor mass from peptide monoisotopic mass
    peptide_prec = (peptide_mass + (charge * proton)) / charge



    intensity_over_time = []

    # Choose the scan numbers from the index (watch out that scan is 1-indexed ..)
    nearby_scans = mzml.scan_idx[np.abs(mzml.rt_idx - peptide_rt.item()) <= rt_tolerance]

    # All nearby ms1 scans:
    nearby_ms1_scans = nearby_scans[mzml.mslvl_idx[nearby_scans] == 1]

    # Loop through each spectrum, check if it is an MS1 spectrum, check if it is within 1 minute of retention time
    for scan in nearby_ms1_scans:

        spec = mzml.msdata[scan-1]  # Going back to 0 index

        # Get the spectrum based on the spectrum number
        for iso in iso_to_do:
            # Set upper and lower bound
            prec_iso_am = peptide_prec + (iso * iso_added_mass / charge)
            delta_mass = prec_iso_am*mass_tolerance/2

            matching_int = np.sum(spec[np.abs(spec[:,0] - prec_iso_am) <= delta_mass, 1])

            intensity_over_time.append([iso,
                                        mzml.rt_idx[mzml.scan_idx == scan].item(),
                                        matching_int,
                                        #prec_iso_am,
            ]
            )

    if not intensity_over_time:
        raise Exception("No intensity profile for peptide {0}".format(
            prec_iso_am
        )
        )

    # self.peak_logger.debug(intensity_over_time)

    print(intensity_over_time)

    intensity_array = np.array(intensity_over_time)

    print(array)
    #print(array[:, 0])
    #print(array[:, 0] == 1)
    #print(array[array[:, 0] == 1])

    if not intensity_over_time:
        print('Empty intensity over time')

    result = [index] + [(id.loc[index, 'pep_id'])] + integrate_isotope_intensity(intensity_array,
                                                                                 iso_to_do=iso_to_do)

    return result


def integrate_isotope_intensity(intensity_over_time,
                                iso_to_do
                                ):
    """
    Given a list of isotopomer intensity over time, give the integrated intensity of each isotopomer

    :return: Integrated intensity of each isotopomer
    """
    # Integrate the individual isotopomers
    iso_intensity = []

    for j in iso_to_do:

        isotopomer_profile = [[rt, I] for rt, iso, I, mz_value in intensity_over_time if iso == j]

        # If there is no isotopomer profile, set area to 0
        if isotopomer_profile:
            iso_df = pd.DataFrame(isotopomer_profile)
            iso_area = scipy.integrate.trapz(iso_df[1], iso_df[0])
            # Remove all negative areas
            iso_area = max(iso_area, 0)
            # Round to 1 digit
            iso_area = np.round(iso_area, 1)

            # nearbyScan_rt, iso, matching_int, peptide_prec_iso_am
            # self.peak_logger.debug('intensity_over_time {0)'.format(np.array(iso_df[1])))

        else:
            iso_area = 0

        iso_intensity.append(iso_area)

    if not iso_intensity:
        raise Exception("No positive numerical value integrated for isotopmer {0}".format(
            intensity_over_time
        )
        )

    return iso_intensity