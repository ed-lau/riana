# -*- coding: utf-8 -*-

""" Functions to integrate isotopomer peaks """

import numpy as np
from riana import params


def integrate_one(index,
                  id,
                  iso_to_do,
                  rt_tolerance,
                  mass_tolerance,
                  mzml,
                  deuterium_mass_defect,
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

    proton = params.proton_mass

    if deuterium_mass_defect:
        iso_added_mass = params.deuterium_mass_diff # 1.003354835
    else:
        iso_added_mass = params.c13_mass_diff

    # get peptide mass, scan number, and charge
    peptide_mass = float(id.loc[index, 'peptide mass'])
    scan_number = int(id.loc[index, 'scan'])
    charge = float(id.loc[index, 'charge'])

    # get retention time from Percolator scan number
    peptide_rt = mzml.rt_idx[np.searchsorted(mzml.scan_idx, scan_number, side='left')]

    assert isinstance(peptide_rt.item(), float), '[error] cannot retrieve retention time from scan number'

    # calculate precursor mass from peptide monoisotopic mass
    peptide_prec = (peptide_mass + (charge * proton)) / charge

    intensity_over_time = []

    # choose the scan numbers from the index (beware that scan number is 1-indexed)
    nearby_ms1_scans = mzml.scan_idx[np.abs(mzml.rt_idx - peptide_rt) <= rt_tolerance]

    # loop through each MS1 spectrum within retention time range
    for scan in nearby_ms1_scans:

        try:
            spec = mzml.msdata[np.array(np.where(mzml.scan_idx == scan)).item()]

        except ValueError or IndexError:
            raise Exception('Scan does not correspond to MS1 data.')

        for iso in iso_to_do:

            # set upper and lower mass tolerance bounds
            prec_iso_am = peptide_prec + (iso * iso_added_mass / charge)
            delta_mass = prec_iso_am*mass_tolerance/2

            # sum intensities within range
            matching_int = np.sum(spec[np.abs(spec[:, 0] - prec_iso_am) <= delta_mass, 1])

            intensity_over_time.append([iso,
                                        mzml.rt_idx[mzml.scan_idx == scan].item(),
                                        matching_int,
            ]
            )

    if not intensity_over_time:
        raise Exception(
            "No intensity profile for peptide {0}".format(prec_iso_am)
        )

    intensity_array = np.array(intensity_over_time)

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

    # integrate the individual isotopomers
    iso_intensity = []

    for j in iso_to_do:

        isotopomer_profile = intensity_over_time[intensity_over_time[:, 0] == j]

        if isotopomer_profile.size > 0:

            # use np.trapz rather than scipy.integrate to integrate
            iso_area = np.trapz(isotopomer_profile[:, 2], x=isotopomer_profile[:, 1])

        # if there is no isotopomer profile, set area to 0
        else:
            iso_area = 0

        iso_intensity.append(iso_area)

    if not iso_intensity:
        raise Exception("No positive numerical value integrated for isotopmer {0}".format(intensity_over_time)
        )

    return iso_intensity
