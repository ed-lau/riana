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
