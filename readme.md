# RIAna Relative Isotope Abundance Analyzer v.0.4.0

RIAna (Relative Isotope Abundance Analyzer) takes in standard mass spectrometry spectra and spectral ID files,
and returns mass isotopomer distributions, e.g., for protein turnover analysis.

## Update v.0.4.0

Updated to use python 3.5+, up-to-date scipy and numpy, pymzml.

Multi-fraction runs are now supported - the mzml files in the directory need to be the same as the order they appear
in the Percolator indices (check the Percolator output log file if unsure). For the most cases this shouldn't present a problem if
the mzml directory contains exactly the same mzml files used for the database search, unless there is a difference
in how the operating system order files on the search computer (e.g., file_10.mzml vs. file_2. mzml) or if one of the
fractions contained no protein ID and Crux/Percolator decided to skip the file in its indexing.

To support multi-fraction analysis, RIAna now takes in Percolator tab delimited files for protein ID rather than mzid.
Mzid support will be added back in in a future version.


## Getting Started

Requirements:

On Linux or OSX

	* Install Python 3.5+ and pip
	See instructions on Python website for specific instructions for your operating system

	* Set up virtual environment with venv at a designed path, e.g., ./venv/ria
		$ python3.5 -m venv ./venv/ria

	* Activate the ria environment
		$ source ./venv/ria/bin/activate

	* Check to ensure that the specific python build inside the venv is used
		$ which python

	* Install the packages in the requirements.txt file
		$ pip install -r requirements.txt

On Windows

	* Install Anaconda from Continuum Analytics

	* Install pymzml via pip

	* Get psi-ms-4.0.1.obo from bioontology


Running
	
	* Launch RIAna (Usage/Help)
		$ python3 riana.py --help

	* Example command: This integrates the 0th and 6th isotopomer, requires one lysine, and requires unique peptides
	For heavy water experiments, replace -i 0,6 with -i 0,1,2,3,4; replace -k 1 with -k 0
		$ python3 riana.py ~/test_mzid ~/test_mzml -u -i 0,6 -q 0.01 -r 0.25 -k 1

	* Deactivate the Virtual Environment upon completion
		$ deactivate


Input files

	* RIAna was tested on the percolator output file from Crux Tide/Percolator or standalone Comet/Percolator.

	* The following workflow has been tested for both amino acid and heavy water labeling data gathered on a QE:

	* Convert raw files to mzML, using pwiz 3.0 msconvert in command line, with the following option:
		** --filter "peakPicking vendor"

	* Download Crux 3.1

	* Run Tide index with the following options:
	    ** --digestion partial-digest
	    ** --missed-cleavages

	* Run Tide search with the following options:
		** --isotope-error 1,2 (for HW) or 6,12 (for AA)
		** --compute-sp T
		** --mz-bin-width 0.02
		** --mz-bin-offset 0.0
		** --precursor-window 20
		** --precursor-window-type ppm

	* Run Percolator with the following options:
		** --protein T
		** --fido-empirical-protein-q T

    * Input to RIAna:

	** Take the paths to the directories containing the mzML files (unzipped!) and the percolator.target.psms.txt file

(Obviously you should be able to adjust most if not all of the above parameters to optimize the database search for
particular data sets.)

### Prerequisites

RIAna requires the following:

```
Python 3.5+
pymzml
scipy
numpy
tqdm

```


## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Versioning

We use [SemVer](http://semver.org/) for versioning.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam] (https://github.com/Maggie-Lam)

See also the list of [contributors](https://github.com/ed-lau/pymzml_integrator/graphs/contributors) who participated in this project.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


## Acknowledgments

* [PurpleBooth](https://github.com/PurpleBooth) for Github Readme template.



