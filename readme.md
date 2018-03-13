
# RIAna Relative Isotope Abundance Analyzer v.0.3.0

RIAna (Relative Isotope Abundance Analyzer) takes in standard mass spectrometry spectra and spectral ID files and returns mass isotopomer distributions for protein turnover analysis.

## Update

Updated to use python 3.6.4, up-to-date scipy and numpy, pymzml 2.0.2.
Removed dependency on docopt.
Added gitignore

## Getting Started

Requirements:

On Linux or OSX

	* Install Python 3.6.4
	See instructions on Python website for specific instructions for your operating system
	Python 3.6 should come with the package manager PIP

	* Set up virtual environment with venv at a designed path, e.g., ~/ria
		$ python3.6 -m venv ~/ria

	* Activate the ria environment
		$ source ~/ria/bin/activate

	* Check to ensure that the specific python build inside the venv is used
		$ which python3.6

	* Install the packages in the requirements.txt file
		$ pip install -r requirements.txt

On Windows

	* Install Anaconda from Continuum Analytics

	* Install pymzml via pip

	* Get psi-ms-4.0.1.obo from bioontology


Running
	
	* Launch RIAna (Usage/Help)
		$ python4 riana.py --help

	* Example command: This integrates the 0th and 6th isotopomer, requires one lysine, and requires unique peptides
	For heavy water experiments, replace -i 0,6 with -i 0,1,2,3,4; replace -k 1 with -k 0
		$ python3 riana.py ~/test.mzid ~/test.mzML -u -i 0,6 -q 0.01 -r 0.25 -k 1

	* Deactivate the Virtual Environment upon completion
		$ deactivate


Input files

	* Riana takes in mzid file output by Crux Tide/Percolator. It may also take in Comet/Percolator (not tested).

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
		** --mzid-output T
		** --fido-empirical-protein-q T

    * Input to Riana:

	** Take the mzML file (unzipped!) and the percolator.target.mzid file inside the percolator output directory.



### Prerequisites

RIAna requires the following:

```
Python 3.6.4

```


## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Versioning

We use [SemVer](http://semver.org/) for versioning.


## Authors

* **Edward Lau, PhD** - *Code* - [ed-lau](https://github.com/ed-lau)


See also the list of [contributors](https://github.com/ed-lau/pymzml_integrator/graphs/contributors) who participated in this project.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


## Acknowledgments

* [PurpleBooth](https://github.com/PurpleBooth) for Github Readme template.



