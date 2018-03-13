
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
	For heavy water experiments, replace iso 0,6 with iso 0,1,2,3,4,5
		$ python3 riana.py ~/test.mzid ~/test.mzML -u -i 0,6 -q 0.005 -r 0.5 -k 1

	* Deactivate the Virtual Environment upon completion
		$ deactivate


Input files

	* Riana takes in mzid file output by Crux Tide/Percolator. It may also take in Comet/Percolator (not tested).

	* The following workflow is recommended for both amino acid and heavy water labeling.
	 (Bullseye - Tide - Percolator) for both amino acid and heavy water analysis. Comet actually works better for
	 heavy water labeling data with the --isotope-error tag that is disabled in Tide, but I think the Crux/Comet distribution
	 automatically performs protein inference for some reason which does not help us remove non-unique peptides

	* Convert raw files to mzML,  ms1, and ms2 files using PWiz (Use PWiz command line and not GUI to ensure correct format)
		$ msconvert /Path/to/your.raw 
		$ msconvert /Path/to/your.raw --ms1
		$ msconvert /Path/to/your.raw --ms2

	* Download Crux 3.0 or above

	* Run Bullseye (helps correct for heavy water labeling to bring the Tide performance closer to Comet)
		$ bin/crux bullseye /Path/to/your.ms1 /Path/to/your.ms2s --output-dir YOUR_OUTPUT_DIR --overwrite T

	* For heavy water labeling data (except day 0), add this flag to Bullseye
		$ --averagine-mod 0.05H1

	* Run Tide index
		$ bin/crux tide-index /Path/to/your.fasta SWISSPROT_MM_DB --overwrite T --digestion partial-digest --missed-cleavage 1

	* Run Tide search (Note that in the absence of the —isotope-error flag I had to widen the precursor window to 1.005 Da for the heavy water experiments. Far from ideal but I think currently the best option until we figure out how to make Comet link each PSM sequence to all entries in the fasta file.)
		$ bin/crux tide-search --compute-sp T --output-dir YOUR_OUTPUT_DIR --overwrite T --percursor-window 1.005 --precursor-window-type mz --mz-bin-offset 0.02 --mz-bin-offset 0.0 --isotope-error 1 YOUR_OUTPUT_DIR/bullseye.pid.ms2 SWISSPROT_MM_DB

	* Run Percolator
		$ bin/crux percolator --protein T --overwrite T --output-dir YOUR_OUTPUT_DIR --fido-empirical-protein-q T --mzid-output T --post-processing-qvality T YOUR_OUTPUT_DIR/tide-search.target.txt

    * Input to Riana:

	* The percolator.target.mzid file inside YOUR_OUTPUT_DIR is the mzIdentML file we need.



### Prerequisites

RIAna requires the following:

```
Python 3.5.2

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



