# RIANA - Relative Isotope Abundance Analyzer

RIANA (Relative Isotope Abundance Analyzer) takes in standard mass spectrometry spectra and spectral ID files,
and returns mass isotopomer distributions, e.g., for protein turnover analysis.



## Getting Started

Installation:

Install Python 3.7+ and pip. See instructions on Python website for specific instructions for your operating system.

Riana can be installed from PyPI via pip. We recommend using a virtual environment.

    $ pip install riana

Launch riana as a module (Usage/Help):
	
	$ python -m riana

Alternatively as a console entry point:

    $ riana
    
To test that the installation can load test data files in tests/data:

    $ pip install tox
    $ tox

To run the riana test dataset (a single fraction bovine serum albumin file from a Q-Exactive) and print the result
to the home directory:

    $ python -m riana tests/data/ -u -i 0,1,2,3,4,5 -q 0.1 -r 0.5 -t 10 -o ~/
    
Notes on the expected input files:

	* RIANA.py was tested on the percolator output file from Crux Tide/Percolator or standalone Comet/Percolator.

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

    * Input to RIANA.py:

    ** Note that the riana argument path should point to the project directory, where each individual sample mzML and
     search result files are placed under a sub-directory (e.g., sample1/)
	

### Prerequisites

RIANA.py is tested in Python 3.7 and 3.8 and uses the following packages:

```
matplotlib==3.2.1
pandas==1.0.4
pymzml==2.4.6
scipy==1.4.1
tqdm==4.46.0
scikit-learn==0.23.1
```


## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)

See also the list of [contributors](https://github.com/ed-lau/pymzml_integrator/graphs/contributors) who participated in this project.


## License

This project is licensed under the MIT License - see the LICENSE.md file for details


