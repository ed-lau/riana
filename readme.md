Riana Relative Isotope Abundance Analyzer v.0.2.0

—-------------
- Requirements:
—-------------

On Linux or OSX

	# Install Python 3.5.2
	# See instructions on Python website for specific instructions for your operating system
	# Python 3.5 should come with the package manager PIP

	# Set up a Virtual Environment for Riana
	# Install Virtualenv and Virtualenvwrapper if you have not done so
		$ pip3 install virtualenv
		$ pip3 install virtualenvwrapper

	# Set up env, go for example to project folder, then set up a
	# Virtual Environment folder under a particular destination, e.g., ~/ria
		$  virtualenv —-python=python3.5 ~/ria

	# Activate the ria environment
		$ source ~/ria/bin/activate

	# Check to ensure that the specific python build inside the VirtualEnv is used
		$ which python

	# Install the packages in the requirements.txt file
		$ pip3 install -r requirements.txt

On Windows
	Install Anaconda from Continuum Analytics
	Install docopt via conda; pymzml via pip
	Go get psi-ms-4.0.1.obo from bioontology

—-------------
- Running:
—-------------
	
	# Launch Riana (Usage/Help)
		$ python3 riana.pyc --help

	# Example command: This integrates the 0th and 6th isotopomer, requires one lysine, and requires unique peptides
	# For heavy water experiments, replace iso 0,6 with iso 0,1,2,3,4,5 
		$ python3 riana.pyc integrate /Path/to/your.mzid /Path/to/your.mzML iso 0,6 -K 1 -u --rt 0.5

	# Deactivate the Virtual Environment upon completion
		$ deactivate

—-------------
- Input files:
—-------------

	# Riana takes in mzIdentML file output by Crux Tide/Percolator. It should also take in Crux Comet/Percolator output.
	# It may take in any Percolator output mzid file or even any mzid file, but I recommend the following workflow
	# (Bullseye - Tide - Percolator) for both amino acid and heavy water analysis. Comet actually works better for
	# heavy water labeling data with the --isotope-error tag that is disabled in Tide, but I think the Crux/Comet distribution
	# automatically performs protein inference for some reason which does not help us remove non-unique peptides

	# Convert raw files to mzML,  ms1, and ms2 files using PWiz (Use PWiz command line and not GUI to ensure correct format)
		$ msconvert /Path/to/your.raw 
		$ msconvert /Path/to/your.raw --ms1
		$ msconvert /Path/to/your.raw --ms2

	# Download Crux 3.0 or above

	# Run Bullseye (helps correct for heavy water labeling to bring the Tide performance closer to Comet)
		$ bin/crux bullseye /Path/to/your.ms1 /Path/to/your.ms2s --output-dir YOUR_OUTPUT_DIR --overwrite T

	# For heavy water labeling data (except day 0), add this flag to Bullseye
		$ --averagine-mod 0.05H1

	# Run Tide index
		$ bin/crux tide-index /Path/to/your.fasta SWISSPROT_MM_DB --overwrite T --digestion partial-digest --missed-cleavage 1

	# Run Tide search (Note that in the absence of the —isotope-error flag I had to widen the precursor window to 1.005 Da for the heavy water experiments. Far from ideal but I think currently the best option until we figure out how to make Comet link each PSM sequence to all entries in the fasta file.)
		$ bin/crux tide-search --compute-sp T --output-dir YOUR_OUTPUT_DIR --overwrite T --percursor-window 1.005 --precursor-window-type mz --mz-bin-offset 0.02 --mz-bin-offset 0.0 --isotope-error 1 YOUR_OUTPUT_DIR/bullseye.pid.ms2 SWISSPROT_MM_DB

	# Run Percolator
		$ bin/crux percolator --protein T --overwrite T --output-dir YOUR_OUTPUT_DIR --fido-empirical-protein-q T --mzid-output T --post-processing-qvality T YOUR_OUTPUT_DIR/tide-search.target.txt

	# Input to Riana: 
	# The percolator.target.mzid file inside YOUR_OUTPUT_DIR is the mzIdentML file we need.



Contact: Edward Lau
lau1@stanford.edu