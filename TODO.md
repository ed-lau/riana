# TODO:

Short term:
- write unit tests for integrate peak class methods with correct assertions
- build in function for kinetic curve fittings - DONE!
- add args for users to specify iso/neutron mass
- add function to calculate peptide accurate mass from scratch and check against percolator
- add support for noncanonical aa
- add support for multiple AAs in the aa tag for SILAC experiment
- validate label free/ibaq against ionquant or maxquant
- reimplement S-G smoothing - DONE!

Mid term:
- write new prep module to do match between run functions and compile new Percolator result
- Probably perform mass error correction there too
- For SILAC, perform isotope incursion correction
- For heavy water, we should calculate the RIA from all isotopomers to reduce effect of outliers/contaminants. 

Long term:
- Graphical user interface

TODO: mass error correction
TODO: isotope correction
TODO: sparse array usage
TODO: mbr
TODO: GUI