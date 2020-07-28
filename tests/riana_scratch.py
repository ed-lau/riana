import os, datetime, re
from riana.read_peptide import ReadPercolator
from riana.read_directory import ReadDirectory
from riana.parse_mzml import Mzml
from functools import partial
from riana import integrate


now = datetime.datetime.now()
directory_to_write = os.path.join('out', 'riana_' + now.strftime('%Y%m%d%H%M%S'))
os.makedirs(directory_to_write, exist_ok=True)

iso_to_do = [0,1,2,3,4,5]

dir_loc = '/Users/edwardlau/Data/Lab_Dox_Turnover/data/'
assert os.path.isdir(dir_loc), '[error] project directory not valid'

# This is the directory that holds the entire project
project = ReadDirectory(dir_loc)

# Get the master peptide ID list
# if input_type == 'Percolator':
mzid = ReadPercolator(project, directory_to_write)
mzid.read_all_project_psms()
mzid.make_master_match_list(lysine_filter=0,
                            peptide_q=0.01,
                            unique_only=True, )


current_sample = project.samples[0]
sample_loc = os.path.join(project.path, current_sample)

mzid.get_current_sample_psms(current_sample=current_sample)
mzid.get_current_sample_mzid_indices()

mzml_files = [f for f in os.listdir(sample_loc) if re.match('^.*.mzML', f)]

# Sort the mzML files by names
# Note this may create a problem if the OS Percolator runs on has natural sorting (xxxx_2 before xxxx_10)
# But we will ignore for now
mzml_files.sort()

idx = 0
mzid.get_current_fraction_psms(idx)
mzid.filter_current_fraction_psms(lysine_filter=0,
                                  protein_q=1,
                                  peptide_q=0.01,
                                  unique_only=True,
                                  use_soft_threshold=True,
                                  match_across_runs=True)

mzml = Mzml(os.path.join(sample_loc, mzml_files[idx]))


# #
# # Read the spectra into dictionary and also create MS1/MS2 indices
# #
mzml.parse_mzml()

integrate_one_partial = partial(integrate.integrate_one,
                                id=mzid.curr_frac_filtered_id_df.copy(),
                                iso_to_do=iso_to_do,
                                mzml=mzml,
                                rt_tolerance=1,
                                mass_tolerance=100e-6,
                                )

integrate_one_partial(0)