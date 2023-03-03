""" setup module """

from setuptools import setup, find_packages
import os.path

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
   name='riana',
   version='0.7.2',
   description='riana integrates the relative abundance of isotopomers in mass spectra',

   long_description=long_description,
   long_description_content_type='text/markdown',

   url='https://github.com/ed-lau/riana',

   author='Edward Lau',
   author_email='edward.lau@cuanschutz.edu',

   classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],

   keywords='scientific turnover kinetics proteomics mass-spectrometry',  # Optional

   packages=find_packages(),

   python_requires='>=3.9, <4',

   install_requires=['pymzml>=2,<3',
                     'scipy>=1,<2',
                     'pandas>=1,<2',
                     'matplotlib>=3,<4',
                     'tqdm>=4,<5',
                     'scikit-learn>=1'], #external packages as dependencies

   entry_points={
           'console_scripts': [
               'riana=riana.main:main',
           ],
       },

   project_urls={
        'Source': 'https://github.com/ed-lau/riana',
        'Edward Lau Lab': 'https://www.laulab.net',
    },

   # data_files=[('tests',
   #               [os.path.join('tests', 'data', 'sample1', 'percolator.target.mzid'),
   #                os.path.join('tests', 'data', 'sample1', 'percolator.target.psms.txt'),
   #                os.path.join('tests', 'data', 'sample1', '20180216_BSA.mzML.gz'),
   #                ]),
   #              ],

)
