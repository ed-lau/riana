# -*- coding: utf-8 -*-

""" Main. """
import logging
import re
import os
import sys
from functools import partial

from riana.project import ReadDirectory
from riana.peptides import ReadPercolator
from riana.spectra import Mzml

from riana import integrate, params, __version__

import tqdm
import pandas as pd


def runfit(args):
    print('Fit function not implemented yet.')
    pass
