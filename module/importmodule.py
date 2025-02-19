# import all need modules...
# Check python version
import sys, argparse, os
ver = sys.version_info # Get Python version
if ver.major != 3:
    sys.exit('Python 3 is required to run phot-d!')
# Astropy -----------------------------------------------------
from astropy.table import Table, Column, hstack, unique, join
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astropy import units as u
# Others  -----------------------------------------------------
import numpy as np
import timeit
import multiprocessing as mp
from functools import partial
from time import localtime, strftime, time, gmtime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# -------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
