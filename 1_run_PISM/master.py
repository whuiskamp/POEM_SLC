# Overarching script for dynamic sea level coupling code
# This script reads data from restarts, creates python data structures and initialises
# the main scipts.

from netCDF4 import Dataset as CDF
import copy as cp
import argparse
from create_data_structs import MOM_vars, SIS_vars
import numpy as np
from chk_water_col import calc_coast
from regrid_lateral import get_param

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    test = True # We'll use different datasets while running tests
    





















