# This script checks the integrated column thickness of ocean cells and determines whether or not they should be changed to land.
# 
# This script requires the following inputs:
#  - MOM6 restart file
#  - Updates bathymetry
#  - Mask of points to check
#
# The script works as follows: We have from perious steps already checked points that have become land/ocean due to changes in ice
# extent as well as due to changes in bathymetry. The final check is to investigate column height, relative to the updated bathymetry.

## Import packages ##

import sys
import os
import numpy as np
import copy
import collections as col
import time
import matplotlib.pyplot as plt
import argparse
from netCDF4 import Dataset as CDF

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"


## Define functions ##

def halo_eta(data,col,row)
	# calculates the mean water column thicnkess in a halo
	# around point of interest (i,j) and saves this to 
	# array at (i,j).
	thk = np.sum()





if __name__ == "__main__":

