# This script re-distributes mass and tracers to/from ocean cells are either created or removed 
# This script requires the following inputs:
#  - MOM6 restart file
#  - SIS2 restart file
#  - Static grid information for MOM
#  - Mask information showing which cells are to be created/removed
#
# This script iterates through the cell change mask and for each cell that changes
# type, the following is done:
# 1: A halo is created around the cell and a check is performed to ensure that
#    the halo is large enough (eg: 10x surface area of target cell).
# 2: Mass and tracers to be created/removed calculated from target cell
# 3: Mass and tracers redistributed to/from halo cells (weighted by cell area)
# 4: Once all cells have been altered, checks are performed to ensure no 
#    unnaturally large gradients in SSH or tracers.
# (5): Where such gradients exist, smooth them.
#

import sys
sys.path.append('../6_check_ocean_cells')
import os
import numpy as np
import copy
import collections as col
import time
import matplotlib.pyplot as plt
import argparse
from netCDF4 import Dataset as CDF
# Custom functions
import chk_water_col

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"

## Define functions ##

def halo_check(col,row,halo_mask):
    
    return

def redist_mass():
    
    return

def chk_grads():
    
    return
    
if __name__ == "__main__":
