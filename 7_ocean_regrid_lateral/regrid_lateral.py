# This script identifies regions in which ocean cells are either created or removed and re-distributes mass and tracers
# on the lateral grid accordingly.
# This script requires the following inputs:
#  - MOM6 restart file
#  - SIS2 restart file
#  - Static grid information for MOM
#  - New and Old topogrpahy files
#
# The script works as follows: First we identify where in the topography an ocean cell has turned into a land cell or vice versa. 
# Next, we look at the integrated column thickness in the ocean restart and determine if is too shallow (< 0.5m) - if so, it is 
# is also flagged to become a land cell. 
#
# An array is created with the indices of all cells that will change (either wet or dry) and relevant attributes are extracted 
# for these cells including h (for ice and ocean), tracer concentrations and SSH. Next we identify a halo around each wettin/drying 
# cell. The initial halo will form a square two gridcells from the cell in question. A check is conducted that enough wet ocean cells 
# fall within this halo (minimum 2). If not, the halo is expanded by one and repeated until this criteria is met.
#   
# The mass to be removed from a drying cell is calculated and divided amongst the number of halo cells. This distribution is weighted
# by halo-cell area.



## Should we allow cells to use the same halo cells? ##

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

if __name__ == "__main__":
