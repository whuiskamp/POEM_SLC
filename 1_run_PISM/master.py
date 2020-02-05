# Over-arching script for dynamic sea level coupling code
# This script reads data from restarts, creates python data structures and initialises
# the main scipts.

#from netCDF4 import Dataset as CDF
#import copy as cp
#import argparse
#import numpy as np
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/1_run_PISM')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/2_check_ocean_cells')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/3_check_channels')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/4_ocean_regrid_lateral')
# Import custom functions
from create_data_structs import init_data_structs 
from chk_water_col import check_water_col
from regrid_lateral import redist_vals

import matplotlib.pyplot as plt
__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2020"
__credits__ = ["Willem Huiskamp", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"

########### FOR TESTING ###########################
# Run after init of data structs
ICE.I_mask[1,82] = 1
ICE.I_mask[79,111] = 1

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    test = True # We'll use different datasets while running tests
    verbose = True # Activates verbose output and error-checking
    # Read in model files and create data structures
    MOM,SIS,OLD,ICE,ETH,FLAGS = init_data_structs('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/',test,verbose)
    # Check to see if cells need to change to/from land/ocean
    check_water_col(MOM,ICE,FLAGS)
    if FLAGS.cont == False:
        # No cells are changing, copy new input files and restart model
        print('No cells change during this coupling step') # Can add to this later? Which coupling step?
    # Implement changes to land sea mask and redistribute relevant tracers
    else:
        # There are cells that need altering, so run applicable scripts.
        redist_vals




# Create new restarts for MOM and SIS - need to find a home for this
#    sp.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc',\
#                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.new.nc'])
#    sp.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc',\
#                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.new.nc'])













