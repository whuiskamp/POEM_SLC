# Over-arching script for dynamic sea level coupling code
# This script reads data from restarts, creates python data structures and initialises
# the main scipts.

from netCDF4 import Dataset as CDF
import copy as cp
import argparse
import numpy as np
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/1_run_PISM')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/2_check_ocean_cells')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/3_check_channels')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/4_ocean_regrid_lateral')
# Import custom functions
from create_data_structs import init_data_structs 
import chk_water_col 
import chk_cells
import SIS2_funcs

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    test = True # We'll use different datasets while running tests
    debugging = True # Activates verbose output and error-checking
    # Read in model files and create data structures
    MOM,SIS,OLD = init_data_structs('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/',test)
    # Check to see if cells need to change to/from land/ocean
    check_water_col
    # Implement changes to land sea mask and redistribute relevant tracers
    redist_vals



















