#!/usr/bin/env python3

# Over-arching script for dynamic sea level & land sea mask tool designed for use
# with the OM4 ocean and sea ice model, with optional additional input from the
# ice sheet model PISM and solid Earth model VILMA.
# This script has been designed to use within a bash job script wrapper, but can
# also be used standalone.
# For more information on the tool as a whole, see the README
# This script requires the use of cdo

import time as t
import sys
import matplotlib.pyplot as plt
import argparse
import subprocess as sp

# Path where the tool is saved
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/common_funcs')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/data_processing')

# Import custom functions
from create_data_structs import init_data_structs 
from chk_water_col import check_water_col
from regrid_lateral import redist_vals
import regrid_restarts

__author__ 	   = "Willem Huiskamp"
__copyright__  = "Copyright 2020"
__credits__    = ["Willem Huiskamp", ""]
__license__    = "GPLv3"
__version__    = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ 	   = "huiskamp@pik-potsdam.de"
__status__     = "Production"


if __name__ == "__main__":

########### Argument parsing ###########
	parser = argparse.ArgumentParser(
                description=
                "Identifies cells that should change from ocean to land or vice versa \
                 and changes these accordingly while conserving mass and tracers",
                epilog=
                "Restart files are read in from MOM6, PISM (optional) and VILMA (optional)  \
                 and ocean cells are checked if they should change from ocean to land based \
                 based on water column thickness, grounded ice data from PISM and 			\
                 shallowness of the bathymetry. Similarly, coastal land cells are assessed  \
                 to see if they should become ocean."
            )

    parser.add_argument('-p', '--path', action="store", dest="exp_path",
                        required=True, 
                        help="Path to the experiment directory")
    parser.add_argument('-i', '--iteration', action="store", dest="iteration",
                        required=True, 
                        help="The current coupling interation")
    parser.add_argument('--PISM', action="store", dest="store",
                        required=True, 
                        help="Check for changes in land ice extent?")
    parser.add_argument('--VILMA', action="store", dest="store",
                        required=True, nargs='+',
                        help="list of variable names not to copy to output \
                                file")
    parser.add_argument('-t', '--time', action="store_true", 
                        help="print script timings")
    parser.add_argument('-v', '--verbose', action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()

########### Setup ###########
    t_master_start = t.time()

    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()

    if args.verbose | args.PISM
        print("Running with PISM")
    if args.verbose | args.VILMA
        print("Running with VILMA")
        print()

########### Regridding PISM/VILMA restarts (Optional) ###########

    if args.PISM or args.VILMA
        regrid_rest(args.path)



































