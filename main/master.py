#!/usr/bin/env python3

# Over-arching script for dynamic sea level & land sea mask tool designed for use
# with the OM4 ocean (MOM6) and sea ice (SIS2) model, with optional additional input from the
# ice sheet model PISM and solid Earth model VILMA (or any comparable models that provide
# grounded ice and relative sea level fields, respectively).
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
import regrid_restarts
from create_data_structs import init_data_structs 
from chk_water_col import check_water_col
from regrid_lateral import redist_vals
import prep_restarts

__author__ 	   = "Willem Huiskamp"
__copyright__  = "Copyright 2023"
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
                        help="Path to the experiment (SLC data) directory")
    parser.add_argument('-i', '--iteration', action="store", dest="iteration",
                        required=True, 
                        help="The current coupling interation")
    parser.add_argument('--PISM', action="store_true", dest="PISM",
                        required=False, 
                        help="Running with coupled PISM (optional)")
    parser.add_argument('--VILMA', action="store_true", dest="VILMA",
                        required=False,
                        help="Running with coupled VILMA (optional)")
    parser.add_argument('-v', '--verbose', action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()

########### Settings ###########
    min_depth = 5  # The depth at which a water column is made land
    min_thk   = 2  # The minimum water column thickness allowable before a cell is made land
    new_depth = 5  # The depth of a newly initialised ocean cell
    iso_depth = 10 # The depth at which isolated cells (essentially inland ocean) are made land
    iso_size  = 5  # The number of isolated cells (in a group) below which, these cells are made land
    def_halo  = 10 # The default halo size when redistributing mass/tracers. Defined as a factor of the target cell area.
                   # E.g., a value of 5 would mean that the halo must be >= 5 times the area of the target cell. 
 
########### Setup ###########
    t_master_start = t.time()

    if args.verbose:
        print(" -> verbose output = True")
        print()
    if args.verbose and args.PISM:
        print("Running with PISM")
    if args.verbose and args.VILMA:
        print("Running with VILMA")
        print()
    if args.verbose:
        print(" Configurations settings:")
        print("min_depth = "+str(min_depth)+"m")
        print("min_thk   = "+str(min_thk)+"m")
        print("new_depth = "+str(new_depth)+"m")
        print("iso_depth = "+str(iso_depth)+"m")
        print("iso_size  = "+str(iso_size)+" cells")
        print("def_halo  = "+str(def_halo)+" x target cell area")
########### Regridding PISM/VILMA restarts (Optional) ###########
    t_data_start = t.time()

    if args.PISM or args.VILMA:
        t_regr_start = t.time()
        if args.verbose:    
            print("Regridding inputs files from PISM/VILMA...")
        regrid_rest(args.path)
        t_regr = t.time() - t_regr_start

########### Read in model files and create data structures ###########
    MOM,SIS,OLD,ICE,OPTS = init_data_structs((str(args.exp_path)),args.PISM,args.VILMA,args.verbose)
    # Apply settings
    OPTS.min_depth = min_depth
    OPTS.new_depth = new_depth
    OPTS.iso_depth = iso_depth
    OPTS.iso_size  = iso_size
    OPTS.def_halo  = def_halo

    t_data = t.time()-t_data_start
########### Check if any cells need to change from land-ocean & vice versa ###########
    t_check_start = t.time()
    check_water_col(MOM,ICE,OPTS)
    t_check = t.time()-t_check_start
    if (OPTS.cont == False) and (OPTS.bathy_chg == False):
        print("No cells are changing, copy new input files and restart model")
        pass
    elif (OPTS.cont == False) and (OPTS.bathy_chg == True):
        t_change_start = t.time()
        update_bathy(MOM)
    else:
        t_change_start = t.time()
        # There are cells that need altering, redistribute mass and tracers
        redist_vals(MOM,SIS,OLD,OPTS)
        # Now we need to update secondary fields in the restart files
        prep_fields(MOM,SIS)
        # Write out new restart files
        write_rest(MOM,SIS,OPTS)
        t_change = t.time() - t_change_start
        t_total = t.time() - t_master_start
    
    # Optional performance information
    if args.verbose:
        print("SLC tool finished. Time taken for each operation:")
        print("Pre-processing: "+str(t_regr)+" secs")
        print("Data importing (total): "+str(t_data)+" secs")
        print("Check cells: "+str(t_check)+" secs")
        if t_change in locals():
            print("Updating restarts: "+str(t_change)+" secs")
        print("Total: "+str(t_total)+" secs")