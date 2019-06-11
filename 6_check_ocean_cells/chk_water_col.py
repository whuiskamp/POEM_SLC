# This script checks the integrated column thickness of ocean cells and determines whether or not they should be changed to land.
# 
# This script requires the following inputs:
#  - MOM6 restart file
#  - Updated bathymetry
#  - Mask of points to check
#
# The script works as follows: We have from perious steps already checked points that have become land/ocean due to changes in ice
# extent as well as due to changes in bathymetry. The final check is to investigate column height, relative to the updated bathymetry.

## Import packages ##

import sys
import os
import numpy as np
import copy as cp
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

def get_halo(data, col, row, size):
    """Creates a halo for a grid cell in a 2D array of radius 'size' 
    """
    grid = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    lat = grid.variables['geolat'][:,:]
    lon = grid.variables['geolon'][:,:]
    # Step 1: Identidy halo cells
    # We are below polar fold, so calculation can be done easy way.
    if lat[row,col] < 65.5:
        for i in range(size):
            # save neighbor indices
            col_E[i] = col+size
            col_W[i] = col-size
            row_N[i] = row+size
            row_S[i] = row-size
 
        ## correction for domain edges (prevents out of bounds)
            # rightmost column
            if col_E[i] >= data.shape[0]-size:
                col_E[i] = size-(data.shape[0]-col);
            # leftmost column
            if col_W[i] <= size-1:
                col_W[i] = data.shape[0] - (size-col);
            # lowermost row
            if row_S[i] >= data.shape[1]-1:
                row_S[i] = row;
    else:
        # We have to account for the North Pole overlap
        
    # Step 2: Check if they are ocean

def halo_eta(eta,cols,rows)
	# calculates the mean ssh in a halo
	# around point of interest (i,j) and returns an value at (i,j).
    # Eta must already have land values set to nan.
    if np.isnan(eta[0,0]) == False:
        raise ValueError(str('eta not properly formatted. Land = '+str(eta[0,0]) \
                             + ', not NaN'))
    
    
        
    
	
def calc_coast(data):
    # calculates whether an ocean cells is a coastal cell and 
    # creates a mask. Coastal cells are 1, all others are 0.




if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    
    # Import all relevant datasets
    new_bathy = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
    MOM6_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/MOM.res.nc','r')
    #SIS2_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/ice_model.res.nc','r')
    #PISM_data = CDF('path/to/PISM/DATA','r') # this should already have been regridded
    # seaice data probably not required in this step.
    Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
    
    # Extract vars.
    depth      = new_bathy.variables['depth'][:,:]; new_bathy.close()
    h          = MOM6_rest.variables['h'][0,:,:,:];
    ave_eta    = MOM6_rest.variables['ave_ssh'][0,:,:];
    eta        = MOM6_rest.variables['sfc'][0,:,:];
    #ice_frac  = PISM_data.variables['ice_frac'][:,:];
    o_mask     = Omask.variables['mask'][:,:];
    o_mask_new = cp.deepcopy(o_mask);
    
    # Variable pre-processing
    h[:,o_mask==0]      = np.nan;
    ave_eta[o_mask==0]  = np.nan;
    thk                 = np.sum(h,0);
    
    # Check 1.1: Has land ice created new land? Make it land and update mask
    depth(ice_frac >= 0.7 & o_mask > 0) = 0;
    o_mask_new(depth==0) = 0;
    
    # Check 1.2: Has land ice created new ocean? 
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if ice_frac[i,j] <=0.3 and o_mask[i,j] == 0 and depth[i,j] >= 5:
                o_mask_new[i,j] = 1;
            else:
                
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    