# This script checks the integrated column thickness of ocean cells and determines whether or not they should be changed to land.
# 
# This script requires the following inputs:
#  - MOM6 restart file
#  - Updated bathymetry
# 
# 

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


################################## Define functions ###################################

def get_halo(omask, col, row, size):
    """Creates a halo for a grid cell in a 2D array of radius 'size' 
    """
    grid = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    lat = grid.variables['geolat'][:,:]
    lon = grid.variables['geolon'][:,:]
    
    # Define variables
    halo_mask = np.zeros(omask.shape, dtype=bool)
    
    # Identify halo cells
    col_W = col+size
    col_E = col-size
    row_N = row+size
    row_S = row-size
    ## correction for domain edges (prevents out of bounds)
    # We have to account for the polar fold at the N Pole.
    if lat[row,col] > 65.5:
       # If northernmost row is less than 0, it means the cells we want
       # sit on the other side of the polar fold.
       if row_N > omask.shape[0]-1:
           fold_row = (omask.shape[0]-1)-(row_N-omask.shape[0])
           # Identify these cells and note them in mask
           PF_cells = fold_cells(col,fold_row,size).astype(int)
           halo_mask[PF_cells[0,:],PF_cells[1,:]] = True;
                
           # Now we have fold cells, so set limit to Northern boundary
           row_N = omask.shape[0]-1
    # Now to obtain the indices for those cells on this side of the fold
    # There is no cyclic i boundary in the dipolar NH grid
    if col_W >= omask.shape[1]:
        col_W = omask.shape[1]
    if col_E <= 0:
        col_E = 0;
            
    # If we are in the regular grid
    # rightmost column (Westmost)
    if col_W >= omask.shape[1]-1:
        col_W = size-(omask.shape[1]-col);
    # leftmost column (Eastmost)
    if col_E <= 0:
        col_E = omask.shape[1]-(size-col);
    # Uppermost row (Southernmost)
    if row_S < 0:
        row_S = 0;
           
    halo_cells = norm_cells(col_E,col_W,row_S,row_N).astype(int) 
    # Add halo cells to mask
    halo_mask[halo_cells[0,:],halo_cells[1,:]] = True; 
    # Finally, remove the cell around which the halo is made 
    # and set all land points to false
    halo_mask[col,row] = False; halo_mask[omask==0] = False    
    
    return halo_mask

def fold_cells(col,row,size):
    # returns indices of halo cells that sit on the other side of the polar fold
    # This assumes a grid 120x80 cells in size
    # Input: col: The column index of the cell around which the halo is created
    #        row: The row index on the other side of the fold 
    #        size: The radius of the halo
    Nrows = 80-row
    rows =  np.arange(row,80,1)
    cols = np.arange(col-size,col+size+1,1)
    # Change the indices to refer to the other side of the polar fold
    cols = 119 - cols
    Ncols = cols.shape[0]
    idx = np.full([2,Nrows*Ncols],np.nan) # Create cell index (rows,cols) 
    for i in range(Nrows):    
        idx[0,i::Nrows] = rows[i]
    for j in range(Ncols):    
         idx[1,j*Nrows:j*Nrows+Nrows] = cols[j]
        
    return idx    

def norm_cells(E_lim,W_lim,S_lim,N_lim):
    # Returns indices of halo cells that sit on the ocean grid (but not across the polar fold)
    # Input: The cell limits in each direction (N,S,E,W).
    Nrows = N_lim - S_lim +1; rows = np.arange(S_lim,N_lim+1,1)
    Ncols = W_lim - E_lim +1; cols = np.arange(E_lim,W_lim+1,1)
    idx = np.full([2,Nrows*Ncols],np.nan) # Create cell index (rows,cols)
    for i in range(Nrows):    
        idx[0,i::Nrows] = rows[i]
    for j in range(Ncols):    
         idx[1,j*Nrows:j*Nrows+Nrows] = cols[j]
           
    return idx

def halo_eta(eta,col,row):
	# calculates the mean ssh in a halo
	# around point of interest (i,j) and returns an value at (i,j).
    # Eta must already have land values set to nan.
    # Assumes o_mask is a global variable and use of numpy
    if np.isnan(eta[0,0]) == False:
        raise ValueError(str('eta not properly formatted. Land = '+str(eta[0,0]) \
                             + ', not NaN'))
    halo = get_halo(o_mask_new,col,row,1)
    eta_halo = eta[halo==True]
    mean_eta = np.nanmean[eta_halo]
    
    return mean_eta
def calc_coast(omask):
    # calculates whether an ocean cells is a coastal cell and 
    # creates a mask. Coastal cells are 1, all others are 0.
    # Input: ocean mask file
    coast_mask = np.full(omask.shape,np.nan)
    for i in range(omask[0]):
        for j in range(omask[1]):
            if omask[i,j] == 1:
                coast_mask[i,j] = check_neighbour(omask,i,j)
    return coast_mask
                
def check_neighbour(data,row,col):
    # checks data at [i,j] and returns true if any adjacent cell is masked
    # (ie, has value 0 in array)
    # We are assuming that the data being input is output from a GFDL MOM model
    # with a tripolar grid.
    col_p1 = col+1
    col_m1 = col-1
    row_p1 = row+1
    row_m1 = row-1
    
    if col_p1 == data.shape[1]:
        col_p1 = 0;
    # leftmost column (Eastmost)
    if col_m1 == -1:
        col_m1 = data.shape[1]-1;
    # Uppermost row (Southernmost)
    if row_m1 < 0:
        row_m1 = 0;
        
    if data[row,col_p1] == 0 or data[row,col_m1] == 0 or data[row_m1,col] == 0:
        land = True
    else: 
        land = False
    
    # If we're in the top row, the row_p1 cell is on the other side of the polar fold.
    # This means it will be the same top row, but in a mirrored column.
    if row_p1 == data.shape[0]:
        col_f = 199 - col
        if data[row,col_f] == 0:
            land = True
        else:
            land = False
    elif data[row_p1,col] == 0:
            land = True
    else:
        land = False
        
    return land

########################################## Main code ###########################################

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    
    # Import datasets
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
    
    # Check 1.2: Has change in ice sheet/ bathymetry created new ocean? 
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if ice_frac[i,j] <= 0.3 and o_mask[i,j] == 0 and depth[i,j] >= 5: # This needs to be updated to include a SSH check
                o_mask_new[i,j] = 1;
            
    # Check 2: Has change in ice sheet/ bathymetry created new land?
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if (ice_frac[i,j] >= 0.7 and o_mask[i,j] == 1) or (depth[i,j] < 5 and thk): # This needs to be updated to include a SSH check
                o_mask_new[i,j] = 0;
    
    
    # Check 3: Have cells become wet/dry due to changes in mass?
     for i in range(thk.shape[0]):
         for j in range(thk.shape[1]):
             eta_mean = halo_eta(eta,i,j)
             # New land?
             if o_mask_new[i,j] == 1 and 
             # Check column thicnkess
             
             # Comapare with new bathy depth
             
             # Compare with surrounding ave SSH
               
        
        
    
    # This script should output a map identifying which cells will become ocean and 
    # which will become land
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    