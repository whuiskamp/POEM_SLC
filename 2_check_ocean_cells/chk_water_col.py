# This script checks the integrated column thickness of ocean cells and determines whether or not they should be changed to land.
# 
# This script requires the following inputs:
#  - MOM6 restart file
#  - Updated bathymetry
# 
# This script outputs the 'change mask', an array indicating which cells will change from
# ocean to land (-1) or land to ocean (1). All other values are NaN.
#
# There are 3 checks for changes in ocean/land. 
# Check 1: Has land ice exntent or change in bathymetry height created new land?
#          This check is unconditional (ie: no matter what sea level is doing, this
#          cell is now land).
# Check 2: Has column thickness decreased enough to dry out the cell?
#          This is conditional on changes in bathymetry and surrounding SSH.
# Check 3: Have cells become ocean due to receding land ice or SLR?          


## Import packages ##

import sys
import os
import numpy as np
import copy as cp
import time
#import matplotlib.pyplot as plt
import argparse
from netCDF4 import Dataset as CDF

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"


############################## Define functions ###############################

def get_halo(omask, row, col, size):
    # Creates a halo for a grid cell in a 2D array of radius 'size' 
    # Requires global variable omask
    grid = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    lat = grid.variables['geolat'][:,:]
    
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
           PF_cells = fold_cells(fold_row,col,size).astype(int)
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
    halo_mask[row,col] = False; halo_mask[omask==0] = False    
    
    return halo_mask

def fold_cells(row,col,grid_x,grid_y,size):
    # This function returns indices of halo cells that sit on the other side of
    # the polar fold.
    # In:    col    - The column index of the cell around which the halo is created
    #        row    - The row index on the other side of the fold
    #        grid_x - Size of the ocean grid (longitude)
    #        grid_y - Size of the ocean grid (latitude)
    #        size   - The radius of the halo (# of gridcells)
    # Out:   idx    - Index of halo cells on the other side of polar fold
    Nrows = grid_y-row
    rows =  np.arange(row,grid_y,1)
    cols = np.arange(col-size,col+size+1,1)
    # Change the indices to refer to the other side of the polar fold
    cols = (grid_x-1) - cols
    Ncols = cols.shape[0]
    idx = np.full([2,Nrows*Ncols],np.nan) # Create cell index (rows,cols) 
    for i in range(Nrows):    
        idx[0,i::Nrows] = rows[i]
    for j in range(Ncols):    
         idx[1,j*Nrows:j*Nrows+Nrows] = cols[j]
        
    return idx    

def norm_cells(E_lim,W_lim,S_lim,N_lim,grid_x):
    # This function returns indices of halo cells that sit on the ocean grid 
    # (but not across the polar fold). Remember that in python the grid is upside 
    # down (South east corner of the map is the 'top left' element of the array)
    # In:  E(W,S,N)_lim - The cell limits in each direction (N,S,E,W).
    #      grid_x       - Size of the ocean grid (longitude)
    # Out: idx          - Index of halo cells
    Nrows = N_lim - S_lim +1; rows = np.arange(S_lim,N_lim+1,1)
    # Account for lon wrap-around
    if W_lim < E_lim:
        E_cells = grid_x - E_lim; 
        Ncols = W_lim - E_lim + (grid_x+1); cols = np.zeros(Ncols);
        cols[0:E_cells] = np.arange(E_lim,grid_x,1)
        cols[E_cells:Ncols] = np.arange(0,W_lim+1,1)
    else:
        Ncols = W_lim - E_lim +1; cols = np.arange(E_lim,W_lim+1,1)
    idx = np.full([2,Nrows*Ncols],np.nan) # Create cell index (rows,cols)
    for i in range(Nrows):    
        idx[0,i::Nrows] = rows[i]
    for j in range(Ncols):    
        idx[1,j*Nrows:j*Nrows+Nrows] = cols[j]
           
    return idx

def halo_eta(eta,row,col):
	# calculates the mean ssh in a halo
	# around point of interest (i,j) and returns an value at (i,j).
    # Eta must already have land values set to nan.
    # Assumes o_mask is a global variable and use of numpy
    if np.isnan(eta[0,0]) == False:
        raise ValueError(str('eta not properly formatted. Land = '+str(eta[0,0]) \
                             + ' not NaN'))
    halo = get_halo(o_mask_new,row,col,1)
    eta_halo = eta[halo==True]
    mean_eta = np.mean(eta_halo)
    
    return mean_eta
def calc_coast(omask):
    # This function calculates whether an ocean cells is a coastal cell or not 
    # and creates a mask. Coastal cells are 1, all others are 0.
    # In:  omask      - ocean mask file
    # Out: coast_mask - mask of coastal cells
    coast_mask = np.full(omask.shape,0)
    for i in range(omask.shape[0]):
        for j in range(omask.shape[1]):
            if omask[i,j] == 0:
                coast_mask[i,j] = check_neighbour(omask,i,j,1)
    return coast_mask
                
def check_neighbour(data,row,col,flag):
    # This function checks data at [i,j] and returns true if any adjacent cell 
    # fulfills a criteria (ie, has value 0 in array). It is assumed that the 
    # data being used is output from a GFDL MOM model with a tripolar grid.
    # In:  data - The ocean data being checked for a certain value
    #      row  - y value of latitude (ie: position in array, not latitude)
    #      col  - x value of longitude
    #      flag - The value we are checking neighbouring cells for
    # Out: criteria - Is the check successful. Returns true or false.
    
    data = np.around(data) # Values in o_mask are inexact 64bit floats
    criteria = False
    
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
        
    if data[row,col_p1] == flag \
        or data[row,col_m1] == flag \
        or data[row_m1,col] == flag:
        criteria = True
    # If we're in the top row, the row_p1 cell is on the other side of the polar fold.
    # This means it will be the same top row, but in a mirrored column.
    if row_p1 == data.shape[0]:
        col_f = (data.shape[1] -1) - col
        if data[row,col_f] == flag:
            criteria = True
    elif data[row_p1,col] == flag:
         criteria = True
    
        
    return criteria

################################# Main Code ###################################

    # Check 1 Have we created new land via changes in ice sheet extent or
    # topography height? Update mask & change mask
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if (ice_frac[i,j] >= 0.7 and o_mask[i,j] > 0) or 0 < depth[i,j] < 5: 
                o_mask_new[i,j] = 0;
                depth[i,j]      = 0; 
                chng_mask[i,j]  = -1;
    
    # Check 2: Has change in SSH/ bathymetry created new land?
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if (ave_eta[i,j] + depth[i,j] < 2) and o_mask[i,j] > 0:
                o_mask_new[i,j] = 0;
                depth[i,j]      = 0;
                chng_mask[i,j]  = -1;
        
    # Check 3: Have cells become ocean due to receding land ice or SLR?    
    for i in range(depth.shape[0]):
        for j in range(depth.shape[1]):
            if ice_frac[i,j] <= 0.3 and o_mask[i,j] == 0 and depth[i,j] >= 5:
                eta_mean = halo_eta(eta,i,j);
                if eta_mean + depth[i,j] > 2: # optionally add: and coast[i,j] == 1: 
                    o_mask_new[i,j] = 1;
                    chng_mask[i,j]  = 1;
     
    # Finally, we need to make sure our new land-sea mask has not created isolated
    # cells or inland seas, if so, updated o_mask and chng_mask where required
    o_mask_new, chng_mask = chk_cells(o_mask_new,chng_mask)               
                    
                    
    # Write change mask to netCDF
    if test:
        id = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc', 'w')
    else:
        id = CDF('dir', 'w')
    # Create dimensions for vars.
    id.createDimension('lonh', lon.shape[1]);
    id.createDimension('lath', lat.shape[0]);
    id.createVariable('lonh', 'f8', ('lonh'));
    id.variables['lonh'].units = 'degrees_east'
    id.variables['lonh'].cartesian_axis = 'X'
    id.variables['lonh'].long_name = 'T-cell longitude'
    id.variables['lonh'][:] = grid.variables['lonh'][:]
    id.createVariable('lath', 'f8', ('lath'));
    id.variables['lath'].units = 'degrees_north'
    id.variables['lath'].cartesian_axis = 'Y'
    id.variables['lath'].long_name = 'T-cell latitude'
    id.variables['lath'][:] = grid.variables['lath'][:]
    # Define variables and fill data
    # Cell change mask
    id.createVariable('chng_mask', 'f8', ('lath','lonh'))
    id.variables['chng_mask'].units = 'none'
    id.variables['chng_mask'].long_name = 'Mask indicating changing ocean/ land cells'
    id.variables['chng_mask'][:] = chng_mask[:]
    # New ocean mask
    id.createVariable('o_mask_new', 'f8', ('lath','lonh'))
    id.variables['o_mask_new'].units = 'none'
    id.variables['o_mask_new'].long_name = 'Updated ocean mask'
    id.variables['o_mask_new'][:] = o_mask_new[:]
    
    id.description = "This file contains a mask which lists cells that should change from \
                        ocean to land (-1) or land to ocean (1) as well as an updated land-sea mask"
    id.history = "Created " + time.ctime(time.time())
    id.source = "created using /p/projects/climber3/huiskamp/POEM/work/slr_tool/6_check_ocean_cells/chk_water_col.py"
    id.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    