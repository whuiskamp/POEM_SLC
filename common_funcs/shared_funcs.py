#!/usr/bin/env python3

# File of shared functions that are part of the sea level change (SLC) tool developed
# for the POEM earth system model.
# !!IMPORTANT!! This script assumes that the model grid is constructed with the transition to di-polar dome
# grid at ~65N as is standard practice with MOM models. If this is (for some weird reason) not the case, this needs
# to be accounted for in 'get_halo'.

import numpy as np
import re

def get_halo(MOM,row,col,size,omask,flag=''):
    # This function creates a halo for a grid cell defined at [row,col]
    # in a 2D array of radius 'size' 
    # In:   MOM      - Data structure containing all ocean restart data
    #       row      - Latitude index of grid cell
    #       col      - Longitude index of grid cell
    #      size      - Radius of halo in gridcells
    #      grid_x    - Size of the ocean grid (longitude)
    #      grid_y    - Size of the ocean grid (latitude)
    #        flag    - Ocean cells only = True, otherwise return whole halo mask
    # Out: halo_mask - Boolean mask indicating halo cells
        
    # Define variables
    halo_mask = np.zeros([MOM.grid_y,MOM.grid_x], dtype=bool)
        
    # Identify halo cells
    col_W = col+size
    col_E = col-size
    row_N = row+size
    row_S = row-size
    ## correction for domain edges (prevents out of bounds)
    # We have to account for the polar fold at the N Pole.
    if MOM.lat[row,col] > 65.5:
       # If northernmost row is less than 0, it means the cells we want
       # sit on the other side of the polar fold.
       if row_N > MOM.grid_y-1:
           fold_row = (MOM.grid_y-1)-(row_N-MOM.grid_y)
           # Identify these cells and note them in mask
           PF_cells = fold_cells(fold_row,col,size,MOM).astype(int)
           halo_mask[PF_cells[0,:],PF_cells[1,:]] = True;
                
           # Now we have fold cells, so set limit to Northern boundary
           row_N = MOM.grid_y-1
       # Now to obtain the indices for those cells on this side of the fold
       # There is no cyclic i boundary in the dipolar NH grid
       if col_W >= MOM.grid_x:
           col_W = MOM.grid_x
       if col_E <= 0:
           col_E = 0;
            
    # If we are in the regular grid
    # rightmost column (Westmost)
    if col_W >= MOM.grid_x-1:
        col_W = size-(MOM.grid_x-col);
    # leftmost column (Eastmost)
    if col_E <= 0:
        col_E = MOM.grid_x-(size-col);
    # Uppermost row (Southernmost)
    if row_S < 0:
        row_S = 0;
           
    halo_cells = norm_cells(col_E,col_W,row_S,row_N,MOM.grid_x).astype(int) 
    # Add halo cells to mask
    halo_mask[halo_cells[0,:],halo_cells[1,:]] = True; 
    # Finally, if we only want ocean cells, remove the cell around which the 
    # halo is made and set all land points to false
    if flag == True:
        halo_mask[row,col] = False; halo_mask[omask==0] = False    
    
    return halo_mask

def fold_cells(row,col,size,MOM):
    # This function returns indices of halo cells that sit on the other side of
    # the polar fold.
    # In:    col    - The column index of the cell around which the halo is created
    #        row    - The row index on the other side of the fold
    #        grid_x - Size of the ocean grid (longitude)
    #        grid_y - Size of the ocean grid (latitude)
    #        size   - The radius of the halo (# of gridcells)
    #        MOM    - Data structure containing all ocean restart data
    # Out:   idx    - Index of halo cells on the other side of polar fold
    Nrows = MOM.grid_y-row
    rows =  np.arange(row,MOM.grid_y,1)
    cols = np.arange(col-size,col+size+1,1)
    # Change the indices to refer to the other side of the polar fold
    cols = (MOM.grid_x-1) - cols
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

def halo_eta(row,col,MOM):
	# This function calculates the mean ssh in a halo
	# around point of interest (i,j) and returns a value at (i,j).
    # Eta must already have land values set to nan.
    # In:   eta        - Sea surface height (m)
    #       row        - Latitude grid index
    #       col        - Longitude grid index
    #       o_mask_new - An ocean mask that is updated during the processing 
    #                    of this code
    #       MOM        - Data structure of MOM restart variables
    # Out:  mean_eta   - The average sea level height calculated in a halo around
    #                    cell (row,col)
    
    if ~np.isnan(MOM.eta[0,0]):
        raise ValueError(str('eta not properly formatted. Land = '+str(MOM.eta[0,0]) \
                             + ' not NaN'))
    halo = get_halo(MOM,row,col,1,MOM.o_mask_new,True)
    eta_halo = MOM.eta[halo==True]
    mean_eta = np.mean(eta_halo)
    
    return mean_eta

def calc_coast(omask):
    # This function calculates whether an ocean cells is a coastal cell or not 
    # and creates a mask. Coastal cells are 1, all others are 0.
    #  In: omask      - ocean mask file
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

def get_param(data,name=''):
    # This function returns values of OM4 model parameters of choice
    # The search through the parameter file is done using the regex 'find',
    # which searches for the string 'name' and returns the value associated 
    # with it.
    # In:  data - the parameter file (from either MOM6 or SIS2)
    #      name - string naming the variable sought after
    # Out: val  - The value of paramter 'name' as a float
    
    find = re.compile(r"^("+str(name)+")\s\=\s(\S*)")
    
    for line in data:
        out = find.match(line)
        if out:       
            _, val = out.groups()
    #out = float(val)
    return float(val)