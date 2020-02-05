#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script checks a given ocean mask for non-advective cells or cells that 
# have been isolated from the ocean. When found, appropriate action is taken to
# either make them land, or leave them alone (so long as they are stable).
# (REMEMBER THAT NORTH IS 'DOWN'when viewing in Spyder)

import numpy as np
import copy as cp

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2020"
__credits__ = ["Willem Huiskamp", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"

# For testing, create dummy o_mask array with isolated cells
# Normal square block
#o_mask = np.full([20,20],1)
#o_mask[0::,0] = 0; o_mask[0::,19] = 0; o_mask[0,0::] = 0; o_mask[19,0::] = 0;
#o_mask[11::,9] = 0; o_mask[11,9::] = 0;  #o_mask[11,10] = 1;
#
## More complex shape
#o_mask = np.full([20,20],1)
#o_mask[0::,0] = 0; o_mask[0::,19] = 0; o_mask[0,0::] = 0; o_mask[19,0::] = 0;
#o_mask[18,9] = 0; o_mask[17,10:16] = 0; o_mask[16,16] = 0; o_mask[16,18] = 0;
#o_mask[15,17] = 0;

def chk_cells(MOM,FLAGS):
    # This function checks the ocean mask for cells or groups of cells isolated
    # from the ocean. For the tracing algorithm to function, we must first 
    # eliminate individual isolated cells.
    # In:  o_mask    - Ocean mask (ocean = 1, land = 0)
    #      chng_mask - A mask indicating which cells are changing from the last 
    #                  simulation to the next (-1 = new land, 1 = new ocean)
    #      MOM       - Data structure containing all ocean restart data
    # Out: o_mask    - Updated ocean mask
    #      chng_mask - Updated change mask
        
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if MOM.chng_mask[i,j] == 1 or MOM.chng_mask[i,j] == -1:
                print('row='+str(i)+'col='+str(j))
                # For algorithm to work, we must search in both 'directions'.
                # Try one, then the other.
                iso_mask = ID_iso_cells(i,j,'right',MOM,FLAGS)
                if iso_mask is None:
                    iso_mask = ID_iso_cells(i,j,'left',MOM,FLAGS)
                # If we have isolated cells, determine what to do with them
                if iso_mask is not None:
                    chng_mask_2, o_mask_2 = fix_iso_cells(iso_mask,MOM)
    return

def ID_iso_cells(row,col,turn1,MOM,FLAGS):
    # This function consists of a square contour-tracing algorithm 
    # (https://tinyurl.com/qrlo5pm) that will identify an isolated group of 
    # ocean cells and return them as a masked array. The algorithm begins from 
    # the cell that been changed from ocean to land, searching out in one of two 
    # different directions (defined in turn1). If cells are isolated, one 
    # direction will lead to the open ocean (for which a stop criteria is 
    # implemented) and the other will identify the isolated cells (if they exist). 
    # In the event that cells are not isolated, the function will not 
    # return anything.
    # In:       row - latitude index
    #           col - logitude index
    #        o_mask - ocean mask
    #         turn1 - The tracing algorithm can turn left or right upon init.
    #                 This means we can end up in either isolated cells or the 
    #                 open ocean. Therefore, if one direction does not work, try 
    #                 the other.
    #          stop - The algorithm's stopping criteria. In this case, if the
    #                 number of cells traced exceeds a certain large value,
    #                 it is assumed we are in the open ocean rather than an isolated
    #                 group of cells, and the function is terminated.
    #           grd - Grid dimensions in form [nrows,ncols]
    # Out: iso_mask - Mask of isolated cells associated with cell (row,col)
    #         
    count       = 0       # Number of isolated ocean cells identified
    steps_taken = 0       # Number of cells the tracker has moved through
    grd = MOM.o_mask_new.shape
    stop = MOM.grid_y/2 # Pick an arbitrary stopping criteria scaled by grid size
    loopy       = stop*10 # If this many steps are taken, you've got an infinite loop 
    # For first step
    row_new, col_new, new_dir = update_orientation('N',row,col,turn1,grd)
    row_1 = cp.deepcopy(row_new); col_1 = cp.deepcopy(col_new);
    dir_1 = cp.deepcopy(new_dir)
    iso_mask = np.full(MOM.o_mask_new.shape,0)
    
    # Extra check is needed if ocean is isolated by land cells only connected
    # by their vertices. The square tracing algorithm cannot normally handle this.
    while count < stop and steps_taken < loopy:
        if MOM.o_mask_new[row_new,col_new] == 0:
            turn = 'right'
            steps_taken += 1
            #print('test 1')
        elif MOM.o_mask_new[row_new,col_new] == 1 and 1 not in iso_mask:
            turn = 'left'
            iso_mask[row_new,col_new] = 1
            count += 1
            steps_taken += 1
            #print('test 2')
        elif MOM.o_mask_new[row_new,col_new] == 1 and bounds_chk(row_new,col_new,iso_mask):
            turn = 'left'
            iso_mask[row_new,col_new] = 1
            count += 1
            steps_taken += 1
            #print('test 3')
        elif MOM.o_mask_new[row_new,col_new] == 1 and not bounds_chk(row_new,col_new,iso_mask):
            turn = 'right'
            steps_taken += 1
            #print('test 4')
        try:    
            row_new, col_new, new_dir = update_orientation(new_dir,row_new,col_new,turn,grd)
        except UnboundLocalError:
            print('Error at row='+str(row)+', col='+str(col))
        
        # If we arrive back at the initial cell in the same manner, we've
        # successfully mapped the isolated region, so break loop
        if row_new == row_1 and col_new == col_1 and new_dir == dir_1:
            break
    # If count gets too big, we're in open ocean - return nothing.
     
    if count == stop: 
        return
    elif 1 in iso_mask:
        if FLAGS.verbose:
            print('Isolated cells found in vicinity of row = '+str(row)+', col = '+str(col))
        return iso_mask
    else:
        return None

def update_orientation(old_dir,row,col,turn,grd):
    # This function computes new indices and directions for the square tracing
    # algorithm using the old direction the tracer is facing as well as whether
    # it now needs to turn 'left' or 'right'. It returns the next grid index
    # for the tracer as well as the new direction it is facing.
    # In:  old_dir - Direction (N,S,E,W) that the tracing algorithm is 'facing'
    #      row     - Latitude index of current grid cell
    #      col     - Longitude index of current grid cell
    #      turn    - The direction the tracer will turn for the next step
    #      grd     - Grid size information (nRows,nCols)
    # Out: row     - Latitude index of next grid cell
    #      col     - Longitude index of next grid cell
    #      new_dir - Direction that the tracing algorithm is 'facing' in next 
    #                grid cell
    
    idx = np.full(2,0)
    idx[0] = row; idx[1] = col
    if turn == 'left':
        if old_dir == 'N':
            idx += [0, 1]
            new_dir = 'W'
        elif old_dir == 'S':
            idx += [0, -1]
            new_dir = 'E'
        elif old_dir == 'E':
            idx += [1, 0]
            new_dir = 'N'
        elif old_dir == 'W':
            idx += [-1, 0]
            new_dir = 'S'
    elif turn == 'right':
        if old_dir == 'N':
            idx += [0, -1]
            new_dir = 'E'
        elif old_dir == 'S':
            idx += [0, 1]
            new_dir = 'W'
        elif old_dir == 'E':
            idx += [-1, 0]
            new_dir = 'S'
        elif old_dir == 'W':
            idx += [1, 0]
            new_dir = 'N'
    
    row = idx[0]; col = idx[1]

    # Check if tracker is crossing prime meridian
    if col == grd[1]:
        col = 0
    if col == -1:
        col = grd[1]-1
    # ... or crossing the polar fold
    if row == grd[0]:
        row = grd[0]-1
        col = (grd[1]-1)-col

     
    return row, col, new_dir

def one_cell_chk(row,col,MOM):
    # This function checks for individual isolated cells. If/when found, the 
    # chng_mask will be updated.
    # In:        row - Latitude index
    #            col - Longitude index
    #            MOM - Ocean data structure
    # Out: chng_mask - Mask specifying which cells will change from ocean>land
    #                  or vice versa
    
    # For cells not interacting with the polar fold
    if row < MOM.o_mask_new.shape[0]-1:
        if col == MOM.o_mask_new.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = MOM.o_mask_new.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        if row == 0:
            row_m1 = 0; row_p1 = 1
        else:
            row_m1 = row - 1; row_p1 = row+1
            
        if MOM.o_mask_new[row_p1,col] == 0 \
            and MOM.o_mask_new[row_m1,col] == 0 \
            and MOM.o_mask_new[row,col_m1] == 0 \
            and MOM.o_mask_new[row,col_p1] == 0:
            MOM.chng_mask[row,col] = -1
            MOM.o_mask_new[row,col] = 0
        else: 
            return
    
    # For cell at the polar fold, row_p1 will lie on other side. ie: same row,
    # different col.    
    elif row == MOM.o_mask_new.shape[0]-1:
        if col == MOM.o_mask_new.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = MOM.o_mask_new.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        row_m1 = row-1
        if MOM.o_mask_new[row,(MOM.o_mask_new.shape[1]-1)-col] == 1 \
           or MOM.o_mask_new[row_m1,col] == 1 \
           or MOM.o_mask_new[row,col_m1] == 1 \
           or MOM.o_mask_new[row,col_p1] == 1:
           MOM.chng_mask[row,col] = -1
           MOM.o_mask_new[row,col] = 0
        else: 
            return
        
def bounds_chk(row,col,data):
    # This checks whether an ocean cell is in fact within the same land-bounded
    # region as other ocean cells, as the tracing algorithm can fail when ocean
    # is bounded by land where only vertices of the cells are touching.
    # Checks/ corrections are included in case prime meridian/ polar fold are 
    # crossed.
    # In:    row - Latitude index
    #        col - Longitude index
    #       data - Mask of cells identified by the ID_iso_cells function. Ocean 
    #              cells are marked as 1, land as 0
    # Out: flag  - When true, these ocean cells can said to be connected.
    #              The opposite is true when false.
    
    # For cells not interacting with the polar fold
    if row < data.shape[0]-1:
        if col == data.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = data.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        if row == 0:
            row_m1 = 0; row_p1 = 1
        else:
            row_m1 = row - 1; row_p1 = row+1
            
        if data[row_p1,col] == 1 or data[row_m1,col] == 1 \
            or data[row,col_m1] == 1 or data[row,col_p1] == 1:
            return True
        else:
            return False
    # For cell at the polar fold, row_p1 will lie on other side. ie: same row,
    # different col.    
    elif row == data.shape[0]-1:
        if col == data.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = data.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        row_m1 = row-1
        
        if data[row,(data.shape[1]-1)-col] == 1 or data[row_m1,col] == 1 \
            or data[row,col_m1] == 1 or data[row,col_p1] == 1:
            return True
        else:
            return False  
    
def chk_sides(row,col,iso_mask,MOM):
    # This function will check to see if isolated cells have been missed by the 
    # tracing algorithm, given that it only identifies the *outline* of the 
    # isolated region. It does so by checking cells neighbouring identified
    # isolated cells and checking to see if they are also ocean.
    # In:        row - Latitude index
    #            col - Longitude index
    #            MOM - Ocean data structure
    #       iso_mask - Mask showing the location of isolated cells associated with
    #                  a single changed gridcell
    #     o_mask_new - Updated ocean mask (includes points that will change for
    #                  next run)
    # Out:  iso_mask - Updated iso_mask
    
    # For cells not interacting with the polar fold
    if row < MOM.o_mask_new.shape[0]-1:
        if col == MOM.o_mask_new.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = MOM.o_mask_new.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        if row == 0:
            row_m1 = 0; row_p1 = 1
        else:
            row_m1 = row - 1; row_p1 = row+1
        
        # cell must be adjacent to isolated cells, but not already identified    
        if (iso_mask[row_p1,col] == 1 \
            or iso_mask[row_m1,col] == 1 \
            or iso_mask[row,col_m1] == 1 \
            or iso_mask[row,col_p1] == 1) \
            and iso_mask[row,col] == 0 \
            and MOM.o_mask_new[row,col] == 1:
                iso_mask[row,col] = 1
        else: 
            return
    # For cell at the polar fold, row_p1 will lie on other side. ie: same row,
    # different col.    
    elif row == MOM.o_mask_new.shape[0]-1:
        if col == MOM.o_mask_new.shape[1]-1:
            col_m1 = col-1; col_p1 = 0
        elif col == 0:
            col_p1 = 1; col_m1 = MOM.o_mask_new.shape[1]-1
        else:
            col_m1 = col-1; col_p1 = col+1
        row_m1 = row-1
        
        # cell must be adjacent to isolated cells, but not already identified
        if (iso_mask[row,(iso_mask.shape[1]-1)-col] == 1 \
           or iso_mask[row_m1,col] == 1 \
           or iso_mask[row,col_m1] == 1 \
           or iso_mask[row,col_p1] == 1) \
           and iso_mask[row,col] == 0 \
           and MOM.o_mask_new[row,col] == 1:
               iso_mask[row,col] = 1
        else: 
            return

def fix_iso_cells(iso_mask,MOM):
    # This function takes a group of isolated ocean cells, examines their 
    # collective properties and determines whether or not to leave them alone 
    # or fill them in.
    # In:        MOM - Ocean data structure
    #       iso_mask - Mask showing the location of isolated cells associated with
    #                  a single changed gridcell
    #      chng_mask - Mask showing which cells will change from land-sea or vice
    #                  versa.
    # Out: chng_mask - Updated change mask.
    #         o_mask - Updated ocean mask
    
    # Firstly, identify any isolated cells the tracing algoritm may have missed
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            chk_sides(i,j,iso_mask,MOM)
    # If cells are less than 5m in depth, fill them in
    tmp = MOM.h_sum[iso_mask==1] 
    if np.mean(tmp) < 5:
        MOM.chng_mask
    # If the group consists of 3 or fewer cells, fill it in
    size = np.sum(iso_mask)
    if size <= 3:
        MOM.chng_mask[iso_mask==1] = -1
    
    # If group consists of 4 or more cells and is deeper than 5m, leave it be
   
    
    
    
    
    
    
    return
#grd = o_mask.shape      

#test = ID_iso_cells(17,10,'right',MOM)        
#if test is None:     
#    test = ID_iso_cells(17,10,o_mask,'left',10,grd)
#
#test = ID_iso_cells(18,9,o_mask,'right',40,grd)        
#if test is None:     
#    test = ID_iso_cells(18,9,o_mask,'left',40,grd)
     
        
        
        
        
        
        
        
        
        
        
        