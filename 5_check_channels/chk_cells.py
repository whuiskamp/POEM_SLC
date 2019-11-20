#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script checks a given ocean mask for non-advective cells or cells that 
# have been isolated from the ocean. When found, appropriate action is taken to
# either make them land, or implement a solution that reconnects them with the 
# ocean

import numpy as np
import copy as cp

# For testing, create dummy o_mask array with isolated cells
o_mask = np.full([20,20],1)
o_mask[0::,0] = 0; o_mask[0::,19] = 0; o_mask[0,0::] = 0; o_mask[19,0::] = 0;
o_mask[11::,9] = 0; o_mask[11,9::] = 0;  #o_mask[11,10] = 1;

def chk_cells(o_mask,chng_mask):
    # This function checks the ocean mask for cells or groups of cells isolated
    # from the ocean. 
    # In:  o_mask    - Ocean mask (ocean = 1, land = 0)
    #      chng_mask - A mask indicating which cells are changing from the last 
    #                  simulation to the next (-1 = new land, 1 = new ocean)
    # Out: o_mask    - Updated ocean mask
    #      chng_mask - Updated change mask
    grd = o_mask.shape
    
    
    
    
    return

def ID_iso_cells(row,col,o_mask,turn1,stop,grd):
    # This function consists of a square contour-tracing algorithm 
    # (https://tinyurl.com/qrlo5pm) that will identify an isolated group of 
    # ocean cells and return them as a masked array. The algorithm begins from 
    # the cell that been changed from ocean to land, searching out in two 
    # different directions (defined in turn1). If cells are isolated, one 
    # direction will lead to the open ocean (for which a stop criteria is 
    # implemented) and the other will identify the isolated cells. 
    # In the event that cells are not 
    # isolated, the function will not return anything.
    # In:       row - latitude index
    #           col - logitude index
    #        o_mask - ocean mask
    # Out: iso_mask - Mask of isolated cells associated with cell (row,col)
    #         
    count = 0
    
    # For first step
    row_new, col_new, new_dir = update_orientation('N',row,col,turn1,grd)
    row_1 = cp.deepcopy(row_new); col_1 = cp.deepcopy(col_new);
    dir_1 = cp.deepcopy(new_dir)
    iso_mask = np.full(o_mask.shape,0)
    
    while count < stop:
        if o_mask[row_new,col_new] == 0:
            turn = 'right'
        else:
            turn = 'left'
            iso_mask[row_new,col_new] = 1
            count += 1
        row_new, col_new, new_dir = update_orientation(new_dir,row_new,col_new,turn,grd)
        
        # If we arrive back at the initial cell in the same manner, we've
        # successfully mapped the isolated region, so break loop
        if row_new == row_1 and col_new == col_1 and new_dir == dir_1:
            break
    # If count gets too big, we're in open ocean - return nothing.    
    if count == stop: 
        return
    else:
        return iso_mask

def update_orientation(old_dir,row,col,turn,grd): # Change to just idx?
    # This function computes new indices and directions for the square tracing
    # algorithm using the old direction the tracer is facing as well as whether
    # it now needs to turn 'left' or 'right'. It returns the next grid index
    # for the tracer as well as the new direction it is facing.
    # In:  old_dir - Direction (N,S,E,W) that the tracing algorithm is 'facing'
    # 
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
    if col == grd[0]:
        col = 0
    if col == -1:
        col = grd[0]
    # ... or crossing the polar fold
    if row == grd[1]:
        row = grd[1]-1
        col =  

     
    return row, col, new_dir
        
        
      
test = ID_iso_cells(11,10,o_mask,'right',40,grd)        
if test is None:     
    test = ID_iso_cells(11,10,o_mask,'left',40,grd)

        
        
        
        
        
        
        
        
        
        
        
        