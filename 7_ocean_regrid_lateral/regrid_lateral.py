# This script re-distributes mass and tracers to/from ocean cells are either created or removed 
# This script requires the following inputs:
#  - MOM6 restart file
#  - SIS2 restart file
#  - Static grid information for MOM
#  - Mask information showing which cells are to be created/removed
#    (this is generated in chk_water_col.py)
#
# This script iterates through the cell change mask and for each cell that changes
# type, the following is done:
# 1: A halo is created around the cell and a check is performed to ensure that
#    the halo is large enough (eg: 10x surface area of target cell).
# 2: Mass and tracers to be created/removed calculated from target cell
# 3: Mass and tracers redistributed to/from halo cells (weighted by cell area)
# 4: Once all cells have been altered, checks are performed to ensure no 
#    unnaturally large gradients in SSH or tracers.
# (5): Where such gradients exist, smooth them.
#
# This script requires the following functions defined in 'chk_water_col.py':
#   - get_halo; halo_eta; calc_coast
# 
# At this point in time, this script will only redistribute mass, energy and salt!
# If additional tracers are required, this will need to be implemented.

import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/6_check_ocean_cells')
import os
import numpy as np
import copy
import time
#import matplotlib.pyplot as plt
import argparse
from netCDF4 import Dataset as CDF
# Custom functions
from chk_water_col import get_halo, fold_cells, norm_cells

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"

############################## Define functions ###############################

def halo_calc(row,col,data,size,operation):
    # calculate the sum or mean of some value in halo cells
    # returns the ave/sum, the halo mask and the values of the data in the halo
    halo = get_halo(o_mask_new,row,col,size) #Must use *new* omask for the halo calc
    dat_halo = data[halo==True]
    if operation == 'sum':
        sum_halo = np.sum(dat_halo);
        return sum_halo, halo, dat_halo
    elif operation == 'mean':
        ave_halo = np.mean(dat_halo);
        return ave_halo, halo, dat_halo
    
def cell_weights(row,col,size):
    # Calculate the cells weights for redistributing mass and tracers.
    # Weights are based on cell area
    # Size input should be the h_size_mask
    halo_sum,_,_ = halo_calc(i,j,cell_area,size,'sum');

def redist_mass(row,col):
    # This function assumes the global variables chng_mask, h_size_mask,
    # cell_area (tcells from ocean) as well as all relevant field from 
    # ice and ocean restarts.
    # h_ice/sno = ice/snow mass/m2; part_size = fraction of cell covered in ice
    # 
    # 1. Are we filling or emptying a cell?
    if chng_mask(row,col) == -1: # Emptying a cell (ocean -> Land)
    # 2. Calculate the mass to remove from the cell
        # Ice model
        ice_mass = h_ice[row,col]*cell_area[row,col]*part_size[row,col]; 
        sno_mass = h_snow[row,col]*cell_area[row,col]*part_size[row,col];
        # Ocean model
        sea_mass = h_sum[row,col]*cell_area[row,col];
        # Total
        tot_mass = ice_mass + sno_mass + sea_mass;
    elif chng_mask(row,col) == 1: # Filling a cell (land -> ocean)
        # Check surrounding SSH
        
        # Subtract from topog depth, then multiply by cell area
        
        # Create new water column with default z levels
    return

def redist_tracers(row,col):
    # Redistributes energy and salinity only at this point.
    
    return
    
def chk_grads():
    
    return
    
################################# Main Code ###################################
    
if __name__ == "__main__":
    test = True # We'll use different datasets while running tests
    # Open data files
    if test == True:
        MOM6_rest = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc','r')
        chng_file = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc','r')
        SIS2_rest = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc','r')
        Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        grid      = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    else:
        # Insert arg. parsing in here later
        
    #### Extract/define variables ####
    chng_mask       = chng_file.variables['chng_mask'][:,:];         # Mask of cells to change
    h_size_mask     = np.zeros(chng_mask.shape,dtype=float);         # Halo size mask
    h_ice           = SIS2_rest.variables['h_ice'][0,:,:,:];         # Ice thickness
    h_sno           = SIS2_rest.variables['h_snow'][0,:,:,:];        # Snow thickness
    s_ice           = SIS2_rest.variables['sal_ice'][0,:,:,:,:];     # Salinity of sea ice in g/kg
    e_ice           = SIS2_rest.variables['enth_ice'][0,:,:,:,:];    # Enthalpy of sea ice in J
    e_snow          = SIS2_rest.variables['enth_snow'][0,:,:,:,:];   # Enthalpy of sea ice in J
    cell_area       = grid.variables['Ah'][:,:];                     # Area of h (tracer) cells
    h_oce           = MOM6_rest.variables['h'][0,:,:,:];             # Ocean layer thickness
    o_mask_new      = Omask.variables['mask'][:,:];                  # Updated ocean mask
    # Variable pre-processing
    
    # 1: Set up and check halos for all change points, making the halos bigger when required
    # Here size is increased by 1 in the while loop, but that means it will be too big when we 
    # find the correct value. Therefore, when we assign it to h_size_mask, we need to subtract 1 again.
    for i in range(chng_mask.shape[0]):
        for j in range(chng_mask.shape[1]):
            if np.isnan(chng_mask[i,j]) == False:
            #if o_mask_new[i,j] > 0: # For testing all possible ocean cells
                halo_sum = 0;
                size = 1;
                while halo_sum < cell_area[i,j]*10:
                    halo_sum,_,_ = halo_calc(i,j,cell_area,size,'sum');
                    size+=1
                h_size_mask[i,j] = size-1; 
                
    # 2: Redistribute mass and tracers in and out of change points
    # 
    
    for i in range(chng_mask.shape[0]):
        for j in range(chng_mask.shape[1]):
            if np.isnan(chng_mask[i,j]) == False:
                redist_mass[i,j]; redist_tracers[i,j]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    