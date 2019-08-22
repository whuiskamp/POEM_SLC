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

import subprocess
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
    dat_masked = np.full(o_mask_new.shape,np.nan)
    dat_masked[halo] = data[halo];
    if operation == 'sum':
        sum_halo = np.sum(dat_halo);
        return sum_halo, halo, dat_masked
    elif operation == 'mean':
        ave_halo = np.mean(dat_halo);
        return ave_halo, halo, dat_masked
    
def cell_weights(row,col,size):
    # Calculate the cells weights for redistributing mass and tracers.
    # Weights are based on cell area
    # Size input should be the h_size_mask
    # Requires global var. cell_area
    c_weight = np.full(cell_area.shape,np.nan)
    halo_sum,halo,dat_halo = halo_calc(row,col,cell_area,size,'sum');
    c_weight[halo] = dat_halo[halo]/halo_sum
    return c_weight, halo

def h2vgrid(delh,h):
    # This function distributes mass evenly over k levels in a single ocean cell
    # when the change in SSH is greater than 1m to avoid over-filling k=1.
    # In:   delh - change in mass (h) for this gridcell 
    #          h - h field from restart
    #     h_lvls - Num. k levels from redist_mass (global var.)
    #     wght   - Weighting for the addition of mass to each cell (global var.)
    # Out: new_h - updated h grid with redistributed mass
    newh = copy.deepcopy(h);
    
    for i in range(delh[0]):
        for j in range(delh[1]):
            if ~np.isnan(delh[i,j]):
                if delh[i,j] <= 1:
                    newh[i,j] = h[0,i,j]+delh[i,j]
                else:
                    for k in h_lvls:
                        new_h[k,i,j] = h[k,i,j]+(delh[i,j] * wght[k,i,j])
    return newh

def redist_mass(row,col):
    # This function assumes the global variables chng_mask, h_size_mask,
    # cell_area (tcells from ocean) as well as all relevant fields from 
    # ice and ocean restarts. In the process of changing h of halo cells, 
    # a correction is applied to tracers.
    # h_ice/sno = ice/snow kg/m2; ice_frac = fraction of cell covered in ice
    # h_oce = layer thickness of ocean cells m, cell_area = area of ocean cells m2
    # 
    size         = h_size_mask[row,col].astype(int);      # Def. halo radius
    c_wgts, halo = cell_weights(row,col,size);            # Get weights for re-distribution
    ice_mass = 0; sno_mass = 0;                           # Initialise remaining vars
    # 1. Are we filling or emptying a cell?
    if chng_mask[row,col] == -1: # Emptying a cell (ocean -> Land)
    # 2. Calculate the mass to remove from the cell
        # Ice model
        for i in range(h_ice.shape[0]):
            ice_mass = ice_mass + h_ice[i,row,col]*cell_area[row,col]*ice_frac[i,row,col]; 
            sno_mass = sno_mass + h_sno[i,row,col]*cell_area[row,col]*ice_frac[i,row,col];
        # Ocean model
        h_cell     = h_oce[:,row,col];
        h_sum      = h_cell.sum(0);                       # All h levels exist at every point, just very small
        sea_mass   = h_sum*cell_area[row,col];
        # Total
        tot_mass   = ice_mass + sno_mass + sea_mass;
    # 3. Redistribute mass to halo cells
        delta_mass = c_wgts*tot_mass
        delta_h    = delta_mass/cell_area
        newh = h2vgrid(delta_h,h_oce);
        
    elif chng_mask(row,col) == 1: # Filling a cell (land -> ocean)
        # Check surrounding SSH
        
        # Subtract from topog depth, then multiply by cell area
        
        # Create new water column with default z levels
    return

def redist_tracers(row,col,tracer):
    # Redistributes energy and salinity only at this point.
    # This function assumes the global variables chng_mask, h_size_mask,
    # cell_area (tcells from ocean) as well as all relevant field from 
    # ice and ocean restarts.
    # h_ice/sno = ice/snow kg/m2; ice_frac = fraction of cell covered in ice
    # h_oce = layer thickness of ocean cells m, cell_area = area of ocean cells m2
    # e_ice/sno
    return
    
def chk_grads():
    
    return
    
################################# Main Code ###################################
    
if __name__ == "__main__":
    # Create new restarts for MOM and SIS
    subprocess.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc',\
                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.new.nc'])
    subprocess.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc',\
                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.new.nc'])
    test = True # We'll use different datasets while running tests
    # Open data files
    if test:
        MOM6_rest     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc','r')
        MOM6_rest_new = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.new.nc','r+')
        chng_file     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc','r')
        SIS2_rest     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc','r')
        SIS2_rest_new = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.new.nc','r+')
        Omask         = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        grid          = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    else:
        # Insert arg. parsing in here later
        
    #### Extract/define variables ####
    chng_mask       = chng_file.variables['chng_mask'][:,:];         # Mask of cells to change
    h_size_mask     = np.zeros(chng_mask.shape,dtype=float);         # Halo size mask
    h_ice           = SIS2_rest.variables['h_ice'][0,:,:,:];         # Ice thickness
    h_sno           = SIS2_rest.variables['h_snow'][0,:,:,:];        # Snow thickness
    ice_frac        = SIS2_rest.variables['part_size'][0,1:6,:,:];   # Ice fraction
    s_ice           = SIS2_rest.variables['sal_ice'][0,0,0,0,0];     # Salinity of sea ice in g/kg - it's a fixed value
    e_ice           = SIS2_rest.variables['enth_ice'][0,:,:,:,:];    # Enthalpy of sea ice in J
    e_sno           = SIS2_rest.variables['enth_snow'][0,:,:,:,:];   # Enthalpy of sea ice in J
    cell_area       = grid.variables['Ah'][:,:];                     # Area of h (tracer) cells
    h_oce           = MOM6_rest.variables['h'][0,:,:,:].data;        # Ocean layer thickness
    o_temp          = MOM6_rest.variables['Temp'][0,:,:,:].data;     # Ocean potential temperature (degC)
    o_salt          = MOM6_rest.variables['Salt'][0,:,:,:].data;     # Ocean salinity (ppt)
    o_mask_new      = Omask.variables['mask'][:,:];                  # Updated ocean mask
    
    # Variable pre-processing
    h_sum           = np.sum(h_oce,0);                               # Depth of water column (NOT depth of topography)
    tmp             = copy.deepcopy(h_oce);                          # Dummy variable
    tmp             = np.where(tmp>1, tmp, np.nan);                  # While all layers always exist, we only want them if they have 
    h_lvls          = np.nansum(tmp/tmp,axis=0).astype(int); del tmp # significant mass (ie: non-0 layers)
    wght            = np.full(h_oce.shape,np.nan)                    # Cell weights for mass redistribution
    for i in range(h_lvls.shape[0]):                                 # we add mass to a layer proportional to its thickness
        for j in range(h_lvls.shape[1]):
            lvls = h_lvls[i,j];
            for k in range(lvls):
                wght[k,i,j]= h_oce[k,i,j]/ \
                h_oce[:h_lvls[i,j],i,j].sum(0)                
    
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
    
    # 3: Write new data to copies of old restarts. Several other variables 
    #    will also need to be ammended.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    