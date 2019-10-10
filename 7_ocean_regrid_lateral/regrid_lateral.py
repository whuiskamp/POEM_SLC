# This script re-distributes mass and tracers to/from ocean cells are either created or removed
#
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
# 3: Mass redistributed to/from halo cells (weighted by cell area)
# 4: Correction applied to halo cells after mass redist. to ensure conservation
# 5: Tracers redistributed (weighted by cell area)
# 4: Once all cells have been altered, checks are performed to ensure no 
#    unnaturally large gradients in SSH or tracers.
# (5): Where such gradients exist, smooth them.
#
# This script requires the following functions defined in 'chk_water_col.py':
#   - get_halo; halo_eta; calc_coast
# 
# At this point in time, this script will only redistribute mass, energy and salt
# If additional tracers are required, this will need to be implemented.

import subprocess
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/6_check_ocean_cells')
import os
import numpy as np
import copy
import time
import math
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
    # This function distributes mass changes evenly over k levels in a single 
    # ocean cell. When the change in SSH is greater than 1m, mass is added to all 
    # h layers proportional to their existing thickness. It assumes that for a 
    # removal of mass, the delh field should be *negative*.
    # In:   delh - change in mass (h) for this gridcell 
    #          h - h field from restart
    #     h_lvls - Num. k levels from redist_mass (global var.)
    #     wght   - Weighting for the addition of mass to each cell (global var.)
    # Out: newh - updated h grid with redistributed mass
    newh = copy.deepcopy(h);
    
    for i in range(delh.shape[0]):
        for j in range(delh.shape[1]):
            if ~np.isnan(delh.data[i,j]):
                if 0 < delh[i,j] <= 1:              # Adding mass
                    newh[i,j] = h[0,i,j]+delh[i,j]
                elif delh[i,j] > 1:                 # Adding lots of mass
                    for k in h_lvls-1:
                        newh[k,i,j] = h[k,i,j]+(delh[i,j] * wght[k,i,j])
                elif -1 < delh[i,j] < 0:            # Removing mass
                    newh[i,j] = h[0,i,j]+delh[i,j]
                elif delh[i,j] < -1:                # Removing lots of mass
                    for k in h_lvls-1:
                        newh[k,i,j] = h[k,i,j]+(delh[i,j] * wght[k,i,j])
    return newh

def tracer_correct(h,newh):
    # This function corrects tracer concentrations in cells where h levels have
    # been changed due to an addition/ removal of mass to avoid non-conservation 
    # of tracers. At the moment, only temp. and salt are considered.
    # In:      h - old h field from before h was altered for this cell
    #       newh - Altered h field, post mass redistribution.
    #     o_temp - ocean pot. temp. (global var)
    #     o_salt - ocean salinity (global var)
    # Out: new_t - corrected temperature field
    #      new_s - corrected salinity field
    m_ratio = h/newh;
    new_s   = m_ratio*o_salt
    new_t   = m_ratio*o_temp
    
    return new_t, new_s

def redist_mass(row,col):
    # This function redistributes mass between a target cell and a series of halo
    # cells during the process of creating or removing ocean cells.
    # 
    # In:  chng_mask (gloabal var)        - Mask indicating if cell is to change 
    #                                       from land>ocean or vice versa (1 or -1)
    #      h_size_mask (gloabal var)      - Halo size for cell to be altered (radius in # of cells)
    #      h_oce (global var)             - Ocean vertical layer thickness (m)             
    #      cell_area (global var)         - Tracer cell area from ocean model (m2)
    #      depth (global var)             - Updated ocean cell depth (m)
    #      h_ice (global var)             - Ice mass (kg/m2)
    #      ice_frac (global var)          - Ocean ice fraction
    #      h_sno (global var)             - Snow mass (kg/m2)
    #      o_temp (global var)            - Ocean potential temp. (degC)
    #      o_salt (global var)            - Ocean salinity (ppt)
    # Out: new_h                          - Updates ocean vertical layer thickness (m)
    #      o_temp                         - Updated temp. field (degC)
    #      o_salt                         - Updated salinity field  (ppt)
    #     
    # In the process of changing h of halo cells, a correction is applied to 
    # tracers to avoid spurious creation/ removal of energy/ salt 
    # (and later, other relevant tracers).
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
        h_cell          = h_oce[:,row,col];
        h_sum           = h_cell.sum(0);                  # All h levels exist at every point, just very small
        sea_mass        = h_sum*cell_area[row,col];
        # Total
        tot_mass        = ice_mass + sno_mass + sea_mass;
    # 3. Redistribute mass to halo cells
        delta_mass      = c_wgts*tot_mass           # mass going to each halo cell
        delta_h         = delta_mass/cell_area      # converted to change in h thickness
        newh            = h2vgrid(delta_h,h_oce);   # Apply del h to each depth level in halo
        
    # 4. Correct tracer concentrations for change in h thickness
        o_temp, o_salt = tracer_correct(h_oce,newh) # This alters exising tracer fields and returns
                                                    # values corrected for change in mass
    # 5. Remove ocean cell from h field
        newh[:,row,col] = 0.001;                    # Remove mass from cell and set layers to be land
    elif chng_mask(row,col) == 1: # Filling a cell (land -> ocean)
    # 2. Check surrounding SSH
        eta_mean = halo_eta(eta,row,col);
    # 3. Subtract from topog depth
        init_h = eta_mean + depth(row,col);
    # 4. Create new water column with default z levels
        newh[:,row,col] = newcell(init_h)
    # 5. How much mass is required to initialise this cell?
        tot_mass        = init_h*cell_area[row,col];
    # 6. Remove mass from halo cells
        delta_mass      = c_wgts*tot_mass           # mass coming from each halo cell
        delta_h         = delta_mass/cell_area      # converted to change in h thickness
        newh            = h2vgrid(-delta_h,h_oce);  # Apply del h to each depth level in halo
    # 7. Correct tracer concentrations for change in h thickness
        o_temp, o_salt = tracer_correct(h_oce,newh) # This alters exising tracer fields and returns 
                                                    # values corrected for change in mass
    return newh, o_salt, o_temp

def newcell(hsum):
    # This function initialises a new ocean cell of depth 'hsum' and discretises
    # it into defualt h layers, as defined by the models v-grid
    # In: hsum     - Total depth of new water column (excluding 'land depth'
    #                where land has h of 0.001m)
    #     vgrid    - The default model vertical grid spacing (m)
    #
    # Out: h       - The newly created ocean cell, vertically discretised using 
    #                model's vgrid.
    data  = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/vgrid.nc','r')
    vgrid = data.variables['dz'][:]; 
    n_lvls = 0;
    tmp = hsum;
    for i in range(vgrid.shape[0]):
        tmp -= vgrid[i]
        n_lvls += 1;
        if tmp <= 0:
            break
        
    h = np.full(vgrid.shape[0], 0.001) # 0.001 is the layer thickness of 'land'
    if n_lvls == 1:
        h = hsum
    else:
        for i in range(0,n_lvls-1):
            h[i] = vgrid[i]
        resid = sum(vgrid[0:n_lvls]) - hsum # some frac of the lowest gridcell is land
        h[n_lvls-1] = vgrid[n_lvls-1] - resid
    
    # Run check to make sure operation worked correctly
    if math.floor(sum(h)) != hsum:
        raise ValueError(str('Water column not properly discretised. Total depth is '+str(sum(h)) \
                             + ', where it should be '+str(hsum)))
    
    return h

def redist_tracers(row,col,tracer=''):
    # This function redistributes tracers between a cell and target halo cells.
    # This function works with energy and salinity only at this point.
    # This function assumes the global variables chng_mask, h_size_mask,
    # cell_area (tcells from ocean) as well as all relevant field from 
    # ice and ocean restarts.
    # In:  row                       - latitude index
    #      col                       - longitude index
    #      tracer                    - The tracer to be redistributed
    #      h_ice (gloabal var)       - Ice mass (kg/m2)
    #      h_sno (gloabal var)       - Snow mass (kg/m2)
    #      ice_frac (global var)     - fraction of cell covered in ice (0-1)
    #      h_oce  (gloabal var)      - Ocean vertical layer thickness (m) 
    #      cell_area (gloabal var)   - Tracer cell area from ocean model (m2)
    #      e_ice (global_var)        - Enthalpy of ice (J)
    #      e_sno (global var)        - Enthalpy of snow (J)
    #      o_temp (global var)       - Ocean potential temp. (degC)
    #      o_salt (global var)       - Ocean salinity (ppt)
    #      s_ice (global var)        - Salinity of sea (g/kg)
    #      chng_mask (gloabal var)   - Mask indicating if cell is to change 
    #                                  from land>ocean or vice versa (1 or -1)
    #      h_size_mask (gloabal var) - Halo size for cell to be altered (radius in # of cells)
    #      h_oce (global var)        - Ocean vertical layer thickness (m)
    #      h_ice (global var)        - Ice mass (kg/m2)
    #      ice_frac (global var)     - Ocean ice fraction
    #      h_sno (global var)        - Snow mass (kg/m2)
    # Out: o_temp_new                - Updated temperature field
    #      o_salt_new                - Updated salinity field
    
    size         = h_size_mask[row,col].astype(int);      # Def. halo radius
    c_wgts, halo = cell_weights(row,col,size);            # Get weights for re-distribution
    sum_enth_ice = 0; sum_enth_sno = 0; I_Nk = 1/e_ice.shape[0];
    if tracer == 'temp':
    # Sum enthalpy in sea ice and snow (the latter only has one layer) - During test, should equal -3.144585077871032E+21 across whole domain
        if np.sum(ice_frac,0)[row,col] > 0:       # If no ice exists, skip it
            for row in range(ice_frac.shape[1]):
                for col in range(ice_frac.shape[2]):
                    for l in range(e_ice.shape[1]):
                        for k in range(e_ice.shape[0]):
                            sum_enth_ice = sum_enth_ice + (cell_area[row,col] \
                            * ((h_ice[l,row,col]*ice_frac[l,row,col])*I_Nk)) * e_ice[k,l,row,col]
                        sum_enth_sno = sum_enth_sno + (cell_area[row,col] \
                        * ((h_sno[l,row,col]*ice_frac[l,row,col]))) * e_sno[0,l,row,col]
    total = sum_enth_ice + sum_enth_sno + 3.144585077871032E+21;
    # We are out by 0.21549873394247008e+21 or 7.3%
    # Other function returns -3.566105016138741e+21
    # Sum energy in sea water
    
    elif tracer == 'salt':
        
    
    
    
    
    
    
    return
    
def chk_ssh():
    # This function checks for sea surface height gradients after redistribution
    # of mass and tracers between ocean cells. 
    
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
    chng_mask    = chng_file.variables['chng_mask'][:,:];                 # Mask of cells to change
    h_size_mask  = np.zeros(chng_mask.shape,dtype=float);                 # Halo size mask
    h_ice        = SIS2_rest.variables['h_ice'][0,:,:,:].data;            # Ice thickness
    h_sno        = SIS2_rest.variables['h_snow'][0,:,:,:].data;           # Snow thickness
    ice_frac     = SIS2_rest.variables['part_size'][0,1:6,:,:].data;      # Ice fraction
    s_ice        = SIS2_rest.variables['sal_ice'][0,0,0,0,0].data;        # Salinity of sea ice in g/kg - it's a fixed value
    e_ice        = SIS2_rest.variables['enth_ice'][0,:,:,:,:].data;       # Enthalpy of sea ice in J
    e_sno        = SIS2_rest.variables['enth_snow'][0,:,:,:,:].data;      # Enthalpy of snow in J
    cell_area    = grid.variables['Ah'][:,:];                             # Area of h (tracer) cells
    h_oce        = MOM6_rest.variables['h'][0,:,:,:].data;                # Ocean layer thickness
    o_temp       = MOM6_rest.variables['Temp'][0,:,:,:].data;             # Ocean potential temperature (deg C)
    o_salt       = MOM6_rest.variables['Salt'][0,:,:,:].data.astype(int); # Ocean salinity (ppt)
    o_mask_new   = Omask.variables['mask'][:,:];                          # Updated ocean mask
    
    # Variable pre-processing
    h_sum        = np.sum(h_oce,0);                               # Depth of water column (NOT depth of topography)
    tmp          = copy.deepcopy(h_oce);                          # Dummy variable
    tmp          = np.where(tmp>.01, tmp, np.nan);                # While all layers always exist, we only want them if they have 
    h_lvls       = np.nansum(tmp/tmp,axis=0).astype(int); del tmp # significant mass (ie: non-0 layers, defined has having h > 0.001)
    wght         = np.full(h_oce.shape,np.nan)                    # Cell weights for mass redistribution
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
                newh, o_salt, o_temp = redist_mass[i,j]; redist_tracers[i,j]
    
    # 3: Write new data to copies of old restarts. Several other variables 
    #    will also need to be ammended.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    