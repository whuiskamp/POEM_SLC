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
# 3: Tracers redistributed to/from halo cells (weighted by cell area)
# 4: Mass redistributed (weighted by cell area)
# 5: Correction applied to halo cells after mass redist. to ensure conservation
# 6: Once all cells have been altered, checks are performed to ensure no 
#    unnaturally large gradients in SSH or tracers.
# (7): Where such gradients exist, smooth them.
#
# This script requires the following functions defined in 'chk_water_col.py':
#   - get_halo; halo_eta; calc_coast
# 
# At this point in time, this script will only redistribute mass, energy and salt
# If additional tracers are required, this will need to be implemented.
# 
# Note that the variable o_mask_new should be the same as the old ocean mask,
# but with points changing to land also masked out. Points changing to ocean should
# remain land in this mask! This is done to ensure a smoother initialisation of new
# ocean points.

import subprocess as sp
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/6_check_ocean_cells')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/7_ocean_regrid_lateral')
import os
import numpy as np
import copy as cp
import time
import math
import re
#import matplotlib.pyplot as plt
import argparse
from netCDF4 import Dataset as CDF
# Custom functions
from chk_water_col import get_halo, fold_cells, norm_cells
from SIS2_funcs import sum_ice_enth

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Prototype"

############################## Define functions ###############################
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
    out = float(val)
    return out

def halo_calc(row,col,data,omask,size,operation):
    # This function calculates the sum or mean of some variable in halo cells
    # Note: this should NOT be used for tracers, as no consideration is made 
    # regarding volume of the gridcell/water column.
    # In:  row       - latitude index
    #      col       - longitude index
    #      data      - Data for which to calculate halo values
    #      omask     - Ocean mask (ocean is 1, land is 0)
    #      size      - size of the halo in # grid cells
    #      operation - Sum or mean of halo cells
    # Out: ave_halo   - Average value for a variable in halo cells
    #      sum_halo   - Sum of values for a variable in halo cells
    #      halo       - Mask of the halo
    #      dat_masked - Returns original data array, but masked with the halo
    halo = get_halo(omask,row,col,size) # Must use *new* omask for the halo calc
    dat_halo = data[halo==True]
    dat_masked = np.full(omask.shape,np.nan)
    dat_masked[halo] = data[halo];
    if operation == 'sum':
        sum_halo = np.sum(dat_halo);
        return sum_halo, halo, dat_masked
    elif operation == 'mean':
        ave_halo = np.mean(dat_halo);
        return ave_halo, halo, dat_masked
    
def cell_weights(row,col,size,omask):
    # This function calculate the cells weights for redistributing mass and tracers.
    # Weights are based on cell area
    # In:  row       - latitude index
    #      col       - longitude index
    #      size      - size of the halo in # grid cells
    #      omask     - Ocean mask (ocean is 1, land is 0)
    # Out: c_weight  - Masked array containing weights (summing to 1)
    #      halo      - The halo for which weights were calculated
    c_weight = np.full(cell_area.shape,np.nan)
    halo_sum,halo,dat_halo = halo_calc(row,col,cell_area,omask,size,'sum');
    c_weight[halo] = dat_halo[halo]/halo_sum
    return c_weight, halo

def h2vgrid(delh,h):
    # This function distributes mass changes evenly over k levels in a single 
    # ocean cell. When the change in SSH is greater than 1m, mass is added to all 
    # h layers proportional to their existing thickness. It assumes that for a 
    # removal of mass, the delh field should be *negative*.
    # In:  delh   - change in mass (h) for this gridcell 
    #      h      - h field from restart
    #      h_lvls - Num. k levels from redist_mass (global var.)
    #      wght   - Weighting for the addition of mass to each cell (global var.)
    # Out: newh   - updated h grid with redistributed mass
    newh = cp.deepcopy(h);
    
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
    # In:  h      - old h field from before h was altered for this cell
    #      newh   - Altered h field, post mass redistribution.
    #      o_temp - ocean pot. temp. (global var)
    #      o_salt - ocean salinity (global var)
    # Out: new_t  - corrected temperature field
    #      new_s  - corrected salinity field
    m_ratio = h/newh;
    new_s   = m_ratio*o_salt
    new_t   = m_ratio*o_temp
    
    return new_t, new_s

def redist_mass(row,col):
    # This function redistributes mass between a target cell and a series of halo
    # cells during the process of creating or removing ocean cells. This function
    # should only be run *after* redist_tracer, as that subroutine relies on old
    # h field.
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
    c_wgts, halo = cell_weights(row,col,size,o_mask_new); # Get weights for re-distribution
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
    # In:  hsum     - Total depth of new water column (excluding 'land depth'
    #                 where land has h of 0.001m)
    #      vgrid    - The default model vertical grid spacing (m)
    #
    # Out: h        - The newly created ocean cell, vertically discretised using 
    #                 model's vgrid.
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

def sum_ocean_enth(row,col, T, h):
    # This function calculates the total energy in seawater in a single cell.
    # The method used is taken from MOM_sum_ouput.F90, line 675. Units are in J
    # In:  h                         - Ocean vertical layer thickness (m) 
    #      T                         - Ocean potential temp. (degC)
    #      cell_area (gloabal var)   - Tracer cell area from ocean model (m2)
    #      C_P     (gloabal var)     - Heat capacity of seawater (J/K/kg)
    #      H_to_kg_m2                - Grid cell thickness to mass conversion factor
    #                                  Default is 1 for Boussinesq numerics
    # Out: total                     - Total ocean energy (J)
        
    total = 0
    for i in range(o_temp.shape[0]):
        total = total + (C_P * T[i,row,col]) * \
                (h[i,row,col]*(H_to_kg_m2 * cell_area[row,col]))
    
    return total
                  
def redist_tracers(row,col,tracer=''):
    # This function redistributes tracers between a cell and target halo cells.
    # It works with energy and salinity only at this point.
    # This function assumes the global variables chng_mask, h_size_mask,
    # cell_area (tcells from ocean) as well as all relevant field from 
    # ice and ocean restarts.
    # Functions from SIS2_funcs.py are utilised.
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
    
    size         = h_size_mask[row,col].astype(int);         # Def. halo radius
    c_wgts, halo = cell_weights(row,col,size,o_mask_redist); # Get weights for re-distribution
    ice_tot = 0;
    
    # 1. Are we filling or emptying a cell?
    if chng_mask[row,col] == -1: # Emptying a cell (ocean -> Land)
        # 1.1. Calculate the tracer to remove from the cell
        if tracer == 'temp':
            # 1.2 Define variables
            global o_temp, e_ice, e_sno, ice_frac, h_ice, h_sno
            # Save old fields for consistency checks
            o_temp_old = cp.deepcopy(o_temp); e_ice_old = cp.deepcopy(e_ice);
            e_sno_old = cp.deepcopy(e_sno)
            # 1.3 Sum enthalpy in sea ice and snow (the latter only has one layer) - During test, should equal -3.144585077871032E+21 across whole domain
            if np.sum(ice_frac[:,row,col],0) > 0:
                ice_tot = sum_ice_enth(row,col,e_ice,e_sno,h_ice,h_sno,ice_frac, \
                               cell_area,s_ice,H_to_kg_m2)
        
            # 1.4 Sum energy in sea water
            ocean_tot = sum_ocean_enth(row,col,o_temp,h_oce)
            tot_heat = ocean_tot + ice_tot
        
            # 1.5 Calculate how much energy each halo cell recieves 
            delta_heat      = c_wgts*tot_heat           # energy going to each halo cell
            o_temp          = heat2cell(delta_heat,o_temp)

                  
    if debugging == True:
        total_old = 0; total_new = 0
        for row in range(cell_area.shape[0]):
            for col in range(cell_area.shape[1]):
                total_old = total_old + sum_ocean_enth(row,col,o_temp_old[:,row,col],
                            h_oce_old[:,row,col]) \
                            + sum_ice_enth(row,col,e_ice_old[:,:,row,col], \
                                          e_sno_old[:,:,row,col],h_ice[:,row,col], \
                                          h_sno_old[:,row,col],ice_frac[:,row,col])
                total_new = total_new + sum_ocean_enth(row,col,o_temp[:,row,col],
                            h_oce[:,row,col]) \
                            + sum_ice_enth(row,col,e_ice[:,:,row,col], \
                                          e_sno[:,:,row,col],h_ice[:,row,col], \
                                          h_sno[:,row,col],ice_frac[:,row,col])           
    elif tracer == 'salt':
        
    
    
    
    
    
    
    return
    
def heat2cell(delta_E,T):
    # This function adds or removes energy to/from a grid cell and updates the 
    # temperature field of the water column. Currently, the approach is to spread
    # the energy change over the entire water column, minimising the change in
    # T in any one layer.
    # In:  delta_E   - The energy being added/removed from a series of halo cells 
    #      T         - The temperature field to be altered
    #      cell_area - Tracer cell area
    #      h_oce     - Ocean layer thickness
    # Out: T_new     - Updated temperature field
    wght = np.full(T.shape,0)
    for i in range(h_oce.shape[0]):
        mass[i,:,:] = h_oce[i,:,:] * cell_area
        wght[i,:,:] = h_oce[i,:,:]/np.sum(h_oce,0)
        T_new = T[i,:,:] + delta_E*wght[i,:,:]/(mass*C_P)  
        
    return T_new

def chk_ssh(h_sum): # Unfinished- double check everything
    # This function checks for sea surface height gradients after redistribution
    # of mass and tracers between ocean cells. For the moment, if neighbouring cells 
    # have a SSH difference of 1m or more, they will get flagged.
    # In:  h_sum       - Total water column height (m)
    #      bathy       - Ocean bathymetry depth (m)
    # Out: lrg_grad    - Field indicating which cells have excessive SSH gradients
    lrg_grd = np.full(h_sum.shape)
    for i in range(h_sum.shape[0]):
        for j in range(h_sum.shape[1]-1): # Cells on other side of zonal split handled later
            if bathy[i,j] != 0 and bathy[i,j+1] != 0:
                diff_x = abs((h_sum[i,j] - bathy[i,j]) - (h_sum[i,j+1] - bathy[i,j+1]))
                if diff_x >= 1:
                    lrg_grd[i,j] = 1
   for i in range(h_sum.shape[0]-1):
        for j in range(h_sum.shape[1]):
            if bathy[i,j] != 0 and bathy[i+1,j] != 0:
                diff_y = abs((h_sum[i,j] - bathy[i,j]) - (h_sum[i+1,j] - bathy[i+1,j]))
                if diff_y >= 1:
                    lrg_grd[i,j] = 1
    return
    
################################# Main Code ###################################
    
if __name__ == "__main__":
    # Create new restarts for MOM and SIS
    sp.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc',\
                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.new.nc'])
    sp.run(['cp','/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc',\
                   '/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.new.nc'])
    test = True # We'll use different datasets while running tests
    debugging = True # Use for consistency checks in script
    # Open data files
    if test:
        MOM6_rest     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc','r')
        MOM6_rest_new = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.new.nc','r+')
        chng_file     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc','r')
        SIS2_rest     = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc','r')
        SIS2_rest_new = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.new.nc','r+')
        Omask         = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        grid          = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
        params_MOM    = open('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM_parameter_doc.all','r').readlines()
        params_SIS    = open('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/SIS_parameter_doc.all','r').readlines()
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
    C_P          = get_param(params_MOM,'C_P');                           # The heat capacity of seawater in MOM6 (J kg-1 K-1)
    H_to_kg_m2   = get_param(params_SIS,'H_TO_KG_M2');                    # grid cell to mass conversion factor (1 by default)
    
    # Variable pre-processing
    h_sum        = np.sum(h_oce,0);                               # Depth of water column (NOT depth of bathymetry)
    tmp          = cp.deepcopy(h_oce);                            # Dummy variable
    tmp          = np.where(tmp>.01, tmp, np.nan);                # While all layers always exist, we only want them if they have 
    h_lvls       = np.nansum(tmp/tmp,axis=0).astype(int); del tmp # significant mass (ie: non-0 layers, defined has having h > 0.001)
    wght         = np.full(h_oce.shape,np.nan)                    # Cell weights for mass redistribution
    for i in range(h_lvls.shape[0]):                              # we add mass to a layer proportional to its thickness
        for j in range(h_lvls.shape[1]):
            lvls = h_lvls[i,j];
            for k in range(lvls):
                wght[k,i,j]= h_oce[k,i,j]/ \
                h_oce[:h_lvls[i,j],i,j].sum(0)                
    o_mask_redist = cp.deepcopy(o_mask_new)                       # New ocean calls masked out. We need a seperate ocean mask for the redistribution  
    o_mask_redist[chng_mask==1] = 0;                              # code, as we cannot allow cells that are yet to be initialised act as halo cells.
    
        
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
                    halo_sum,_,_ = halo_calc(i,j,cell_area,o_mask_new,size,'sum');
                    size+=1
                h_size_mask[i,j] = size-1; 
                
    # 2: Redistribute mass and tracers in and out of change points
    # 
    
    for i in range(chng_mask.shape[0]):
        for j in range(chng_mask.shape[1]):
            if np.isnan(chng_mask[i,j]) == False:
                redist_tracers[i,j]; newh, o_salt, o_temp = redist_mass[i,j]
    if debugging == True:
        err_mass, err_T, err_S = chk_consrv()
        print('Redistribution of mass and tracers complete.' \
              'Error in mass   = '+str(err_mass) \
              'Error in energy = '+str(err_T) \
              'Error in salt   = '+str(err_S))
    # 3: Write new data to copies of old restarts. Several other variables 
    #    will also need to be ammended.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    