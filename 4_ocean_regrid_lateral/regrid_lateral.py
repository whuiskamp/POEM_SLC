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
# 5: Correction applied to halo cells after mass redist. to ensure tracer conservation
# 6: Once all cells have been altered, checks are performed to ensure no 
#    unnaturally large gradients in SSH.
# (7): Where such gradients exist, smooth them. (currently not implemented)
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

#import subprocess as sp
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/2_check_ocean_cells')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/1_run_PISM')
#import os
import numpy as np
import copy as cp
import time
import math
#import matplotlib.pyplot as plt
#from netCDF4 import Dataset as CDF
# Custom functions
from shared_funcs import get_halo
from shared_funcs import halo_eta
from SIS2_funcs import sum_ice_enth
from SIS2_funcs import sum_ice_sal

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2020"
__credits__ = ["Willem Huiskamp", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Alpha"

############################## Define functions ###############################
def halo_calc(row,col,data,MOM,size,operation):
    # This function calculates the sum or mean of some variable in halo cells
    # Note: this should NOT be used for tracers with the 'sum' flag, 
    # as no consideration is made regarding volume of the gridcell/water column.
    # In:  row               - latitude index
    #      col               - longitude index
    #      data              - Data for which to calculate halo values (2D field)
    #      o_mask_redist     - Ocean mask (ocean is 1, land is 0) new ocean cells still marked as land
    #      size              - size of the halo in # grid cells
    #      operation         - Sum or mean of halo cells
    # Out: ave_halo          - Average value for a variable in halo cells
    #      sum_halo          - Sum of values for a variable in halo cells
    #      halo              - Mask of the halo
    #      dat_masked        - Returns original data array, but masked with the halo
    halo = get_halo(MOM,row,col,size,MOM.o_mask_redist,True) # Must use *new* omask for the halo calc
    dat_halo = data[halo==True]
    dat_masked = np.full(MOM.o_mask.shape,np.nan)
    dat_masked[halo] = data[halo];
    if operation == 'sum':
        sum_halo = np.sum(dat_halo);
        return sum_halo, halo, dat_masked
    elif operation == 'mean':
        ave_halo = np.mean(dat_halo);
        return ave_halo, halo, dat_masked
        
    
def cell_weights(row,col,size,MOM):
    # This function calculate the cells weights for redistributing mass and tracers.
    # Weights are based on cell area
    # In:  row       - latitude index
    #      col       - longitude index
    #      size      - size of the halo in # grid cells
    #      MOM       - Ocean data structure
    # Out: c_weight  - Masked array containing weights (summing to 1)
    #      halo      - The halo for which weights were calculated
    c_weight = np.full([MOM.grid_y, MOM.grid_x],np.nan) # For some stupid reason, this function only works if this is nan
    halo_sum,halo,dat_halo = halo_calc(row,col,MOM.cell_area,MOM,size,'sum');
    c_weight[halo] = dat_halo[halo]/halo_sum
    c_weight[np.isnan(c_weight)] = 0
    return c_weight, halo

def h2vgrid(delh,h,MOM):
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
    lvls = MOM.h_lvls.astype(int)
    
    for i in range(delh.shape[0]):
        for j in range(delh.shape[1]):
            if delh.data[i,j] != 0:
                if 0 < delh[i,j] <= 1:              # Adding mass
                    newh[0,i,j] = np.sum([h[0,i,j],delh[i,j]])
                elif delh[i,j] > 1:                 # Adding lots of mass
                    for k in range(lvls[i,j]):
                        newh[k,i,j] = np.sum([h[k,i,j],np.multiply(delh[i,j],MOM.wght[k,i,j])])
                elif -1 < delh[i,j] < 0:            # Removing mass
                    newh[0,i,j] = np.sum([h[0,i,j],delh[i,j]])
                elif delh[i,j] < -1:                # Removing lots of mass
                    for k in range(lvls[i,j]):
                        newh[k,i,j] = np.sum([h[k,i,j],np.multiply(delh[i,j],MOM.wght[k,i,j])])
    if np.any(np.less(newh,0)): # check if we've accidentally created cells of -tive thickness
        raise ValueError(str('Oops! You have negative layer thicknesses'))
    return newh

def tracer_correct(h,newh,MOM):
    # This function corrects tracer concentrations in cells where h levels have
    # been changed due to an addition/ removal of mass to avoid non-conservation 
    # of tracers. At the moment, only temp. and salt are considered.
    # In:  h      - old h field from before h was altered for this cell
    #      newh   - Altered h field, post mass redistribution.
    #      o_temp - ocean pot. temp. (global var)
    #      o_salt - ocean salinity (global var)
    # Out: new_t  - corrected temperature field
    #      new_s  - corrected salinity field
    m_ratio = np.divide(h,newh);
    new_s   = np.multiply(m_ratio,MOM.o_salt)
    new_t   = np.multiply(m_ratio,MOM.o_temp)
    
    return new_t, new_s

def redist_mass(MOM,SIS,row,col):
    # This function redistributes mass between a target cell and a series of halo
    # cells during the process of creating or removing ocean cells. This function
    # should only be run *after* redist_tracer, as that subroutine relies on old
    # h field.
    # 
    # In:  chng_mask         - Mask indicating if cell is to change 
    #                          from land>ocean or vice versa (1 or -1)
    #      h_size_mask       - Halo size for cell to be altered (radius in # of cells)
    #      h_oce             - Ocean vertical layer thickness (m)             
    #      cell_area         - Tracer cell area from ocean model (m2)
    #      depth             - Updated ocean cell depth (m)
    #      h_ice             - Ice mass (kg/m2)
    #      ice_frac          - Ocean ice fraction
    #      h_sno             - Snow mass (kg/m2)
    #      o_temp            - Ocean potential temp. (degC)
    #      o_salt            - Ocean salinity (ppt)
    # Out: new_h             - Updates ocean vertical layer thickness (m)
    #      o_temp            - Updated temp. field (degC)
    #      o_salt            - Updated salinity field  (ppt)
    #     
    # In the process of changing h of halo cells, a correction is applied to 
    # tracers to avoid spurious creation/ removal of energy/ salt 
    # (and later, other relevant tracers).
    # 
    
    size         = MOM.h_size_mask[row,col].astype(int) # Def. halo radius
    c_wgts, halo = cell_weights(row,col,size,MOM)       # Get weights for re-distribution
    ice_mass = 0.0; sno_mass = 0.0                      # Initialise remaining vars
    # 1. Are we filling or emptying a cell?
    if MOM.chng_mask[row,col] == -1: # Emptying a cell (ocean -> Land)
        # 2. Calculate the mass to remove from the cell
        # Ice model
        for i in range(SIS.h_ice.shape[0]):
            ice_mass = ice_mass + SIS.h_ice[i,row,col]*MOM.cell_area[row,col]*SIS.ice_frac[i,row,col]; 
            sno_mass = sno_mass + SIS.h_sno[i,row,col]*MOM.cell_area[row,col]*SIS.ice_frac[i,row,col];
        # Ocean model
        h_cell          = MOM.h_oce[:,row,col];
        h_sum           = math.fsum(h_cell);             # All h levels exist at every point, just very small
        sea_mass        = np.multiply(h_sum,MOM.cell_area[row,col])
        # Total
        tot_mass        = math.fsum([ice_mass,sno_mass,sea_mass])
        # 3. Redistribute mass to halo cells
        delta_mass      = np.multiply(c_wgts,tot_mass)        # mass going to each halo cell
        delta_h         = np.divide(delta_mass,MOM.cell_area) # converted to change in h thickness
        newh            = h2vgrid(delta_h,MOM.h_oce,MOM)      # Apply del h to each depth level in halo
        
        # 4. Correct tracer concentrations for change in h thickness
        o_temp, o_salt = tracer_correct(MOM.h_oce,newh,MOM) # This alters exising tracer fields and returns
                                                     # values corrected for change in mass
        # 5. Remove ocean cell from h field
        newh[:,row,col] = 0                          # Remove mass from cell and set layers to be land
        SIS.ice_frac[:,row,col] = 0                  # Remove ice and snow from cell
    elif MOM.chng_mask(row,col) == 1: # Filling a cell (land -> ocean)
        # 2. Check surrounding SSH
        eta_mean = halo_eta(MOM.eta,row,col);
        # 3. Subtract from topog depth
        init_h = eta_mean + MOM.depth[row,col];
        # 4. Create new water column with default z levels
        newh[:,row,col] = newcell(init_h)
        # 5. How much mass is required to initialise this cell?
        tot_mass        = init_h*MOM.cell_area[row,col];
        # 6. Remove mass from halo cells
        delta_mass      = c_wgts*tot_mass                 # mass coming from each halo cell
        delta_h         = delta_mass/MOM.cell_area        # converted to change in h thickness
        newh            = h2vgrid(-delta_h,MOM.h_oce,MOM) # Apply del h to each depth level in halo
        # 7. Correct tracer concentrations for change in h thickness
        o_temp, o_salt = tracer_correct(MOM.h_oce,newh) # This alters exising tracer fields and returns 
                                                        # values corrected for change in mass
    return newh, o_salt, o_temp

def newcell(MOM,hsum):
    # This function initialises a new ocean cell of depth 'hsum' and discretises
    # it into defualt h layers, as defined by the models v-grid
    # In:  MOM      - Ocean data structure
    #      hsum     - Total depth of new water column (excluding 'land depth'
    #                 where land has h of 0.001m)
    #      vgrid    - The default model vertical grid spacing (m)
    #
    # Out: h        - The newly created ocean cell, vertically discretised using 
    #                 model's vgrid.
     
    n_lvls = 0;
    tmp = hsum;
    for i in range(MOM.vgrid.shape[0]):
        tmp -= MOM.vgrid[i]
        n_lvls += 1;
        if tmp <= 0:
            break
        
    h = np.full(MOM.vgrid.shape[0], 0.001) # 0.001 is the layer thickness of 'land'
    if n_lvls == 1:
        h = hsum
    else:
        for i in range(0,n_lvls-1):
            h[i] = MOM.vgrid[i]
        resid = sum(MOM.vgrid[0:n_lvls]) - hsum # some frac of the lowest gridcell is land
        h[n_lvls-1] = MOM.vgrid[n_lvls-1] - resid
    
    # Run check to make sure operation worked correctly
    if math.floor(sum(h)) != hsum:
        raise ValueError(str('Water column not properly discretised. Total depth is '+str(sum(h)) \
                             + ', where it should be '+str(hsum)))
    
    return h

def sum_ocean_enth(row,col,T,h,MOM,flag):
    # This function calculates the total energy in seawater in a single cell.
    # The method used is taken from MOM_sum_ouput.F90, line 675. Units are in J
    # Option included for a single cell used during tracer_redist
    # In:  MOM                       - Ocean data structure 
    #      h                         - Ocean vertical layer thickness (m) 
    #      T                         - Ocean potential temp. (degC)
    #      cell_area (gloabal var)   - Tracer cell area from ocean model (m2)
    #      C_P     (gloabal var)     - Heat capacity of seawater (J/K/kg)
    #      H_to_m                    - Grid cell thickness to m conversion factor
    #                                  Default is 1 for Boussinesq numerics
    #      flag                      - Indicates whether or not operation is for
    #                                  a single cell (true) or the entire domain (false)
    # Out: total                     - Total ocean energy (J)
        
    total = 0
    if flag:
        for i in range(MOM.grid_z):
            total = total + (MOM.C_P * T[i,row,col]) * \
                    (h[i,row,col]*(MOM.H_to_m * MOM.cell_area[row,col]))    
    else:
        for i in range(MOM.grid_z):
            total = total + np.sum((MOM.C_P * T[i,:,:]) * \
                    (h[i,:,:]*(MOM.H_to_m * MOM.cell_area[row,col])))
    return total

def sum_ocean_salt(row,col,S,h,MOM,flag):
    # This function calculates the total amount of salt in a single grid cell.
    # The method used is taken from MOM_sum_ouput.F90, line 673. Units are in Kg
    # In:  MOM                       - Ocean data structure           
    #      h                         - Ocean vertical layer thickness (m) 
    #      S                         - Ocean salinity (ppt)
    #      H_to_kg_m2                - Grid cell thickness to mass conversion factor
    #                                  Default is 1 for Boussinesq numerics
    #      flag                      - Indicates whether or not operation is for
    #                                  a single cell (true) or the entire domain (false)
    # Out: total                     - Total salinity (kg)
    
    total = 0
    if flag:
        for i in range(MOM.grid_z):
            total = total + S[i,row,col] * \
                    (h[i,row,col]*(MOM.H_to_m * MOM.cell_area[row,col]))
    else:    
        for i in range(MOM.grid_z):
            total = total + np.sum(S[i,:,:] * \
                    (h[i,:,:]*(MOM.H_to_m * MOM.cell_area)))      
    return total*0.001

def redist_tracers(MOM,SIS,row,col,tracer=''):
    # This function redistributes tracers between a cell and target halo cells.
    # It works with energy and salinity only at this point.
    # All data are imported via the MOM and SIS data structures
    # Note that when creating new ocean cells, it is a requirement that we 
    # first initialise its vertical grid structure, complicating matters when
    # we eventually redistribute mass.
    #
    # Functions from SIS2_funcs.py are utilised.
    # In:  row         - latitude index
    #      col         - longitude index
    #      tracer      - The tracer to be redistributed
    #      h_ice       - Ice mass (kg/m2)
    #      h_sno       - Snow mass (kg/m2)
    #      ice_frac    - fraction of cell covered in ice (0-1)
    #      h_oce       - Ocean vertical layer thickness (m) 
    #      cell_area   - Tracer cell area from ocean model (m2)
    #      e_ice       - Enthalpy of ice (J)
    #      e_sno       - Enthalpy of snow (J)
    #      o_temp      - Ocean potential temp. (degC)
    #      o_salt      - Ocean salinity (ppt)
    #      s_ice       - Salinity of sea ice (g/kg)
    #      chng_mask   - Mask indicating if cell is to change 
    #                    from land>ocean or vice versa (1 or -1)
    #      h_size_mask - Halo size for cell to be altered (radius in # of cells)
    #      h_oce       - Ocean vertical layer thickness (m)
    #      h_ice       - Ice mass (kg/m2)
    #      ice_frac    - Ocean ice fraction
    #      h_sno       - Snow mass (kg/m2)
    # Out: o_temp_new  - Updated temperature field
    #      o_salt_new  - Updated salinity field
    
    size         = MOM.h_size_mask[row,col].astype(int); # Def. halo radius
    c_wgts, halo = cell_weights(row,col,size,MOM);       # Get weights for re-distribution
    ice_tot      = 0;
    
    # Is ocean becoming land?
    if MOM.chng_mask[row,col] == -1: # ocean -> Land
        # Which tracer
        if tracer == 'temp':
            # 1.1 Sum enthalpy in sea ice and snow (the latter only has one layer) - During test, should equal -3.144585077871032E+21 across whole domain
            if np.sum(SIS.ice_frac[:,row,col],0) > 0:
                ice_tot = sum_ice_enth(row,col,SIS.e_ice,SIS.e_sno,SIS.h_ice, \
                                       SIS.h_sno,SIS.ice_frac,MOM.cell_area, \
                                       SIS.s_ice,MOM.H_to_m)
        
            # 1.2 Sum energy in sea water
            oce_tot = sum_ocean_enth(row,col,MOM.o_temp,MOM.h_oce,MOM,"False")
            tot_heat = oce_tot + ice_tot
        
            # 1.3 Calculate how much energy each halo cell recieves 
            delta_heat        = c_wgts*tot_heat          # energy going to each halo cell
            o_temp            = heat2cell(delta_heat,MOM)# Put energy in halo cells
            o_temp[:,row,col] = 0                        # This cell is now land
            new_tracer        = o_temp
        elif tracer == 'salt':
            # 1.4 Sum salinity in sea ice (snow has no salt content)
            if np.sum(SIS.ice_frac[:,row,col],0) > 0:
                ice_tot = sum_ice_sal(row,col,SIS.ice_frac,MOM.cell_area,SIS.h_ice, \
                                      SIS.s_ice,SIS.H_to_kg_m2,SIS.nk_ice)
            
            # 1.5 Sum salt in ocean water
            oce_tot = sum_ocean_salt(row,col,MOM.o_salt,MOM.h_oce,MOM,True)
            tot_salt = oce_tot + ice_tot
            
            # 1.6 Calculate how much salt each halo cell recieves
            delta_salt        = c_wgts*(tot_salt*1000) # salt going to each halo cell
            o_salt            = tracer2cell(delta_salt,MOM.o_salt,MOM,'salt')
            o_salt[:,row,col] = 0                        # This cell is now land
            new_tracer        = o_salt
    # Or is land becoming ocean?
    elif MOM.chng_mask[row,col] == 1: # land -> ocean
        # Before we move tracers around, we need to initialise the cell to have
        # h levels and mass information. This follows the same steps as in the 
        # mass redistribution code.
        # 2.1 Check surrounding SSH
        eta_mean = halo_eta(MOM.eta,row,col)
        # 2.2 Subtract from topog depth
        init_h = eta_mean + MOM.depth[row,col]
        # 2.3 Create new water column with default z levels
        newh = newcell(init_h)
        
        # Which tracer?
        if tracer == 'temp':
            # 2.5 Create new temperature profile based on surrounding cells
            T = np.full(MOM.grid_z,0)
            for i in range(MOM.grid_z):
                T[i],_,_ = halo_calc(row,col,MOM.o_temp[i,:,:],MOM,1,'mean')
            # 2.6 Calculate total energy for new cell
            T_tot = sum_ocean_enth(row,col,T,newh,MOM,"True")
            # 2.7 Calculate how much energy must be removed from each halo cell
            delta_heat      = c_wgts*T_tot
            o_temp          = heat2cell(-delta_heat,MOM.o_temp,MOM)
            new_tracer      = o_temp
        elif tracer == 'salt':
            # 2.8 Create new salinity profile based on surrounding cells
            S = np.full(MOM.grid_z,0)
            for i in range(MOM.grid_z):
                S[i],_,_ = halo_calc(row,col,MOM.o_salt[i,:,:],MOM,1,'mean')
            # 2.9 Calculate total salt for new cell
            S_tot = sum_ocean_salt(row,col,S,newh,MOM,True)
            # 2.10 Calculate how much salt must be removed from each halo cell
            delta_salt      = c_wgts*S_tot          # salt going to each halo cell
            o_salt          = tracer2cell(-delta_salt,MOM.o_salt,MOM,'salt')
            new_tracer      = o_salt
    elif MOM.chng_mask[row,col] == 0:
        print('Whoops, you called this function on cells that are not changing!')
        return
    return new_tracer
    
def heat2cell(delta_E,MOM):
    # This function adds or removes energy to/from a grid cell and updates the 
    # temperature field. The approach used is designed to spread out the change
    # in T equally throughout the water column (so, proportional to h) to try and
    # maintain the vertical temperature structure of the grid cell.
    # In:  delta_E   - The energy being added/removed from a series of halo cells 
    #      o_temp    - The temperature field to be altered
    #      cell_area - Tracer cell area
    #      h_oce     - Ocean layer thickness
    # Out: T_new     - Updated temperature field
    wght = np.full(MOM.o_temp.shape,np.nan)
    mass = np.full(MOM.o_temp.shape,np.nan)
    T_new = np.full(MOM.o_temp.shape,np.nan)
    h_oce = cp.deepcopy(MOM.h_oce); h_oce[:,MOM.o_mask==0] = np.nan
    for i in range(MOM.grid_z):
        mass[i,:,:] = h_oce[i,:,:] * MOM.cell_area
        wght[i,:,:] = h_oce[i,:,:]/np.sum(h_oce,0)
        T_new[i,:,:] = MOM.o_temp[i,:,:] + delta_E*wght[i,:,:]/(mass[i,:,:]*MOM.C_P)  
    T_new[np.isnan(T_new)] = 0    
    return T_new

def tracer2cell(delta_t,tracer,MOM,var=''):
    # This function adds or removes a tracer to/from a grid cell and updates
    # that tracer's field. Changes in the tracer are spread evenly over the water
    # column in an attempt to maintain the vertical profiles of the tracer 
    # concentration. This may not be appropriate for all tracers, but for now,
    # only code for salt is implemented.
    # In:  delta_t   - Change in some tracer (kg for salt)
    #      tracer    - The tracer field to be altered
    #      MOM       - The ocean data structure
    #      var       - Name of the tracer being altered (string)
    # Out: new       - Updated tracer field (ppt for salt)
    h_oce = cp.deepcopy(MOM.h_oce); h_oce[:,MOM.o_mask==0] = np.nan
    wght = np.full([MOM.grid_z,MOM.grid_y,MOM.grid_x],np.nan)
    mass = np.full(wght.shape,np.nan); new = np.full(wght.shape,np.nan)
    for i in range(MOM.grid_z):
        mass[i,:,:] = np.multiply(h_oce[i,:,:],MOM.cell_area)
        wght[i,:,:] = np.divide(h_oce[i,:,:],np.sum(h_oce,0))
        if var == 'salt':
            new[i,:,:] = tracer[i,:,:] + (np.divide(np.multiply(delta_t,wght[i,:,:]),mass[i,:,:]))
    return new

def chk_ssh(h_sum,MOM): # Unfinished- double check everything
    # This function checks for sea surface height gradients after redistribution
    # of mass and tracers between ocean cells. For the moment, if neighbouring cells 
    # have a SSH difference of 1m or more, they will get flagged.
    # In:  h_sum       - Total water column height (m)
    #      bathy       - Ocean bathymetry depth (m)
    # Out: lrg_grad    - Field indicating which cells have excessive SSH gradients
    lrg_grd = np.full(h_sum.shape)
    for i in range(h_sum.shape[0]):
        for j in range(h_sum.shape[1]-1): # Cells on other side of zonal split handled later
            if MOM.depth_new[i,j] != 0 and MOM.depth_new[i,j+1] != 0:
                diff_x = abs((h_sum[i,j] - MOM.depth_new[i,j]) - (h_sum[i,j+1] - MOM.depth_new[i,j+1]))
                if diff_x >= 1:
                    lrg_grd[i,j] = 1
    for i in range(h_sum.shape[0]-1):
        for j in range(h_sum.shape[1]):
            if MOM.depth_new[i,j] != 0 and MOM.depth_new[i+1,j] != 0:
                diff_y = abs((h_sum[i,j] - MOM.depth_new[i,j]) - (h_sum[i+1,j] - MOM.depth_new[i+1,j]))
                if diff_y >= 1:
                    lrg_grd[i,j] = 1
    return lrg_grd
    
def chk_conserv(OLD,SIS,MOM,data=""):
    # This function checks conservation for a given quantity/tracer (currently
    # only mass, energy, and salt). Function returns the error between new and 
    # old fields and whether or not this is of an acceptable size.
    # In:  h_oce (new/old)    - Fields of ocean column thickess prior to and after
    #                           redistribution (m)
    #      h_ice (new/old)    - Fields of ice thickess prior to and after
    #                           redistribution (m)
    #      h_sno (new/old)    - Fields of snow thickess prior to and after
    #                           redistribution (m)
    #      o_temp (new/old)   - Fields of ocean temperature prior to and after
    #                           redistribution (deg C)
    #      o_salt (new/old)   - Fields of ocean salinity prior to and after
    #                           redistribution (ppt)
    #      s_ice (new/old)    - Fields of sea ice salinity prior to and after
    #                           redistribution (ppt)
    #      e_ice (new/old)    - Fields of sea ice enthalpy prior to and after
    #                           redistribution (J)
    #      e_sno (new/old)    - Fields of snow enthalpy prior to and after
    #                           redistribution (J)
    #      ice_frac (new/old) - Fraction of a cell with sea ice prior to and after
    #                           redistribution (J)
    #      MOM                - Ocean data structure (contains all updated values)
    #      SIS                - Sea ice data structure (contains all updated values)
    #      OLD                - Data structure with values prior to redistribution
    #      data               - Name of the variable we are checking (mass,temp,salt)
    total_old = 0; total_new = 0
    err_tol = 1e-14 # Error tolerance for conservation
    if data == 'temp':
        # Calculate total ocean enthalpy
        total_old += sum_ocean_enth(None,None,OLD.o_temp,OLD.h_oce,MOM,False)
        total_new += sum_ocean_enth(None,None,MOM.o_temp,MOM.h_oce,MOM,False)
        # Calculate total ice enthalpy
        for row in range(MOM.grid_y):
            for col in range(MOM.grid_x):
                total_old += sum_ice_enth(row,col,OLD.e_ice,OLD.e_sno,OLD.h_ice, \
                                   OLD.h_sno,OLD.ice_frac,MOM.cell_area,SIS.s_ice,SIS.H_to_kg_m2)
                total_new += sum_ice_enth(row,col,SIS.e_ice,SIS.e_sno,SIS.h_ice, \
                                   SIS.h_sno,SIS.ice_frac,MOM.cell_area,SIS.s_ice,SIS.H_to_kg_m2)
    elif data == 'salt':
        total_old = sum_ocean_salt(None,None,OLD.o_salt,OLD.h_oce,MOM,False)
        total_new = sum_ocean_salt(None,None,MOM.o_salt,MOM.h_oce,MOM,False)
        for row in range(MOM.grid_y):
            for col in range(MOM.grid_x):
                total_old += sum_ice_sal(row,col,OLD.ice_frac, \
                                           MOM.cell_area,OLD.h_ice,SIS.s_ice, \
                                           SIS.H_to_kg_m2,SIS.nk_ice)
                total_new += sum_ice_sal(row,col,SIS.ice_frac, \
                                           MOM.cell_area,SIS.h_ice,SIS.s_ice, \
                                           SIS.H_to_kg_m2,SIS.nk_ice)
        
    elif data == 'mass':
        ice_mass_old = 0; sno_mass_old = 0; ice_mass_new = 0; sno_mass_new = 0
        # Ice model
        for i in range(SIS.h_ice.shape[0]):
            ice_mass_old += np.sum(OLD.h_ice[i,:,:]*MOM.cell_area[:,:]*OLD.ice_frac[i,:,:]) 
            sno_mass_old += np.sum(OLD.h_sno[i,:,:]*MOM.cell_area[:,:]*OLD.ice_frac[i,:,:])
            ice_mass_new += np.sum(SIS.h_ice[i,:,:]*MOM.cell_area[:,:]*SIS.ice_frac[i,:,:])
            sno_mass_new += np.sum(SIS.h_sno[i,:,:]*MOM.cell_area[:,:]*SIS.ice_frac[i,:,:])
        # Ocean model
        sea_mass_old = np.sum(np.sum(OLD.h_oce,0)*MOM.cell_area)
        sea_mass_new = np.nansum(np.sum(MOM.h_oce,0)*MOM.cell_area)
        # Total
        total_old    = ice_mass_old + sno_mass_old + sea_mass_old
        total_new    = ice_mass_new + sno_mass_new + sea_mass_new
        
    if math.isclose(total_old,total_new,abs_tol=err_tol): # Past 1e-14, choice of summing algorith begins to impact
        print(data+' is conserving within a tolerance of '+str(err_tol))
    else:
        tot_diff = total_old - total_new
        raise ValueError(str(data+' is not conserving. Total_old - Total_new = '+str(tot_diff)))
    return
################################# Main Code ###################################
    
def redist_vals(MOM,SIS,OLD,FLAGS):
# This function uses the chng_mask variable to either initialise or remove 
# ocean cells from the mondel domain. The process is currently designed to 
# conservatively redistribute mass and tracers of energy (temperature) and
# salt between the cell that requires modification and a series of surrounding
# halo cells. 
# 
#  In: chng_mask    - Mask of cells to change
#      h_size_mask  - Halo size mask
#      h_ice        - Ice thickness
#      h_sno        - Snow thickness
#      ice_frac     - Ice fraction
#      s_ice        - Salinity of sea ice in g/kg - it's a fixed value
#      e_ice        - Enthalpy of sea ice in J
#      e_sno        - Enthalpy of snow in J
#      cell_area    - Area of h (tracer) cells
#      h_oce        - Ocean layer thickness
#      o_temp       - Ocean potential temperature (deg C)
#      o_salt       - Ocean salinity (ppt)
#      o_mask_redist- Updated ocean mask (cells becoming ocean are marked as land)
#      C_P          - The heat capacity of seawater in MOM6 (J kg-1 K-1)
#      H_to_kg_m2   - grid cell to mass conversion factor (1 by default)
    
    t_start = time.time()
    # Variable pre-processing
    scale_fac    = 20                                             # A scaling factor determining how many times the surface are of the target cell is required for redist.  
    tmp          = cp.deepcopy(MOM.h_oce);                        # Dummy variable
    tmp          = np.where(tmp>0.01, tmp, np.nan);               # While all layers always exist, we only want them if they have 
    MOM.h_lvls   = np.nansum(tmp/tmp,axis=0).astype(int); del tmp # significant mass (ie: non-0 layers, defined has having h > 0.001)
    MOM.wght     = np.full(MOM.h_oce.shape,0).astype(float)       # Cell weights for mass redistribution
    for i in range(MOM.h_lvls.shape[0]):                          # we add mass to a layer proportional to its thickness
        for j in range(MOM.h_lvls.shape[1]):
            lvls = MOM.h_lvls[i,j];
            for k in range(lvls):
                MOM.wght[k,i,j]= MOM.h_oce[k,i,j]/ \
                MOM.h_oce[:MOM.h_lvls[i,j],i,j].sum(0)                
    MOM.o_mask_redist = cp.deepcopy(MOM.o_mask_new)               # New ocean cells masked out. We need a seperate ocean mask for the redistribution  
    MOM.o_mask_redist[MOM.chng_mask==1] = 0;                      # code, as we cannot allow cells that are yet to be initialised act as halo cells.
    
        
    # 1: Set up and check halos for all change points, making the halos bigger when required
    # Here size is increased by 1 in the while loop, but that means it will be too big when we 
    # find the correct value. Therefore, when we assign it to h_size_mask, we need to subtract 1 again.
    for i in range(MOM.chng_mask.shape[0]):
        for j in range(MOM.chng_mask.shape[1]):
            if MOM.chng_mask[i,j] != 0:
            #if o_mask_new[i,j] > 0: # For testing all possible ocean cells
                halo_sum = 0;
                size = 1;
                while halo_sum < MOM.cell_area[i,j]*scale_fac:
                    halo_sum,_,_ = halo_calc(i,j,MOM.cell_area,MOM,size,'sum');
                    size+=1
                MOM.h_size_mask[i,j] = size-1; del halo_sum
    
   # 2: Redistribute mass and tracers in and out of change points
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if MOM.chng_mask[i,j] != 0:
                MOM.o_temp = redist_tracers(MOM,SIS,i,j,'temp')
                MOM.o_salt = redist_tracers(MOM,SIS,i,j,'salt') 
                MOM.h_oce, MOM.o_salt, MOM.o_temp = redist_mass(MOM,SIS,i,j)
                # Change NaN vals in T and S fields back to 0's
                MOM.o_temp[np.isnan(MOM.o_temp)] = 0
                MOM.o_salt[np.isnan(MOM.o_salt)] = 0
    if FLAGS.verbose:
        print('Checking for conservation of mass...')
        chk_conserv(OLD,SIS,MOM,'mass')
        print('Checking for conservation of energy...')
        chk_conserv(OLD,SIS,MOM,'temp')
        print('Checking for conservation of salt...')
        chk_conserv(OLD,SIS,MOM,'salt')
#        print('Redistribution of mass and tracers complete.' \
#              '\n Error in mass   = '+str(err_mass) + \
#              '\n Error in energy = '+str(err_T) +\
#              '\n Error in salt   = '+str(err_S))
    FLAGS.t_redist = time.time() - t_start
    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    