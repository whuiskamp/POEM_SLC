#!/usr/bin/env python3

# This is part of the SLC tool
# Written by Willem Huiskamp
# 2023

# These scripts update any fields derived from others we may have altered 
# while running the SLC tool. Updated fields are then saved to the restart files.
# Additionally, a 'change mask' which records changes to the land-sea mask is also
# saved and updated (appended to) if this mask was altered during this iteration
# of the SLC tool.

import time
from netCDF4 import Dataset as CDF
import numpy as np

def update_bathy(MOM):
    # This function saves the new bathymetry field to the new topog.nc topography file for MOM6.
    # This function assumes no changes were made to the land sea mask, but that the bathymetry was modified due
    # to forcing from the solid earth model. 
    new_bathy  = CDF(work_dir + '/topog.nc','r+')
    new_bathy.variables['depth'][:,:] = MOM.depth_new[:,:]
    new_bathy.close()


def prep_fields(MOM,SIS):
    # This function alters relevant secondary fields in the MOM/SIS restart fields that should change due to changes in interface thickness, 
    # temperature and salinity changes.
    # In:  MOM   - The ocean data structure containing all relevant restart variables
    #      SIS   - The sea ice data structure containing all relevant restart variables
    # Out: Modified fields in data structures

    # 1. The following fields always need to be changed in the restart file
    MOM.h2      = np.square(MOM.h_oce)         # Auxiliary ocean layer thickness (m)
    MOM.uh      = np.multiply(MOM.u,MOM.h_oce) # Zonal thickness flux (m3 s-1)
    MOM.vh      = np.multiply(MOM.v,MOM.h_oce) # Meridional thickness flux (m3 s-1)

    # First element of each variable is essentially a land-sea mask. Update it here.
    SIS.rough_mom[0,:,:]   = MOM.o_mask_new
    SIS.rough_heat[0,:,:]  = MOM.o_mask_new
    SIS.rough_moist[0,:,:] = MOM.o_mask_new

    # Change NaN values back to 0 as the model requires...


    # 2. These fields either do not need to be changed, or can be optionally changed. Uncomment to use. 
    # MOM.u        # Zonal velocity (m s-1)
    # MOM.v        # Meridional velocity (m s-1)
    # MOM.u2       # Auxiliary zonal velocity (m s-1)
    # MOM.v2       # Auxiliary meridional velocity (m s-1)


def write_rest(MOM,SIS,FLAGS):
    # Load new restart files
    # Out: Restart files ice_model.res.nc and MOM.res.nc; bathymetry file topog
    #    Change all ocean/ ice fields that need to be changed

    # Save new fields to netCDF
    # Write change mask to netCDF
    if os.path.isfile(str(path)+'history/change_mask.nc') == False: # Does the change mask file exist yet? If not, make one
        chg_msk = CDF(str(path)+'history/change_mask', 'w')
        # Create dimensions for vars.
        chg_msk.createDimension('lonh', MOM.lon.shape[1])
        chg_msk.createDimension('lath', MOM.lat.shape[0])
        chg_msk.createDimension('time', None)
        chg_msk.createVariable('lonh', 'f8', ('lonh'));
        chg_msk.variables['lonh'].units = 'degrees_east'
        chg_msk.variables['lonh'].cartesian_axis = 'X'
        chg_msk.variables['lonh'].long_name = 'T-cell longitude'
        chg_msk.variables['lonh'][:] = grid.variables['lonh'][:]
        chg_msk.createVariable('lath', 'f8', ('lath'));
        chg_msk.variables['lath'].units = 'degrees_north'
        chg_msk.variables['lath'].cartesian_axis = 'Y'
        chg_msk.variables['lath'].long_name = 'T-cell latitude'
        chg_msk.variables['lath'][:] = grid.variables['lath'][:]
        # chg_msk.createVariable('time') #Not sure we want this?? Can we carry through a time var?
        # Define variables and fill data
        # Cell change mask
        chg_msk.createVariable('chng_mask', 'f8', ('time','lath','lonh'))
        chg_msk.variables['chng_mask'].units = 'none'
        chg_msk.variables['chng_mask'].long_name = 'Mask indicating changing ocean/ land cells'
        chg_msk.variables['chng_mask'][-1,:,:] = MOM.chng_mask[:,:]
        # New ocean mask
        chg_msk.createVariable('o_mask_new', 'f8', ('time','lath','lonh'))
        chg_msk.variables['o_mask_new'].units = 'none'
        chg_msk.variables['o_mask_new'].long_name = 'Updated ocean mask'
        chg_msk.variables['o_mask_new'][-1,:,:] = MOM.o_mask_new[:,:]
        
        chg_msk.description = "This file contains a mask which lists cells that should change from \
                        ocean to land (-1) or land to ocean (1) as well as an updated land-sea mask"
        chg_msk.history = "Created " + time.ctime(time.time())
        chg_msk.source = "Created using https://github.com/whuiskamp/POEM_SLC"
        chg_msk.close()
    else:
        # If the file already exists...
        chg_msk = CDF(str(path)+'history/change_mask.nc','r+');
        chg_msk.variables['chng_mask'][-1,:,:] = MOM.chng_mask[:,:]
        chg_msk.variables['o_mask_new'][-1,:,:] = MOM.o_mask_new[:,:]
        chg_msk.close()
    # Write appropriate diagnostics to restart's metadata
    
    # Open files to write to        
    new_bathy  = CDF(work_dir + '/topog.nc','r+')        
    MOM6_rest  = CDF(work_dir + '/MOM.res.nc','r+')
    SIS2_rest  = CDF(work_dir + '/ice_model.res.nc','r+')
    Omask      = CDF(work_dir + 'ocean_mask.nc','r+')
            
            
            
            
            
            
            
            
            
            