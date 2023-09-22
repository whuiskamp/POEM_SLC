#!/usr/bin/env python3

# This is part of the SLC tool
# Written by Willem Huiskamp
# 2023

# These scripts update any fields derived from others we may have altered 
# while running the SLC tool. Updated fields are then saved to the restart files.
# Additionally, diagnostics are also saved including a 'change mask' which records 
# changes to the land-sea mask is also saved and updated (appended to) if this mask 
# was altered during this iteration of the SLC tool.

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
    if os.path.isfile(str(path)+'history/slc_diag.nc') == False: # Does the slc tool diagnostic file exist yet? If not, make one
        slc = CDF(str(path)+'history/slc_diag.nc', 'w')
        # Create dimensions for vars.
        slc.createDimension('lonh', MOM.lon.shape[1])
        slc.createDimension('lath', MOM.lat.shape[0])
        slc.createDimension('time', None)
        slc.createVariable('lonh', 'f8', ('lonh'));
        slc.variables['lonh'].units = 'degrees_east'
        slc.variables['lonh'].cartesian_axis = 'X'
        slc.variables['lonh'].long_name = 'T-cell longitude'
        slc.variables['lonh'][:] = grid.variables['lonh'][:]
        slc.createVariable('lath', 'f8', ('lath'));
        slc.variables['lath'].units = 'degrees_north'
        slc.variables['lath'].cartesian_axis = 'Y'
        slc.variables['lath'].long_name = 'T-cell latitude'
        slc.variables['lath'][:] = grid.variables['lath'][:]
        # slc.createVariable('time') #Not sure we want this?? Can we carry through a time var?
        ########## Define variables and fill data ##########
        # Cell change mask
        slc.createVariable('chng_mask', 'f8', ('time','lath','lonh'))
        slc.variables['chng_mask'].units = 'none'
        slc.variables['chng_mask'].long_name = 'Mask indicating changing ocean/ land cells'
        slc.variables['chng_mask'][:,:,:] = MOM.chng_mask[:,:]
        # New ocean mask
        slc.createVariable('o_mask_new', 'f8', ('time','lath','lonh'))
        slc.variables['o_mask_new'].units = 'none'
        slc.variables['o_mask_new'].long_name = 'Updated ocean mask'
        slc.variables['o_mask_new'][:,:,:] = MOM.o_mask_new[:,:]
        # Conservation errors
        # Enthalpy
        slc.createVariable('err_enth','f8','time')
        slc.variables['err_enth'].units = 'Joules'
        slc.variables['err_enth'].long_name = 'Error in total ocean-ice enthalpy after redistribution'
        slc.variables['err_enth'][:] = MOM.err_enth
        # Salt
        slc.createVariable('err_salt','f8','time')
        slc.variables['err_salt'].units = 'kg'
        slc.variables['err_salt'].long_name = 'Error in total ocean-ice salt after redistribution'
        slc.variables['err_salt'][:] = MOM.err_salt
        # Mass
        slc.createVariable('err_mass','f8','time')
        slc.variables['err_mass'].units = 'kg'
        slc.variables['err_mass'].long_name = 'Error in total ocean-ice mass after redistribution'
        slc.variables['err_mass'][:] = MOM.err_mass
        
        slc.description = "This file contains diagnostics from the SLC tool including: a mask which  \
                        lists cells that should change from ocean to land (-1) or land to ocean (1), \
                        an updated land-sea mask and error terms for conservation checking."
        slc.history = "Created " + time.ctime(time.time())
        slc.source = "Created using https://github.com/whuiskamp/POEM_SLC"
        slc.close()
    else:
        # If the file already exists...
        slc = CDF(str(path)+'history/slc_diag.nc','a');
        slc.variables['chng_mask'][:,:,:] = MOM.chng_mask[:,:]
        slc.variables['o_mask_new'][:,:,:] = MOM.o_mask_new[:,:]
        slc.variables['err_enth'][:] = MOM.err_enth
        slc.variables['err_salt'][:] = MOM.err_salt
        slc.variables['err_mass'][:] = MOM.err_mass
        slc.close()
    
    # Write appropriate diagnostics to restart's metadata
    # Open files to write to        
    new_bathy  = CDF(work_dir + '/topog.nc','r+')        
    MOM6_rest  = CDF(work_dir + '/MOM.res.nc','r+')
    SIS2_rest  = CDF(work_dir + '/ice_model.res.nc','r+')
    Omask      = CDF(work_dir + 'ocean_mask.nc','r+')

    # Save new bathymetry/ mask data        
    new_bathy.variables['depth'][:,:] = MOM.depth_new
    Omask.variables['mask'][:,:].data = MOM.o_mask_new

    # Update MOM6 restart
    MOM6_rest.variables[''] 

    # Update SIS2 restart
    SIS2_rest.variables
            
            
            
            
            
            