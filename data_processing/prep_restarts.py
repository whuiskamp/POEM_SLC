# This script updates restart files with new variables that have been generated 
# in previous scripts. In addition, secondary fields that need to be altered
# based on redistribution of values/ changes in the land sea mask are also 
# altered here before being saved to appropriate restarts.

import time
from netCDF4 import Dataset as CDF
import numpy as np

__author__     = "Willem Huiskamp"
__copyright__  = "Copyright 2020"
__credits__    = ["Willem Huiskamp", ""]
__license__    = "GPLv3"
__version__    = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__      = "huiskamp@pik-potsdam.de"
__status__     = "Prototype"

def prep_fields(MOM,SIS,FLAGS):
    # This function alters relevant secondary fields in the MOM/SIS restart fields that should change due to changes in interface thickness, 
    # temperature and salinity changes.
    # In:  MOM   - The ocean data structure containing all relevant restart variables
    #      SIS   - The sea ice data structure containing all relevant restart variables
    #      FLAGS - The option flags and directory information data structure
    # Out: Restart files ...

    # 1. The following fields always need to be changed in the restart file
    MOM.h_oce   =  # Ocean layer thickness (m)
    MOM.o_temp  =  # Ocean potential temperature (deg C)
    MOM.o_salt  =  # Ocean salinity (ppt)
    MOM.ave_eta =  # Time-average sea surface height (m)
    MOM.eta     =  # Sea surface height (m)    
    MOM.h2      = np.square(MOM.h_oce) # Auxiliary ocean layer thickness (m)
    MOM.uh      = np.multiply(MOM.u,MOM.h_oce) # Zonal thickness flux (m3 s-1)
    MOM.vh      = np.multiply(MOM.v,MOM.h_oce) # Meridional thickness flux (m3 s-1)

    # 2. These fields either do not need to be changed, or can be optionally changed (this is still under development)
    #    Uncomment to use. 
    MOM.u        # Zonal velocity (m s-1)
    MOM.v        # Meridional velocity (m s-1)
    MOM.u2       # Auxiliary zonal velocity (m s-1)
    MOM.v2       # Auxiliary meridional velocity (m s-1)


def write_rest(MOM,SIS,FLAGS):
    # Load new restart files
    # Out: Restart files ice_model.res.nc and MOM.res.nc; bathymetry file topog
    #    Change all ocean/ ice fields that need to be changed

    # Save new fields to netCDF
    # Write change mask to netCDF
    if FLAGS.test:
        chg_msk = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc', 'w')
    else:
        if os.path.isfile(str(path)+'history/change_mask.nc') == False: # Does the change mask file exist yet? If not, make one
            chg_msk = CDF(str(path)+'history/change_mask', 'w')
            # Create dimensions for vars.
            chg_msk.createDimension('lonh', lon.shape[1])
            chg_msk.createDimension('lath', lat.shape[0])
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
            chg_msk.createVariable('time')
            # Define variables and fill data
            # Cell change mask
            chg_msk.createVariable('chng_mask', 'f8', ('lath','lonh'))
            chg_msk.variables['chng_mask'].units = 'none'
            chg_msk.variables['chng_mask'].long_name = 'Mask indicating changing ocean/ land cells'
            chg_msk.variables['chng_mask'][:] = chng_mask[:]
            # New ocean mask
            chg_msk.createVariable('o_mask_new', 'f8', ('lath','lonh'))
            chg_msk.variables['o_mask_new'].units = 'none'
            chg_msk.variables['o_mask_new'].long_name = 'Updated ocean mask'
            chg_msk.variables['o_mask_new'][:] = o_mask_new[:]
    
            chg_msk.description = "This file contains a mask which lists cells that should change from \
                        ocean to land (-1) or land to ocean (1) as well as an updated land-sea mask"
            chg_msk.history = "Created " + time.ctime(time.time())
            chg_msk.source = "created using /p/projects/climber3/huiskamp/POEM/work/slr_tool/5_prep_restarts/prep_restarts.py"
            chg_msk.close()
        else:
            chg_msk = CDF(str(path)+'history/change_mask.nc','r+');
# Write appropriate diagnostics to restart's metadata
            
            
            
            
            
            
            
            
            
            
            
            
            