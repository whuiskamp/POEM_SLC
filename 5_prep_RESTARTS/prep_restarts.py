# This script updates restart files with new variables that have been generated 
# in previous scripts. In addition, secondary fields that need to be altered
# based on redistribution of values/ changes in the land sea mask are also 
# altered here before being saved to appropriate restarts.

import time
from netCDF4 import Dataset as CDF

def write_rests(MOM,SIS,FLAGS):
# Load new restart files

# Change all ocean/ ice fields that need to be changed

# Save new fields to netCDF
# Write change mask to netCDF
    if FLAGS.test:
        id = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/change_mask.nc', 'w')
    else:
        id = CDF('dir', 'w')
    # Create dimensions for vars.
    id.createDimension('lonh', lon.shape[1]);
    id.createDimension('lath', lat.shape[0]);
    id.createVariable('lonh', 'f8', ('lonh'));
    id.variables['lonh'].units = 'degrees_east'
    id.variables['lonh'].cartesian_axis = 'X'
    id.variables['lonh'].long_name = 'T-cell longitude'
    id.variables['lonh'][:] = grid.variables['lonh'][:]
    id.createVariable('lath', 'f8', ('lath'));
    id.variables['lath'].units = 'degrees_north'
    id.variables['lath'].cartesian_axis = 'Y'
    id.variables['lath'].long_name = 'T-cell latitude'
    id.variables['lath'][:] = grid.variables['lath'][:]
    # Define variables and fill data
    # Cell change mask
    id.createVariable('chng_mask', 'f8', ('lath','lonh'))
    id.variables['chng_mask'].units = 'none'
    id.variables['chng_mask'].long_name = 'Mask indicating changing ocean/ land cells'
    id.variables['chng_mask'][:] = chng_mask[:]
    # New ocean mask
    id.createVariable('o_mask_new', 'f8', ('lath','lonh'))
    id.variables['o_mask_new'].units = 'none'
    id.variables['o_mask_new'].long_name = 'Updated ocean mask'
    id.variables['o_mask_new'][:] = o_mask_new[:]
    
    id.description = "This file contains a mask which lists cells that should change from \
                        ocean to land (-1) or land to ocean (1) as well as an updated land-sea mask"
    id.history = "Created " + time.ctime(time.time())
    id.source = "created using /p/projects/climber3/huiskamp/POEM/work/slr_tool/6_check_ocean_cells/chk_water_col.py"
    id.close()
# Write appropriate diagnostics to restart's metadata