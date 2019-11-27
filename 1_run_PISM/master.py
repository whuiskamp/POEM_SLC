# Overarching script for dynamic sea level coupling code
# This script reads data from restarts, creates python data structures and initialises
# the main scipts.

from netCDF4 import Dataset as CDF
import copy as cp
import argparse
from create_data_structs import MOM_vars

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    test = True # We'll use different datasets while running tests
    if test:
        new_bathy = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
        ctrl_bathy= CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
        MOM6_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/MOM.res.nc','r')
        PISM_data = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        grid      = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
    else:
        #old_bathy = CDF('')
        new_bathy = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
        MOM6_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/MOM.res.nc','r')
        #SIS2_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/ice_model.res.nc','r')
        #PISM_data = CDF('path/to/PISM/DATA','r') # this should already have been regridded
        # seaice data probably not required in this step.
        Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
    
    # Extract/ define variables
    depth      = new_bathy.variables['depth'][:,:]; new_bathy.close()
    ctrl_depth = ctrl_bathy.variables['depth'][:,:]; ctrl_bathy.close()
    h          = MOM6_rest.variables['h'][0,:,:,:]; 
    ave_eta    = MOM6_rest.variables['ave_ssh'][0,:,:];
    eta        = MOM6_rest.variables['sfc'][0,:,:];
    lat        = grid.variables['geolat'][:,:];
    lon        = grid.variables['geolon'][:,:];
    grid_x     = lon.shape[1];
    grid_y     = lat.shape[0];
    
    if test:
        ice_frac  = PISM_data.variables['mask'][:,:];
        ice_frac[ice_frac==0] = 2; ice_frac[ice_frac<2] = 0;
    else:
        ice_frac  = PISM_data.variables['ice_frac'][:,:];
    o_mask     = Omask.variables['mask'][:,:];
    o_mask_new = cp.deepcopy(o_mask);
    chng_mask  = np.full(o_mask.shape, np.nan);
    
    # Variable pre-processing
    h[:,o_mask==0]      = np.nan;  # Change land to NaN
    ave_eta[o_mask==0]  = np.nan;  # Change land to NaN
    eta[o_mask==0]      = np.nan   # Change land to NaN
    
    # Identify coastal cells
    coast = calc_coast(o_mask)
###############################################################################    
    # Initialise ocean data strucure
    ocn_vars = MOM_vars()
    # Initialise sea ice data structure
    ice_vars = SIS_vars()





















