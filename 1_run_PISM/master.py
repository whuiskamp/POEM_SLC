# Overarching script for dynamic sea level coupling code
# This script reads data from restarts, creates python data structures and initialises
# the main scipts.

from netCDF4 import Dataset as CDF
import copy as cp
import argparse
from create_data_structs import MOM_vars, SIS_vars
import numpy as np
from chk_water_col import calc_coast
from regrid_lateral import get_param

if __name__ == "__main__":
# For now, we ignore argument parsing - this will be implemented once the test script works
    test = True # We'll use different datasets while running tests
    if test:
        new_bathy = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
        ctrl_bathy= CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/topog.nc','r')
        MOM6_rest = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM.res.nc','r')
        SIS2_rest = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ice_model.res.nc','r')
#        PISM_data = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
        grid      = CDF('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/ocean_geometry.nc','r')
        params_MOM    = open('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/MOM_parameter_doc.all','r').readlines()
        params_SIS    = open('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/SIS_parameter_doc.all','r').readlines()
    else:
        #old_bathy = CDF('')
        #new_bathy = CDF('','r')
        #MOM6_rest = CDF('','r')
        #SIS2_rest = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/history/MOM6_2019_04_25_11_05_29/RESTART/ice_model.res.nc','r')
        #PISM_data = CDF('path/to/PISM/DATA','r') # this should already have been regridded
        # seaice data probably not required in this step.
        #Omask     = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/INPUT/ocean_mask.nc','r')
    
    # Extract/ define variables
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
    
    
    depth      = new_bathy.variables['depth'][:,:]; new_bathy.close()     # Ocean depth of current simulation (m)
    ctrl_depth = ctrl_bathy.variables['depth'][:,:]; ctrl_bathy.close()   # The control, pre-industrial bathymetry (m)
    h          = MOM6_rest.variables['h'][0,:,:,:];                       # Ocean layer thickness (m)
    ave_eta    = MOM6_rest.variables['ave_ssh'][0,:,:];                   # Time-average sea surface height (m)
    eta        = MOM6_rest.variables['sfc'][0,:,:];                       # Sea surface height (m)
    lat        = grid.variables['geolat'][:,:];                           # Latitude of tracer cell-centres
    lon        = grid.variables['geolon'][:,:];                           # Longitude of tracer cell-centres
    grid_x     = lon.shape[1];                                            # Size of ocean longitude domain 
    grid_y     = lat.shape[0];                                            # Size of ocean latitude domain
    
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
    MOM = MOM_vars()
    # Initialise sea ice data structure
    SIS = SIS_vars()





















