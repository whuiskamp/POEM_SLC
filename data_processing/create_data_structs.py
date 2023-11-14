#!/usr/bin/env python3

# This script creates classes to handle the restart variables
# from the climate model output, in addition to requiring the topography file
# from MOM6 and the ocean and ice parameter files created at runtime. 
# The script assumes that these files are moved to a working directory at the end
# of a simulation. 
# Note that for Boussinesq numerics, H_to_m is used instead of H_to_kg_m2.
# In any event, this should not matter either way because they default to 1

from netCDF4 import Dataset as CDF
import numpy as np
import copy as cp
import time as t
import sys
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/common_funcs')
from shared_funcs import get_param
from shared_funcs import calc_coast

# Class containing options/ settings
class options:
    pass
# This class contains all variables from the MOM6 restart file
class MOM_vars:
    pass
# This class contains all variables from the SIS2 restart file    
class SIS_vars:
    pass
# This class contains all relevant variables from the PISM restart file
class PISM_vars:
    pass
# This class contains all relevant variables from the VILMA restart file
#class VILMA_vars:
#    pass
# This class contains all 'old' fields of variables that will require stock-
# checking for conservation.    
class old_vars:
    pass         

def init_data_structs(work_dir,ICE,EARTH,verbose):
    # Identidy restart, parameter, and gridspec files
    old_bathy  = CDF(work_dir + '../INPUT/topog.nc','r')
    new_bathy  = CDF(work_dir + '/topog.nc','r')
    ctrl_bathy = CDF(work_dir + '../INPUT/topog_ctrl','r')
    MOM6_rest  = CDF(work_dir + '/MOM.res.nc','r')
    SIS2_rest  = CDF(work_dir + '/ice_model.res.nc','r')
    vgrid      = CDF(work_dir + '../INPUT/vgrid.nc','r')
    params_MOM = open(work_dir + '/MOM_parameter_doc.all','r').readlines()
    params_SIS = open(work_dir + '/SIS_parameter_doc.all','r').readlines()
    Omask     = CDF(work_dir + '../INPUT/ocean_mask.nc','r')
    grid      = CDF(work_dir + '/ocean_static.nc','r')

    if ICE:
        PISM_data = CDF(work_dir + 'ice-data','r') # This should already have been regridded to the ocean grid
    if EARTH:
        VILMA_data = CDF(work_dir + 'earth_data','r') # This should already have been regridded to the ocean grid

    if verbose:
        print('Relevant restart and prarameter files found. Extracting data...')
        
    # Extract/ define variables
    # Ocean
    h_oce        = MOM6_rest.variables['h'][0,:,:,:].data                # Ocean layer thickness (m)
    o_temp       = MOM6_rest.variables['Temp'][0,:,:,:].data             # Ocean potential temperature (deg C)
    o_salt       = MOM6_rest.variables['Salt'][0,:,:,:].data             # Ocean salinity (ppt)
    ave_eta      = MOM6_rest.variables['ave_ssh'][0,:,:].data            # Time-average sea surface height (m)
    eta          = MOM6_rest.variables['sfc'][0,:,:].data                # Sea surface height (m)
    u            = MOM6_rest.variables['u']                              # Zonal velocity (m s-1)
    v            = MOM6_rest.variables['v']                              # Meridional velocity (m s-1)
    u2           = MOM6_rest.variables['u2']                             # Auxiliary zonal velocity (m s-1)
    v2           = MOM6_rest.variables['v2']                             # Auxiliary meridional velocity (m s-1)
    h2           = MOM6_rest.variables['h2']                             # Auxiliary ocean layer thickness (m)
    uh           = MOM6_rest.variables['uh']                             # Zonal thickness flux (m3 s-1)
    vh           = MOM6_rest.variables['vh']                             # Meridional thickness flux (m3 s-1)
    diffu        = MOM6_rest.variables['diffu']                          # Zonal horizontal viscous acceleration (m s-2)
    diffv        = MOM6_rest.variables['diffv']                          # Meridional horizontal viscous acceleration (m s-2)
    ubtav        = MOM6_rest.variables['ubtav']                          # Time mean baratropic zonal velocity (m s-1)
    vbtav        = MOM6_rest.variables['vbtav']                          # Time mean baratropic meridional velocity (m s-1)
    ubt_IC       = MOM6_rest.variables['ubt_IC']                         # Next init. cond. for baratropic zonal velocity (m s-1)
    vbt_IC       = MOM6_rest.variables['vbt_IC']                         # Next init. cond. for baratropic meridional velocity (m s-1)
    uhbt_IC      = MOM6_rest.variables['uhbt_IC']                        # Next init. cond. for baratropic zonal transport (m3 s-1)
    vhbt_IC      = MOM6_rest.variables['uhbt_IC']                        # Next init. cond. for baratropic meridional transport (m3 s-1)
    Kd_shear     = MOM6_rest.variables['Kd_shear']                       # Shear-driven turbulent diffusivity at interfaces (m2 s-1)
    Kv_shear     = MOM6_rest.variables['Kv_shear']                       # Shear-driven turbulent viscosity at interfaces (m2 s-1)
    TKE_turb     = MOM6_rest.variables['TKE_turb']                       # Turbulent kinetic energy per unit mass at interfaces (m2 s-1)
    Kv_shear_Bu  = MOM6_rest.variables['Kv_shear_Bu']                    # Shear-driven turbulent viscosity at vertex interfaces (m2 s-1)
    Kv_slow      = MOM6_rest.variables['Kv_slow']                        # Vertical turbulent viscosity at interfaces due to slow processes (m2 s-1)
    time_mom     = MOM6_rest.variables['Time']                           # Model time in days since simulation start (days)
    
    ######################### Optional ocean vars ########################
    # This section should include all optional tracers such as biogeochemical fields
    # Not all of these will require alteration, they are just included here for completeness
    try:
        age      = MOM6_rest.variables['age1']                           # Passive water mass age tracer (yrs)
        print('Age tracer in use!')
    except(KeyError):
        print('Age tracer not in use')

    ######################### End Optional ocean vars #####################

    if verbose:
        print('Ocean vars successfully extracted')
#    chng_mask    = chng_file.variables['chng_mask'][:,:];                # Mask of cells to change
    # Sea Ice
    # We do not import the incoming shortwave radiation as we have no reason to alter these fields
    h_ice        = SIS2_rest.variables['h_ice'][0,:,:,:].data;           # Ice thickness (m)
    h_sno        = SIS2_rest.variables['h_snow'][0,:,:,:].data;          # Snow thickness (m)
    ice_frac     = SIS2_rest.variables['part_size'][0,1:6,:,:].data;     # Ice fraction (0-1)
    s_ice        = SIS2_rest.variables['sal_ice'][0,0,0,0,0].data;       # Salinity of sea ice in g/kg - it's a fixed value
    e_ice        = SIS2_rest.variables['enth_ice'][0,:,:,:,:].data;      # Enthalpy of sea ice in J
    e_sno        = SIS2_rest.variables['enth_snow'][0,:,:,:,:].data;     # Enthalpy of snow in J
    flux_u       = SIS2_rest.variables['flux_u'][0,:,:].data;            # The flux of x-momentum into the ocean (Pa)
    flux_v       = SIS2_rest.variables['flux_v'][0,:,:].data;            # The flux of y-momentum into the ocean (Pa)
    flux_t       = SIS2_rest.variables['flux_t'][0,:,:].data;            # The flux of sensible heat out of the ocean (W m-2)
    flux_q       = SIS2_rest.variables['flux_q'][0,:,:].data;            # The evaporative moisture flux out of the ocean (kg m-2 s-1)
    flux_salt    = SIS2_rest.variables['flux_salt'][0,:,:].data;         # The flux of salt out of the ocean (kg m-2)
    flux_lw      = SIS2_rest.variables['flux_lw'][0,:,:].data;           # The longwave flux out of the ocean (W m-2)
    lprec        = SIS2_rest.variables['lprec'][0,:,:].data;             # The liquid precipitation flux into the ocean (kg m-2)
    fprec        = SIS2_rest.variables['fprec'][0,:,:].data;             # The frozen precipitation flux into the ocean (kg m-2)
    runoff       = SIS2_rest.variables['runoff'][0,:,:].data;            # Liquid water runoff into the ocean (kg m-2)
    calving      = SIS2_rest.variables['calving'][0,:,:].data;           # Calving of ice or runoff of frozen fresh water into the ocean (kg m-2)
    runoff_hflx  = SIS2_rest.variables['runoff_hflx'][0,:,:].data;       # Runoff heatflux. This variable should not be required normally. May be used when coupling with PICO-PISM
    calving_hflx = SIS2_rest.variables['calving_hflx'][0,:,:].data;      # The heat flux associated with calving, based on the temperature difference relative to a reference temperature (no units given in model code)
    p_surf       = SIS2_rest.variables['p_surf'][0,:,:].data;            # The pressure at the ocean surface (Pa) - may or may not include atmospheric pressure (from SIS2 code)
    t_surf_ice   = SIS2_rest.variables['t_surf_ice'][0,:,:].data;        # Surface temperature of sea ice (deg K)
    h_pond       = SIS2_rest.variables['h_pond'][0,:,:,:].data;          # Pond water thickness (m)
    u_ice        = SIS2_rest.variables['u_ice'][0,:,:].data;             # Zonal ice velocity (ms-1)
    v_ice        = SIS2_rest.variables['v_ice'][0,:,:].data;             # Meridional ice velocity (ms-1)
    sig11        = SIS2_rest.variables['sig11'][0,:,:].data;             # The xx component of the stress tensor in Pa m (or N m-1)
    sig12        = SIS2_rest.variables['sig12'][0,:,:].data;             # The xy and yx component of the stress tensor in Pa m (or N m-1)
    sig22        = SIS2_rest.variables['sig22'][0,:,:].data;             # The yy component of the stress tensor in Pa m (or N m-1)
    rough_mom    = SIS2_rest.variables['rough_mom'][0,:,:,:].data;       # The roughness for momentum at the ocean surface, as provided by ocean_rough_mod, apparently (?!) in m
    rough_heat = rough_mom; rough_moist = rough_mom                      # These are identical land-sea masks
    #coszen                                                               # Cosine of solar zenith angle
    #T_skin                                                               # The sea ice surface skin temperature (deg C)
    if verbose:
        print('Sea ice vars successfully extracted')

    # Parameters and new vars
    o_mask       = np.rint(Omask.variables['mask'][:,:].data);           # Old ocean mask (boolean)
    chng_mask    = np.full(o_mask.shape, 0);                             # Mask indicating where ocean/land cells are changing (0=no change, 1=land>ocean, -1=ocean>land)
    h_size_mask  = np.zeros(chng_mask.shape,dtype=float);                # Halo size mask
    o_mask_new   = np.rint(Omask.variables['mask'][:,:].data);           # Updated ocean mask (boolean)
    C_P          = get_param(params_MOM,'C_P');                          # The heat capacity of seawater in MOM6 (J kg-1 K-1)
    H_to_m       = get_param(params_MOM,'H_TO_M');
    H_to_kg_m2   = get_param(params_SIS,'H_TO_KG_M2');                   # Grid cell to mass conversion factor (1 by default)
    err_mass     = 0                                                     # Error in ocean mass after re-distribution between cells
    err_enth     = 0                                                     # Error in ocean temp after re-distribution between cells 
    err_salt     = 0                                                     # Error in ocean salt after re-distribution between cells
        
    # Geography
    depth_new  = new_bathy.variables['depth'][:,:];                      # Ocean depth of current simulation (m)
    depth_old  = old_bathy.variables['depth'][:,:];                      # Ocean depth of previous simulation (m)
    depth_ctrl = ctrl_bathy.variables['depth'][:,:];                     # The control, pre-industrial bathymetry (m)
    
    if verbose:
        print('Masks and topography vars successfully extracted')

    # Grid variables
    lat        = grid.variables['geolat'][:,:];                          # Latitude of tracer cell-centres
    lon        = grid.variables['geolon'][:,:];                          # Longitude of tracer cell-centres
    grid_x     = lon.shape[1];                                           # Size of ocean longitude domain 
    grid_y     = lat.shape[0];                                           # Size of ocean latitude domain
    grid_z     = int(get_param(params_MOM,'NK'));                        # Number of vertical levels in ocean
    cell_area  = grid.variables['Ah'][:,:];                              # Area of h (tracer) cells
    nk_ice     = get_param(params_SIS,'NK_ICE');                         # Number of z-levels in sea ice model (default=4)
    grid_dz    = vgrid.variables['dz'][:];                               # Default vertical grid spacing (m) (at the moment this is not used)
    
    if verbose:
        print('Grid vars successfully extracted')

    # Solid Earth and Ice sheet fields
    if EARTH:
        topo_chng = VILMA_data.variables['rsl'][-1,:,:]                  # Relative sea level in m referenced to the start of the simulation (we only want the most recent time-sclice)
        print('RSL field from VILMA successfully extracted')
        depth_new = depth_ctrl[:,:] + topo_chng[:,:]                     # Update topography with RSL fields from VILMA
        print('Topography updated with RSL field from VILMA')
        bathy_chg = True
    else:
        print('Model running without solid earth. Bypassing RSL field')
        bathy_chg = False
    if ICE:
        ice_mask  = PISM_data.variables['mask'][-1,:,:]                  # Land ice mask (1=ice-free bedrock, 2=grounded ice, 3=floating ice shelf, 4=open ocean) (we only want the most recent time-sclice)
        print('Grounded ice field from PISM successfully extracted')
    else:
        ice_mask = np.zeros((grid_y,grid_x))                             # With no ice data, we just init an array of 0's 
        print('Model running without land ice. Bypassing grounded ice field')

    # Variable pre-processing
    h_oce[:,o_mask==0]  = 0                                              # Change land to 0
    ave_eta[o_mask==0]  = np.nan                                         # Change land to NaN
    eta[o_mask==0]      = np.nan                                         # Change land to NaN
    h_sum               = np.sum(h_oce,0);                               # Thickness of water column (NOT depth of bathymetry)
    
    # Create copies of original fields for conservation checks 
    o_temp_old = cp.deepcopy(o_temp); e_ice_old    = cp.deepcopy(e_ice);
    e_sno_old  = cp.deepcopy(e_sno);  h_oce_old    = cp.deepcopy(h_oce);
    h_ice_old  = cp.deepcopy(h_ice);  o_salt_old   = cp.deepcopy(o_salt);
    h_sno_old  = cp.deepcopy(h_sno);  ice_frac_old = cp.deepcopy(ice_frac);

    if verbose:
        print('Preprocessing complete. Calculating coast...')
        
    # Identify coastal cells
    coast = calc_coast(o_mask)

    if verbose:
        print('Calculation of coast complete')
    
    # Close all input files (update with PISM/VILMA files later)
    new_bathy.close(); old_bathy.close(); ctrl_bathy.close(); MOM6_rest.close()
    SIS2_rest.close(); Omask.close(); grid.close(); vgrid.close()
    
###############################################################################    
    # Initialise options class
    OPTS = options()
    OPTS.verbose      = verbose
    OPTS.ICE          = ICE
    OPTS.EARTH        = EARTH
    OPTS.w_dir        = work_dir
    OPTS.bathy_chg    = bathy_chg

    # Initialise ocean data strucure
    MOM = MOM_vars()
    MOM.o_mask         = o_mask
    MOM.o_mask_new     = o_mask_new
    MOM.depth_new      = depth_new
    MOM.depth_old      = depth_old
    MOM.depth_ctrl     = depth_ctrl
    MOM.h_oce          = h_oce
    MOM.h_sum          = h_sum
    MOM.ave_eta        = ave_eta
    MOM.eta            = eta
    MOM.o_salt         = o_salt
    MOM.o_temp         = o_temp
    MOM.u              = u
    MOM.v              = v
    MOM.u2             = u2
    MOM.v2             = v2
    MOM.h2             = h2
    MOM.uh             = uh
    MOM.vh             = vh
    MOM.diffu          = diffu
    MOM.diffv          = diffv
    MOM.ubtav          = ubtav
    MOM.vbtav          = vbtav
    MOM.ubt_IC         = ubt_IC
    MOM.vbt_IC         = vbt_IC
    MOM.uhbt_IC        = uhbt_IC
    MOM.vhbt_IC        = vhbt_IC
    MOM.Kd_shear       = Kd_shear
    MOM.Kv_shear       = Kv_shear
    MOM.TKE_turb       = TKE_turb
    MOM.Kv_shear_Bu    = Kv_shear_Bu
    MOM.Kv_slow        = Kv_slow
    MOM.coast          = coast
    MOM.chng_mask      = chng_mask
    MOM.h_size_mask    = h_size_mask
    MOM.lat            = lat
    MOM.lon            = lon
    MOM.grid_x         = grid_x
    MOM.grid_y         = grid_y
    MOM.grid_z         = grid_z
    MOM.grid_dz        = grid_dz
    MOM.cell_area      = cell_area
    MOM.C_P            = C_P
    MOM.H_to_m         = H_to_m
    MOM.err_mass       = err_mass
    MOM.err_enth       = err_enth
    MOM.err_salt       = err_salt
    MOM.time           = time_mom
    # Optional fields
    try:
        MOM.age        = age
    except:
        pass

    # Initialise sea ice data structure
    SIS = SIS_vars()
    SIS.ice_frac       = ice_frac
    SIS.h_ice          = h_ice
    SIS.h_sno          = h_sno
    SIS.s_ice          = s_ice
    SIS.e_ice          = e_ice
    SIS.e_sno          = e_sno
    SIS.flux_u         = flux_u
    SIS.flux_v         = flux_v
    SIS.flux_t         = flux_t
    SIS.flux_q         = flux_q
    SIS.flux_salt      = flux_salt
    SIS.flux_lw        = flux_lw
    SIS.lprec          = lprec
    SIS.fprec          = fprec
    SIS.runoff         = runoff
    SIS.calving        = calving
    SIS.calving_hflx   = calving_hflx
    SIS.runoff_hflx    = runoff_hflx
    SIS.p_surf         = p_surf
    SIS.t_surf_ice     = t_surf_ice
    SIS.h_pond         = h_pond
    SIS.u_ice          = u_ice
    SIS.v_ice          = v_ice
    SIS.sig11          = sig11
    SIS.sig12          = sig12
    SIS.sig22          = sig22
    SIS.rough_mom      = rough_mom
    SIS.rough_moist    = rough_moist
    SIS.rough_heat     = rough_heat
    SIS.nk_ice         = nk_ice
    SIS.H_to_kg_m2     = H_to_kg_m2
    
    # Initialise ice sheet data structure
    ICE = PISM_vars()
    ICE.I_mask         = ice_mask
    
    # Initialise stock checking data structure
    OLD = old_vars()
    OLD.o_temp         = o_temp_old
    OLD.h_oce          = h_oce_old
    OLD.o_salt         = o_salt_old
    OLD.e_sno          = e_sno_old
    OLD.e_ice          = e_ice_old
    OLD.h_ice          = h_ice_old
    OLD.h_sno          = h_sno_old
    OLD.ice_frac       = ice_frac_old
      
    return MOM, SIS, OLD, ICE, OPTS 
    