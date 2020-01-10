# This script creates classes to handle the restart variables
# from the climate model output

from netCDF4 import Dataset as CDF
import numpy as np
import copy as cp

class MOM_vars:
    # This class contains all variables from the MOM6 restart file
    def __init__ (self):
        self.o_mask         = o_mask
        self.o_mask_new     = o_mask_new
        self.depth          = depth
        self.ctrl_depth     = ctrl_depth
        self.h_oce          = h_oce
        self.h_sum          = h_sum
        self.ave_eta        = ave_eta
        self.eta            = eta
        self.u              = u
        self.v              = v
        self.u2             = u2
        self.v2             = v2
        self.h2             = h2
        self.uh             = uh
        self.vh             = vh
        self.diffu          = diffu
        self.diffv          = diffv
        self.ubtav          = ubtav
        self.vbtav          = vbtav
        self.ubt_IC         = ubt_IC
        self.vbt_IC         = vbt_IC
        self.uhbt_IC        = uhbt_IC
        self.vhbt_IC        = vhbt_IC
        try:
            self.age        = age
        except:
            print('Age tracer not used')
        self.Kd_shear       = Kd_shear
        self.Kv_shear       = Kv_shear
        self.TKE_turb       = TKE_turb
        self.Kv_shear_Bu    = Kv_shear_Bu
        self.Kv_slow        = Kv_slow
        self.coast          = coast
        self.chng_mask      = chng_mask
        self.lat            = lat
        self.lon            = lon
        self.grid_x         = grid_x
        self.grid_y         = grid_y
        self.C_P            = C_P
        self.H_to_kg_m2     = H_to_kg_m2
        
        
class SIS_vars:
    # This class contains all variables from the SIS2 restart file
    def __init__ (self):
        self.ice_frac       = ice_frac
        self.h_ice          = h_ice
        self.h_sno          = h_sno
        self.s_ice          = s_ice
        self.e_ice          = e_ice
        self.e_sno          = e_sno
        self.flux_u         = flux_u
        self.flux_v         = flux_v
        self.flux_t         = flux_t
        self.flux_q         = flux_q
        self.flux_salt      = flux_salt
        self.flux_lw        = flux_lw
        self.lprec          = lprec
        self.fprec          = fprec
        self.runoff         = runoff
        self.calving        = calving
        self.runoff_hflx    = runoff_hflx
        self.p_surf         = p_surf
        self.t_surf_ice     = t_surf_ice
        self.h_pond         = h_pond
        self.u_ice          = u_ice
        self.v_ice          = v_ice
        self.sig11          = sig11
        self.sig12          = sig12
        self.sig22          = sig22
        self.rough_mom      = rough_mom
        self.rough_moist    = rough_moist
        self.rough_heat     = rough_heat
        
def init_data_structs(work_dir,test):
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
#        old_bathy = CDF(work_dir+'INPUT/topog.nc,'r')
        new_bathy = CDF(work_dir + '/INPUT/topog.nc','r')
        MOM6_rest = CDF(work_dir + '/RESTART/MOM.res.nc,'r')
        SIS2_rest = CDF(work_dir + '/RESTART/ice_model.res.nc','r')
        #PISM_data = CDF('path/to/PISM/DATA','r') # this should already have been regridded
        #Omask     = CDF(work_dir + '/INPUT/ocean_mask.nc','r')
    
    # Extract/ define variables
    # Ocean
    h_oce        = MOM6_rest.variables['h'][0,:,:,:].data                # Ocean layer thickness (m)
    o_temp       = MOM6_rest.variables['Temp'][0,:,:,:].data             # Ocean potential temperature (deg C)
    o_salt       = MOM6_rest.variables['Salt'][0,:,:,:].data #.astype(int) # Ocean salinity (ppt)
    ave_eta      = MOM6_rest.variables['ave_ssh'][0,:,:]                 # Time-average sea surface height (m)
    eta          = MOM6_rest.variables['sfc'][0,:,:]                     # Sea surface height (m)
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
    
    # Optional ocean vars
    try:
        age      = MOM6_rest.variables['age1']                           # Passive water mass age tracer (yrs)
    except(KeyError):
        print('Age tracer not in use')
#    chng_mask    = chng_file.variables['chng_mask'][:,:];                # Mask of cells to change
    # Sea Ice
    # We do not import the incoming shortwave radiation as we have no reason to alter these fields
    h_ice        = SIS2_rest.variables['h_ice'][0,:,:,:].data;           # Ice thickness (m)
    h_sno        = SIS2_rest.variables['h_snow'][0,:,:,:].data;          # Snow thickness (m)
    ice_frac     = SIS2_rest.variables['part_size'][0,1:6,:,:].data;     # Ice fraction (0-1)
    s_ice        = SIS2_rest.variables['sal_ice'][0,0,0,0,0].data;       # Salinity of sea ice in g/kg - it's a fixed value
    e_ice        = SIS2_rest.variables['enth_ice'][0,:,:,:,:].data;      # Enthalpy of sea ice in J
    e_sno        = SIS2_rest.variables['enth_snow'][0,:,:,:,:].data;     # Enthalpy of snow in J
    flux_u       = SIS2_rest.variables['flux_u'][0,:,:].data;            # 
    flux_v       = SIS2_rest.variables['flux_v'][0,:,:].data;            # 
    flux_t       = SIS2_rest.variables['flux_t'][0,:,:].data;            #
    flux_q       = SIS2_rest.variables['flux_q'][0,:,:].data;            #
    flux_salt    = SIS2_rest.variables['flux_salt'][0,:,:].data;         #
    flux_lw      = SIS2_rest.variables['flux_lq'][0,:,:].data;           #
    lprec        = SIS2_rest.variables['lprec'][0,:,:].data;             #
    fprec        = SIS2_rest.variables['fprec'][0,:,:].data;             #
    runoff       = SIS2_rest.variables['runoff'][0,:,:].data;            # Liquid water runoff 
    calving      = SIS2_rest.variables['calving'][0,:,:].data;           # 
    runoff_hflx  = SIS2_rest.variables['runoff_hflx'][0,:,:].data;       # Runoff heatflux. This variable should not be required normally. May be used when coupling with PICO-PISM
    calving_hflx = SIS2_rest.variables['calving_hflx'][0,:,:].data;      #
    p_surf       = SIS2_rest.variables['p_surf'][0,:,:].data;            #
    t_surf_ice   = SIS2_rest.variables['t_surf_ice'][0,:,:].data;        # Surface temperature of sea ice (deg K)
    h_pond       = SIS2_rest.variables['h_pond'][0,:,:,:].data;          # Pond water thickness (m)
    u_ice        = SIS2_rest.variables['u_ice'][0,:,:].data;             # 
    v_ice        = SIS2_rest.variables['v_ice'][0,:,:].data;             #
    sig11        = SIS2_rest.variables['sig11'][0,:,:].data;             #
    sig12        = SIS2_rest.variables['sig12'][0,:,:].data;             #
    sig22        = SIS2_rest.variables['sig22'][0,:,:].data;             #
    rough_mom    = SIS2_rest.variables['rough_mom'][0,:,:,:].data;       #
    rough_heat = rough_mom; rough_moist = rough_mom                      # These are identical land-sea masks
    coszen # ??
    T_skin # ??
    
    # Parameters and new vars
    o_mask       = Omask.variables['mask'][:,:]
    chng_mask    = np.full(o_mask.shape, np.nan);                        # Mask indicating where ocean/land cells are changing
    h_size_mask  = np.zeros(chng_mask.shape,dtype=float);                # Halo size mask
    o_mask_new   = Omask.variables['mask'][:,:];                         # Updated ocean mask
    C_P          = get_param(params_MOM,'C_P');                          # The heat capacity of seawater in MOM6 (J kg-1 K-1)
    H_to_kg_m2   = get_param(params_SIS,'H_TO_KG_M2');                   # grid cell to mass conversion factor (1 by default)
    
    # Geography
    depth      = new_bathy.variables['depth'][:,:]; new_bathy.close()    # Ocean depth of current simulation (m)
    ctrl_depth = ctrl_bathy.variables['depth'][:,:]; ctrl_bathy.close()  # The control, pre-industrial bathymetry (m)
    
    # Grid variables
    lat        = grid.variables['geolat'][:,:];                          # Latitude of tracer cell-centres
    lon        = grid.variables['geolon'][:,:];                          # Longitude of tracer cell-centres
    grid_x     = lon.shape[1];                                           # Size of ocean longitude domain 
    grid_y     = lat.shape[0];                                           # Size of ocean latitude domain
    grid_z     = get_param(params_MOM,'NK');                             # Number of vertical levels in ocean
    cell_area  = grid.variables['Ah'][:,:];                              # Area of h (tracer) cells
    
    if test:
        ice_frac  = PISM_data.variables['mask'][:,:];
        ice_frac[ice_frac==0] = 2; ice_frac[ice_frac<2] = 0 
    else:
        ice_frac  = PISM_data.variables['ice_frac'][:,:] 
     
      
    # Variable pre-processing
    h_oce[:,o_mask==0]  = np.nan                                         # Change land to NaN
    ave_eta[o_mask==0]  = np.nan                                         # Change land to NaN
    eta[o_mask==0]      = np.nan                                         # Change land to NaN
    h_sum               = np.sum(h_oce,0);
    
    # Identify coastal cells
    coast = calc_coast(o_mask)
    
###############################################################################    
    # Initialise ocean data strucure
    MOM = MOM_vars()
    # Initialise sea ice data structure
    SIS = SIS_vars()
    
    return MOM, SIS
    