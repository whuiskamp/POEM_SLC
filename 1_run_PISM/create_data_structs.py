# This script creates classes to handle the restart variables
# from the climate model output

class MOM_vars:
    # This class contains all variables from the MOM6 restart file
    def __init__ (self):
        self.o_mask         = o_mask
        self.o_mask_new     = o_mask_new
        self.depth          = depth
        self.ctrl_depth     = ctrl_depth
        self.h              = h
        self.ave_eta        = ave_eta
        self.eta            = eta
        self.coast          = coast
        self.chng_mask      = chng_mask
        self.lat            = lat
        self.lon            = lon
        self.grid_x         = grid_x
        self.grid_y         = grid_y
        self.chng_mask      = chng_mask
        
class SIS_vars:
    # This class contains all variables from the SIS2 restart file
    def __init__ (self):
        self.ice_frac       = ice_frac
        self.h_ice          = h_ice
        self.h_sno          = h_sno
        self.s_ice          = s_ice
        self.e_ice          = e_ice
        self.e_sno          = e_sno
        
        