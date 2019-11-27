# This script creates classes to handle the restart variables
# from the climate model output

import os
import sys

class MOM_vars:
    # This class contains all variables from the MOM6 restart file
    def __init__ (self):
        self.o_mask         = o_mask
        self.o_mask_new     = o_mask_new
        self.depth          = depth
        self.h              = h
        self.ave_eta        = ave_eta
        self.eta            = eta
        self.coast          = coast
        self.chng_mask      = chng_mask
        self.lat            = lat
        self.lon            = lon
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
class SIS_vars:
    # This class contains all variables from the SIS2 restart file
    def __init__ (self):