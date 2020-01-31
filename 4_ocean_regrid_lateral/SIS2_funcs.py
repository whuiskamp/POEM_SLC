#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script contains functions for calculating enthalpy and salinity
# in the sea ice model SIS2. These functions are taken directly from the model.
# There appear to be many redundant variables here, such as enth_unit,
# but these are retained to match the original functions and because they can
# be altered from default values when required (however unlikely).

import numpy as np

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2020"
__credits__ = ["Willem Huiskamp", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Production"

def enthalpy_liquid_freeze(S):
    # enthalpy_liquid_freeze returns the enthalpy of liquid water at the freezing
    # point for a given salinity. Taken from SIS2_ice_thm.F90
    # In: S          - The ice bulk salinity in g/kg
    #     enth_unit  - A conversion factor for enthalpy from Joules kg-1
    #     Cp_water   - Heat capacity of seawater in SIS2
    #     dTf_dS     - The derivative of the freezing temperature with salinity
    #     ENTH_LIQ_0 - The enthalpy of liquid fresh water at 0 C.
    enth_unit  = 1
    Cp_water   = 4200
    dTf_dS     = -0.054
    ENTH_LIQ_0 = 0
    
    enthalpy_liquid_freeze = enth_unit * \
    (Cp_water*(dTf_dS*S) + ENTH_LIQ_0)

    return enthalpy_liquid_freeze

def energy_melt_enthS(En,S):
    # energy_melt_enthS returns the energy needed to melt a given snow/ice
    # configuration, in J kg-1. Taken from SIS2_ice_thm.F90
    # In: S          - The ice bulk salinity in g/kg
    #     enth_unit  - A conversion factor for enthalpy from Joules kg-1
    #     En         - The ice enthalpy, in enthalpy units (often J/kg)
    enth_unit  = 1
    
    e_to_melt = enth_unit * (enthalpy_liquid_freeze(S) - En)
    
    return e_to_melt

def glob_sum_ice_enth(e_ice,e_sno,ice_frac,cell_area,h_sno,h_ice,s_ice,kg_H,nk_ice):
    # Calculates the global sum of enthalpy for ice and snow in sea ice model
    # This function is taken from the equivalent function ice_stock_pe from 
    # ice_type.F90
    # In:  e_ice      - Enthalpy of ice, in enthalpy units (often J/kg)
    #      e_sno      - Enthalpy of snow, in enthalpy units (often J/kg)
    #      ice_frac   - Fraction of grid cell covered in ice (0-1)
    #      cell_area  - Tracer grid cell area, from ocean grid
    #      h_sno      - Snow mass (kg/m2)
    #      h_ice      - Ice mass (kg/m2)
    #      s_ice      - Salinity of ice (g/kg). Default is 5.0
    # Out: total      - total energy in cell (J)
    
    kg_H_Nk = kg_H/nk_ice; total = 0
    
    for row in range(ice_frac.shape[1]):
        for col in range(ice_frac.shape[2]):
            if np.sum(ice_frac,0)[row,col] > 0:       # If no ice exists, skip it
                for cat in range(e_ice.shape[1]):
                    part_wt = cell_area[row,col]*ice_frac[cat,row,col]
                    total = total - (part_wt * (kg_H * h_sno[cat,row,col])) * \
                        energy_melt_enthS(e_sno[0,cat,row,col], 0)
                    for lvl in range(e_ice.shape[0]):
                        total = total - (part_wt * (kg_H_Nk * h_ice[cat,row,col])) * \
                        energy_melt_enthS(e_ice[lvl,cat,row,col], s_ice)
   
    return total

def sum_ice_enth(row,col,e_ice,e_sno,h_ice,h_sno,ice_frac,cell_area,s_ice,kg_H):
    # Calculates the sum of enthalpy for ice and snow in sea ice model in
    # a single grid cell.
    # This function is taken from the equivalent function ice_stock_pe from 
    # ice_type.F90
    # In:  e_ice                  - Enthalpy of ice, in enthalpy units (often J/kg)
    #      e_sno                  - Enthalpy of snow, in enthalpy units (often J/kg)
    #      ice_frac               - Fraction of grid cell covered in ice (0-1)
    #      cell_area              - Tracer grid cell area, from ocean grid
    #      h_sno                  - Snow mass (kg/m2)
    #      h_ice                  - Ice mass (kg/m2)
    #      S                      - Salinity of ice (g/kg). Default is 5.0
    # Out: total                  - total energy in cell (J)
    
    kg_H_Nk = kg_H/e_ice.shape[0]; total = 0
    
    for cat in range(e_ice.shape[1]):
        part_wt = cell_area[row,col]*ice_frac[cat,row,col]
        total = total - (part_wt * (kg_H * h_sno[cat,row,col])) * \
            energy_melt_enthS(e_sno[0,cat,row,col], 0)
        for lvl in range(e_ice.shape[0]):
            total = total - (part_wt * (kg_H_Nk * h_ice[cat,row,col])) * \
            energy_melt_enthS(e_ice[lvl,cat,row,col], s_ice)
   
    return total

def glob_sum_ice_sal(h_ice,ice_frac,S,kg_H,nk_ice,cell_area):
    # Calculates the sum of salinity for ice in sea ice model. Note there is 
    # no calculation for snow, as it has a salinity of 0.
    # This function is taken from the equivalent function ice_stock_pe from 
    # ice_type.F90
    # In: ice_frac   - Fraction of grid cell covered in ice (0-1)
    #     cell_area  - Tracer grid cell area, from ocean grid
    #     h_ice      - Ice mass (kg/m2)
    #     S          - Salinity of ice (g/kg). Default is 5.0
    #     kg_H       - kg to h grid conversion factor. Default is 1.
    #     nk_ice     - Number of vertical levels in sea ice field
    # Out: total     - Total salinity in sea ice
    kg_H_Nk = kg_H/nk_ice; total = 0
    sal_ice = np.full([4,5,80,120], S)
    for row in range(ice_frac.shape[1]):
        for col in range(ice_frac.shape[2]):
            for cat in range(ice_frac.shape[0]):
                for lvl in range(sal_ice.shape[0]):
                    total = total + (ice_frac[cat,row,col] * cell_area[row,col]) * \
                        (0.001*(kg_H_Nk*h_ice[cat,row,col])) * sal_ice[lvl,cat,row,col]
                        
    return total

def sum_ice_sal(row,col,ice_frac,cell_area,h_ice,S,kg_H,nk_ice):
    # Calculates the sum of salinity for ice in sea ice model at grid cell 
    # [row,col]. Note there is no calculation for snow, as it has a salinity of 0.
    # This function is taken from the equivalent function ice_stock_pe from 
    # ice_type.F90
    # In:  ice_frac   - Fraction of grid cell covered in ice (0-1)
    #      cell_area  - Tracer grid cell area, from ocean grid
    #      h_ice      - Ice mass (kg/m2)
    #      S          - Salinity of ice (g/kg). Default is 5.0
    #      kg_H       - kg to h grid conversion factor. Default is 1.
    #      nk_ice     - Number of vertical levels in sea ice field
    # Out: total      - Total salinity in sea ice
    kg_H_Nk = kg_H/nk_ice; total = 0
    sal_ice = np.full([4,5,80,120], S)
    
    for cat in range(ice_frac.shape[0]):
        for lvl in range(sal_ice.shape[0]):
            total = total + (ice_frac[cat,row,col] * cell_area[row,col]) * \
                (0.001*(kg_H_Nk*h_ice[cat,row,col])) * sal_ice[lvl,cat,row,col]
                
    return total








