#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script contains functions for calculating enthalpy in the sea ice model 
# SIS2. These functions are taken directly from the model.
# There appear to be many redundant variables here, such as enth_unit,
# but these are retained to match the original functions and because they can
# be altered from default values when required (however unlikely).

import numpy as np

def enthalpy_liquid_freeze(S):
    # enthalpy_liquid_freeze returns the enthalpy of liquid water at the freezing
    # point for a given salinity. Taken from SIS2_ice_thm.F90
    # In: S          - The ice bulk salinity in g/kg
    #     enth_unit  - A conversion factor for enthalpy from Joules kg-1
    #     Cp_water   - Heat capacity of seawater
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

def sum_ice_enth(e_ice,e_sno):
    # Calculates the global sum of enthalpy for ice and snow in sea ice model
    # This function is taken from the equivalent function ice_stock_pe from 
    # ice_type.F90
    # In: e_ice      - Enthalpy of ice, in enthalpy units (often J/kg)
    #     e_sno      - Enthalpy of snow, in enthalpy units (often J/kg)
    #     ice_frac   - Fraction of grid cell covered in ice (0-1)
    #     cell_area  - Tracer grid cell area, from ocean grid
    #     h_oce      - Ocean vertical layer thickness (m)
    #     h_ice      - Ice mass (kg/m2)
    #     S          - Salinity of ice (g/kg). Default is 5.0
    kg_H = 1; kg_H_Nk = kg_H / e_ice.shape[0]; value = 0
    
    for row in range(ice_frac.shape[1]):
        for col in range(ice_frac.shape[2]):
            if np.sum(ice_frac,0)[row,col] > 0:       # If no ice exists, skip it
                for cat in range(e_ice.shape[1]):
                    part_wt = cell_area[row,col]*ice_frac[cat,row,col]
                    value = value - (part_wt * (kg_H * h_sno[cat,row,col])) * \
                        energy_melt_enthS(e_sno[0,cat,row,col], 0)
                    for lvl in range(e_ice.shape[0]):
                        value = value - (part_wt * (kg_H_Nk * h_ice[cat,row,col])) * \
                        energy_melt_enthS(e_ice[lvl,cat,row,col], s_ice)
    
    
    
    
    return total









