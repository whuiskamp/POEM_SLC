#!/bin/bash
# Use second order conservative regridding to remap fields to MOM grid
# Regrid VILMA output to the MOM ocean grid...
cdo remapcon2,MOM.res.nc rsl.nc vilma2mom.nc

# Do the same for PISM output
cdo remapcon2,MOM.res.nc pism_paleo06_6165_125ka.nc pism2mom.nc


