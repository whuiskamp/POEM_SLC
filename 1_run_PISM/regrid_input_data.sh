#!/bin/bash

# Regrid VILMA output to the MOM ocean grid...
cdo remapbil,MOM.res.nc rsl.nc vilma2mom.nc

# Do the same for PISM output
cdo remapbil,MOM.res.nc pism_paleo06_6165_125ka.nc pism2mom.nc


