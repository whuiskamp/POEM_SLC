#!/bin/bash

# delete PISM history
ncatted -a history_of_appended_files,global,d,c, pism_paleo06_6165_125ka.nc

# remap from PISM (stere) grid to VILMA (gaussian) grid 
cdo remapbil,rsl.nc pism_paleo06_6165_125ka.nc pism2vilma_cdo.nc

# remove missing and fill values, such that python can replace them
ncatted -a _FillValue,thk,d,, pism2vilma_cdo.nc
ncatted -a missing_value,thk,d,, pism2vilma_cdo.nc

# make a copy and remoce first time slice to start from -123 kyr
#cp pism2vilma_cdo.nc pism2vilma.nc
cdo seltimestep,2/125 pism2vilma_cdo.nc pism2vilma.nc

# replace northern hemisohere values with ice5g
python replace_ant_from_pism.py 

# volker prefers this definition to avoid cdo issues
ncrename -v time,epoch -d time,epoch -v thk,Ice pism2vilma.nc
ncatted -O -a units,epoch,o,c,'ka BP' -a long_name,epoch,o,c,'Epoch' -a calendar,epoch,o,c,'proleptic_gregorian' pism2vilma.nc

