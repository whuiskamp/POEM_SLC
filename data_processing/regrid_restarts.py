# This script imports the restart files from VILMA (when available) and PISM and regrids them to the ocean grid using cdo.
# The VILMA restart is then imported and changes in relative sea level are applied to the ocean model's default bathymetry.
# This new altered topography file is saved as topog_new.nc and is required for the main coupling script.

# Note that VILMA writes consecutive time-slices to the same restart file, so when extracting the field for the current time-step, 
# we always access the most recent slice.

import os.path
import subprocess as sp

def regrid_rest(path):
	# This reads in the experiment path and regrids the PISM (pism_res.nc) and, when necessary, the VILMA restart (rsl.nc) 
	# file to the ocean grid. Regridding is done using cdo and employs second-order conservative method. 
	# This script requires a master topography file to be stored in the 'path' directory (typically named SLC). 
	# 
	#  In: path 		- the path of the experiment directory's coupling data folder
	# Out: vilma2mom.nc - VILMA restart file on ocean grid
	#	   pism2mom.nc  - PISM restart file on ocean grid
	#      topog_new.nc - Default topography updated with VILMA data

	# First, check if this is a VILMA coupling step
	if os.path.isfile(str(path)+'rsl.nc'):
		print('Relative sea level file found. Regridding to ocean grid ...')
		sp.run(['cdo', 'remapcon2,'+str(path)+'MOM.res.nc', str(path)+'rsl.nc', str(path)+'vilma2mom.nc'])
	else:
		print('Relative sea level file not found. Proceeding to ice sheet input ...')
    # Read PISM restart file    
	if os.path.isfile(str(path)+'pism_res.nc'):
		print('PISM restart file found. Regridding to ocean grid ...')
		sp.run(['cdo', 'remapcon2,'+str(path)+'MOM.res.nc', str(path)+'pism_res.nc', str(path)+'pism2mom.nc'])
	else:
		print('PISM restart file not found.')
	


    	



































