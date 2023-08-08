# This script imports the restart files from VILMA (when available) and PISM and regrids them to the ocean grid using cdo.
# The VILMA restart is then imported and changes in relative sea level are applied to the ocean model's default bathymetry.
# This new altered topography file is saved as topog_new.nc and is required for the main coupling script.

# Note that VILMA writes consecutive time-slices to the same restart file, so when extracting the field for the current time-step, 
# we always access the most recent slice.

import time
import os.path
from netCDF4 import Dataset as CDF
import glob
import subprocess as sp

def regrid_rest(path):
	# This function reads in the experiment path and regrids the PISM and, when necessary, the VILMA restart file to the ocean grid
	# Regridding is done using cdo and employs second-order conservative method. This script requires a master topography file
	# to be stored in the 'path' directory (typically named SLC). It should include terrestrial topography as well as bathymetry
	# as it exists at the beginning of the simulation.
	# Next, if this is a VILMA coupling step, new topography field is calculated and saved to netCDF
	#  In: path 		- the path of the experiment directory's coupling data folder
	# Out: vilma2mom.nc - VILMA restart file on ocean grid
	#	   	pism2mom.nc - PISM restart file on ocean grid
	#      topog_new.nc - Default topography updated with VILMA data

	# First, check if this is a VILMA coupling step
    if os.path.isfile(str(path)+'rsl.nc'): # Change this path later for production runs (should not link to 'test_data')
        print('Relative sea level file found. Regridding to ocean grid ...')
        sp.run(['cdo', 'remapcon2,MOM.res.nc', 'rsl.nc', 'vilma2mom.nc'])
        VILMA_couple = True
	else:
		print('Relative sea level file not found. Proceeding to ice sheet input ...')
    # Read PISM restart file    
    try:
        sp.run(['cdo', 'remapcon2,'+str(path)+'MOM.res.nc', '+str(path)+'pism_res.nc, +str(path)+'vilma2mom.nc'])
    except:
        s = ("PISM extra file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_extra_file) )
	

    if VILMA_couple:
    	sp.run(['cp', str(path)+'../INPUT/topog.nc', str(path)+'topog_new.nc'])
    	master_topo = CDF(str(path)+'master_topo.nc','r',) # Import master topog. on MOM grid
    	VILMA_data  = CDF(str(path)+'vilma2mom.nc','r',) # Import VILMA restart on MOM grid
    	topo_m = master_topo.variables['depth'][:,:]; master_topo.close()  
    	slr    = VILMA_data.variables['slr'][-1,:,:]; VILMA_data.close() # Note that VILMA appends restarts, so we always select the most
    	# recent value
    	depth_new = topo_m - slr

    	topog_new = CDF(str(path)+'topog_new.nc','r+')   # Import what will be the new topogrpahy file
    	depth     = topog_new.variables['depth'][:,:]
    	depth[:,:]= depth_new
    	topog_new.close()

    else:
    	sp.run(['cp', str(path)+'INPUT/topog.nc', str(path)+'topog_new.nc']) # Doesn't need to be updated, just copy the file as-is.

    # Finally, we create a .nc file which tracks the history of the model topography & bathymetry at each coupling time-step
    # This file is updated at the end of the SLC operations (changes due to column thickness and ice sheets must be accounted for)
    if os.path.isfile(str(path)+'../history/topog_history.nc') == False: # At the start of a run, we need to creat this file as it doesn't exist yet
    	topog = CDF(str(path)+'INPUT/topog.nc','r');                  # The current topog file will represent the first time-step
    	depth = topo.variables['depth'][:,:]; topog.close()

    	hist  = CDF.Dataset('topog_history.nc', 'w')				  # Create the topog_history file and add metadata
		hist.createDimension('nx', depth.shape[1])
		hist.createDimension('ny', depth.shape[0])
		hist.createDimension('time', None)

		hist.createVariable('topog',f8,('time','ny','nx'))            
		hist.topog[0,:,:] = depth[:,:]
		hist.close()


    	



































