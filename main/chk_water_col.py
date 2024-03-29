# This script 1) Updates topography depth with a uniform correction of +-1m depending on 
# changes in global mean sea level and 2) checks the integrated column thickness of ocean 
# cells and determines whether or not they should be changed to land.
# 
# This script requires the following inputs:
#  - MOM6 restart file
#  - Updated bathymetry
# 
# This script outputs the 'change mask', an array indicating which cells will change from
# ocean to land (-1) or land to ocean (1). All other values are NaN.
#
# There are 3 checks for changes in ocean/land. 
# Check 1: Has land ice exntent or change in bathymetry height created new land?
#          This check is unconditional (ie: no matter what sea level is doing, this
#          cell is now land).
# Check 2: Has column thickness decreased enough to dry out the cell?
#          This is conditional on changes in bathymetry and surrounding SSH.
# Check 3: Have cells become ocean due to receding land ice or SLR?          


## Import packages ##

import numpy as np
import sys
import copy as cp
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/1_run_PISM')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/3_check_channels')
from shared_funcs import halo_eta
from chk_cells import check_cells
import matplotlib.pyplot as plt

################################# Main Code ###################################
def check_water_col(MOM,ICE,OPTS):
    # Before checking individual cells, check if an an adjustment needs to be made for 
    # global mean sea level. If so, re-scale the topography and SSH fields.
    
    glob_ave_ssh = np.nanmean(MOM.ave_eta)
    if OPTS.verbose:
        print('Checking global mean seal level for changes...')
        print('Global mean sea level = '+str(glob_ave_ssh)+'m')
    if glob_ave_ssh >= 1:
        MOM.depth_new += 1
        MOM.eta       -= 1
        MOM.ave_eta   -= 1  
        OPTS.bathy_chg = True
        if OPTS.verbose:
            print('Sea level has increased; adjusting topography by +1m...')
    elif glob_ave_ssh <= -1:
        MOM.depth_new -= 1
        MOM.eta       += 1
        MOM.ave_eta   += 1
        OPTS.bathy_chg = True
        if OPTS.verbose:
            print('Sea level has decreased; adjusting topography by -1m...')
    elif not abs(glob_ave_ssh) >= 1 and OPTS.verbose:
        print('GMS has not changed by 1m or more, moving on...')

    # Check 1: Have we created new land via changes in ice sheet extent or
    # topography height? Update mask & change mask
    # Remember, in PISM 1 = Ice free bedrock, 2 = Grounded ice, 3 = Floating ice shelf
    # 4 = Open ocean.
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if (ICE.I_mask[i,j] == 2 and MOM.o_mask[i,j] > 0) or (MOM.depth_new[i,j] < OPTS.min_depth and MOM.o_mask[i,j] > 0):
                MOM.o_mask_new[i,j] = 0;
                MOM.depth_new[i,j]  = 0; 
                MOM.chng_mask[i,j]  = -1;
                if OPTS.verbose:
                    print('cell i='+str(i)+', j='+str(j)+' is becoming land \
                          due to changes in ice mask/ bathymetry')
    
    # Check 2: Has a change in SSH created new land?
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if  MOM.h_sum[i,j] < OPTS.min_thk and MOM.o_mask[i,j] > 0:
                MOM.o_mask_new[i,j] = 0;
                MOM.depth_new[i,j]  = 0;
                MOM.chng_mask[i,j]  = -1;
                if OPTS.verbose:
                    print('cell i='+str(i)+', j='+str(j)+' is becoming land \
                          due to a decrease in SSH')
    # Check 3: Have cells become ocean due to receding land ice or SLR?    
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if ICE.I_mask[i,j] == 4 and MOM.o_mask[i,j] == 0 and MOM.coast == 1 and MOM.depth_new[i,j] >= OPTS.new_depth:
                eta_mean = halo_eta(MOM.h_sum,i,j);
                if eta_mean >= OPTS.new_depth: 
                    MOM.chng_mask[i,j]  = 1;
                    MOM.o_mask_new[i,j] = 1;
                    if OPTS.verbose:
                        print('cell i='+str(i)+', j='+str(j)+' is becoming ocean')

    # Having completed these checks/ changes, set all land values to 0.
    MOM.depth_new[MOM.o_mask <= 0] = 0   

    # Finally, we need to make sure our new land-sea mask has not created isolated
    # cells or inland seas, if so, update o_mask and chng_mask where required
    
    if np.any(MOM.chng_mask): # We only need to perform this test if cells change
        check_cells(MOM,OPTS)               
        cells = np.count_nonzero(MOM.chng_mask)
        if OPTS.verbose:
            print('Number of cells changing this time-step is = '+ str(cells))
        OPTS.cont = True
        # Print figures showing new land mask and which cells are changing
        # Create dummy variable where changing cells are added to an o_mask
        dummy = cp.deepcopy(MOM.o_mask)
        dummy[MOM.chng_mask==1] = 3; dummy[MOM.chng_mask==-1] = -3;
        extent = (0, MOM.grid_x, MOM.grid_y, 0)
        _, ax = plt.subplots()
        ax.imshow(np.flipud(dummy[:,:]), extent=extent,cmap='Set3');
        # Set minor ticks
        ax.set_xticks(np.arange(0, 120, 1), minor=True);
        ax.set_yticks(np.arange(0, 80, 1), minor=True);
        # Add lines to make cells clearer
        ax.grid(which='minor', color='w', linewidth=0.5)
        #ax.set_frame_on(False)
        plt.savefig(OPTS.w_dir + 'chng_mask.pdf', dpi=150)
    else:
        # If no cells require changing, return a flag indicating that the rest
        # of the tool is not required for this coupling step
        OPTS.cont = False      
    