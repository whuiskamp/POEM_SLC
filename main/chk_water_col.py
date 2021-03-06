# This script checks the integrated column thickness of ocean cells and determines whether or not they should be changed to land.
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
import time
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/1_run_PISM')
sys.path.append('/p/projects/climber3/huiskamp/POEM/work/slr_tool/3_check_channels')
from shared_funcs import halo_eta
from chk_cells import check_cells
import matplotlib.pyplot as plt

__author__ = "Willem Huiskamp"
__copyright__ = "Copyright 2020"
__credits__ = ["Willem Huiskamp", ""]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Willem Huiskamp"
__email__ = "huiskamp@pik-potsdam.de"
__status__ = "Alpha"

################################# Main Code ###################################
def check_water_col(MOM,ICE,FLAGS):
    # Check 1 Have we created new land via changes in ice sheet extent or
    # topography height? Update mask & change mask
    t_start = time.time()
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if (ICE.I_mask[i,j] >= 0.7 and MOM.o_mask[i,j] > 0) or 0 < MOM.depth_new[i,j] < 5: 
                MOM.o_mask_new[i,j] = 0;
                MOM.depth_new[i,j]  = 0; 
                MOM.chng_mask[i,j]  = -1;
    
    # Check 2: Has change in SSH/ bathymetry created new land?
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if (MOM.ave_eta[i,j] + MOM.depth_new[i,j] < 2) and MOM.o_mask[i,j] > 0:
                MOM.o_mask_new[i,j] = 0;
                MOM.depth_new[i,j]  = 0;
                MOM.chng_mask[i,j]  = -1;
        
    # Check 3: Have cells become ocean due to receding land ice or SLR?    
    for i in range(MOM.grid_y):
        for j in range(MOM.grid_x):
            if ICE.I_mask[i,j] <= 0.3 and MOM.o_mask[i,j] == 0 and MOM.depth_new[i,j] >= 5:
                eta_mean = halo_eta(MOM.eta,i,j);
                if eta_mean + MOM.depth_new[i,j] > 2: # optionally add: and coast[i,j] == 1: 
                    MOM.chng_mask[i,j]  = 1;
                    MOM.o_mask_new[i,j] = 1;
    # Finally, we need to make sure our new land-sea mask has not created isolated
    # cells or inland seas, if so, updated o_mask and chng_mask where required
    
    if np.any(MOM.chng_mask): # We only need to perform this test if cells change
        check_cells(MOM,FLAGS)               
        cells = np.count_nonzero(MOM.chng_mask)
        if FLAGS.verbose:
            print('Number of cells changing this time-step is = '+ str(cells))
        FLAGS.cont = True
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
        plt.savefig('/p/projects/climber3/huiskamp/POEM/work/slr_tool/test_data/chng_mask.pdf', dpi=150)
        t_end = time.time()
        FLAGS.t_chk_cells = t_end - t_start
    else:
        # If no cells require changing, return a flag indicating that the rest
        # of the tool is not required for this coupling step
        FLAGS.cont = False
        t_end = time.time()
        FLAGS.t_chk_cells = t_end - t_start
    