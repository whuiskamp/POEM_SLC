#! /bin/sh
#
# This file is for running MOM6-SIS2 in a standalone configuration with a changing land-sea mask
# due to changes in ocean mass or externally forced changes in bathymetry depth.

#SBATCH --job-name=MOM6_SLC

#SBATCH --tasks=32

#SBATCH --tasks-per-node=16
##SBATCH --partition=haswell

# Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
#                   "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#SBATCH --qos=medium
#SBATCH --time=7-00:00:00
#SBATCH --constraint=tasksmax

#SBATCH --output=sbatch.%j.out
##SBATCH --error=sbatch.%j.err

#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT,TIME_LIMIT_90,TIME_LIMIT

# ----------- Define paths ----------

SLC_work = $PWD/SLC
SLR_tool = /p/projects/climber3/huiskamp/POEM/work/slr_tool
runoff_regrid = /p/projects/climber3/huiskamp/regrid_runoff

# ----------- Run parameters ----------

CPL_TIMESTEP   = 10 # Coupling time-step in years (must be >=1)
CPL_ITERATIONS = 50 # Number of coupling iterations to simulate (must be >=1)

PISM_forcing   = False  # Will this simulation include external forcing from PISM
VILMA_forcing  = False  # Will this simulation include external forcing from VILMA 




# ----------- Run the model ----------


for RUN in `seq 1 $CPL_ITERATIONS`
do
    SIM_START_TIME=$(echo "$CPL_TIMESTEP * ($RUN -1)" | bc )
    SIM_END_TIME=$(echo "$CPL_TIMESTEP * $RUN" | bc )
    echo
    echo ' >> CPL_ITERATION=' $RUN

    # Run poem for coupling timestep
    poem_run
    
    # Run sea level change tool
    slc_run
    

done

















