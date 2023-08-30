#! /bin/sh
#
# This file is for running MOM6-SIS2 in a standalone configuration with a changing land-sea mask
# due to changes in ocean mass or externally forced changes in bathymetry depth.
# It assumes the INPUT directory contains a master topography file 'topog_ctrl.nc' with both topography and bathymetry in the same field.

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
diag_dir      = history
SLC_work      = mkdir $PWD/SLC
SLR_tool      = /p/projects/climber3/huiskamp/POEM/work/slr_tool
runoff_regrid = /p/projects/climber3/huiskamp/regrid_runoff # Maybe incorporate this inside the tool??

# ----------- Run parameters ----------
START          = 1  # For a new run, start at 1. For a restart, set this to previous end + 1 (duh..)

CPL_TIMESTEP   = 10 # Coupling time-step in years (must be >=1)
CPL_ITERATIONS = 50 # Number of coupling iterations to simulate (must be >=1)

if [ $START != 1 ]; then
    CPL_ITERATIONS += $START # We want numbering consistency for a run upon restart. 
fi

DO_PISM        = False  # Will this simulation include external forcing from PISM
DO_VILMA       = False  # Will this simulation include external forcing from VILMA 

OM4_run(){
    # Run MOM6-SIS2
    echo "Running MOM6-SIS2..."
    srun --propagate=ALL ./MOM6_SIS2 > MOM_log 2>&1
    mkdir -p "${diag_dir}/MOM6_run${RUN}"
    
    mv *.nc sbatch.62* MOM_parameter_doc* time_stamp.out CPU_stats ocean.stats logfile* MOM_log available_diags.* V_velocity_truncations U_velocity_truncations SIS_parameter_doc.layout SIS.available_diags SIS_fast.available_diags SIS_parameter_doc.debugging SIS_parameter_doc.short SIS_parameter_doc.all MOM_IC.nc seaice.stats "${diag_dir}/MOM6_run${RUN}"
    
    cp -a MOM_input input.nml MOM_override diag_table field_table SIS2_memory.h MOM_memory.h SIS_input SIS_override static_input.nml README "${diag_dir}/MOM6_run${RUN}"
}

slc_run(){
    # Run sea level change tool...
    echo "Running sea level change tool..."
    $SLR_tool/main/master.py
        --path      SLC_work
        --iteration RUN
        --PISM      DO_PISM
        --VILMA     DO_VILMA
        --time
        --verbose
}


# ----------- Run the model ----------


for RUN in `seq 1 $CPL_ITERATIONS`
do
    SIM_START_TIME=$(echo "$CPL_TIMESTEP * ($RUN -1)" | bc )
    SIM_END_TIME=$(echo "$CPL_TIMESTEP * $RUN" | bc )
    echo
    echo ' >> CPL_ITERATION=' $RUN
    time_tot = 0
    # Run poem for coupling timestep
    OM4_run
    if [$RUN == 1]; then
       cp ${diag_dir}/MOM6_run${RUN}/{MOM_parameter_doc.all,SIS_parameter_doc.all} $SLC_work
    # Run sea level change tool
    slc_run

    # If change mask of altered cells has been generated, save figure properly.
    if [-f chng_mask.pdf]; then
        mv chng_mask.pdf ${diag_dir}/MOM6_run${RUN}/chng_mask_${RUN}.pdf
    
    
    # 
done

















