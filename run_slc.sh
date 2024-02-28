#!/bin/sh
#
# This file is for running MOM6-SIS2 in a standalone configuration with a changing land-sea mask
# due to changes in ocean mass or externally forced changes in bathymetry depth.
# It assumes the INPUT directory contains a master topography file 'topog_ctrl.nc' with both topography and bathymetry in the same field.

#SBATCH --job-name=MOM6_SLC

#SBATCH --tasks=32

##SBATCH --tasks-per-node=16
##SBATCH --partition=haswell

# Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
#                   "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
##SBATCH --qos=priority
##SBATCH --time=2:00:00
#SBATCH --qos=medium
#SBATCH --time=7-00:00:00
##SBATCH --constraint=tasksmax

#SBATCH --output=sbatch.%j.out
##SBATCH --error=sbatch.%j.err

#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT,TIME_LIMIT_90,TIME_LIMIT

# ----------- Load python env. ----------
source activate /home/huiskamp/.conda/envs/mom6

# ----------- Define paths ----------
diag_dir=history
SLC_work=$PWD/SLC
SLR_tool=/p/projects/climber3/huiskamp/POEM/work/slr_tool
runoff_regrid=/p/projects/climber3/huiskamp/regrid_runoff # Maybe incorporate this inside the tool??

# ----------- Run parameters ----------
START=1  # For a new run, start at 1. For a restart, set this to previous end + 1 (duh..)

CPL_TIMESTEP=10 # Coupling time-step in years (must be >=1)
CPL_ITERATIONS=50 # Number of coupling iterations to simulate (must be >=1)

if [ $START != 1 ]; then
    CPL_ITERATIONS+=$START # We want numbering consistency for a run upon restart. 
fi

DO_PISM=false  # Will this simulation include external forcing from PISM
DO_VILMA=false  # Will this simulation include external forcing from VILMA

mkdir -p $SLC_work

OM4_run(){
    # Run MOM6-SIS2
    echo "Running MOM6-SIS2..."
    srun --propagate=ALL ./MOM6 > MOM_log 2>&1
    result=$?
    if [ $result != 0 ]; then
        exit $result
    fi

    mkdir -p "${diag_dir}/MOM6_run_${RUN_long}"
    
    # Copy param. files to working dir...
    cp -a MOM_parameter_doc.all SIS_parameter_doc.all $SLC_work
    # Copy mask file and old topography to history dir.
    cp -a INPUT/{topog.nc,ocean_mask.nc} "${diag_dir}/MOM6_run_${RUN_long}"

    # Move diagnostics and log files to history dir.
    mv *.nc sbatch.62* MOM_parameter_doc* time_stamp.out CPU_stats ocean.stats logfile* MOM_log available_diags.* \
    V_velocity_truncations U_velocity_truncations SIS_parameter_doc.layout SIS.available_diags SIS_fast.available_diags \
    SIS_parameter_doc.debugging SIS_parameter_doc.short SIS_parameter_doc.all MOM_IC.nc seaice.stats "${diag_dir}/MOM6_run_${RUN_long}"
    
    cp -a MOM_input input.nml MOM_override diag_table field_table SIS2_memory.h MOM_memory.h SIS_input SIS_override \
    static_input.nml README "${diag_dir}/MOM6_run_${RUN_long}"

    # Move restarts to working dir (also create a backup)...
    cp -ar RESTART/  "${diag_dir}/MOM6_run_${RUN_long}"
    cp -a RESTART/{MOM.res.nc,ice_model.res.nc} $SLC_work
    cp -a INPUT/{topog.nc,ocean_mask.nc} $SLC_work
}

slc_run(){
    # Run sea level change tool...
    if [[ $DO_PISM = false && $DO_VILMA = false ]]; then
        echo "Running sea level change tool without solid earth and ice sheets..."
        $SLR_tool/main/master.py  \
            --path      $SLC_work  \
            --iteration $RUN       \
            --verbose
        result=$?
    elif [[ $DO_PISM = true && $DO_VILMA = false ]]; then
        echo "Running sea level change tool with ice sheets..."
        $SLR_tool/main/master.py  \
            --path      $SLC_work  \
            --iteration $RUN       \
            --PISM                 \
            --verbose
        result=$?
    elif [[ $DO_PISM = false && $DO_VILMA = true ]]; then
        echo "Running sea level change tool with solid earth..."
        $SLR_tool/main/master.py  \
            --path      $SLC_work  \
            --iteration $RUN       \
            --VILMA                 \
            --verbose
        result=$?   
    fi

    if [ $result != 0 ]; then
        exit $result
    fi
}   


# ----------- Run the model ----------


for RUN in `seq $START $CPL_ITERATIONS`
do
    SIM_START_TIME=$(echo "$CPL_TIMESTEP * ($RUN -1)" | bc )
    SIM_END_TIME=$(echo "$CPL_TIMESTEP * $RUN" | bc )
    echo
    echo ' >> CPL_ITERATION=' $RUN
    RUN_long=$(printf "%04d" "$RUN") # Leading 0's make organising directories easier
    # Run poem for coupling timestep
    OM4_run
    if [ $RUN == 1 ]; then
       cp ${diag_dir}/MOM6_run_${RUN_long}/{MOM_parameter_doc.all,SIS_parameter_doc.all} $SLC_work
       cp ${diag_dir}/MOM6_run_${RUN_long}/*.ocean_static.nc $SLC_work/ocean_static.nc
    fi
    # Run sea level change tool
    slc_run

    # If change mask of altered cells has been generated, save figure properly.
    if [ -f chng_mask.pdf ]; then
        mv chng_mask.pdf ${diag_dir}/MOM6_run_${RUN_long}/chng_mask_${RUN_long}.pdf
    fi
    # In any case, copy restarts and topography back to INPUT for restart.
    cp ${SLC_work}/{MOM.res.nc,ice_model.res.nc,coupler.res,topog.nc} INPUT
done

















