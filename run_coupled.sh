#! /bin/sh

#  Copyright (C) 2019, 2020 PISM-MOM_coupling authors, see AUTHORS file
#
#  This file is part of PISM-MOM_coupling
#
#  PISM-MOM_coupling is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PISM-MOM_coupling is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PISM-MOM_coupling.  If not, see <https://www.gnu.org/licenses/>.
#
# This is a test file, created 20/05/21 for integrating MOM6 and the SLC tool with PISM-PICO and VILMA


#SBATCH --job-name={{settings.experiment}}

#SBATCH --tasks={{settings.slurm_tasks}}
{% if settings.slurm_exclusive==True -%}
#SBATCH --exclusive
{% endif %}

{% if settings.slurm_partition=="broadwell" -%}
#SBATCH --tasks-per-node=32
#SBATCH --partition=broadwell
{% else %}
#SBATCH --tasks-per-node=16
##SBATCH --partition=haswell
{% endif %}

# Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
#                   "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#SBATCH --qos={{settings.slurm_qos}}
#SBATCH --time={{settings.slurm_time}}
{{settings.slurm_add_directives}}

#SBATCH --output=sbatch.%j.out
##SBATCH --error=sbatch.%j.err

#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT,TIME_LIMIT_90,TIME_LIMIT


TIME_START_SCRIPT=$(date +%s.%N)

# ---------------------------- define project paths ----------------------------
ROOT_WORK_DIR=$PWD
POEM_WORK_DIR=$ROOT_WORK_DIR/POEM
PISM_WORK_DIR=$ROOT_WORK_DIR/PISM

POEM_TOOLS_DIR={{settings.poem_proj_dir}}/bin


# -------------------------------- set parameters ------------------------------
# time & coupling
CPL_TIMESTEP={{settings.coupling_timestep}}            # in years, must be greater or equal 1
MAX_CPL_ITERATION={{settings.max_cpl_iteration}}       # number of coupling iterations

# initial PISM input files
PISM_RESTART_FILE={{settings.pism_exp_dir}}/initdata/{{settings.pism_infile}}

# when continuing from a previous coupled run:
# set filename of last PISM output file processed by PISM-to-MOM_processing.py
#   file must be in $ROOT_WORK_DIR/x_PISM-to-MOM
#   -> only set if restarting from a coupled run
{% if settings.coupled_restart==True %}
PISM_TO_MOM_FLUX_RESTART_FILE={{settings.pism_to_mom_flux_restart_file}}
{% else %}
#PISM_TO_MOM_FLUX_RESTART_FILE=
{% endif %}



# ------------------------------ define functions-------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
{% raw %}
arbitrary() { 
    # creates list of SLURM nodes for specified number of tasks
    #   -> parameters: list of numbers of tasks to assigned to each node
    [[ SHELLOPTS =~ $xtrace ]] && minusx=-x # save the flag setting
    set +x
    declare -a tasks nodes
    tasks=($*)
    nodes=(`scontrol show hostnames $SLURM_NODELIST`)
    node_cnt=${#nodes[*]}
    task_cnt=${#tasks[*]}

    if [ $node_cnt -lt $task_cnt ]
    then
	echo "ERROR: You only have $node_cnt nodes, but requested layout on $task_cnt nodes." >&2
	task_cnt=$node_cnt
    fi

    cnt=0
    layout=""
    #echo tasks ${tasks[*]}
    #echo nodes ${nodes[*]}
    for cnt in `seq 0 $(( $task_cnt - 1 ))`
    do
	task=${tasks[$cnt]}
	node=${nodes[$cnt]}
	for i in `seq 1 $task`
	do
	    [ "" != "$layout" ] && layout="${layout},"
	    layout="${layout}${node}"
	done
    done
    echo "$layout"
    set +x $minusx # restore the flag setting
}
{% endraw %}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
poem_run() {
    echo '  > run POEM'
    echo " ... running POEM from $SIM_START_TIME to $SIM_END_TIME of $SIMULATION_TIME model years"
    START_POEM_RUN=$(date +%s.%N)
    cd $POEM_WORK_DIR

    # decision between parallel and serial case
    #   -> concatenate the values of both variables, one might be the empty string
    if [ 0"$SLURM_NTASKS" -gt 1 -o 0"$SLURM_NNODES" -gt 1 ]
    then
        export MPD_CON_EXT=Slurm_Job_$SLURM_JOBID
    
        # 16 tasks atm + 32 tasks ocn = 48 tasks total
        #local arb="`arbitrary 16 16 16`"
        local arb="`arbitrary 16 16`"
        
        time srun --propagate=ALL -m arbitrary -n $SLURM_NTASKS -w "$arb" \
                -o fms.out-$SLURM_JOBID-$SLURM_NNODES-$SLURM_NTASKS-'%02t' \
                ./fms_MOM_SIS.x
        RESULT=$?
        #mpdallexit

        # merge single task output to one file
        ulimit -c 0 # we dont want coredumps from mppnccombine
        local bases=`ls -1 *.nc.[0-9][0-9][0-9][0-9]* | \
                        sed -e "s/.nc.[0-9][0-9][0-9][0-9]\+/.nc/" | sort -u`
        time for b in $bases ; do 
            echo mppnccombine $b 
            $POEM_TOOLS_DIR/mppnccombine.pik-hlrs2015-ifort -r $b $b.[0-9][0-9][0-9][0-9]* 
        done
    else
       unset SLURM_PMI_KVS_DUP_KEYS
       unset SLURM_JOB_ID SLURM_STEPID SLURM_NPROCS SLURM_PROCID SLURM_GTIDS
       time      ./fms_MOM_SIS.x > fms.out-1-$SLURM_JOBID 2>&1
       RESULT=$?
    fi

    local enddate=`$POEM_TOOLS_DIR/time_stamp.csh -ef digital`
    return_check $RESULT "POEM_run.$enddate"

    END_POEM_RUN=$(date +%s.%N)
    TIME_POEM_RUN=$(echo "$END_POEM_RUN - $START_POEM_RUN" | bc)
    TIME_POEM=$(echo "$TIME_POEM + $TIME_POEM_RUN" | bc)
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
poem_prerun() {
    echo 'poem_prerun func'
    cd $POEM_WORK_DIR

    # create pre-run directory
    if [ ! -d $POEM_WORK_DIR/prerun/out ] ; then 
        mkdir -p $POEM_WORK_DIR/prerun/out
    fi

    # ----------- create prerun version of  ----------------
    # ------- diag_table, data_table, input.nml ------------

    DIAG_TABLE_STORE=$(find diag_table -maxdepth 0 -printf %l)
    DATA_TABLE_STORE=$(find data_table -maxdepth 0 -printf %l)
    INPUT_NML_STORE=$(find input.nml -maxdepth 0 -printf %l)

    # -------- modify input.nml -----------------
    cp input.nml-regular input.nml-prerun
    # comment out all lines with months, days, hours
    sed -e '/^\s*months/ s/^./!&/' -i input.nml-prerun
    sed -e '/^\s*days/   s/^./!&/' -i input.nml-prerun
    sed -e '/^\s*hours/  s/^./!&/' -i input.nml-prerun
    # add new line: days = 1
    sed -e '/&coupler_nml/a\' -e '\tdays = 1,' -i input.nml-prerun
    
    # -------- modify data_table -----------------
    cp data_table-dummy data_table-prerun
    # comment out all lines with FILE_NAME_DUMMY
    #   -> no field overwrite in prerun
    sed -e '/.*FILE_NAME_DUMMY/ s/^./#&/' -i data_table-prerun

    # direct symlinks to prerun versions of lists
    ln -sf diag_table-prerun diag_table
    ln -sf input.nml-prerun input.nml
    ln -sf data_table-prerun data_table     # for pre-run and CPL_ITERATION==1
    
    # ------------- run POEM/FMS --------------------------
    # decision between parallel and serial case
    #   -> concatenate the values of both variables, one might be the empty string
    if [ 0"$SLURM_NTASKS" -gt 1 -o 0"$SLURM_NNODES" -gt 1 ]
    then
        export MPD_CON_EXT=Slurm_Job_$SLURM_JOBID
    
        # 16 tasks atm + 32 tasks ocn = 48 tasks total
        #local arb="`arbitrary 16 16 16`"
        local arb="`arbitrary 16 16`"
        
        time srun --propagate=ALL -m arbitrary -n $SLURM_NTASKS -w "$arb" \
                -o fms.out-$SLURM_JOBID-$SLURM_NNODES-$SLURM_NTASKS-'%02t' \
                ./fms_MOM_SIS.x
        RESULT=$?
        #mpdallexit

        # merge single task output to one file
        ulimit -c 0 # we dont want coredumps from mppnccombine
        local bases=`ls -1 *.nc.[0-9][0-9][0-9][0-9]* | \
                        sed -e "s/.nc.[0-9][0-9][0-9][0-9]\+/.nc/" | sort -u`
        time for b in $bases ; do 
            echo mppnccombine $b 
            $POEM_TOOLS_DIR/mppnccombine.pik-hlrs2015-ifort -r $b $b.[0-9][0-9][0-9][0-9]* 
        done
    else
       unset SLURM_PMI_KVS_DUP_KEYS
       unset SLURM_JOB_ID SLURM_STEPID SLURM_NPROCS SLURM_PROCID SLURM_GTIDS
       time      ./fms_MOM_SIS.x > fms.out-1-$SLURM_JOBID 2>&1
        RESULT=$?
    fi

    # revert links to previous location 
    ln -sf $DIAG_TABLE_STORE diag_table
    #ln -sf $DATA_TABLE_STORE data_table
    ln -sf data_table-prerun data_table
    ln -sf $INPUT_NML_STORE input.nml

    return_check $RESULT "POEM_pre-run"
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
{% raw %}
poem_run_postprocess() {
    echo '  > run POEM postprocessing'
    START_POEM_POSTP=$(date +%s.%N)

    # extract begin/end dates
    local begindate=`$POEM_TOOLS_DIR/time_stamp.csh -bf digital`
    if [ "$begindate" == "" ] ; then 
        begindate=tmp`date +%Y%j%H%M%S`
    fi
    local enddate=`$POEM_TOOLS_DIR/time_stamp.csh -ef digital`
    if [ "$enddate" == "" ] ; then 
        enddate=tmp`date +%Y%j%H%M%S`
    fi
    rm -f time_stamp.out
    # write to global variables
    POEM_TIME_BEGIN=$begindate
    POEM_TIME_END=$enddate
    echo POEM TIME BEGIN $POEM_TIME_BEGIN
    echo POEM TIME END $POEM_TIME_END

    # rename fms.out file
    for f in fms.out-$SLURM_JOBID* ; do
        mv $f $enddate.$(echo $f | sed -e s/$SLURM_JOBID-// )
    done

    # put files to history dir
    for i in *.nc 
      do mv $i history/$enddate.$i
    done
    time tar cvjf history/$enddate.out.tar.bz2 $enddate.fms.out-* logfile*out diag_field_log*out diag_integral.out stocks.out
    
    mkdir -p fms.out/$SLURM_JOBID
    mv $enddate.fms.out* fms.out/$SLURM_JOBID

    # copy files to RESTART dir
    cp -p input.nml data_table diag_table field_table RESTART/.
    ( cd RESTART
      ulimit -c 0 # we dont want coredumps from mppnccombine
      # Concatenate blobs restart files. mppnccombine would not work on them.
      if [ -f ocean_blobs.res.nc.0000 ]
      then
         ncecat ocean_blobs.res.nc.???? ocean_blobs.res.nc
         rm -f ocean_blobs.res.nc.????
      fi
    
      # Concatenate iceberg restarts
      if [ -f ocean_blobs.res.nc.0000 ]
      then
         ncrcat icebergs.res.nc.???? icebergs.res.nc
         rm icebergs.res.nc.????
      fi
    
      # Land restarts need to be combined with  combine-ncc
      # More simply just tar them up in this version
      land_files="cana glac lake land snow soil vegn1 vegn2"
      for file in $land_files
      do
         input_files="`/bin/ls ${file}.res.nc.????`"
         if [ ${#input_files} -gt 0 ]
         then
             tar czf ${file}.res.nc.tar $input_files
             if [ $? != 0 ]
            then
                 echo "ERROR: in creating land restarts tarfile"
                 exit 1
             fi
             rm $input_files
         fi
      done
    
      # now combine all the remaining restart files
      bases=`ls -1 *.nc.[0-9][0-9][0-9][0-9]* | sed -e "s/.nc.[0-9][0-9][0-9][0-9]\+/.nc/" | sort -u`
      [ -n "$bases" ] && time for b in $bases ; do echo mppnccombine $b ; rm -f $b; $POEM_TOOLS_DIR/mppnccombine.pik-hlrs2015-ifort -64 -r $b $b.[0-9][0-9][0-9][0-9]* ; done
    )

    # tar files for restart and put to history
    mv RESTART $enddate.RESTART
    time tar cvjf history/$enddate.RESTART.tar.bz2 $enddate.RESTART/.
    mv $enddate.RESTART RESTART
    
    # copy RESTART files to INPUT
    cd RESTART
    for i in *.res*; do
        rm ../INPUT/$i
        cp $i ../INPUT/. 
    done
    cd $OLDPWD

    END_POEM_POSTP=$(date +%s.%N)
    TIME_POEM_POSTP=$(echo "$END_POEM_POSTP - $START_POEM_POSTP" | bc)
    TIME_POEM_POSTPROC=$(echo "$TIME_POEM_POSTPROC + $TIME_POEM_POSTP" | bc)
}
{% endraw %}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
poem_prerun_postprocess() {

    # rename fms.out file
    for f in fms.out-$SLURM_JOBID* ; do
        mv $f prerun.$(echo $f | sed -e s/$SLURM_JOBID-// )
    done

    # put files to pre-run dir
    for i in *.nc 
      do mv $i prerun/prerun.$i
    done
    mv prerun.fms.out-* logfile*out diag_field_log*out diag_integral.out \
        stocks.out time_stamp.out $POEM_WORK_DIR/prerun/out/
    cp -p input.nml data_table diag_table field_table $POEM_WORK_DIR/prerun
    
    cd $OLDPWD
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pism_run() {
    # $1 : pism runtime (years)
    echo '  > run PISM'
    echo " ... running PISM from $SIM_START_TIME to $SIM_END_TIME of $SIMULATION_TIME model years"
    local __TIMESTEP=$1

    START_PISM_RUN=$(date +%s.%N)
    cd $PISM_WORK_DIR

    # extract model time from input file
    PISM_TIME_SEC=$(ncdump -v time $PISM_RESTART_FILE | grep "time = [0-9]\+ ;" | \
                        awk '{print $3}')
    PISM_TIME_YR=$(echo $PISM_TIME_SEC / '( 365 * 24 * 60 * 60 )' | bc )

    PISM_TIME_BEGIN=$PISM_TIME_YR
    PISM_TIME_END=$(echo $PISM_TIME_BEGIN + $__TIMESTEP | bc )
    PISM_TIME_END_PAD=$(printf "%06d\n" $PISM_TIME_END)
    echo PISM_TIME_BEGIN $PISM_TIME_BEGIN
    echo PISM_TIME_END $PISM_TIME_END

    # save first and last year of coupled simulation output
    if [ $CPL_ITERATION -eq 2 ]; then 
        SIM_START_YEAR=$(printf "%06d\n" $PISM_TIME_BEGIN)
    fi
    SIM_END_YEAR=$(printf "%06d\n" $PISM_TIME_END)

    # export all bash variables required by pism_run.sh
    #export SLURM_NTASKS ROOT_WORK_DIR PISM_PROJ_DIR PISM_WORK_DIR PISM_RESTART_FILE __TIMESTEP  PISM_TIME_BEGIN PISM_TIME_END POEM_TIME_END
    export PISM_RESTART_FILE PISM_TIME_BEGIN PISM_TIME_END PISM_TIME_END_PAD POEM_TIME_END
    # run pism executable
    time ./pism_run_script.sh
    RESULT=$?
    return_check $RESULT "PISM_run.$PISM_TIME_END_PAD"

    # use PISM output file for next restart
    PISM_RESTART_FILE=$(readlink -f results/$PISM_TIME_END_PAD.pism_out.nc)
    echo PISM_RESTART_FILE: $PISM_RESTART_FILE


    END_PISM_RUN=$(date +%s.%N)
    TIME_PISM_RUN=$(echo "$END_PISM_RUN - $START_PISM_RUN" | bc)
    TIME_PISM=$(echo "$TIME_PISM + $TIME_PISM_RUN" | bc)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pism_prerun() {
    # $1 : timestep for pism run 
    echo 'pism_pre_run func'
    cd $PISM_WORK_DIR
    local __TIMESTEP=$1
    
    # extract model time from input file
    PISM_TIME_SEC=$(ncdump -v time $PISM_RESTART_FILE | grep "time = [0-9]\+ ;" | \
                        awk '{print $3}')
    PISM_TIME_YR=$(echo $PISM_TIME_SEC / '( 365 * 24 * 60 * 60 )' | bc )

    PISM_TIME_BEGIN=$PISM_TIME_YR
    PISM_TIME_END=$(echo $PISM_TIME_BEGIN + $__TIMESTEP | bc )
    PISM_TIME_END_PAD=$(printf "%06d\n" $PISM_TIME_END)
    echo PISM_TIME_END $PISM_TIME_END

    if [ ! -d $PISM_WORK_DIR/prerun ] ; then mkdir $PISM_WORK_DIR/prerun ; fi

    # export all bash variables required by pism_run.sh
    export PISM_RESTART_FILE __TIMESTEP PISM_TIME_BEGIN PISM_TIME_END 
    #export SLURM_NTASKS ROOT_WORK_DIR PISM_PROJ_DIR PISM_WORK_DIR PISM_RESTART_FILE PISM_OCEAN_START_FILE __TIMESTEP  PISM_TIME_BEGIN PISM_TIME_END 
    #export P_Mx P_My P_Lz P_Lbz P_Mz P_Mbz 

    # run pism executable
    time ./pism_prerun_script.sh
    RESULT=$?
    return_check $RESULT "PISM_pre-run"
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
process_mom_to_pism() {
    echo " >> process_mom-to-pism function"
    echo
    #echo POEM TIME BEGIN $POEM_TIME_BEGIN
    #echo POEM TIME END $POEM_TIME_END
    
    START_MOM_to_PISM_PROCESS=$(date +%s.%N)
    cd $ROOT_WORK_DIR/inter-model-processing

    echo "  > bilinear regridding MOM output to PISM grid"
    OCEAN_OUT=$POEM_WORK_DIR/history/$POEM_TIME_END.ocean-yearly.nc
    {% if settings.ocean_to_ice_timeseries==True -%}
    # don't average ocean output over coupling time step before regridding
    # but extract temp, salt for regridding
    REGRID_IN=$ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.ocean-yearly.sparse.nc
    ncks -O -v temp,salt,time_bounds $OCEAN_OUT $REGRID_IN
    {% else %}
    # average ocean output over coupling time step before regridding
    REGRID_IN=$ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.ocean-yearly.mean.nc
    ncra -O -v temp,salt -F -d time,1,,1 $OCEAN_OUT $REGRID_IN
    #ncks -O -v  temp,salt -d time,-1 $OCEAN_OUT $REGRID_IN     # only last timestep 
    {% endif %}
    REGRID_OUT=$ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.regrid.MOM-to-PISM.bil.cdo.nc
    WEIGHTS_PATH=$ROOT_WORK_DIR/pre-processing/$WEIGHTS
    cdo -b F64 -f nc4c remap,$PISM_PRE_OUT,$WEIGHTS_PATH $REGRID_IN $REGRID_OUT
    RESULT=$?
    return_check $RESULT "cdo_regridding"

    echo "  > regriddedMOM-to-PISM_processing script"
    ./regriddedMOM-to-PISM_processing.py                                    \
        -i $REGRID_OUT                                                      \
        -b $PISM_WORK_DIR/prerun/prerun.pism_extra.nc                       \
        -e $ROOT_WORK_DIR/pre-processing/pism_edges.nc                      \
        -f temp salt                                                        \
        -d $BASIN_SHELF_DEPTH_FILE                                          \
        -o $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc     \
        -v
    RESULT=$?
    return_check $RESULT "regriddedMOM-to-PISM_processing.py"
        
    {% if settings.do_ocean_anomaly==True -%}
    ### anomaly approach: use pism_ocn_file (e.g. Schmidtko data) for 
    #       PISM-PICO forcing and add difference of 
    #       regriddedMOM-to-PISM_processing.py output from current to first 
    #       coupling time step
    echo OCEAN_ANOMALY_REFERENCE_FILE=$OCEAN_ANOMALY_REFERENCE_FILE
    echo OCEAN_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_ANOMALY_REFERENCE_FILE_NOTIME
    # difference from current ocean output to reference
    ncbo -O --op_typ=subtract $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $OCEAN_ANOMALY_REFERENCE_FILE_NOTIME \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # add difference to given ocean forcing file
    ncbo -O --op_typ=add $ROOT_WORK_DIR/PISM/initdata/{{settings.pism_ocn_file}} \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.anomaly.nc \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time variable to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncks"
        {% if settings.ocean_to_ice_timeseries==True -%}
    ncks -A -v time_bounds \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncks"
        {% endif %}

    {% else %}
    # don't use ocean anomaly approach:
    #       make PISM input file the same as regriddedMOM-to-PISM_processing.py output
    cd $ROOT_WORK_DIR/x_MOM-to-PISM
    ln -sf $POEM_TIME_END.processed_MOM.nc $POEM_TIME_END.PISM_input.nc
    cd $ROOT_WORK_DIR/inter-model-processing
    {% endif %}

    END_MOM_to_PISM_PROCESS=$(date +%s.%N)
    TIME_MOM_to_PISM_PROCESS_T=$(echo "$END_MOM_to_PISM_PROCESS - $START_MOM_to_PISM_PROCESS" | bc)
    TIME_MOM_to_PISM_PROCESS=$(echo "$TIME_MOM_to_PISM_PROCESS + $TIME_MOM_to_PISM_PROCESS_T" | bc)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
process_pism_to_mom(){
    echo " >> process_pism-to-mom function"
    START_PISM_to_MOM_PROCESS=$(date +%s.%N)
    cd $ROOT_WORK_DIR/inter-model-processing

    ./PISM-to-MOM_processing.py                                             \
        -o $PISM_PRE_OUT                                                    \
        -e $PISM_WORK_DIR/results/$PISM_TIME_END_PAD.pism_extra.nc          \
        -m $ROOT_WORK_DIR/pre-processing/PISMbasin-to-MOMcell_mapping.nc    \
        -a $MOM_PRE_OUT                                                     \
        -f $ROOT_WORK_DIR/x_PISM-to-MOM/$PISM_TIME_END_PAD.fluxes.nc        \
        -d $ROOT_WORK_DIR/x_PISM-to-MOM/$PISM_TIME_END_PAD.basin_shelf_depth.nc \
        -v
    RESULT=$?

    {% if not settings.use_prescribed_basin_shelf_depth==True %}
    # save file path for usage in next process_mom_to_pism call
    BASIN_SHELF_DEPTH_FILE=$ROOT_WORK_DIR/x_PISM-to-MOM/$PISM_TIME_END_PAD.basin_shelf_depth.nc
    {% endif %}

    # modify data table to include PISM output fluxes
    cd $ROOT_WORK_DIR/POEM
    cp data_table-dummy data_table-pism_in 
    sed "s/FILE_NAME_DUMMY/$PISM_TIME_END_PAD.fluxes.nc/" -i data_table-pism_in
    ln -sf data_table-pism_in data_table

    return_check $RESULT "PISM-to-MOM_processing.py"

    END_PISM_to_MOM_PROCESS=$(date +%s.%N)
    TIME_PISM_to_MOM_PROCESS_T=$(echo "$END_PISM_to_MOM_PROCESS - $START_PISM_to_MOM_PROCESS" | bc)
    TIME_PISM_to_MOM_PROCESS=$(echo "$TIME_PISM_to_MOM_PROCESS + $TIME_PISM_to_MOM_PROCESS_T" | bc)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
concat_output_files(){
START_CONCAT_FILES=$(date +%s.%N)

    # concatenate MOM output
    cd $POEM_WORK_DIR/history
    INPUT_FILES=$(echo `seq -f "%06g0101.ice-monthly.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-monthly.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ice-yearly.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-yearly.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ice-decadal.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-decadal.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ocean-scalar.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-scalar.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ocean-monthly.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-monthly.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ocean-yearly.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-yearly.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

    INPUT_FILES=$(echo `seq -f "%06g0101.ocean-decadal.nc" \
        $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-decadal.nc" \
        $SIM_START_YEAR $SIM_END_YEAR`)
    ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE
    cd $ROOT_WORK_DIR

    # concatenate PISM output
    cd $PISM_WORK_DIR/results
    ncrcat --overwrite $(echo `seq -f "%06g.pism_extra.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.pism_extra.nc
    ncrcat --overwrite $(echo `seq -f "%06g.pism_snap.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.pism_snap.nc
    ncrcat --overwrite $(echo `seq -f "%06g.pism_ts.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.pism_ts.nc
    cd $ROOT_WORK_DIR
    
    # concatenate inter-model-processing files
    cd $ROOT_WORK_DIR/x_MOM-to-PISM
    ncrcat --overwrite $(echo `seq -f "%06g0101.processed_MOM.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.processed_MOM.nc
    ncrcat --overwrite $(echo `seq -f "%06g0101.processed_MOM.anomaly.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.processed_MOM.anomaly.nc
    ncrcat --overwrite $(echo `seq -f "%06g0101.PISM_input.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.PISM_input.nc
    cd $ROOT_WORK_DIR/x_PISM-to-MOM
    ncrcat --overwrite $(echo `seq -f "%06g.basin_shelf_depth.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.basin_shelf_depth.nc
    ncrcat --overwrite $(echo `seq -f "%06g.fluxes.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
        $SIM_START_YEAR-$SIM_END_YEAR.fluxes.nc
    cd $ROOT_WORK_DIR

END_CONCAT_FILES=$(date +%s.%N)
TIME_CONCAT_FILES=$(echo "$END_CONCAT_FILES - $START_CONCAT_FILES" | bc)

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
return_check() {
    # $1 : RETURN_VAL 
    # $2 : FUNC_STRING
    # usage: return_check $RETURN_VAL "func xy" 

    if [ 0 = "$1" ] ; then
        echo "'"$2"'" ended successfully
    else
        echo "'"$2"'" exited with code $1
        exit_script $1
        #cd $ROOT_WORK_DIR
        #rm -f .workdir_locked_by_batchjob
        #exit $1
    fi
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
exit_script() {
    # $1 : RETURN_VAL 
    # usage: exit_script $RESULT
    cd $ROOT_WORK_DIR
    rm -f .workdir_locked_by_batchjob
    exit $1
}


# ------------------------------- begin script ---------------------------------

export LANG=C

echo $0
date
echo

# ---- initialize time measuring variables ---
TIME_POEM=0
TIME_PISM=0
TIME_POEM_POSTPROC=0
TIME_MOM_to_PISM_PROCESS=0
TIME_PISM_to_MOM_PROCESS=0

# --- load packages ---
module purge
module load poem/2017.4
module load pism/stable1.0
#echo $PETSC_DIR
#PETSC_DIR=/home/kreuzer/software/petsc-3.9.1_db0/
#echo $PETSC_DIR

# load anaconda python environment for pre & inter-model processing scripts
module load anaconda
source activate py3_netcdf_xarray
conda info --envs
conda list

module list

# print some environment variables
set +x
env | egrep SBATCH\|SLURM\|MPI\|OMP\|KMP | sort
set -x

if [ "" = "$I_MPI_PMI_LIBRARY" ]
then
    export I_MPI_PMI_LIBRARY=/p/system/slurm/lib/libpmi.so
fi


# check if executable is currently running to avoid interfering runs
set -x
ln -s $SLURM_JOBID .workdir_locked_by_batchjob
if [ $? -ne 0 ]
then
    # sigh. it might be a restart attempt initated by loadleveler after a vacate order
    # in that case the job ID is still the same that created the existing lock file.
    # And if so, we just start anew silently.
    oldjobid="`ls -l .workdir_locked_by_batchjob | awk '{print $11}'`"
    if [ "$oldjobid" != "$SLURM_JOBID" ]
    then
        echo Error: Directory `pwd` is locked by another batch job: 
        ls -l .workdir_locked_by_batchjob
        exit 1
    else
        echo own lockfile already there - looks like a slurm restart attempt >> fms.$SLURM_JOBID.restart-log
        echo `date` `hostname` >> fms.$SLURM_JOBID.restart-log
    fi
fi


# unset some limits
#ulimit -c unlimited             # unlimit maximum size of core files created
ulimit -s unlimited             # unlimit maximum stack size
ulimit -d unlimited             # unlimit maximum size of a process's data segment
#ulimit -m unlimited             # unlimit maximum resident set size
#ulimit -v unlimited             # unlimit maximum amount of virtual memory for shell
#ulimit -f unlimited             # unlimit maximum size of files written
ulimit -t unlimited             # unlimit maximum amount of cpu time (sec)

echo Limits:
ulimit -a
echo

# --------------- check input parameters -------------------------------

if ! [[ $CPL_TIMESTEP =~ ^[0-9]+$ ]] ; then
    echo "ERROR: CPL_TIMESTEP is $CPL_TIMESTEP  -> must be a positive integer"
    exit -1
elif [ "$CPL_TIMESTEP" -lt "1" ]; then
    echo "ERROR: CPL_TIMESTEP is $CPL_TIMESTEP  -> must be greater or equal 1"
    exit -2
fi
if ! [[ $MAX_CPL_ITERATION =~ ^[0-9]+$ ]] ; then
    echo "ERROR: MAX_CPL_ITERATION is $MAX_CPL_ITERATION  -> must be a positive integer"
    exit -1
elif [ "$MAX_CPL_ITERATION" -lt "1" ]; then
    echo "WARNING: MAX_CPL_ITERATION is $MAX_CPL_ITERATION  -> no iterations will be computed"
fi

SIMULATION_TIME=$(echo "$CPL_TIMESTEP * $MAX_CPL_ITERATION" | bc )
echo 'CPL_TIMESTEP (yrs):    ' $CPL_TIMESTEP
echo 'MAX_CPL_ITERATION:     ' $MAX_CPL_ITERATION
echo 'SIMULATION_TIME (yrs): ' $SIMULATION_TIME




# prepare POEM directories
cd $POEM_WORK_DIR
mkdir -p RESTART history
rm -f *.nc *.nc.[0-9][0-9][0-9][0-9]*
cd $ROOT_WORK_DIR


# -------- modify input.nml -----------------
## comment out all lines with months, days, hours
sed -e '/^\s*months/ s/^./!&/' -i POEM/input.nml-regular
sed -e '/^\s*days/   s/^./!&/' -i POEM/input.nml-regular
sed -e '/^\s*hours/  s/^./!&/' -i POEM/input.nml-regular
 

## add new line: months = $CPL_TIMESTEP * 12
CPL_TIMESTEP_MONTHS=$(echo "$CPL_TIMESTEP * 12" | bc)
# delete all similar occurances first
sed "/months\s*=\s*$CPL_TIMESTEP_MONTHS\s*,/d"  -i POEM/input.nml-regular
# add new line
sed -e '/&coupler_nml/a\' -e "\tmonths = $CPL_TIMESTEP_MONTHS," \
    -i POEM/input.nml-regular


### --------------------- synchronise timestamp of models ---------------------

cd $ROOT_WORK_DIR/pre-processing
PISM_SHIFTED_RESTART_PATH_FILE=$ROOT_WORK_DIR/pre-processing/PISM_restart_timeshift_path.txt

#
echo " >> synchronise timestamps of PISM and MOM/POEM"
./PISM_timeshift.py                         \
    -p $PISM_RESTART_FILE                   \
    -c $POEM_WORK_DIR/INPUT/coupler.res     \
    -s $PISM_SHIFTED_RESTART_PATH_FILE      \
    -v
RESULT=$?
return_check $RESULT "PISM_timeshift.py"

# renew PISM restart file pointer to shifted restart file
PISM_RESTART_FILE=$(cat $PISM_SHIFTED_RESTART_PATH_FILE | tail -1)

cd $ROOT_WORK_DIR

### --------------------------------- preruns ---------------------------------
START_PRERUNS=$(date +%s.%N)
echo
echo ">>> POEM pre-run"
echo
poem_prerun 
poem_prerun_postprocess

echo
echo ">>> PISM pre-run"
echo
pism_prerun 2

END_PRERUNS=$(date +%s.%N)
TIME_PRERUNS=$(echo "$END_PRERUNS - $START_PRERUNS" | bc)
echo
echo ">>> finished preruns"
echo

### ---------------------------- preprocessing --------------------------------
echo
echo ">>> preprocessing"
echo
echo
START_PREPROC=$(date +%s.%N)
cd $ROOT_WORK_DIR/pre-processing

echo " >> bilinear interpolation"
echo "  > remove timestamp variable which confuses CDO "
ncks -x -v timestamp $PISM_WORK_DIR/prerun/prerun.pism_out.nc \
    $PISM_WORK_DIR/prerun/prerun.pism_out-no_timestamp.nc -O
RESULT=$?
return_check $RESULT "nco-operation"


echo "  > create weights to regrid from MOM to PISM grid"
PISM_PRE_OUT=$PISM_WORK_DIR/prerun/prerun.pism_out-no_timestamp.nc
MOM_PRE_OUT=$POEM_WORK_DIR/prerun/prerun.ocean-prerun.nc
WEIGHTS=weights.MOM-to-PISM.bil.cdo.nc
export REMAP_EXTRAPOLATE=off
cdo -b F64 -f nc4c genbil,$PISM_PRE_OUT $MOM_PRE_OUT $WEIGHTS 
RESULT=$?
return_check $RESULT "cdo-weight-generation"

echo "  > regrid prerun output from MOM to PISM grid "
REGRID_IN=$MOM_PRE_OUT
REGRID_OUT=$ROOT_WORK_DIR/x_MOM-to-PISM/prerun.regrid.MOM-to-PISM.bil.cdo.nc
cdo -b F64 -f nc4c remap,$PISM_PRE_OUT,$WEIGHTS $REGRID_IN $REGRID_OUT
RESULT=$?
return_check $RESULT "cdo_regridding"

# identify edges of missing/non-missing values of regridded ocean output 
#   on PISM grid
echo " >> identify edges on PISM grid"
./find_edge.py                      \
    -i $REGRID_OUT                  \
    -o pism_edges.nc                \
    -f temp                         \
    -e temp salt area_t average_DT  \
    -v
RESULT=$?
return_check $RESULT "find_edge.py"


# create mapping between PISM basin and southern MOM cells
echo " >> calculate mapping of PISM basins to MOM cells"
./PISMbasin-to-MOMcell_mapping.py       \
    -b $PISM_PRE_OUT                    \
    -m $MOM_PRE_OUT                     \
    -o PISMbasin-to-MOMcell_mapping.nc  \
    -v
RESULT=$?
return_check $RESULT "PISMbasin-to-MOMcell_mapping.py"


{% if not settings.use_prescribed_basin_shelf_depth==True %}
# create basin_shelf_depth_file
echo " >> create basin_shelf_depth_file from prerun output"
$ROOT_WORK_DIR/inter-model-processing/PISM-to-MOM_processing.py \
    -o $PISM_PRE_OUT                                            \
    -e $PISM_WORK_DIR/prerun/prerun.pism_extra.nc               \
    -m PISMbasin-to-MOMcell_mapping.nc                          \
    -a $MOM_PRE_OUT                                             \
    -f $ROOT_WORK_DIR/x_PISM-to-MOM/prerun.fluxes.nc            \
    -d $ROOT_WORK_DIR/x_PISM-to-MOM/prerun.basin_shelf_depth.nc \
    -v
RESULT=$?
return_check $RESULT "PISM-to-MOM_processing.py"
# save file path for usage in next process_mom_to_pism call
BASIN_SHELF_DEPTH_FILE=$ROOT_WORK_DIR/x_PISM-to-MOM/prerun.basin_shelf_depth.nc
{% else %}
# use prescribed basin shelf depths from file
BASIN_SHELF_DEPTH_FILE={{settings.prescribed_basin_shelf_depth_path}}
{% endif %}


# set PISM-to-MOM when restarting from coupled run
if [ -z ${PISM_TO_MOM_FLUX_RESTART_FILE+x} ]; then 
    echo "PISM_TO_MOM_FLUX_RESTART_FILE is unset"; 
    echo " -> no restart from coupled run"
else 
    echo "PISM_TO_MOM_FLUX_RESTART_FILE is set to '$PISM_TO_MOM_FLUX_RESTART_FILE'"; 
    echo " -> restart from coupled run"
    echo " -> modifying data_table for first coupling iteration"

    # modify data table to include PISM output fluxes for first restart iteration
    cd $ROOT_WORK_DIR/POEM
    cp data_table-dummy data_table-pism_in 
    sed "s/FILE_NAME_DUMMY/$PISM_TO_MOM_FLUX_RESTART_FILE/" -i data_table-pism_in
    ln -sf data_table-pism_in data_table    # for CPL_ITERATION==1
    cd $ROOT_WORK_DIR
fi

### ocean anomaly reference file
# ____________________________________________________________________________ #
{% if settings.do_ocean_anomaly==True -%}
echo
echo "ocean anomaly approach enabled"
    # ........................................................................ #
    {% if settings.use_ocean_anomaly_from_prev_run==True %}
echo "using ocean anomaly reference file {{settings.ocean_anomaly_reference_path}}"
OCEAN_ANOMALY_REFERENCE_FILE=$ROOT_WORK_DIR/x_MOM-to-PISM/{{settings.ocean_anomaly_reference_file}}
echo OCEAN_ANOMALY_REFERENCE_FILE=$OCEAN_ANOMALY_REFERENCE_FILE

    # ........................................................................ #
    {% else %} # compute ocean anomaly reference state
echo "computing ocean anomaly reference state from " \
    {{settings.calc_ocn_anomaly['path']}}

cd {{settings.calc_ocn_anomaly['path']}}
echo "  > [ocean anomaly reference] average yearly model output (possibly "\
    "spread across multile files)"
INPUT_FILES=$(echo `seq -f "{{settings.calc_ocn_anomaly['name_format_in']}}" \
                           "{{settings.calc_ocn_anomaly['yr_start']}}"       \
                           "{{settings.calc_ocn_anomaly['yr_step']}}"        \
                           "{{settings.calc_ocn_anomaly['yr_end']}}" `)
AVG_OUT_FILE=$(echo `printf "{{settings.calc_ocn_anomaly['name_format_out']}}" \
                           "{{settings.calc_ocn_anomaly['yr_start']}}"         \
                           "{{settings.calc_ocn_anomaly['yr_end']}}" `)
ncra --overwrite -v temp,salt $INPUT_FILES $AVG_OUT_FILE
RESULT=$?
return_check $RESULT "ncra-operation"
cd $ROOT_WORK_DIR

echo "  > [ocean anomaly reference] regrid averaged output from MOM to PISM grid "
REGRID_IN={{settings.calc_ocn_anomaly['path']}}/$AVG_OUT_FILE
REGRID_OUT=$ROOT_WORK_DIR/x_MOM-to-PISM/ocean_mean.regrid.MOM-to-PISM.bil.cdo.nc
WEIGHTS_PATH=$ROOT_WORK_DIR/pre-processing/$WEIGHTS
cdo -b F64 -f nc4c remap,$PISM_PRE_OUT,$WEIGHTS_PATH $REGRID_IN $REGRID_OUT
RESULT=$?
return_check $RESULT "cdo_regridding"

echo "  > [ocean anomaly reference] regriddedMOM-to-PISM_processing script"
cd $ROOT_WORK_DIR/inter-model-processing

OUTFILE=$(echo `printf "%06g-%06g.mean.processed_MOM.nc"              \
                       "{{settings.calc_ocn_anomaly['yr_start']}}"  \
                       "{{settings.calc_ocn_anomaly['yr_end']}}" `)
OCEAN_ANOMALY_REFERENCE_FILE=$ROOT_WORK_DIR/x_MOM-to-PISM/$OUTFILE
./regriddedMOM-to-PISM_processing.py                                    \
    -i $REGRID_OUT                                                      \
    -b $PISM_WORK_DIR/prerun/prerun.pism_extra.nc                       \
    -e $ROOT_WORK_DIR/pre-processing/pism_edges.nc                      \
    -f temp salt                                                        \
    -d $BASIN_SHELF_DEPTH_FILE                                          \
    -o $OCEAN_ANOMALY_REFERENCE_FILE                                    \
    -v
RESULT=$?
return_check $RESULT "regriddedMOM-to-PISM_processing.py"
cd $ROOT_WORK_DIR

    {% endif %}
    # ........................................................................ #

# removing time dimension from ocean anomaly reference file 
#  -> needed for clean broadcasting with ncbo
FILE_EXTENSION="${OCEAN_ANOMALY_REFERENCE_FILE##*.}"
FILE_PATH_NO_EXTENSION="${OCEAN_ANOMALY_REFERENCE_FILE%.*}"
OCEAN_ANOMALY_REFERENCE_FILE_NOTIME=$FILE_PATH_NO_EXTENSION.no-time.$FILE_EXTENSION

echo OCEAN_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_ANOMALY_REFERENCE_FILE_NOTIME
ncwa -O -a time $OCEAN_ANOMALY_REFERENCE_FILE $OCEAN_ANOMALY_REFERENCE_FILE_NOTIME


# ____________________________________________________________________________ #
{% else %}
echo
echo "ocean anomaly approach disabled"
{% endif %}
# ____________________________________________________________________________ #


END_PREPROC=$(date +%s.%N)
TIME_PREPROC=$(echo "$END_PREPROC - $START_PREPROC" | bc)
echo
echo ">>> finished preprocessing"
echo 
echo

### -------------------------- coupling iterations ----------------------------
echo 
echo ">>> starting coupling iterations"
echo 
echo

for CPL_ITERATION in `seq 1 $MAX_CPL_ITERATION`
do
    SIM_START_TIME=$(echo "$CPL_TIMESTEP * ($CPL_ITERATION -1)" | bc )
    SIM_END_TIME=$(echo "$CPL_TIMESTEP * $CPL_ITERATION" | bc )
    echo
    echo ' >> CPL_ITERATION=' $CPL_ITERATION

    # run poem for coupling timestep
    poem_run
    poem_run_postprocess

    # process MOM output for PISM/PICO input
    process_mom_to_pism

    # run PISM for coupling timestep
    pism_run $CPL_TIMESTEP 

    # process PISM output for MOM input
    process_pism_to_mom
        
done
echo 
echo ">>> finished coupling iterations"
echo 
echo 

### ---------------------------- postprocessing --------------------------------
echo
echo ">>> postprocessing"
echo 
echo 
concat_output_files


### --------------------------- runtime statistics -----------------------------
set +x

print_stat() {
    # $1 : subroutine name
    # $2 : subroutine time (sec)
    # $3 : total script time (sec)

    format="  %-30s %10.2f \t %6.2f\n"
    TIME_PERCENT=$(echo "$2 / $3 * 100" | bc -l)

    printf "$format" $1 $2  $TIME_PERCENT
}

TIME_END_SCRIPT=$(date +%s.%N)
TIME_SCRIPT=$(echo "$TIME_END_SCRIPT - $TIME_START_SCRIPT" | bc)


str="-------------"
sep_str=$str$str$str$str$str
header="\n  %-30s %10s %11s \n"
width=54

echo
echo 
echo
echo
echo " ---------------- runtime statistics ----------------- "
echo
echo " model time (yrs)          " $SIMULATION_TIME
echo " coupling timestep (yrs)   " $CPL_TIMESTEP
echo
echo
printf "$header" "ROUTINE" "TIME(s)" "ratio(%)"
printf "%$width.${width}s\n" " $sep_str"
print_stat total $TIME_SCRIPT $TIME_SCRIPT
print_stat preruns $TIME_PRERUNS $TIME_SCRIPT
print_stat preprocessing $TIME_PREPROC $TIME_SCRIPT
print_stat POEM_runs $TIME_POEM $TIME_SCRIPT
print_stat POEM_postprocessing $TIME_POEM_POSTPROC $TIME_SCRIPT
print_stat PISM_runs $TIME_PISM $TIME_SCRIPT
print_stat MOM-to-PISM_processing $TIME_MOM_to_PISM_PROCESS $TIME_SCRIPT
print_stat PISM-to-MOM_processing $TIME_PISM_to_MOM_PROCESS $TIME_SCRIPT
print_stat concat_output_files $TIME_CONCAT_FILES $TIME_SCRIPT
printf "%$width.${width}s\n" " $sep_str"

echo
echo


set -x

exit_script $RESULT