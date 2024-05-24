#!/bin/csh
#########################################################################
# Script: driver.csh
#
# Purpose: Top-level driver script to run MPAS forecasts
#
# Spring 2024: modified from Craig Schwartz's scripts for application to 
#        IC blending experiments on Hera
#########################################################################
#
#SBATCH -J driver
#SBATCH -o driver.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 6:00:00 
#SBATCH -A fv3lam


#atj: from env.rocky8.build:
module purge
module load cmake/3.28.1
module load gnu
module load intel/2023.2.0
module load impi/2023.2.0
module load pnetcdf/1.12.3
module load szip
module load hdf5parallel/1.10.5
module load netcdf-hdf5parallel/4.7.0
setenv PNETCDF /apps/pnetcdf/1.12.3/intel_2023.2.0-impi
setenv CMAKE_C_COMPILER mpiicc
setenv CMAKE_CXX_COMPILER mpiicpc
setenv CMAKE_Fortran_COMPILER mpiifort




# Basic system information
setenv batch_system         SBATCH # Type of submission system. 
setenv num_procs_per_node   40  # Number of processors per node on the machine
setenv run_cmd              "mpirun" # Command to run MPI executables #"mpirun.lsf"

# Total number of processors to use (not the number of nodes)
setenv NUM_PROCS_MPAS_ATM   1200 
setenv NUM_PROCS_MPAS_INIT  1 #360 
setenv num_mpas_procs_per_node  40 #$num_procs_per_node # Number of processors per node you want to use for MPAS forecasts (not intialization)

setenv mpas_init_walltime 240      #Run-time (minutes) for MPAS initialization
setenv mpas_fcst_walltime 240      #Run-time (minutes) for MPAS forecasts
setenv mpas_account   "fv3lam"  # core-hour account
setenv mpas_queue     "hera"   # system queue

# Decide which stages to run (run if true; lowercase):
setenv RUN_UNGRIB              false   # (true, false )
setenv RUN_MPAS_INITIALIZE     false
setenv RUN_MPAS_FORECAST       true

######################################
#Directories pointing to source code #
######################################

setenv   SCRIPT_DIR           /scratch2/BMC/fv3lam/ajohns/simple_mpas_scripts_dev1  # Location of all these .csh scripts

# Path to MPAS initialization code
setenv   MPAS_INIT_CODE_DIR   /scratch2/BMC/fv3lam/HWT/code/MPAS-Model-NSSL_develop
#setenv   MPAS_INIT_CODE_DIR   /scratch2/BMC/fv3lam/ajohns/test/MPAS-Model

# Path to MPAS model code (could be different from initialization code)
setenv   MPAS_CODE_DIR        $MPAS_INIT_CODE_DIR

# Path to WPS; need ungrib.exe for MPAS initialization
setenv   WPS_DIR              /scratch2/BMC/fv3lam/ajohns/WPS

setenv   TOOL_DIR             /scratch2/BMC/fv3lam/ajohns/WRFDA/var/build/  # Location of ./da_advance_time.exe, built from WRFDA
setenv   VTABLE_DIR           $SCRIPT_DIR # Location of where ungrib.exe VTables are
setenv   WPS_GEOG_DIR         /scratch2/BMC/fv3lam/HWT/WPS_GEOG  #Directory with WPS geogrid information, needed for MPAS_INIT

############################################################
# This controls the top-level directory of the experiment  #
############################################################

setenv MESH      wofs_mpas  # The MPAS mesh. Can really name whatever you want
#setenv EXPT      dev1      # Generate and use own static file
#setenv EXPT      dev2      # use static file from Michelle
setenv EXPT      retroenkf      #production runs

###########################
# Time/Experiment control #
###########################
 
#setenv start_init    2022050100  # starting and ending forecast initialization times
#setenv start_init    2022050200  # starting and ending forecast initialization times
#setenv start_init    2022050300  # starting and ending forecast initialization times
#setenv start_init    2022050400  # starting and ending forecast initialization times
#setenv start_init    2022050500  # starting and ending forecast initialization times
#setenv start_init    2022050600  # starting and ending forecast initialization times
#setenv start_init    2022050700  # starting and ending forecast initialization times
#setenv start_init    2022050800  # starting and ending forecast initialization times
#setenv start_init    2022050900  # starting and ending forecast initialization times
#setenv start_init    2022051000  # starting and ending forecast initialization times
#setenv start_init    2022051100  # starting and ending forecast initialization times
setenv start_init    2022051200  # starting and ending forecast initialization times
setenv end_init      $start_init
setenv inc_init      24  #Time (hours) between forecast initializations

setenv FCST_RANGE              36    #Length of MPAS forecasts (hours)
setenv diag_output_interval    1     # Diagnostic file output frequency (hours)
setenv restart_output_interval 24000 # Restart file output frequency (hours)

setenv update_sst  .false.    # Use SST from a different dataset, other than the IC files?

setenv MPAS_REGIONAL .true.   # If .true., we are doing ungrib and mpas_init for a regional run, and running MPAS for a regional run
setenv LBC_FREQ 3          # LBC frequency for a regional run (hours)

#This is the model providing initial conditions for MPAS. Mostly needed to tell the MPAS initialization
#  how many vertical levels to expect in the GRIB files.  See run_mpas_init.csh
#atj: separate this into ic model, lbc ctl model, lbc pert model; and add "RRFS" option
setenv  COLD_START_INITIAL_CONDITIONS_MODEL RRFS # (GFS_FV3, GFS_new, GFS,GEFS, RRFS)
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_CTL GFS #atj: new variable
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_PERT GEFS #atj: new variable
 
setenv ENS_SIZE             9 # Ensemble size for the forecasts. Ensemble forecasts for all members are run all at once.
set    ie        =          0 # What ensemble member to start with?  Usually set to 1. Members ${ie} - ${ENS_SIZE} will be run

#atj: only need to generate static for one member
#setenv ENS_SIZE             0
#set    ie        =          0

####################################################
# These control input into the sequence of scripts #
####################################################



setenv      NAMELIST_TEMPLATE        ${SCRIPT_DIR}/namelist_template.csh  #MPAS namelist tempalte. Much is "hard-wired" here, some filled on-the-fly.
setenv      STREAMS_TEMPLATE         ${SCRIPT_DIR}/streams_template.csh # Template for MPAS streams. Much is hard-wired, but some filled on-the-fly.

  # example directory for input data "by ensemble member and date" : /glade/scratch/schwartz/GEFS_grib_0.5deg/ens_3/2023052300
  #  this directory needs to be there, with the data. For the input directories, the date needs to be last (after the member).
setenv      GRIB_INPUT_DIR_MODEL     /scratch2/BMC/fv3lam/ajohns/allgribdata # Location of global GRIB files (sub-directories by ensemble member and initialization time)
setenv      GRIB_INPUT_DIR_SST       $GRIB_INPUT_DIR_MODEL  # Where you could put GRIB files for SST

setenv      ungrib_prefx_model   "FILE"  # Probably never need to change, but there might be some need to in the future
setenv      ungrib_prefx_sst     "FILE"  #"SST"  # Usually "SST"...but can set to "FILE" (or $ungrib_prefx_model) if wanting to use SST from GFS

####################################################
# Where input files are stored
####################################################




   # example of directory "by ensemble member and date": /glade/scratch/schwartz/MPAS/ungrib_met/2023052300/ens_3
   #   these will be created

setenv   UNGRIB_OUTPUT_DIR_MODEL      /scratch2/BMC/fv3lam/ajohns/mpasruns/ungrib_met  #Top-level WPS/Ungrib output directory for model data (sub-dirs by ensemble member and date)
setenv   UNGRIB_OUTPUT_DIR_SST        /scratch2/BMC/fv3lam/ajohns/mpasruns/ungrib_sst  #Top-level WPS/Ungrib output directory for SST data    (sub-dirs by ensemble member and date)
setenv   MPAS_INIT_OUTPUT_DIR_TOP     /scratch2/BMC/fv3lam/ajohns/mpasruns/${MESH}/mpas_init   #Top-level directory holding MPAS initialization files (sub-dirs by ens member date)
setenv   EXP_DIR_TOP                  /scratch2/BMC/fv3lam/ajohns/mpasruns/${MESH}/${EXPT}  #Directory where MPAS forecasts are run

########################
# MPAS mesh settings
########################

# -------------------------------------------------------
# Settings for regional or global forecasts
# -------------------------------------------------------
setenv    num_mpas_cells    1792182
#setenv    num_mpas_cells    1894063 

setenv    MPAS_GRID_INFO_DIR      /scratch2/BMC/fv3lam/HWT/CONUS_grid_static # Directory containing MPAS grid files, must be there
setenv    graph_info_prefx        conus.graph.info.part. #x1.${num_mpas_cells}.graph.info.part.     # Should be located in $MPAS_GRID_INFO_DIR, files must be there
setenv    grid_file_netcdf        ${MPAS_GRID_INFO_DIR}/conus.grid.nc #${MPAS_GRID_INFO_DIR}/x1.${num_mpas_cells}.grid.nc  # Needs to be there before doing anything for global...not for regional though b/c you'll start with static.nc
#setenv    mpas_static_data_file   /scratch2/BMC/fv3lam/ajohns/conus.static.nc      # Will be created during initialization, just needs to be created once.
setenv    mpas_static_data_file   ${MPAS_GRID_INFO_DIR}/conus.static.nc      # Will be created during initialization, just needs to be created once.





# --------------------------------------------------------------------------------------------
# Vertical grid dimensions, same for both the ensemble and high-res determinsitic forecasts
# --------------------------------------------------------------------------------------------
setenv    num_mpas_vert_levels   59      # Number of vertical levels ( mass levels )
setenv    num_mpas_soil_levels   9       # Number of soil levels
setenv    num_soilcat            16     #atj: added to be consistent with NSSL
setenv    z_top_meters           25878.712      # atj: set directly in meters since not integer number of km
#setenv    z_top_km               20      # MPAS model top (km)

############################################################################
# MPAS model settings; can also hardcode in $NAMELIST_TEMPLATE if you want
############################################################################
setenv time_step            20.0   # Seconds. Typically should be 4-6*dx; use closer to 4 for cycling DA
setenv radiation_frequency  30     # Minutes. Typically the same as dx (for dx = 15 km, 15 minutes)
setenv config_len_disp      3000. # Meters, diffusion length scale, which should be finest resolution in mesh (not needed in MPASv8.0+)
setenv soundings_file       dum #${SCRIPT_DIR}/sounding_locations.txt   # set to a dummy to disable soundings
setenv physics_suite        "convection_permitting" #"mesoscale_reference" #"convection_permitting_wrf390" 
setenv deep_conv_param      "off" #"cu_ntiedtke" #"cu_grell_freitas" # "tiedtke"

# You can override physics suite individual parameterizations with these
#  Make sure you put these options in the MPAS namelist
#setenv microphysics_param        "mp_thompson"
#setenv lsm_param                 "noah"
#setenv pbl_param                 "bl_mynn"
#setenv longwave_rad_param        "rrtmg_lw"
#setenv shortwavae_rad_param      "rrtmg_sw"
#setenv sfc_layer_param           "sf_mynn"

#########################################################
#
# NOTHING BELOW HERE SHOULD NEED CHANGING (HOPEFULLY...)
#
#########################################################

set ff = ( $NAMELIST_TEMPLATE  $STREAMS_TEMPLATE )
foreach f ( $ff ) 
   if ( ! -e $f ) then
     echo "$f doesn't exist. Abort!"
     exit
   endif
end

#-------------------------

setenv DATE $start_init
while ( $DATE <= $end_init )

   echo "Processing $DATE"

   setenv FCST_RUN_DIR   ${EXP_DIR_TOP}/${DATE}/fc # Where MPAS forecasts are run

   if ( $RUN_UNGRIB == true ) then
#     bsub -J "ungrib_${DATE}" -q geyser -n 1 < ${SCRIPT_DIR}/run_ungrib.csh
      foreach mem ( `seq $ie 1 $ENS_SIZE` )
         ${SCRIPT_DIR}/run_ungrib.csh $mem # send the ensemble member into run_ungrib.csh
      end
   endif

   if ( $RUN_MPAS_INITIALIZE == true ) then
      set queue_opts = "-q economy -n $NUM_PROCS_MPAS_INIT -P $mpas_account -W 30"
      set this_num_procs_per_node = $num_procs_per_node
      set this_num_needed_nodes = `echo "$NUM_PROCS_MPAS_INIT / $this_num_procs_per_node" | bc`
      if ( $this_num_needed_nodes == 0 ) then # Mostly useful for testing in serial mode to avoid lots of log files
	set this_num_needed_nodes = 1
	set this_num_procs_per_node = 1
      endif
      set last_member = $ENS_SIZE
      if ( $batch_system == LSF ) then
	 bsub -J "mpas_init_${DATE}" $queue_opts < ${SCRIPT_DIR}/run_mpas_init.csh
      else if ( $batch_system == PBS ) then
	 set pp_init = `qsub -N "mpas_init_${DATE}" -A "$mpas_account" -q "$mpas_queue" \
	                     -l walltime=${mpas_init_walltime}:00 -J ${ie}-${last_member} \
	                     -l "select=${this_num_needed_nodes}:ncpus=${this_num_procs_per_node}:mpiprocs=${this_num_procs_per_node}" \
			      ${SCRIPT_DIR}/run_mpas_init.csh`
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $ie 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_mpas_init.slm  $mem #atj
          #${SCRIPT_DIR}/run_mpas_init.csh  $mem #atj: running on same single processor as driver for now
        end
        #sbatch ${SCRIPT_DIR}/run_mpas_init.csh 
      else if ( $batch_system == none ) then
         ${SCRIPT_DIR}/run_mpas_init.csh 1 # really just for testing
      endif
   endif

   if ( $RUN_MPAS_FORECAST == true ) then
      if ( $batch_system == LSF ) then
	 set pp_fcst = `bsub -J "run_mpas_${DATE}[${ie}-${ENS_SIZE}]" -n $NUM_PROCS_MPAS_ATM -W $mpas_fcst_walltime -q regular < ${SCRIPT_DIR}/run_mpas.csh`
      else if ( $batch_system == PBS ) then
	 set num_needed_nodes = `echo "$NUM_PROCS_MPAS_ATM / $num_mpas_procs_per_node" | bc`
	     # Use the next line if you DON'T want a dependency condition based on completion of run_mpas_init.csh
	#set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${ie}-${ENS_SIZE} \
	 set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${ie}-${ENS_SIZE} -W depend=afterok:${pp_init} \
		     -q "$mpas_queue" \
		     -l "select=${num_needed_nodes}:ncpus=${num_mpas_procs_per_node}:mpiprocs=${num_mpas_procs_per_node}" \
		     -l walltime=${mpas_fcst_walltime}:00 ${SCRIPT_DIR}/run_mpas.csh`
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $ie 1 $ENS_SIZE` )
        #foreach mem ( `seq 0 1 0` )
          ${SCRIPT_DIR}/run_mpas.csh  $mem #atj: just create the batch file to be submitted separately for now
        end
      else if ( $batch_system == none ) then
         ${SCRIPT_DIR}/run_mpas.csh 1 # really just for testing
      endif
   endif

   # Done with this initialization; go to next one
   setenv DATE `$TOOL_DIR/da_advance_time.exe $DATE $inc_init`

end # loop over time/initializations

exit 0
