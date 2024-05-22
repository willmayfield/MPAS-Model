#!/bin/csh
#
##BSUB -n 16
##BSUB -J MPAS_init
##BSUB -o mpas_init.out
##BSUB -e mpas_init.err
##BSUB -q regular
##BSUB -P NMMM0001
##BSUB -W 10
##BSUB -R "span[ptile=16]"

##PBS -S /bin/csh
##PBS -N init
##PBS -A NMMM0021
##PBS -l walltime=8:00   
##PBS -q economy
##PBS -o init.out
##PBS -j oe 
##PBS -l select=16:ncpus=36:mpiprocs=36
##PBS -k eod 
##PBS -V 

#
#SBATCH -J logs/mpas_init
#SBATCH -o logs/mpas_init.o%j
#SBATCH -N 30
#SBATCH -n 1200
##SBATCH -p hera
#SBATCH -t 01:00:00
#SBATCH -A fv3lam

#
# This script runs MPAS initialization
#
set DATE = $DATE   # From driver.csh
set START_DATE_MPAS = `${TOOL_DIR}/da_advance_time.exe $DATE 0 -w`
set END_DATE_MPAS   = `${TOOL_DIR}/da_advance_time.exe $DATE $FCST_RANGE -w` # Used for regional MPAS for LBCs
set START_DATE_LBC = `${TOOL_DIR}/da_advance_time.exe $DATE $LBC_FREQ -w`
echo "Start LBC ${START_DATE_LBC}"

# Deal with batch stuff and specifying ensemble member
if ( $batch_system == LSF ) then
   set mem = $LSB_JOBINDEX
else if ( $batch_system == PBS ) then
   set mem = $PBS_ARRAY_INDEX #PBS with "-J" flag (PBS pro)
else if ( $batch_system == none ) then
   set mem = $1
else if ( $batch_system == SBATCH ) then
   set mem = $1
   echo "running MPAS_init on driver.csh node for member $mem"
else
   echo "what is your batch sytem? exit"
   exit 2
endif

########

set sdate = `echo "$START_DATE_MPAS" | cut -c 1-13`
set sdate_lbc = `echo "$START_DATE_LBC" | cut -c 1-13`
#setenv z_top_meters `expr $z_top_km \* 1000` # needed for namelist
setenv z_top_meters 

# MH Took out GFS and GFS_new -- this is for GFS FV3
if ( $COLD_START_BOUNDARY_CONDITIONS_MODEL == GFS ) then # GFS FV3
   setenv num_ungrib_vertical_levels_lbc  51 # Could be as many as 50
else if ( $COLD_START_BOUNDARY_CONDITIONS_MODEL == GEFS ) then
   setenv num_ungrib_vertical_levels_lbc 32 # 27 32
else
   echo "$COLD_START_BOUNDARY_CONDITIONS_MODEL not valid.  Specify a valid one. exit"
   exit
endif

# MH Took out GFS and GFS_new -- this is for GFS FV3
# MH Added HRRR.pressure
if ( $COLD_START_INITIAL_CONDITIONS_MODEL == GFS ) then # GFS FV3
   setenv num_ungrib_vertical_levels  34 # Could be as many as 50
else if ( $COLD_START_INITIAL_CONDITIONS_MODEL == GEFS ) then
   setenv num_ungrib_vertical_levels  32 # 27
else if ( $COLD_START_INITIAL_CONDITIONS_MODEL == RRFS ) then
   setenv num_ungrib_vertical_levels  41
else if ( $COLD_START_INITIAL_CONDITIONS_MODEL == HRRR.pressure ) then
   setenv num_ungrib_vertical_levels  46
else
   echo "$COLD_START_INITIAL_CONDITIONS_MODEL not valid.  Specify a valid one. exit"
   exit
endif

# MH changes for HRRR and GFS
setenv num_ungrib_soil_levels  9
setenv num_ungrib_soil_levels_lbc  4

########

set rundir = ${MPAS_INIT_OUTPUT_DIR_TOP}/${DATE}/ens_${mem} # make and go to the run directory
mkdir -p $rundir
cd $rundir
rm -f ./*.err # clean-up any error files from last time

ln -sf ${MPAS_INIT_CODE_DIR}/init_atmosphere_model . # link executable

#
# First iteration  : Interpolate static data onto the domain--only need to do once per domain/mesh
# Second iteration : Interpolate SST update data onto the domain, if $update_sst = .true.
# Third  iteration : Interpolate meteorological data onto domain with proper number of vertical levels
# Fourth iteration : Make boundary conditions, using f00 HRRR -- GFS only has LBC data from 6-h onward
# Fifth iteration : Make boundary conditions
#
# Changed to 2 3 4 5 --> already have static file
#foreach iter ( 1 2 3 4 5 )
foreach iter ( 3 4 5 )

   if ( $iter == 1 ) then

      if ( -e $mpas_static_data_file ) continue  # Only need to do this once per mesh

      setenv case_number   7
      setenv config_stop_time   $START_DATE_MPAS

      setenv config_static_interp   .true.
      setenv config_vertical_grid   .false.
      setenv config_met_interp      .false.
      setenv config_input_sst       .false.
      setenv config_frac_seaice     .false.

      setenv config_extrap_airtemp  linear

      setenv config_input_name    $grid_file_netcdf       # Full path, from driver.csh # MH needed ro change to lower case
      setenv config_output_name   `basename $mpas_static_data_file` # Full path, from driver.csh

      set suffix = "static_data"
      
   else if ( $iter == 2 ) then
      if ( $update_sst != .true. && $update_sst != true ) continue
      if ( ! -e ${UNGRIB_OUTPUT_DIR_SST}/${DATE}/ens_${mem}/${ungrib_prefx_sst}:${sdate} ) continue

      ln -sf ${UNGRIB_OUTPUT_DIR_SST}/${DATE}/ens_${mem}/${ungrib_prefx_sst}* .

      setenv case_number   8
      setenv config_stop_time   $START_DATE_MPAS

      setenv config_static_interp   .false.
      setenv config_vertical_grid   .false.
      setenv config_met_interp      .false.
      setenv config_input_sst       .true.
      setenv config_frac_seaice     .true.

      setenv config_extrap_airtemp  linear

      setenv config_input_name    $mpas_static_data_file
      setenv config_output_name   "./none"   # Output is really sfc_update.nc

      set suffix = "sst_data"

   else if ( $iter == 3 ) then
      echo "In iter 3"
      if ( ! -e ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/ens_${mem}/${ungrib_prefx_model}:${sdate} ) continue
      ln -sf ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/ens_${mem}/${ungrib_prefx_model}:${sdate} .

      echo "Linked iter 3 files."

      setenv case_number   7
      setenv config_start_time   $START_DATE_MPAS
      setenv config_stop_time   $START_DATE_MPAS

      setenv config_static_interp   .false.
      setenv config_vertical_grid   .true.
      setenv config_met_interp      .true.
      setenv config_input_sst       .false.
      setenv config_frac_seaice     .true.

      setenv config_extrap_airtemp  linear

      ln -sf $mpas_static_data_file  ./static.nc
      ln -sf $MPAS_GRID_INFO_DIR/${graph_info_prefx}* .
      setenv config_input_name    "./static.nc"
      setenv config_output_name   "./init.nc"
      setenv ic_output_interval   "initial_only"
      setenv lbc_output_interval   "none"

      setenv this_ungrib_soil_levels $num_ungrib_soil_levels #atj added
      setenv this_ungrib_vertical_levels $num_ungrib_vertical_levels #atj added

      set suffix = "met_data"

   else if ( $iter == 4 ) then

      if ( ! -e ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/ens_${mem}/${ungrib_prefx_model}:${sdate} ) continue
      ln -sf ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/ens_${mem}/${ungrib_prefx_model}:${sdate} .

      if ( $MPAS_REGIONAL == true || $MPAS_REGIONAL == .true. ) then # only do this if regional
	 setenv case_number   9
         setenv config_start_time  $START_DATE_MPAS
	 setenv config_stop_time   $START_DATE_MPAS

	 setenv config_static_interp   .false. # set these to false, but probably doesn't matter
	 setenv config_vertical_grid   .false.
	 setenv config_met_interp      .false.
	 setenv config_input_sst       .false.
	 setenv config_frac_seaice     .false.

         setenv config_extrap_airtemp  linear

	 ln -sf $MPAS_GRID_INFO_DIR/${graph_info_prefx}* .

	 setenv config_input_name    ./init.nc  # use the initial conditions file as input
	 if ( ! -e $config_input_name ) then
	    echo "MPAS initialization failed for iteration ${iter}."
	    echo "The file $config_input_name is missing."
	    exit 1
	 endif
	 setenv config_output_name   "./none"
	 setenv ic_output_interval   "initial_only"
	 setenv lbc_output_interval  "${LBC_FREQ}:00:00"

         setenv this_ungrib_soil_levels $num_ungrib_soil_levels #atj added
         setenv this_ungrib_vertical_levels $num_ungrib_vertical_levels #atj added

	 set suffix = "met_data_lbcs_f00"
      endif

   else if ( $iter == 5 ) then

      if ( ! -e ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/lbc_${mem}/${ungrib_prefx_model}:${sdate_lbc} ) continue
      ln -sf ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/lbc_${mem}/${ungrib_prefx_model}* .

      if ( $MPAS_REGIONAL == true || $MPAS_REGIONAL == .true. ) then # only do this if regional
         echo ${START_DATE_LBC}

         setenv case_number   9
	 setenv config_start_time  ${START_DATE_LBC}
         setenv config_stop_time   $END_DATE_MPAS

         setenv config_static_interp   .false. # set these to false, but probably doesn't matter
         setenv config_vertical_grid   .false.
         setenv config_met_interp      .false.
         setenv config_input_sst       .false.
         setenv config_frac_seaice     .false.

         setenv config_extrap_airtemp  linear

	 ln -sf $MPAS_GRID_INFO_DIR/${graph_info_prefx}* .

         setenv config_input_name    ./init.nc  # use the initial conditions file as input
         if ( ! -e $config_input_name ) then
            echo "MPAS initialization failed for iteration ${iter}."
            echo "The file $config_input_name is missing."
            exit 1
         endif
         setenv config_output_name   "./none"
         setenv ic_output_interval   "initial_only"
         setenv lbc_output_interval  "${LBC_FREQ}:00:00"

         setenv this_ungrib_soil_levels $num_ungrib_soil_levels_lbc #atj added
         setenv this_ungrib_vertical_levels $num_ungrib_vertical_levels_lbc #atj added

         set suffix = "met_data_lbcs"
      else
         continue # move on (b/c $iter == 4 here, this will exit the loop)
      endif

    echo "LEAVING ITER 5" 

   else
      echo "Iter = $iter is wrong.  Should be <= 5.  Big bug somehwere.  Exit"
      exit 0
   endif


   echo "Made it here."
   # ---------------------------
   # Fill namelist and streams
   # ---------------------------
   $NAMELIST_TEMPLATE  mpas_init 
   set namelist_file = namelist.init_atmosphere # Output of above command

   $STREAMS_TEMPLATE   mpas_init  
   set stream_file = streams.init_atmosphere # Output of above command

   # -----------------------
   # Run MPAS initialization
   # -----------------------
   mpirun -np 1200 ./init_atmosphere_model

   # Error check
   if ( $status != 0 ) then
      echo "MPAS initialization failed for iteration ${iter}."
      exit 1
   endif

   # Cleanup
   if ( $iter == 1 ) then
      set dd = `dirname $mpas_static_data_file`  # Get directory where this file is located
      mkdir -p $dd # make sure the directory is there; it might not be at the very beginning
      mv $config_output_name $dd  # Once moved, this file is $mpas_static_data_file
      cp ./${namelist_file}  ${dd}/${namelist_file}_${suffix}
      cp ./${stream_file}    ${dd}/${stream_file}_${suffix}
   endif

   if ( $iter == 2 ) rm -f $config_output_name  # This file is created, even though we don't want it to be. It is the same as $mpas_static_data_file

   mv ./${namelist_file}  ./${namelist_file}_${suffix}  # Save a namelist for each step
   mv ./${stream_file}    ./${stream_file}_${suffix}  # Save a stream for each step
   mv ./log.0000.err      ./log.0000.err_${suffix}
   mv ./log.0000.out      ./log.0000.out_${suffix}

end # Iter

exit 0
