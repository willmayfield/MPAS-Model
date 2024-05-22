#!/bin/csh

#
# This script runs UNGRIB, which is a serial code.
#

set DATE = $DATE   # From driver
set START_DATE_WRF = `${TOOL_DIR}/da_advance_time.exe $DATE 0 -w`
set date_metgrid   = `${TOOL_DIR}/da_advance_time.exe $DATE 0 -f ccyy-mm-dd_hh`
set mem = $1 # ensemble member from input

#----------------------
# Zero iteration (added for HWT): IC data
# First iteration: MODEL fields, LBC data
# Second iteration : SST, if $update_sst = .true., assuming a 24-hr interval
#
foreach iter ( 0 1 2 )
   if ( $MPAS_REGIONAL == .true. || $MPAS_REGIONAL == true ) then # $MPAS_REGIONAL from driver
      set END_DATE     =  `${TOOL_DIR}/da_advance_time.exe $DATE $FCST_RANGE`  # $FCST_RANGE from driver
      @ LBC_FREQ_SEC = `expr $LBC_FREQ \* 3600`
   else
      set END_DATE     =  $DATE # For global run, only process SST and model at initial time...no boundary conditions needed
      @ LBC_FREQ_SEC = 86400
   endif

   if ( $iter == 0 ) then
      set output_dir = ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/ens_${mem}
      set output_prefx = $ungrib_prefx_model # From driver
      set VTABLE_TYPE = $COLD_START_INITIAL_CONDITIONS_MODEL # From driver
      set UNGRIB_INPUT_DIR = ${GRIB_INPUT_DIR_MODEL}/ens_${mem} # From driver
      set END_DATE     =  $DATE # atj: IC only needed at initial time.

   else if ( $iter == 1 ) then
      set output_dir = ${UNGRIB_OUTPUT_DIR_MODEL}/${DATE}/lbc_${mem}
      set output_prefx = $ungrib_prefx_model # From driver
      set VTABLE_TYPE = $COLD_START_BOUNDARY_CONDITIONS_MODEL
      set UNGRIB_INPUT_DIR = ${GRIB_INPUT_DIR_MODEL}/lbc_${mem} # From driver

   else if ( $iter == 2 ) then
      if ( $update_sst != .true. && $update_sst != true ) continue

      set output_dir = ${UNGRIB_OUTPUT_DIR_SST}/${DATE}/ens_${mem}
      set output_prefx = $ungrib_prefx_sst # From driver
      set VTABLE_TYPE = "SST" 
      set UNGRIB_INPUT_DIR = ${GRIB_INPUT_DIR_SST}/ens_${mem} # From driver

      set END_DATE     =  $DATE # SST only needed at initial time. Even if MPAS regional, SST not needed for LBCs

   endif

   set END_DATE_WRF =  `${TOOL_DIR}/da_advance_time.exe $END_DATE 0 -w`

   #---- Make and go to working directory ---
   mkdir -p $output_dir
   cd $output_dir

   # Fill WPS namelist
   rm -f namelist.wps
   cat > namelist.wps << EOF
&share
 wrf_core = 'ARW',
 max_dom = 1,
 start_date = '$START_DATE_WRF',
 end_date   = '$END_DATE_WRF',
 interval_seconds = ${LBC_FREQ_SEC},
 io_form_geogrid = 2,
 debug_level = 0
/

&ungrib
 out_format = 'WPS'
 prefix     = '${output_prefx}',
/
EOF

   # Link VTable
   if ( -e ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} ) then
      ln -sf ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} ./Vtable # $VTABLE_DIR from driver
   else
      echo "Vtable ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} is missing."
      exit
   endif

   # Link GRIB files
   set files = "${UNGRIB_INPUT_DIR}/${DATE}/*"
   rm -f ./GRIBFILE*
   echo "$files" #atj debug
   ${WPS_DIR}/link_grib.csh $files
   if ( ! -e ./GRIBFILE.AAA ) then
      echo "All files missing.  Exit"
      exit
   endif

   # Link executables and run
   ln -sf ${WPS_DIR}/ungrib/ungrib.exe .
   ln -sf ${WPS_DIR}/util/rd_intermediate.exe .  # For querying the output file only, otherwise not used

   ./ungrib.exe # Run ungrib.exe

   if ( $status != 0) then
      echo "Ungrib.exe failed"
      exit 1
   endif

   rm -f ./PFILE:* # Clean temporary files

end # end loop over "iter"

exit 0
