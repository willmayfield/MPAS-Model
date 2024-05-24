#!/bin/csh
#


#atj: cut out a lot of batch-job stuff from Craig's script to clarify the parts that are actually used in this script

set mem = $1


set DATE = $DATE   # From driver
set lbc_dir = "dummy"  # will be reset below if doing regional MPAS


   # ----------------------------------
   # Make and go to working directory
   # ----------------------------------
set rundir = ${FCST_RUN_DIR}/ens_${mem}
mkdir -p $rundir
cd $rundir
   # ----------------------------------------------
   # Set boundary condition directory if regional
   # ----------------------------------------------
if ( $MPAS_REGIONAL == .true. || $MPAS_REGIONAL == true ) then 
   set lbc_dir = ${MPAS_INIT_OUTPUT_DIR_TOP}/${DATE}/ens_${mem}
endif

echo "lbc_dir: $lbc_dir"

   # ----------------------------------------------
   # Get the initial conditions
   # ---------------------------------------------- 
set mpas_ics = ${MPAS_INIT_OUTPUT_DIR_TOP}/${DATE}/ens_${mem}/init.nc # This is the only input file needed

set pp = ( $mpas_ics ) # array incase you want to add other files to look for
foreach p ( $pp )
   if ( ! -e $p ) then
       echo "$p does not exist" > ./MISSING_FILE
       exit 3
   endif
end

ln -sf $mpas_ics ./init.nc # link the initial conditions to the working directory

   # If true, we want MPAS to look for a surface stream that reads a file with external sst/xice
if ( $update_sst == .true. || $update_sst == true ) then
   set sst_fname = ${MPAS_INIT_OUTPUT_DIR_TOP}/${DATE}/ens_${mem}/sfc_update.nc
   if ( -e $sst_fname ) then
      ln -sf $sst_fname  ./sfc_update.nc
   else
      echo "SST file ${sst_fname} not there. Exit." > MISSING_SST
      exit 1
   endif
   setenv update_sst_interval "99_00:00:00" # setting to "initial_only" doesn't work. give it an absurd length; will still read at initial time
else 
   setenv update_sst_interval  none # namelist.template.csh uses this to set the MPAS namelist. also used in streams_template.csh
endif

   #---------------------------------------------------------------------
   # Link necessary input files and code, and fill namelist and streams
   #---------------------------------------------------------------------


setenv NUM_NODES_MPAS_ATM `expr $NUM_PROCS_MPAS_ATM \/ $num_mpas_procs_per_node` #atj added

ln -sf ${MPAS_CODE_DIR}/atmosphere_model .
ln -sf ${MPAS_CODE_DIR}/*.TBL .
ln -sf ${MPAS_CODE_DIR}/*.DBL .
ln -sf ${MPAS_CODE_DIR}/*DATA .

ln -sf ${MPAS_GRID_INFO_DIR}/${graph_info_prefx}* . # graph files
if ( -e $soundings_file ) ln -sf $soundings_file ./sounding_locations.txt

$NAMELIST_TEMPLATE mpas ${NUM_NODES_MPAS_ATM} ${num_mpas_procs_per_node} # Output is ./namelist.atmosphere #atj modified
 
$STREAMS_TEMPLATE mpas # Output is ./streams.atmosphere

if ( $MPAS_REGIONAL == .true. || $MPAS_REGIONAL == true ) ln -sf ${lbc_dir}/lbc*.nc .

   #----------------------------------------------------
   # Run MPAS
   #----------------------------------------------------
rm -f ./restart* # remove these just in case they were there from a previous failure
rm -f ./log*.abort
rm -f ./core*
rm -f ./*.err


##atj: create a batch job to separately submit this member
#sed "s/DDDD/debug/g" $SCRIPT_DIR/mpas_slurm_template.slm > ./run.slm #atj added
sed "s+DDDD+$rundir+g" $SCRIPT_DIR/mpas_slurm_template.slm > ./run.slm #atj added
sed -i "s/nnnn/$NUM_PROCS_MPAS_ATM/g" ./run.slm #atj added
sed -i "s/NNNN/$NUM_NODES_MPAS_ATM/g" ./run.slm #atj added
sed -i "s/MMMM/$mem/g" ./run.slm #atj added

#   $run_cmd ./atmosphere_model # Run the model!

#   if ( $status != 0 ) then
#      echo "MPAS failed. Exit." >> ./FAIL
#      exit 6
#   endif

   # Cleanup to save space
#   rm -f ./output*.nc # not used

   #--------------------
   # Go to next member
   #--------------------

exit 0
