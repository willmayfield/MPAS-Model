#!/bin/csh
#
##BSUB -n 256
##BSUB -J MPAS_fcst
##BSUB -o output.fc
##BSUB -e output.fc
##BSUB -q regular
##BSUB -P NMMM0032
##BSUB -W 60
##BSUB -R "span[ptile=16]"

##PBS -S /bin/csh
##PBS -N forecast
##PBS -A NMMM0021
##PBS -l walltime=10:00
##PBS -q regular
##PBS -o forecast.out
##PBS -j oe 
##PBS -k oed
##PBS -l select=6:ncpus=32:mpiprocs=32
##PBS -m n
##PBS -M schwartz@ucar.edu
##PBS -V

#SBATCH -J mpas
#SBATCH -o logs/mpas.%j
#SBATCH -e logs/mpas.%j
#SBATCH -n 1200
#SBATCH --exclusive
#SBATCH --partition=hera
#SBATCH -t 04:00:00
#SBATCH -A fv3lam

#
# This script runs MPAS atmosphere
#

# MH Load modules and set envars for MPAS on Hera (updated for Rocky 8):
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

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/scratch2/BMC/fv3lam/HWT/code/jasper/miniconda3_RL/lib

# -----------------------------------------------------
# Deal with batch stuff and specifying ensemble member
# -----------------------------------------------------
if ( $batch_system == LSF ) then
   set mem = $LSB_JOBINDEX
else if ( $batch_system == PBS ) then
   set mem = $PBS_ARRAY_INDEX #PBS with "-J" flag (PBS pro)
else if ( $batch_system == none ) then
   set mem = $1
else if ( $batch_system == SBATCH ) then
   set mem = $1
   echo "running MPAS_atmosphere on driver.csh node for member $mem"
else
   echo "what is your batch sytem? exit"
   exit 2
endif

set DATE = $DATE   # From driver
echo  "DATE is $DATE"

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

# ----------------------------------------------
# Get the initial conditions
# ---------------------------------------------- 
#set mpas_ics = ${MPAS_INIT_OUTPUT_DIR_TOP}/${DATE}/ens_${mem}/init.nc # This is the only input file needed
set mpas_ics = ${BLEND_OUTPUT_DIR_TOP}/${DATE}/ens_${mem}/init_blend_combined.nc # This is the only input file needed

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

$run_cmd ./atmosphere_model # Run the model!

if ( $status != 0 ) then
   echo "MPAS failed. Exit." >> ./FAIL
   exit 6
endif

# Cleanup to save space
#   rm -f ./output*.nc # not used

exit 0
