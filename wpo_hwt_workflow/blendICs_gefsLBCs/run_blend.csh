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

#SBATCH -J blend
#SBATCH -o logs/blend.%j
#SBATCH -e logs/blend.%j
#SBATCH -N 30
#SBATCH -n 270
#SBATCH --exclusive
#SBATCH --partition=hera
#SBATCH -t 03:00:00
#SBATCH -A hmtb

#
# This script runs BLEND-ing
#

module --force purge

module use /scratch2/BMC/fv3lam/mayfield/wpo/blending/ufs-srweather-app/modulefiles
module load build_hera_intel
#module list

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
   echo "running BLEND on driver.csh node for member $mem"
else
   echo "what is your batch sytem? exit"
   exit 2
endif

set start_init = $start_init   # From driver
set DATE = $DATE   # From driver
set diag_output_interval = $diag_output_interval   # From driver
#set end_time = `$TOOL_DIR/da_advance_time.exe $DATE $FCST_RANGE`

# -------------------------------------
# Get current date into proper format
# -------------------------------------
set date_file_format =  `${TOOL_DIR}/da_advance_time.exe $DATE 0 -f ccyy-mm-dd_hh.nn.ss`
echo "Date in format is: ${date_file_format}"
#set vhr = `echo ${date_file_format}`

# ----------------------------------
# Make and go to working directory
# ----------------------------------
set rundir = $BLEND_OUTPUT_DIR_TOP/$start_init/ens_${mem}
mkdir -p $rundir
cd $rundir

#-----------------------------------------------------------
# Link necessary input files and code and fill namelist
#-----------------------------------------------------------
ln -sf ${BLEND_CODE_DIR}/mpas_blending .
ln -sf ${SCRIPT_DIR}/blend_files/* .
ln -sf ${BLEND_GRID_DIR}/*.grid.nc .

#setenv grid_file $MPAS_INIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/init.nc
#setenv hist_file $EXP_DIR_TOP/$start_init/ens_${mem}/history.${date_file_format}.nc
#setenv diag_file $EXP_DIR_TOP/$start_init/ens_${mem}/diag.${date_file_format}.nc
#setenv output_file $MPASSIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/proc.${date_file_format}.nc

#enkf NOTE: these ensemble members go from 0 to 9!
setenv small_file /scratch2/BMC/fv3lam/ajohns/enkf_forec_upponly/$start_init/ens_$((mem - 1))/init.nc
#gefs out of maps_init
setenv large_file $MPAS_INIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/init.nc
#blend
setenv blend_file $BLEND_OUTPUT_DIR_TOP/$start_init/ens_${mem}/init_blended.nc

$NAMELIST_TEMPLATE blend

#----------------------------------------------------
# Run blending
#----------------------------------------------------

set run_cmd = "mpiexec -n 270 -ppn 18" # MPI run command
#mpirun -np 400 ./mpassit namelist.input
rm -f PET*
$run_cmd ./mpas_blending
rm -f PET*

if ( $status != 0 ) then
  echo "BLEND failed. Exit." >> ./FAIL
  exit 6
endif

# recombine the small file with the blended fields
module load gnu
module load intel/2023.2.0
module load netcdf/4.7.0
module load nco

ncks -h -A $small_file init_blend_combined.nc
ncks -h -A $blend_file init_blend_combined.nc

exit 0
