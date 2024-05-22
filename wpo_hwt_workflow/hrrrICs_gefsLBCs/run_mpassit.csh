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

#SBATCH -J mpassit
#SBATCH -o mpassit.%j
#SBATCH -e mpassit.%j
#SBATCH -n 1200
#SBATCH --exclusive
#SBATCH --partition=hera
#SBATCH -t 02:00:00
#SBATCH -A fv3lam

#
# This script runs MPASSIT
#

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
   echo "running MPASSIT on driver.csh node for member $mem"
else
   echo "what is your batch sytem? exit"
   exit 2
endif

set start_init = $DATE   # From driver
set DATE = $DATE   # From driver
set diag_output_interval = $diag_output_interval   # From driver
set end_time = `$TOOL_DIR/da_advance_time.exe $DATE $FCST_RANGE`

while ( $DATE <= $end_time)

   # -------------------------------------
   # Get current date into proper format
   # -------------------------------------
   set date_file_format =  `${TOOL_DIR}/da_advance_time.exe $DATE 0 -f ccyy-mm-dd_hh.nn.ss`
   echo "Date in format is: ${date_file_format}"
   set fhr = `echo ${date_file_format} | cut -c12-13`
   echo "Forecast hour is: ${fhr}"

   # ----------------------------------
   # Make and go to working directory
   # ----------------------------------
   set rundir = $MPASSIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/f${fhr}
   mkdir -p $rundir
   cd $rundir

   #-----------------------------------------------------------
   # Link necessary input files and code and fill namelist
   #-----------------------------------------------------------
   ln -sf ${MPASSIT_CODE_DIR}/bin/mpassit .
   ln -sf ${SCRIPT_DIR}/mpassit_files/* .

   setenv grid_file $MPAS_INIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/init.nc
   setenv hist_file $EXP_DIR_TOP/$start_init/ens_${mem}/history.${date_file_format}.nc
   setenv diag_file $EXP_DIR_TOP/$start_init/ens_${mem}/diag.${date_file_format}.nc
   setenv output_file $MPASSIT_OUTPUT_DIR_TOP/$start_init/ens_${mem}/proc.${date_file_format}.nc

   $NAMELIST_TEMPLATE mpassit

   #----------------------------------------------------
   # Run MPASSIT
   #----------------------------------------------------
   rm -f ./*.log
   rm -f ./core*
   rm -f ./*.err

   #   $run_cmd ./mpassit namelist.input # Run MPASSIT!
   mpirun -np 400 ./mpassit namelist.input

   if ( $status != 0 ) then
      echo "MPASSIT failed. Exit." >> ./FAIL
      exit 6
   endif

   # Modify dx in MPASSIT output (it is wonky)
   module load intel/2022.1.2
   module use /contrib/METplus/modulefiles
   module load metplus/5.1.0
   echo "python changedx.py ${output_file}"
   /scratch1/BMC/dtc/miniconda/miniconda3/envs/metplus_v5.1_py3.10/bin/python changedx.py ${output_file} 

   # Done with this foreecast hour; go to next one
   set DATE = `$TOOL_DIR/da_advance_time.exe ${DATE} $diag_output_interval`
end # loop over time/initializations

exit 0
