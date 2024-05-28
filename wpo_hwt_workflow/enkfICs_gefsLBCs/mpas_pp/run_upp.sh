#!/bin/bash

member=${3}

mpasdir=/scratch2/BMC/fv3lam/ajohns/mpasruns/wofs_mpas/retroenkf/2022051200/fc/ens_${member}

ftime=${1} #202205122000
indir=${mpasdir}/mpassit_${ftime}
rundir=${mpasdir}/upp_${ftime}
outdir=${2}

templatesdir=/scratch2/BMC/fv3lam/ajohns/mpas_pp/templatesupp
 
#####################################
#

module purge
module load cmake/3.28.1
module load gnu
module load intel/2023.2.0
module load impi/2023.2.0
module load pnetcdf/1.12.3
module load szip
module load hdf5parallel/1.10.5
module load netcdf-hdf5parallel/4.7.0


mkdir ${rundir}
infile=${indir}/out_hist_diag.nc

#python ${templatesdir}/../templatesmpassit/changedx.py ${infile}

cp templatesupp/* ${rundir}/


iyear=${ftime:0:4}
imonth=${ftime:4:2}
iday=${ftime:6:2}
ihour=${ftime:8:2}

sed "s!FFFF!${infile}!g" ${templatesdir}/itag_templ > ${rundir}/itag
sed -i "s!YYYY!${iyear}!g" ${rundir}/itag
sed -i "s!MMMM!${imonth}!g" ${rundir}/itag
sed -i "s!DDDD!${iday}!g" ${rundir}/itag
sed -i "s!HHHH!${ihour}!g" ${rundir}/itag



sed "s!DDDD!${rundir}!g" ${templatesdir}/run_templ.slm > ${rundir}/run.slm
sed -i "s!EEEE!${outdir}!g" ${rundir}/run.slm


cd ${rundir}


sbatch run.slm
exit




