#!/bin/bash

member=${2}

mpasdir=/scratch2/BMC/fv3lam/ajohns/mpasruns/wofs_mpas/retroenkf/2022051200/fc/ens_${member}

ftime=${1} #202205122000
rundir=${mpasdir}/mpassit_${ftime}

################################

 
timestr=${ftime:0:4}-${ftime:4:2}-${ftime:6:2}_${ftime:8:2}.${ftime:10:2}.00

mkdir ${rundir}

sed "s!DDDD!${rundir}!g" templatesmpassit/run_mpassit.slm > ${rundir}/run.slm


sed "s!DDDD!${mpasdir}!g" templatesmpassit/mpassit_namelist.input > ${rundir}/namelist.input
sed -i "s!TTTT!${timestr}!g" ${rundir}/namelist.input
sed -i "s!RRRR!${rundir}!g" ${rundir}/namelist.input

cp templatesmpassit/diaglist ${rundir}/
cp templatesmpassit/histlist_soil ${rundir}/
cp templatesmpassit/histlist_2d ${rundir}/
cp templatesmpassit/histlist_3d ${rundir}/

cd ${rundir}
sbatch run.slm
exit 


