#!/bin/ksh --login

#SBATCH -A fv3lam
#SBATCH -p service
#SBATCH -J hpss_MPAS_Ensemble_2
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -o /scratch2/BMC/fv3lam/ajohns/tohpss/hpss_MPAS_Ensemble_2.log
#SBATCH -e /scratch2/BMC/fv3lam/ajohns/tohpss/hpss_MPAS_Ensemble_2.log


mpas_init_path=/scratch2/BMC/fv3lam/ajohns/mpasruns/wofs_mpas/mpas_init
mpas_run_path=/scratch2/BMC/fv3lam/ajohns/mpasruns/wofs_mpas/retroenkf

uppdir=uppfiles #name of directory where the grib2 files from upp end up for each member

DIR="/ESRL/BMC/fv3lam/3year/HIWT/MPAS_Ensemble_Spring_2" #_2 for enkf ICs, _1 for hrrr+gefs ICs

module load hpss

hsi mkdir -p ${DIR}

date=("2022050200") # "2021051200" "2021051300" "2021051400")

for cycle in "${date[@]}"
do
  cd ${mpas_init_path}/${cycle}/
  hsi mkdir -p ${DIR}/${cycle}/
  for mem in `ls -d ens*`
  do
    cd ${mpas_init_path}/${cycle}/${mem}
    hsi mkdir -p ${DIR}/${cycle}/${mem}

#get stuff from mpas_init directory
    cd ${mpas_init_path}/${cycle}/${mem}
    htar -cvf ${DIR}/${cycle}/${mem}/iclbc.tar init.nc lbc*nc
    htar -cvf ${DIR}/${cycle}/${mem}/ungribdata.tar FILE*
    htar -cvf ${DIR}/${cycle}/${mem}/initnamelists.tar streams* namelist*
    htar -cvf ${DIR}/${cycle}/${mem}/initlogfiles.tar log*

#get stuff from mpas run directory
    cd ${mpas_run_path}/${cycle}/fc/${mem}   #not sure if both experiments use same directory structure here
    htar -cvf ${DIR}/${cycle}/${mem}/runnamelists.tar streams* namelist*
    htar -cvf ${DIR}/${cycle}/${mem}/runlogfiles.tar log*

    for file in `ls diag*nc`
    do
      hsi put ${file} : ${DIR}/${cycle}/${mem}/${file}
    done

    for file in `ls history*nc`
    do
      hsi put ${file} : ${DIR}/${cycle}/${mem}/${file}
    done

#get mpassit files. will probably need modification for other experiment
    for dir in `ls -d mpassit_*`
    do
       mv ${dir}/out_hist_diag.nc ./${dir}.nc
       hsi put ${dir}.nc : ${DIR}/${cycle}/${mem}/${dir}.nc
    done

#get upp files. will probably need modification for other experiment
    htar -cvf ${DIR}/${cycle}/${mem}/uppoutputs.tar ${uppdir}/*

  done




done



  

