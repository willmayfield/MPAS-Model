#!/bin/bash

member=$1

outdir=/scratch2/BMC/fv3lam/ajohns/mpasruns/wofs_mpas/retroenkf/2022051200/fc/ens_${member}/uppfiles
mkdir ${outdir}

dateone=20220512
datetwo=20220513

ftimes[1]=${dateone}0000
ftimes[2]=${dateone}0100
ftimes[3]=${dateone}0200
ftimes[4]=${dateone}0300
ftimes[5]=${dateone}0400
ftimes[6]=${dateone}0500
ftimes[7]=${dateone}0600
ftimes[8]=${dateone}0700
ftimes[9]=${dateone}0800
ftimes[10]=${dateone}0900
ftimes[11]=${dateone}1000
ftimes[12]=${dateone}1100
ftimes[13]=${dateone}1200
ftimes[14]=${dateone}1300
ftimes[15]=${dateone}1400
ftimes[16]=${dateone}1500
ftimes[17]=${dateone}1600
ftimes[18]=${dateone}1700
ftimes[19]=${dateone}1800
ftimes[20]=${dateone}1900
ftimes[21]=${dateone}2000
ftimes[22]=${dateone}2100
ftimes[23]=${dateone}2200
ftimes[24]=${dateone}2300
ftimes[25]=${datetwo}0000
ftimes[26]=${datetwo}0100
ftimes[27]=${datetwo}0200
ftimes[28]=${datetwo}0300
ftimes[29]=${datetwo}0400
ftimes[30]=${datetwo}0500
ftimes[31]=${datetwo}0600
ftimes[32]=${datetwo}0700
ftimes[33]=${datetwo}0800
ftimes[34]=${datetwo}0900
ftimes[35]=${datetwo}1000
ftimes[36]=${datetwo}1100
ftimes[37]=${datetwo}1200

itime=1
while [ ${itime} -le 37 ];do #37

echo "${itime}"

./run_upp.sh ${ftimes[${itime}]} ${outdir} ${member}

itime=$((itime+1))
done


