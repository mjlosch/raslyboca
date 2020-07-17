#!/bin/bash
#
# Loop over years for increasing (CORE2) runoff over the Arctic:
#
# Elena Gerwing, Nov 2018

module load cdo/1.7.2

set -ex 

FIRSTYEAR=1979      # first year to be used
LASTYEAR=2011       # last year to be used
#
#Loop over years
YY=${FIRSTYEAR}
while [ ${YY} -le ${LASTYEAR} ]
do

INFILE=/work/ollie/egerwing/datasets/ERAinterim/t2m_ERAi_6hourly_${YY}
OUTFILE=t2m_ERAi_increased_${YY}

python increase_forcing_over_arctic.py ${INFILE} ${OUTFILE}


YY=`expr ${YY} + 1`

done


exit
