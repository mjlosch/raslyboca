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

INFILE=/work/ollie/egerwing/datasets/CORE2_IAF/runoff_CORE2_iaf_monthly_${YY}
OUTFILE=runoff_CORE2_iaf_increased_${YY}
GRIDPATH=/work/ollie/egerwing/MITgcm/llc90/grid/

python increase_runoff_over_arctic.py ${INFILE} ${OUTFILE} ${GRIDPATH} --latitude 60.0
#python increase_runoff_over_arctic.py ${INFILE} ${OUTFILE} ${GRIDPATH}

YY=`expr ${YY} + 1`

done


exit
