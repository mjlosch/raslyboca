#!/bin/bash

set -ex 

FIRSTYEAR=1979      # first year to be used
LASTYEAR=2009       # last year to be used
#
#Loop over years
YY=${FIRSTYEAR}
while [ ${YY} -le ${LASTYEAR} ]
do
    FILE=http://data1.gfdl.noaa.gov/~nnz/mom4/COREv2/data_IAF/CORRECTED/calendar_years/q_10.${YY}.05APR2010.nc
    wget ${FILE}

YY=`expr ${YY} + 1`
done

# rm ${INPUT} ${OUTPUT}

exit
