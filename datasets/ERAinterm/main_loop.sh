#!/bin/bash
#
# Loop over years for processing of ERA-Interim netcdf data for usage in MITgcm:
# - calculate specitic humidity
# - transform into binary format
# - convert units of tp, ssrd, strd
#
# Elena Gerwing, March 2017

module load cdo/1.7.2

set -ex 

FIRSTYEAR=1979      # first year to be used
LASTYEAR=1979       # last year to be used
#
#Loop over years
YY=${FIRSTYEAR}
while [ ${YY} -le ${LASTYEAR} ]
do

#cdo merge ERAi_instantaneous_6hourly_${YY}.nc ERAi_accumulated_6hourly_${YY}.nc ERA_allparams_${YY}.nc

INFILE=ERA_allparams_${YY}.nc
OUTFILE=_ERAi_6hourly_${YY}
TP_OUTFILE=tp_corrected_ERAi_6hourly_${YY}


# Python program for calculating specific humidity
python specHum_from_T_and_P.py ${INFILE}
# Convert to binary and convert units of ssrd, strd and tp
python prepare_forcing_for_MITgcm.py ${INFILE} ${OUTFILE} --convert --flip_y
# Correct precipitation forcing
python correct_precipitation.py ${INFILE} ${TP_OUTFILE} --convert --flip_y

YY=`expr ${YY} + 1`
done

# rm ${INPUT} ${OUTPUT}

exit
