#!/bin/bash
#
# Loop over years for processing of CORE2 (runoff) data for usage in MITgcm:
# - interpolate to llc grid
# - transform into binary format
#
# Elena Gerwing, March 2017

module load cdo/1.7.2

set -ex 

FIRSTYEAR=1979      # first year to be used
LASTYEAR=2011       # last year to be used
#
#Loop over years
YY=${FIRSTYEAR}
while [ ${YY} -le ${LASTYEAR} ]
do

INFILE=/work/ollie/egerwing/datasets/CORE2_IAF/netcdf/runoff_diatren_iaf_${YY}.nc
OUTFILE=runoff_CORE2_iaf_monthly_${YY}
GRIDPATH=/work/ollie/egerwing/MITgcm/IAF/final_run/
TMP=tmp

python move_runoff_to_MITgcm_llc_grid_points.py ${INFILE} ${TMP} ${GRIDPATH}

python move_points_to_ocean.py ${TMP} ${OUTFILE} ${GRIDPATH}

YY=`expr ${YY} + 1`
rm ${TMP}

done


exit
