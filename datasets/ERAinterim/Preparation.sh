#!/bin/bash

#cdo -b B -f srv copy ERA_1314_allparams.nc ERA_1314_allparams.srv
#cdo -b B -f srv copy ERA_1314_allparams.nc ERA_1314_allparams.srv
#cdo mergetime ERA_1314_allparams.srv ERA_1516_allparams.srv ERA_1316_allparams.srv

cdo splitvar ERA_1316_allparams.srv ERA_1316_

cdo splityear ERA_1316_d2m.srv d2m_12hourly-step12_
cdo splityear ERA_1316_sp.srv sp_12hourly-step12_
cdo splityear ERA_1316_ssrd.srv ssrd_12hourly-step12_
cdo splityear ERA_1316_strd.srv strd_12hourly-step12_
cdo splityear ERA_1316_t2m.srv t2m_12hourly-step12_
cdo splityear ERA_1316_tp.srv tp_12hourly-step12_
cdo splityear ERA_1316_u10.srv u10_12hourly-step12_
cdo splityear ERA_1316_v10.srv v10_12hourly-step12_




