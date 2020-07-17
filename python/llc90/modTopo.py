#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, readfield, writefield

topofile='../input/bathy_eccollc_90x50_min2pts_mod.bin'
newtopofile='../input/bathy_eccollc_90x50_min2pts_mod3.bin'
dfile='/uv/home1/jhauck/MITgcm/verification/global_oce_input_fields'
tfile= os.path.join(dfile,'T_OWPv1_M_eccollc_90x50.bin')
sfile= os.path.join(dfile,'S_OWPv1_M_eccollc_90x50.bin')
sssfile= os.path.join(dfile,'SSS_WPv1_M_eccollc_90x50_pm05atl.bin')
tnewfile='../input/T_OWPv3_M_eccollc_90x50.bin'
snewfile='../input/S_OWPv3_M_eccollc_90x50.bin'
sssnewfile= '../input/SSS_WPv3_M_eccollc_90x50_pm05atl.bin'
topo = readfield(topofile,[1170,90],'float32')
t = readfield(tfile,[12,50,1170,90],'float32')
s = readfield(sfile,[12,50,1170,90],'float32')
sss = readfield(sssfile,[12,1170,90],'float32')

tn = np.copy(topo)
tt = np.copy(t)
ss = np.copy(s)
ssss = np.copy(sss)
tf = llc.flat(topo)
# # open Nares Strait
tn[604:606,37]     =  topo[604:606,36]
tt[:,:,604:606,37] = t[:,:,604:606,36]
ss[:,:,604:606,37] = s[:,:,604:606,36]
ssss[:,604:606,37] = sss[:,604:606,36]
tn[606:608,38]     =  topo[606:608,37]
tt[:,:,606:608,38] = t[:,:,606:608,37]
ss[:,:,606:608,38] = s[:,:,606:608,37]
ssss[:,606:608,38] = sss[:,606:608,37]
tn[609,37]     =  topo[609,38]
tt[:,:,609,37] = t[:,:,609,38]
ss[:,:,609,37] = s[:,:,609,38]
ssss[:,609,37] = sss[:,609,38]
#
tn[607:610,36]     =     tn[607:610,37]
tt[:,:,607:610,36] = tt[:,:,607:610,37]
ss[:,:,607:610,36] = ss[:,:,607:610,37]
ssss[:,607:610,36] = ssss[:,607:610,37]
# open Barrow Strait
tn[623,50]     =  topo[623,51]
tn[625,50]     =  topo[625,49]
tt[:,:,623,50] = t[:,:,623,51]
tt[:,:,625,50] = t[:,:,625,49]
ss[:,:,623,50] = s[:,:,623,51]
ss[:,:,625,50] = s[:,:,625,49]
# close bays in CAA
tn[620:622,43:49] = 0.
tn[615:617,47:50] = 0.


writefield(newtopofile,tn)
writefield(tnewfile,tt)
writefield(snewfile,ss)
writefield(sssnewfile,ssss)
