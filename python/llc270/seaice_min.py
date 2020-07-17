import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol
import datetime
import os, sys, re

myrun = 'spinup_era_interim'
#myrun = 'runz00'
bdir = '/work/ollie/mlosch/raslyboca/llc270'
gdir = '/home/ollie/mlosch/MITgcm/MITgcm/llc270/grid'
xg,yg=rdmds(os.path.join(gdir,'XG')),rdmds(os.path.join(gdir,'YG'))
xc,yc=rdmds(os.path.join(gdir,'XC')),rdmds(os.path.join(gdir,'YC'))

dxg=rdmds(os.path.join(gdir,'DXG'))
dyg=rdmds(os.path.join(gdir,'DYG'))
drf=rdmds(os.path.join(gdir,'DRF'))
rac=rdmds(os.path.join(gdir,'RAC'))
drf=rdmds(os.path.join(gdir,'DRF'))

n = dxg.shape[-1]

#rf=rdmds(os.path.join(gdir,'RF'))[:,0,0]
rf=rdmds(os.path.join(bdir,myrun,'RF'))[:,0,0]

iters=[]
for file in os.listdir(os.path.join(bdir,myrun)):
    m = re.search('SIdiags.*.meta',file)
    if m: iters.append(int(m.group().split('.')[1]))

iters = np.unique(np.asarray(iters))

icevol  = np.zeros(len(iters),)
icearea = np.zeros(len(iters),)
icevolsept = []
iceareasept = []
xtimesept = []

refdate = datetime.datetime(1979,2,1,0,0)

timeday = np.asarray(iters)*1800./86400. - 15.
xtime = np.array([refdate + datetime.timedelta(days=i) for i in timeday])

for k, myiter in enumerate(iters):
    print(k,myiter,xtime[k])
    heff = rdmds('SIdiags',itrs=myiter,rec=1,usememmap=True)
    area = rdmds('SIdiags',itrs=myiter,rec=0,usememmap=True)
    # northern hemisphere
    icearea[k] = ((area*rac)[yc>60.]).sum()
    icevol[k]  = ((heff*rac)[yc>60.]).sum()
#    icevol[k] = (heff[6*n:7*n,:]*rac[6*n:7*n,:]).sum()
    if xtime[k].month==9:
        xtimesept.append(xtime[k])
        icevolsept.append(icevol[k])
        iceareasept.append(icearea[k])


fig, ax = plt.subplots(2,1)
# variable SIatmFW is positive upwards, leading sea level decrease, hence the sign
ax[0].plot(xtime,icevol*1e-12,label='sea ice volume')
ax[0].plot(xtimesept,np.asarray(icevolsept)*1e-12,label='September minimum')
ax[1].plot(xtime,icearea*1e-12,label='sea ice area')
ax[1].plot(xtimesept,np.asarray(iceareasept)*1e-12,label='September minimum')

ax[0].set_ylabel('10$^3$ km$^3$')
ax[1].set_ylabel('10$^3$ km$^2$')
for xx in ax:
    xx.grid()
    xx.legend()

fig.show()
