import sys, os
from getopt import gnu_getopt as getopt
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from MITgcmutils import rdmds
import datetime

bdir = '/work/ollie/mlosch/raslyboca/llc270'
zdir = os.path.join(bdir,'runz00/20years')
pdir = os.path.join(bdir,'runp00/20years')

from glob import glob
zfiles=glob(os.path.join(zdir,'stdout.*'))
pfiles=glob(os.path.join(pdir,'stdout.*'))
#
def getKey(item):
    return item[0]

def get_parms (fname):
    rhoIce = 910.
    with open(fname) as f:
        for line in f:
            if '/* Monitor output interval ( s ). */' in line:
                ll = next(f).strip().split()
                mondt = float(ll[-1].replace('E','e'))
            elif 'mass2rUnit' in line:
                ll = next(f).strip().split()
                m2rUnit = float(ll[-1].replace('E','e'))
            elif '/* Reference density (Boussinesq)  ( kg/m^3 ) */' in line:
                ll = next(f).strip().split()
                rhoConst = float(ll[-1].replace('E','e'))
            elif 'startDate_1' in line:
                startDate = line.split('=')[-1]
            elif '/* Gravitational acceleration ( m/s^2 ) */' in line:
                ll = next(f).strip().split()
                gravity = float(ll[-1].replace('E','e'))
            elif 'gravity orientation relative to vertical coordinate' in line:
                ll = next(f).strip().split()
                gravitySign = float(ll[-1].replace('E','e'))

    return mondt, m2rUnit, rhoConst, gravity*gravitySign, startDate

def get_output (fnames, mystring):
    """parse fname and get some numbers out"""
    timev = []
    myvar = []
    if mystring[:3]=='exf':
        timename = 'exf_time_sec'
    elif mystring[:3]=='sea':
        timename = 'seaice_time_sec'
    else:
        timename = 'time_secondsf'
    for fname in fnames:
        try:
            f = open(fname)
        except:
            print(fname + " does not exist, continuing")
        else:
            for line in f:
                if timename in line:
                    ll = line.split()
                    timev.append(float(ll[-1].replace('D','e')))
                elif mystring in line:
                    ll = line.split()
                    myvar.append(float(ll[-1].replace('D','e')))

            f.close()

    # reverse order
    timevs=np.asarray(timev[::-1])
    myvars=np.asarray(myvar[::-1])
    # timevs=np.asarray(timev)
    # myvars=np.asarray(myvar)
    # This sorts again in ascending order and returns the index of
    # the first occurrence of duplicates. Because we have reverted the order
    # before, in this way we use the values at the beginning of a pickup run
    # rather than the overlapping values of the previous (potentially crashed)
    # run
    timevs, isort = np.unique(timevs,return_index=True)
    myvars=myvars[isort]

    return timevs, myvars
# done

def correct_jumps(x):
    ii = np.where(np.abs(np.diff(x)) > 0.1*x.std())[0]
    y = np.copy(x)
    for i in ii:
#        print(i,y[i],y[i+1],2*y[i+1]-y[i])
        y[i+1:] = y[i+1:]-y[i+1]+y[i]

    return y

mondt, m2zUnits, rhoConst, gravityz, startdate = get_parms(zfiles[0])
mondt, m2pUnits, rhoConst, gravityp, startdate = get_parms(pfiles[0])
refdate = datetime.datetime(int(startdate[:4]),
                            int(startdate[4:6]),int(startdate[6:8]),0,0)
timesecz, etaz = get_output(zfiles, 'dynstat_eta_mean')
timesecp, etap = get_output(pfiles, 'dynstat_eta_mean')

rac = rdmds(os.path.join(zdir,'RAC'))
hfc = rdmds(os.path.join(zdir,'hFacC'),lev=[0])
globalArea = (rac*hfc).sum()

grac = rac/globalArea
files=glob(os.path.join(zdir,'diags2D.*.meta'))
itrs=[]
for f in files:
    itrs.append(int(f.split('.')[-2]))

itrs = np.sort(itrs)

zdiag = os.path.join(zdir,'diags2D')
pdiag = os.path.join(pdir,'diags2D')
etanz=[]
etanp=[]
pbotz=[]
pbotp=[]
for itr in itrs:
    print('myiter: %10.10i'%int(itr))
    fld = rdmds(zdiag,itrs=itr,rec=[0,2])
    etanz.append((fld[0,:,:]*grac).sum())
    pbotz.append((fld[1,:,:]*grac).sum())
    fld = rdmds(pdiag,itrs=itr,rec=[0,2])
    etanp.append((fld[0,:,:]*grac).sum())
    pbotp.append((fld[1,:,:]*grac).sum())

xdays = np.array([refdate + datetime.timedelta(days=int(i/48)) for i in itrs])

# correct the jumps
etanz=correct_jumps(np.asarray(etanz))
pbotz=correct_jumps(np.asarray(pbotz))
etanp=correct_jumps(np.asarray(etanp))
pbotp=correct_jumps(np.asarray(pbotp))
#
etaz=correct_jumps(etaz)
etap=correct_jumps(etap)

timedayz = np.asarray(timesecz)/86400.
timedayp = np.asarray(timesecp)/86400.
zdays = np.array([refdate + datetime.timedelta(days=i) for i in timedayz])
pdays = np.array([refdate + datetime.timedelta(days=i) for i in timedayp])

fig, ax = plt.subplots(1,1)
# ax.plot(zdays,etaz/m2zUnits,label='etaz')
# ax.plot(pdays,etap/m2pUnits,label='etap')
ax.plot(xdays,etanz/m2zUnits,label='etaNz')
ax.plot(xdays,etanp/m2pUnits,label='etaNp')
ax.plot(xdays,(pbotz-pbotz[0])*rhoConst/gravity,label='pbotz')
ax.plot(xdays,(pbotp-pbotp[0])*rhoConst/gravity,label='pbotp')
ax.set_title('mass anomaly per area (kg$\,$m$^{-2}$)')
ax.grid()
ax.legend()

fig.show()

fig, ax = plt.subplots(1,1)
ax.plot(zdays,etaz/m2zUnits,label='etaz')
ax.plot(pdays,etap/m2pUnits,label='etap')
ax.set_title('mass anomaly per area (kg$\,$m$^{-2}$)')
ax.grid()
ax.legend()
fig.savefig('massanomaly')
