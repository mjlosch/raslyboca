#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
######################## -*- coding: utf-8 -*-
"""Usage: extract_rm_mean.py INPUTFILE(S)
"""
import numpy as np
from MITgcmutils import rdmds, rdmnc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import datetime
import os, sys
from getopt import gnu_getopt as getopt
from myutils import *

# parse command-line arguments
try:
    optlist,args = getopt(sys.argv[1:], ':', ['verbose'])
    assert len(args) > 0
except (AssertionError):
    sys.exit(__doc__)

from glob import glob
files=[]
if len(args)<2:
    for infile in glob(args[0]):
        files.append(infile)
else:
    files=args

diagfiles=[]
bdir=os.getcwd()
for infile in glob(os.path.join(bdir,'dynStDiag.*.nc')):
    diagfiles.append(infile)

def get_parms (fname):
    rhoIce = 910.
    with open(fname) as f:
        for line in f:
            if 'startDate_1' in line:
                startDate = line.split('=')[-1]
            elif '/* Model clock timestep ( s ) */' in line:
                ll = next(f).strip().split()
                deltaT = float(ll[-1].replace('E','e'))
            elif '/* Reference density (Boussinesq)  ( kg/m^3 ) */' in line:
                ll = next(f).strip().split()
                rhoConst = float(ll[-1].replace('E','e'))
            elif '/* Fresh-water reference density ( kg/m^3 ) */' in line:
                ll = next(f).strip().split()
                rhoFresh = float(ll[-1].replace('E','e'))
            elif '/* density of sea ice (kg/m3) */' in line:
                ll = next(f).strip().split()
                rhoIce = float(ll[-1].replace('E','e'))


    return rhoFresh, rhoIce/rhoConst, deltaT, startDate

def get_output (fnames, mystring):
    """parse fname and get some numbers out"""
    timev = []
    myvar = []
    for fname in fnames:
        try:
            f = open(fname)
        except:
            print(fname + " does not exist, continuing")
        else:
            for line in f:
                if mystring in line:
                    ll = line.split()
                    timev.append(float(ll[-2].replace('E','e')))
                    myvar.append(float(ll[-1].replace('E','e')))

            f.close()

    # reverse order
    timevs=np.asarray(timev[::-1])
    myvars=np.asarray(myvar[::-1])
    # This sorts again in ascending order and returns the index of
    # the first occurrence of duplicates. Because we have reverted the order
    # before, in this way we use the values at the beginning of a pickup run
    # rather than the overlapping values of the previous (potentially crashed)
    # run
    timevs, isort = np.unique(timevs,return_index=True)
    myvars=myvars[isort]

    return timevs, myvars
# done

def get_diagst(fnames):
    """parse fname and get some numbers out"""
    timev = []
    etav  = []
    heffv = []
    for fname in fnames:
        diags = rdmnc(fname)
        timev = np.append(timev,diags['T'][:])
        etav  = np.append(etav, diags['ETAN_ave'][:,0,0])
        heffv = np.append(heffv,rhoIce*diags['SIheff_ave'][:,0,0])

    # reverse order
    timevs=np.asarray(timev[::-1])
    eta   =np.asarray(etav[::-1])
    heff  =np.asarray(heffv[::-1])
    # This sorts again in ascending order and returns the index of
    # the first occurrence of duplicates. Because we have reverted the order
    # before, in this way we use the values at the beginning of a pickup run
    # rather than the overlapping values of the previous (potentially crashed)
    # run
    timevs, isort = np.unique(timevs,return_index=True)

    return timevs, eta[isort], heff[isort]
# done

timesec, gm = get_output(files,'rm Global mean of SIatmFW')
rhoFresh, rhoIce, deltaT, startdate = get_parms(files[0])
#
timediags, eta, hff = get_diagst(diagfiles)

refdate = datetime.datetime(int(startdate[:4]),
                            int(startdate[4:6]),int(startdate[6:8]),0,0)
timeday = np.asarray(timesec)/86400.
xdays = np.array([refdate + datetime.timedelta(days=i) for i in timeday])
timeday=timediags/86400.
ddays = np.array([refdate + datetime.timedelta(days=i) for i in timeday])

fig, ax = plt.subplots(1,1)
# variable SIatmFW is positive upwards, leading sea level decrease, hence the sign
ax.plot(xdays,-np.cumsum(gm)*deltaT/rhoFresh,label='- cumsum(global mean SIatmFW)')

hff = rhoIce*hff
ax.plot(ddays,eta,label='ETAN_ave')
ax.plot(ddays,hff,label='SIheff_ave')
ax.plot(ddays,eta+hff,label='sum')

ax.grid()
ax.legend()
fig.show()

if len(gm)>0:
    # compute yearly correction from gm ([gm]=kg/m^2/s)
    years=range(xdays[0].year,xdays[-1].year+1)
    gmyearly=[]
    for y in years:
        gmyearly.append(0.)
        kyear = 0
        for k, day in enumerate(xdays):
            if day.year == y:
                kyear = kyear+1
                gmyearly[-1]=gmyearly[-1]+gm[k]

        if kyear>0: gmyearly[-1]=gmyearly[-1]/kyear

    figb, axb = plt.subplots(1,1)
    axb.bar(years,np.asarray(gmyearly)*1e6)
    axb.grid()
    axb.set_title('annual correction (mg m$^{-2}$s$^{-1}$ $\Leftrightarrow$ 10$^{-9}$ m s$^{-1}$))')
    figb.show()

    fname = 'yearly_correction'
    with open(fname+'.txt', 'w') as f:
        for item in gmyearly:
            f.write("%s\n" % item)

    import pickle
    with open(fname, 'wb') as fp:
        pickle.dump([years,gmyearly], fp)

    with open(fname, 'rb') as fp:
        gmyearly_saved = pickle.load(fp)
