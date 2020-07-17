import numpy as np
from MITgcmutils import rdmds
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import datetime
import os, sys
from getopt import gnu_getopt as getopt
from myutils import *
import pickle

# estimate of yearly imbalance of uncorrected model
fname = '/work/ollie/mlosch/raslyboca/llc270/spinup_era_interim/yearly_correction'
with open(fname, 'rb') as fp:
    gmyearly_saved = pickle.load(fp)

years = gmyearly_saved[0]
gmyearly = np.asarray(gmyearly_saved[1])

# grid cell area of model
#rac = rdmds('/work/ollie/mlosch/raslyboca/llc270/spinup_era_interim/RAC')

# observed sea level rise
ds = Dataset('/work/ollie/mlosch/raslyboca/datasets/Eta/sea_surface_elevation.nc')
tm = ds['time'][:]
eta= ds['eta'][:]
tmi = np.zeros(tm.shape,dtype=int16)
for k in range(len(tm.data)):
    tmi[k] = int(tm[k])

etam=[]
for k, year in enumerate(years):
        etam.append(eta[tmi==year].mean())
        print(eta[tmi==year].mean())


fig,ax=plt.subplots(1,1)
ax.plot(ds['time'][:],ds['eta'][:])
ax.plot(years,-cumsum(gmyearly)*86400*365/1000)
ax.plot(asarray(years)+.5,etam)
ax.grid()
fig.show()
