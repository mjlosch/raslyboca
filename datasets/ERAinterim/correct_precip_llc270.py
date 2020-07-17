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

rhoFresh = 1000.
# estimate of yearly imbalance of uncorrected model
fname = '/work/ollie/mlosch/raslyboca/llc270/spinup_era_interim/yearly_correction'
fname1 = '/work/ollie/mlosch/raslyboca/datasets/Eta/sea_surface_elevation.nc'
fname2 = '/work/ollie/mlosch/raslyboca/datasets/Eta/dangendorf-etal2019_ncc.txt'
with open(fname, 'rb') as fp:
    gmyearly_saved = pickle.load(fp)

# years = gmyearly_saved[0]
# gmyearly = np.asarray(gmyearly_saved[1])
years = gmyearly_saved[0]
gmyearly = np.asarray(gmyearly_saved[1])

# grid cell area of model
#rac = rdmds('/work/ollie/mlosch/raslyboca/llc270/spinup_era_interim/RAC')

# observed sea level rise
ds = Dataset(fname1,'r')
tm = ds['time'][:]
eta= ds['eta'][:]
tmi = np.zeros(tm.shape,dtype=np.int16)
for k in range(len(tm.data)):
    tmi[k] = int(tm[k])

# remove sea level anomaly of 1979-01-01
etac = eta-eta[0]

etam=[]
for k, year in enumerate(years):
        etam.append(etac[tmi==year].mean())
        print(etac[tmi==year].mean())


tt=np.loadtxt(fname2,skiprows=1)
# remove sea level anomaly of 1979-01
ttref=tt[tt[:,0]==1979,1]
ttc=tt[:,1]-ttref

dt = 86400.*365 * np.ones((len(years),))
# leap years:
import calendar as cal
for k,year in enumerate(years):
    if cal.isleap(year): dt[k]=86400.*366

fig,ax=plt.subplots(2,1,sharex=True)
ax[0].plot(years,-np.cumsum(gmyearly)*dt/rhoFresh,label='llc270')
ax[0].plot(np.asarray(years)+.5,etam,label='yearly average of CSIRO')
ax[0].plot(ds['time'][:],etac,label='CSIRO sea level estimate')
#ax[0].plot(tt[:,0],ttc*1e-3,label='Dangendorf et al. (2019)')
ax[0].set_ylabel('(m)')
ax[0].grid()
ax[0].legend()

etam.insert(0,0.)
etad = np.diff(np.asarray(etam))/dt
etad = np.where(np.isnan(etad),np.nanmean(etad),etad)

ax[1].bar(years,gmyearly/rhoFresh * 1e9,label='model bias')
ax[1].bar(years,etad*1e9,alpha=0.7,label='implied observed bias')
ax[1].set_ylabel('(10$^{-9}$ m s$^{-1}$)')
ax[1].grid()
ax[1].legend()
fig.show()

yearly_precip_corr = gmyearly/rhoFresh + etad

def readfield(fname,dims):
    import sys
    """Call signatures::

    readfield(filename, dims, numpy.datatype)

    Read unblocked binary data with dimentions "dims".
    """

    try:
        fid = open(fname,"rb")
    except:
        sys.exit( fname+": no such file or directory")
    else:
        v   = np.fromfile(fid,'>f4')
        fid.close()

        if   len(v) == np.prod(dims):     v = v.reshape(dims)
        elif len(v) == np.prod(dims[1:]): v = v.reshape(dims[1:])
        else:
            errstr = (  "dimensions do not match: \n len(data) = " + str(len(v))
                        + ", but prod(dims) = " + str(np.prod(dims)) )
            raise RuntimeError(errstr)

    return v

def writefield(fname,data):
    import sys
    print('writing '+fname)
    if False:
        pass
    else:
        fid = open(fname,"wb")
        data.astype('>f4').tofile(fid)
        fid.close()

for k, year in enumerate(years):
#for k, year in enumerate(range(1979,1981)):
    # load data
    kt = 1460
    if cal.isleap(year): kt = 1464
    if year==2019: kt = 976
    fld = readfield('tp_ERAi_6hourly_'+str(year),[kt,241,480])
    writefield('tp_ERAi_corr_llc270_6hourly_'+str(year),
               fld+yearly_precip_corr[k])
