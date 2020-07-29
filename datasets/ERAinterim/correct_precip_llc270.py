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
rhoConst = 1035.
# estimate of yearly imbalance of uncorrected model
d1 = '/work/ollie/mlosch/raslyboca/llc270/spinup_era_interim'
d2 = '/work/ollie/mlosch/raslyboca/datasets/Eta'
#d1 = '/Users/mlosch/Downloads/Eta'; d2 = d1
fname = os.path.join(d1,'yearly_correction')
fname1 = os.path.join(d2,'sea_surface_elevation.nc')
fname2 = os.path.join(d2,'dangendorf-etal2019_ncc.txt')
with open(fname, 'rb') as fp:
    gmyearly_saved = pickle.load(fp)

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
ax[1].bar(years,etad*1e9*rhoConst/rhoFresh,alpha=0.7,label='implied observed bias')
#ax[1].bar(years,etad*1e9,alpha=0.5,label='wrong implied observed bias')
ax[1].set_ylabel('(10$^{-9}$ m s$^{-1}$)')
ax[1].grid()
ax[1].legend()
fig.show()

yearly_precip_corr = gmyearly/rhoFresh + etad*rhoConst/rhoFresh

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

# erai grid
nx,ny,dx0=480,241,0.75
lon = np.linspace(0.,360.,nx+1)
lat = np.linspace(-90.,90.,ny)
dlat = dx0*np.pi/180.*np.ones(nx,)
dlon = dx0*np.pi/180.*np.cos(lat*np.pi/180.)
# dlon[0],dlon[-1]=0.5*dlon[1],0.5*dlon[-2]
# extrapolation to find dlon at poles
dlon[0],dlon[-1]=2.*dlon[2]-dlon[1],2.*dlon[-2]-dlon[-1]
# area of pole element is pi*(0.5*dlat)**2, and circumference nx*dlon = 2*pi*0.5*dlat
dlon[0], dlon[-1] = np.pi*dlat[0]/nx, np.pi*dlat[-1]/nx
dy,dx=np.meshgrid(dlat,dlon)
rac=dx*dy
lon0,lat0 = np.meshgrid(lon[:-1],lat)
lon1=np.hstack((lon0[:,nx//2:]-360,lon0[:,:nx//2+1]))
lat1=np.hstack((lat0[:,nx//2:],lat0[:,:nx//2+1]))

# LLC270 grid
grid = '/home/ollie/mlosch/MITgcm/MITgcm/llc270/grid'
hf = rdmds(os.path.join(grid,'hFacC'),lev=[0])
xg = rdmds(os.path.join(grid,'XG'))
yg = rdmds(os.path.join(grid,'YG'))
xc = rdmds(os.path.join(grid,'XC'))
yc = rdmds(os.path.join(grid,'YC'))
ra = rdmds(os.path.join(grid,'RAC'))
# interpolate land sea mask hf to erai grid
# add wrap around

eramsk = np.zeros(lon1.shape)
for ix in range(nx+1):
    print(lon1[0,ix])
    for iy in range(ny):
        index = np.logical_and(
            np.logical_and(xc <= lon1[iy,ix]+.5*dx0,xc >= lon1[iy,ix]-.5*dx0),
            np.logical_and(yc <= lat1[iy,ix]+.5*dx0,yc >= lat1[iy,ix]-.5*dx0)
            )
        if index.sum() > 0:
            eramsk[iy,ix]= hf[index].mean()
        elif lat1[iy,ix]>80.:
            eramsk[iy,ix]= 1.
        else:
            eramsk[iy,ix]= -999.999

emsk = np.where(eramsk>.5,1.,0.)
# reshape to lon = [0-360]
emsk = np.hstack((emsk[:,nx//2:-1],emsk[:,:nx//2]))
racmskd = rac*emsk

if True:
    for k, year in enumerate(years):
    #for k, year in enumerate(range(1979,1981)):
        # load data
        kt = 1460
        if cal.isleap(year): kt = 1464
        if year==2019: kt = 976
        fld = readfield('tp_ERAi_6hourly_'+str(year),[kt,ny,nx])
        fldc = fld + yearly_precip_corr[k]
        if yearly_precip_corr[k] < 0.:
            # correct for negative precip:
            ract=np.tile(rac.reshape((1,rac.shape[0],rac.shape[1])),(kt,1,1))
            mskt=np.tile(emsk.reshape((1,emsk.shape[0],emsk.shape[1])),(kt,1,1))
            ntot = len(fldc.ravel())
            # mask the fld
            fldm = fldc*ract*mskt
            # 1. find all negative entries
            ineg = fldm<0
            nneg = ract[ineg].sum()
            sumneg = np.sum(fldm[ineg])
            # 2. find all entries over the ocean at least 3 times larger then
            #    the mean of all negative entries and add the mean of
            #    the negatives
            ipos = fldc*mskt > -3*sumneg/nneg
            npos = ract[ipos].sum()
            fldc[fldc<0.] = 0.
            if npos > 0.: fldc = np.where(ipos,fldc + sumneg/npos ,fldc)


        writefield('tp_ERAi_corr_llc270_6hourly_'+str(year),fldc)
