# read "raw" netcdf files downloaded from ECMWF/MARS
# these are "fluxes" that need to be changed to rates

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os

# there 4 files with data between 2017-01-01 to 2019-08-31 (3892 records)
# ds0 = Dataset('dlwrd.nc','r')
# ds1 = Dataset('dswrd.nc','r')
# ds2 = Dataset('snowfall.nc','r')
# ds3 = Dataset('total-pre.nc','r')

def check_flds(fld):

    print(fld.shape)
    print(fld.dtype)
    fmt='mean: %12.6e std: %12.6e min: %12.6e max: %12.6e'
    print(fmt%(fld.mean(),fld.std(),fld.min(),fld.max()))

def writefield(fname,data):
    import sys
    print('writing '+fname)
    if False:
        pass
    else:
        fid = open(fname,"wb")
        data.astype('>f4').tofile(fid)
        fid.close()

def convertflux2rate(fldf):
    recip_dt = 1./(6*3600)
    fldr = np.copy(fldf)
    # substract accumulation after 6hrs from accumulation after 12hrs to
    # get 6hrly values
    fldr[1::2,:,:]=fldf[1::2,:,:]-fldf[0::2,:,:]
    # en passant correct for negative values
    return np.where(fldr<0.,0,fldr)*recip_dt

def savefield(fname,fldr):
    # convert to single precision and save
    writefield(fname+'_2017',fldr[:1460,:,:])
    writefield(fname+'_2018',fldr[1460:1460*2,:,:])
    writefield(fname+'_2019',fldr[1460*2:,:,:])

    return

datadir = '/work/ollie/cwekerle/erai/6_12/'
# radiative fluxes are simpleb
vars = ['strd','ssrd']
for v in vars:
    if v=='strd': ds = Dataset(os.path.join(datadir,'dlwrd.nc'),'r')
    if v=='ssrd': ds = Dataset(os.path.join(datadir,'dswrd.nc'),'r')
    fldf = ds[v][:,::-1,:]
    fldr = convertflux2rate(fldf)

    print(v)
    check_flds(fldr)

    savefield(v+'_ERAi_6hourly',fldr)

dssf = Dataset(os.path.join(datadir,'snowfall.nc'),'r')
dstp = Dataset(os.path.join(datadir,'total-pre.nc'),'r')

# add snow to total precipitation
fldf = dssf['sf'][:,::-1,:] + dstp['tp'][:,::-1,:]
fldr = convertflux2rate(fldf)

print('precip')
check_flds(fldr)

savefield('tp_ERAi_6hourly',fldr)
