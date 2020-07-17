# convert ERAiterim data from (Claudia Wekerle's) netcdf files to ieee-be
# compute specific humidity from dew point temperature
# only for wind, temperature, specific humidity (requires sea level pressure)
# and only for 2017, 2018, 2019 (up to Aug31, when ERAinterim ends)

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os



ncdir = '/work/ollie/clidyn/forcing/erai/raw'
input_variables = ['tdew', 't2m', 'u10', 'v10']
#input_variables = ['t2m']

def check_flds(fld):

    print(fld.shape)
    print(fld.dtype)
    fmt='mean: %12.6e std: %12.6e min: %12.6e max: %12.6e'
    print(fmt%(fld.mean(),fld.std(),fld.min(),fld.max()))

def specific_humidity(Td, P):
    '''calculate specific humidity from
    dew-point temperature and surface pressure
    following ECMWF Equation:
    http://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf,equations 7.4/7.5
    data: data from NetCDF file
    a1,a3,a4: Parameters according to Buck 1981
    Rd,Rv: gas constants of dry air and water vapor
    T0: reference temperature
    '''
    # Parameters according to Buck 1981 for saturation over water
    a1      = 611.21    # (Pa) Pascal
    a3      = 17.502    #
    a4      = 32.19     # (K) Kelvin

    # Gas constants
    Rd      = 287.06    # (J/(kg*K)) dry air
    Rv      = 461.53    # (J/(kg*K)) water vapor

    T0      = 273.16    # (K) reference temperature

    R = Rd/Rv
    # (Pa) saturation water vapor pressure
    e_sat   = a1*np.exp( a3*( (Td -T0)/(Td -a4) ) )

    q_sat   = (R*e_sat) / (P - (1.0-R)*e_sat)

    return q_sat

def writefield(fname,data):
    import sys
    print('writing '+fname)
    if False:
        pass
    else:
        fid = open(fname,"wb")
        data.astype('>f4').tofile(fid)
        fid.close()

def savefield(fname,fldr):
    # convert to single precision and save
    writefield(fname+'_2017',fldr[:1460,:,:])
    writefield(fname+'_2018',fldr[1460:1460*2,:,:])
    writefield(fname+'_2019',fldr[1460*2:,:,:])

    return


for v in input_variables:
    ds = Dataset(os.path.join(ncdir,v+'.nc'),'r')

    if v=='tdew':
        vnew='q'
        dslp = Dataset(os.path.join(ncdir,'slp.nc'),'r')
        tdew = ds['d2m'][:,::-1,:]
        slp = dslp['msl'][:,::-1,:]
        fld = specific_humidity(tdew, slp)
    else:
        fld = ds[v][:,::-1,:]
        vnew=v

    print(vnew)
    check_flds(fld)
    savefield(vnew+'_ERAi_6hourly',fld.data)
