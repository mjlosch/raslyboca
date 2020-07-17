# convert ERAiterim data from (Claudia Wekerle's) netcdf files to ieee-be
# compute specific humidity from dew point temperature
# unfortunately precipitation and downward radiation are only available
# as daily averages, I don't know why

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os



ncdir = '/work/ollie/clidyn/forcing/erai'
input_variables = ['precip', 'tdew', 'rad', 't_02', 'u_10', 'v_10']

years = range(2017, 2019)

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

def writefield(fname,arr):
    import sys
    print('writing '+fname)
    if True:
        pass
    else:
        data = arr.data
        fid = open(fname,"wb")
        data.astype('>f4').tofile(fid)
        fid.close()

for y in years:
    print(y)
    for invar in input_variables:
        fname=os.path.join(ncdir,'erai.%s.%i.nc'%(invar,y))
        if os.path.isfile(fname): print(fname+' exists and is a file')
        ds = Dataset(fname,'r')
        # flip direction of y-axis because the ERA convention is to have
        # indices (0,0) at the top left corner of the field
        if invar=='t_02':
            outfld = ds['T_2_MOD'][:,::-1,:]
            bfile = 't2m_ERAi_6hourly_'+str(y)
            check_flds(outfld)
            writefield(bfile,outfld)
        elif invar=='u_10':
            outfld = ds['U_10_MOD'][:,::-1,:]
            bfile = 'u10_ERAi_6hourly_'+str(y)
            check_flds(outfld)
            writefield(bfile,outfld)
        elif invar=='v_10':
            outfld = ds['V_10_MOD'][:,::-1,:]
            bfile = 'v10_ERAi_6hourly_'+str(y)
            check_flds(outfld)
            writefield(bfile,outfld)
        elif invar=='tdew':
            tdew = ds['d2m'][:,::-1,:]
            slpname=os.path.join(ncdir,'erai.slp.%i.nc'%(y))
            if os.path.isfile(slpname):
                print(slpname+' exists and is a file')
                dslp = Dataset(slpname,'r')
                slp = dslp['SLP'][:,::-1,:]
            else:
                print(slpname+' does not exist, using slp = 1 bar')
                slp = np.zeros(tdew.shape,dtype='float32') + 1e5

            outfld = specific_humidity(tdew, slp)
            bfile='q_ERAi_6hourly_'+str(y)
            check_flds(outfld)
            writefield(bfile,outfld)
        elif invar=='precip':
            rain = ds['RAIN'][:,::-1,:]
            snow = ds['SNOW'][:,::-1,:]
            outfld = rain+snow
            bfile='tp_ERAi_6hourly_'+str(y)
            check_flds(outfld)
            writefield(bfile,outfld)


        if invar=='rad':
            swdw = ds['SWDW'][:,::-1,:]
            bfile='ssrd_ERAi_6hourly_'+str(y)
            check_flds(swdw)
            writefield(bfile,swdw)
            lwdw = ds['LWDW'][:,::-1,:]
            bfile='strd_ERAi_6hourly_'+str(y)
            check_flds(lwdw)
            writefield(bfile,lwdw)
