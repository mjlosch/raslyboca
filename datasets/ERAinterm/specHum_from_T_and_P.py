# /usr/bin/python
# coding: utf-8

# Python program to calculate specific humidity from dew-point Temperature and surface pressure 
# Theory: http://www.ecmwf.int/en/does-era-40-dataset-contain-near-surface-humidity-data &
# http://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf equations 7.4 and 7.5

# Elena Gerwing, Mar 2017                                                                                                
# Import modules #{{{1
import sys, os, glob

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.ticker as mtick
#from matplotlib import colors
from netCDF4 import Dataset
#import datetime
import argparse #}}}1

# Read command line arguments #{{{1
if __name__=='__main__':
     
    parser = argparse.ArgumentParser(description='Calculate specific humidity', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input NetCDF file containing dew-point temperature and surface pressure, output is written in the same file')                    
#    parser.add_argument('--labels', nargs='+', default=None, help='enter labels displayed in the figure,     if set same number as input files must be defined')
#    parser.add_argument('--startyears',type=int,nargs='*',default=[1979],help='starting year of simulations [yyyy]; if specified, one value or some number as input files required')
    args = parser.parse_args()

#    print 'Years', args.startyears
    print 'Files:',args.infile #}}}1

# Functions #{{{1

def specific_humidity(data,a1,a3,a4,Rd,Rv,T0):
    '''calculate specific humidity from 
    dew-point temperature and surface pressure
    following ECMWF Equation:
    http://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf,equations 7.4/7.5
    data: data from NetCDF file
    a1,a3,a4: Parameters according to Buck 1981
    Rd,Rv: gas constants of dry air and water vapor
    T0: reference temperature
    '''
    
    Td = data.variables['d2m'][:]   # (K) dew-point temperature 
    P = data.variables['sp'][:]     # (Pa) surface pressure

    R = Rd/Rv
    e_sat   = a1*np.exp( a3*( (Td -T0)/(Td -a4) ) ) # (Pa) saturation water vapor pressure

    q_sat   = (R*e_sat) / (P - (1.0-R)*e_sat)

    return q_sat

def specific_humidity_2(data,Rv,T0):
    '''calculate specific humidity from 
    dew-point temperature and surface pressure
    following:http://www.faqs.org/faqs/meteorology/temp-dewpoint/
    data: data from NetCDF file
    Rv: gas constants of dry air and water vapor
    T0: reference temperature
    '''
    Pa2hPa = 1e-2;
    es0 = 6.11   # hPa, reference saturation vapor pressure (es at 0 deg C)
    lv  = 2.5e6  # J/kg, latent heat of vaporization of water 
   
    Td = data.variables['d2m'][:]   # (K) dew-point temperature 
    P = data.variables['sp'][:]     # (Pa) surface pressure

    e_sat = es0 * np.exp( lv/Rv * (1/T0 - 1/Td) )
    q_sat   = 0.622*e_sat/(P*Pa2hPa + 0.378*e_sat)

    return q_sat


def write_q_to_file(data,specHum):
    '''include specific humidity in NetCDF file'''
    print 'write q to file'
    q = data.createVariable('q','f8',('time','latitude','longitude')) # 
    q.standard_name = 'specific_humidity'
    q.long_name = 'saturation specific humidity over water'
    q.units = '1'      
    
    q[:] = specHum
    
#    print 'IN WRITE ROUTINE'
#    print 'MAX',np.max(q)
#    print 'MIN',np.min(q1)
#    print 'MEAN',np.mean(q1)

#}}}1

# Set constants #{{{1

# Parameters according to Buck 1981 for saturation over water
a1      = 611.21    # (Pa) Pascal
a3      = 17.502    #
a4      = 32.19     # (K) Kelvin

# Gas constants
Rd      = 287.06    # (J/(kg*K)) dry air
Rv      = 461.53    # (J/(kg*K)) water vapor

T0      = 273.16    # (K) reference temperature

#}}}1


if __name__=='__main__':
    data =  Dataset(args.infile,'r+')
    q1 = specific_humidity(data,a1,a3,a4,Rd,Rv,T0)
#    q2 = specific_humidity_2(data,Rv,T0)
#    print('q1',q1[0,100,200:206],'q2',q2[0,100,200:206])
#    q1 = np.short(q1)

#    print q1.dtype
#    print 'MAX',np.max(q1)
#    print 'MIN',np.min(q1)
#    print 'MEAN',np.mean(q1)

    write_q_to_file(data,q1)
    data.close() 

#    output = Dataset(args.infile,'r')
#    q2 = output.variables['q'][:] 
#    print "Read from output file"
#    print q2.dtype
#    print 'MAX',np.max(q2)
#    print 'MIN',np.min(q2)
#    print 'MEAN',np.mean(q2)
#    output.close()
