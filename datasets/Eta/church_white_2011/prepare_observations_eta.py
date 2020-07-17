#! /usr/bin/python



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from netCDF4 import Dataset



data = np.fromfile('CSIRO_Recons_gmsl_mo_2015.txt',sep=' ')
print data.shape, data[0:7]

data = data.reshape(1608,3)
print data.shape,data[0:7,0]

# use data from 1979 on 
time = data[1188:,0]
eta = data[1188:,1]
error= data[1188:,2]

# convert units from mm to m
eta = eta/1000.
error = error/1000.



data = Dataset('sea_surface_elevation.nc','w')
# create dimensions
t = data.createDimension('time',420)

var = data.createVariable('eta','f8',('time')) # 
var.standard_name = 'sea_surface_height'
var.long_name = 'sea surface height above/below sea level of 1990'
var.units = 'm'     
var.axis = 'T'

var[:] = eta[:]

err = data.createVariable('eta_err','f8',('time')) # 
err.standard_name = 'sea_surface_height_error'
err.long_name = 'measurement uncertainty of sea surface height'
err.units = 'm' 
err.axis = 'T'

err[:] = error[:]

t = data.createVariable('time','f8',('time')) # 
t.standard_name = 'time'
t.long_name = 'time in years'
t.units = 'years'     
t.axis = 'T'

t[:] = time[:]





data.close()
