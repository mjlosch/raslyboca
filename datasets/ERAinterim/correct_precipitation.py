#! /usr/bin/python
# coding: utf-8

# Python program to correct the precipitation forcing
# by subtractiong the trend of the mean sea surface elevation
# and adding the mean sea surface elevation from observations
#
# Elena Gerwing, May 2017
# Import modules #{{{1
from __future__ import division    #floating point division
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import xarray as xr
from scipy.io import FortranFile
import matplotlib.pyplot as plt
from MITgcmutils import rdmds,wrmds,llc

import argparse #}}}1

# Read command line arguments #{{{1
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Correct the precipitation forcing in order to let the model fit the observed sea level change', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input NetCDF file')
    parser.add_argument('outfile', help='Enter name for binary output file (to append to parameter name)')
    parser.add_argument('--convert',action='store_true',default=False,help='If set, converts units of precipitation tp, long- and shortwave radiation strd and ssrd by dividing by time interval (43200 s)')
    parser.add_argument('--change_shape',action='store_true',default=False,help='Change shape from (730,241,480) to (730,240,480) by cutting first line in y-direction')
    parser.add_argument('--flip_y',action='store_true',default=False,help='Flip the array in y-direction in order to start with -90 degree')
    args = parser.parse_args()

    print ('Prepare_forcing, read file :',args.infile) #}}}1

# Functions #{{{1

def writefield(fname,data):  #{{{2
    '''Write data to binary file
       change from little to big endian
       fname: filename to which to save the data
       data: numpy-array to write to file
    '''

    import sys
    print ('write to file: '+fname)
    if sys.byteorder == 'little':
#        print ('swap bytes')
        data.byteswap(True)
    fid = open(fname,"wb")
    data = data.data       # in case array is masked, extract the data
    data.tofile(fid)
    fid.close()

    if sys.byteorder == 'little': data.byteswap(True)      # Swap bag to former endianness
#}}}2

def calc_correction_term(model,obs,timestep=(30.41*4*21600.)): #{{{2

    diff = obs - model

    # since the sea level values are accumulated over the time period
    # I need to substract the preceding value from every entry
    tmp = np.roll(diff,1) # shift all values by one entry
    tmp[0] = 0. # from the first entry nothing needs to be subtracted

    correction = diff - tmp

    #plt.figure(2)
    #plt.plot(np.arange(0,len(model)),diff,".")
    #plt.plot(np.arange(0,len(model)),correction,".")
#   # plt.plot(np.arange(0,len(obs)),obs)
    #plt.show()


    # devide by the timestep in order to get the correction term per second
    correction = correction/timestep


    # np.save('precip_correction2.npy',correction)
    return correction

    #}}}2

def correct_yearly(data,parameters,year,correction,outputfile,convert=True,change_shape=False,flip_y=True): #{{{2
    '''Read precipitation from NetCDF file, correct it with observations
       (subtract difference between model and observational yearly means)
       and write to binary file
       data: NetCDF dataset (data = Dataset('file.nc','r'))
       parameter: list of parameters
       year: current year of file
       correction: numpy-array with the corrected mean values
       outputfile: name of outputfile (appended to parameter name), i.e '_ERAinterim_2011'
       convert: if true, converts units of tp, ssrd, strd by dividing by time interval
       change_shape: if true, change shape from (730,241,480) to (730,240,480)
    '''

    for param in parameters:
        print ('=== YEARLY CORRECTION ===')
        print ('Processing parameter:',param)

        outfile = str(outputfile)

        var = data.variables[param][:]


        if change_shape == True:
#            print ('transform shape')
#           var = var[:,:-1,:]
            var = var[:,1:,:]

        if flip_y == True:
            var = var[:,::-1,:]

        var = np.float32(var)

        if convert == True:
            if param == 'tp' or param == 'ssrd' or param == 'strd':
                # since fields of step 12 are accumulated over 12 hours
                #but I need an avarage value 6-hourly
                # I need to subtract the first 6 hours by subtracting the field from step 6
                # (see README or ECMWF site for further explanation)
                var[1::2,:,:] = var[1::2,:,:] - var[0::2,:,:]
                var = np.maximum(var,0.)
                # divide by the time interval to get mean value
                var = var/21600.
                print ('Converted unit:',param)


        i = (year - 1979)
        if year > 2010: # IAF runoff ends 2011
            i = (2010 - 1979)


        var = var + correction[i]
        # Account for negative precipitation
        sum_neg = np.sum(var[var < 0.0])
        num = np.sum(var < 0.0)
        neg_correction = sum_neg/(var.size)
        var = var + neg_correction
        var[var <= 0.0] = 0.0

        writefield(outfile,var)
        print("============= Finished year",year," ===========")

    data.close()
#}}}2

#}}}1


if __name__=='__main__':
    #{{{1 Calculate correction term

#ML    obs_ds = Dataset('/work/ollie/egerwing/datasets/Eta/sea_surface_elevation.nc')
    obs_ds = Dataset('/work/ollie/mlosch/raslyboca/datasets/Eta/sea_surface_elevation.nc')
#ML    model_ds = Dataset('/work/ollie/egerwing/MITgcm/llc90/restart_uncorrected_nyf/Eta.0000280512.t001.nc')
    model_ds = Dataset('/work/ollie/mlosch/raslyboca/MITgcm_elena/llc90/uncorrected_nyf/Eta.0000280512.t001.nc')

    model_m = model_ds.variables['ETAN_ave'][:420,0,0] # only data till 2013 cause the observational data ends then
    obs_m = obs_ds.variables['eta_shifted'][:]  #len(model_m)]

    correction = calc_correction_term(model_m,obs_m)
    # calc yearly correction from monthly
    correction_yr = np.zeros(int(np.rint(len(correction)/12.)))
    for year in np.arange(0,int(np.rint(len(correction)/12.))):
        correction_yr[year] = np.mean(correction[year*12:(year+1)*12])

    #}}}1

    data =  Dataset(args.infile,'r')              # read from NetCDF file

    year = int(args.infile.split('_')[-1].split('.')[-2])   # get year out of file name (file needs to end with '_YEAR.nc')
    print ("%%%%%%%%%%%%%%%% Starting processing of year",year," %%%%%%%%%%%%%%%%%%%%")

    parameter = ['tp']

    correct_yearly(data,parameter,year,correction_yr,args.outfile,convert=args.convert,change_shape=args.change_shape,flip_y=args.flip_y)


#    data.close()
