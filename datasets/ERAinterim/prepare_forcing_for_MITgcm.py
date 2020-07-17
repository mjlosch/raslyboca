#! /usr/bin/python
# coding: utf-8

# Python program to convert (ERA-Interim) NetCDF files into binary files (IEEE?)                                      
# and eventually convert units of short- and longwave radiation from J/m^2 to W/m^2 and precipitation from m to m/s
#
# Elena Gerwing, Mar 2017                                                                                                
# Import modules #{{{1
from __future__ import division    #floating point division
import numpy as np
from netCDF4 import Dataset
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import argparse #}}}1

# Read command line arguments #{{{1
if __name__=='__main__':
     
    parser = argparse.ArgumentParser(description='Prepare forcing input data for use in MITgcm: write data from netcdf to binary, eventually convert shape and convert units of tp, ssrd, strd', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input NetCDF file')                    
    parser.add_argument('outfile', help='Enter name for binary output file (to append to parameter name)')  
    parser.add_argument('--convert',action='store_true',default=False,help='If set, converts units of precipitation tp, long- and shortwave radiation strd and ssrd by dividing by time interval (43200 s)')
    parser.add_argument('--change_shape',action='store_true',default=False,help='Change shape from (730,241,480) to (730,240,480) by cutting first line in y-direction')
    parser.add_argument('--flip_y',action='store_true',default=False,help='Flip the array in y-direction in order to start with -90 degree')
    args = parser.parse_args()
    
    print 'Prepare_forcing, read file :',args.infile #}}}1

# Functions #{{{1

def writefield(fname,data):  #{{{2
    '''Write data to binary file
       change from little to big endian
       fname: filename to which to save the data
       data: numpy-array to write to file
    '''

    import sys
    print 'write to file: '+fname
    if sys.byteorder == 'little': 
#        print 'swap bytes'
        data.byteswap(True)
    fid = open(fname,"wb")
    data.tofile(fid)
    fid.close()

    if sys.byteorder == 'little': data.byteswap(True)      # Swap bag to former endianness
#}}}2

def loop_for_conversion(data,parameters,outputfile,convert=True,change_shape=False,flip_y=True): #{{{2
    '''Read parameters from NetCDF file and write to binary file
       data: NetCDF dataset (data = Dataset('file.nc','r'))
       parameter: list of parameters
       outputfile: name of outputfile (appended to parameter name), i.e '_ERAinterim_2011'
       convert: if true, converts units of tp, ssrd, strd by dividing by time interval
       change_shape: if true, change shape from (730,241,480) to (730,240,480)
    '''

    for param in parameters:
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print 'Processing parameter:',param

        outfile = str(param) + str(outputfile)

        var = data.variables[param][:]

        
        if change_shape == True:
#            print 'transform shape'
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
                print 'Converted unit:',param


        writefield(outfile,var)


    data.close()
#}}}2

#}}}1


if __name__=='__main__':

    data =  Dataset(args.infile,'r')              # read from NetCDF file

    parameters = ['q','ssrd','strd','u10','v10','t2m']

    loop_for_conversion(data,parameters,args.outfile,convert=args.convert,change_shape=args.change_shape,flip_y=args.flip_y)



    
