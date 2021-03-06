#! /usr/bin/python
# coding: utf-8

#Increase runoff in the Arctic (>= 55 degree N) by 30 % 
#At the moment only working for llc90 cause the dimensions are set

# Elena Gerwing, September 2018                                                                                               
# Import modules #{{{1
import sys, os, glob

import numpy as np
from scipy import interpolate
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
#import matplotlib.ticker as mtick
from matplotlib import colors
from netCDF4 import Dataset
#import datetime
import argparse

from MITgcmutils import rdmds,wrmds,llc
from myutils import sq,pcol
#}}}1

# Read command line arguments #{{{1
if __name__=='__main__':
     
    parser = argparse.ArgumentParser(description='Increase the runoff into the Arctic ocean by a given percentage (default = 30 %)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input data file (runoff on MITgcm grid)')                    
    parser.add_argument('outfile', help='Enter name for output file')  
    parser.add_argument('gridpath', help='Enter directory containing MITgcm grid (XC,YC)') 
    parser.add_argument('--increasing_factor',default=1.3,type=float,help='factor by what to increase the runoff')
    parser.add_argument('--latitude',default=55.0,type=float,help='latitude north above which to increase the runoff')

    args = parser.parse_args()
    
    print ('Data file :',args.infile) 
    print ('Grid directory:',args.gridpath) 
    print ('Increasing factor',args.increasing_factor)
    print ('Increase above latitude:',args.latitude)
    
    #}}}1

# Functions #{{{1

def readfield(fname,dims,datatype): #{{{2
    """Call signatures::

    readfield(filename, dims, numpy.datatype)
    
    Read unblocked binary data with dimensions "dims".
    """
#    print dims
    fid = open(fname,"rb")
    v   = np.fromfile(fid, datatype)
    fid.close()

#    if sys.byteorder == 'little': v.byteswap(True)

    v.byteswap(True)

    if   len(v) == np.prod(dims):     v = v.reshape(dims)
    elif len(v) == np.prod(dims[1:]): v = v.reshape(dims[1:])
    else:
        errstr = (  "dimensions do not match: \n len(data) = " + str(len(v)) 
                  + ", but prod(dims) = " + str(np.prod(dims)) )
        raise RuntimeError(errstr)

    return v 
#}}}2

def writefield(fname,data):  #{{{2
    '''Write data to binary file
       change from little to big endian
       fname: filename to which to save the data
       data: numpy-array to write to file
    '''

    import sys
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'write to file: '+fname

    data = np.float32(data)

    if sys.byteorder == 'little': 
#        print 'swap bytes'
        data.byteswap(True)
    fid = open(fname,"wb")
    data.tofile(fid)
    fid.close()

    if sys.byteorder == 'little': data.byteswap(True)      # Swap bag to former endianness
#}}}2

def increase_runoff(infile,grid_path,increasing_factor,latitude):

    infile = readfield(infile,(12,1170,90),np.float32).reshape((12,105300))
    
    x_mit = rdmds(grid_path+'XC',astype=None).flatten()
    y_mit = rdmds(grid_path+'YC',astype=None).flatten()

    outfile = np.zeros_like(infile)

    for i,point in enumerate(y_mit):

        if y_mit[i] >= latitude:
            outfile[:,i] = infile[:,i]*increasing_factor
        else:
            outfile[:,i] = infile[:,i]

    outfile = np.reshape(outfile,(12,1170,90))

    return outfile



#}}}1



if __name__=='__main__':
 
    v_new = increase_runoff(args.infile,args.gridpath,args.increasing_factor,args.latitude)
    
#    outname = args.parameter + args.outfile
    outname = args.outfile
    writefield(outname,v_new)

#    fig = plt.figure(figsize=(15,6))
#    fig.clf()
#    pcol(sq(llc.flat(v_new[5,:,:])),cmap='viridis')
##    pcol(llc.flat(t_new2[6,:,:]),cmap='viridis')
#    plt.xlabel('Latitude [$\degree$]', fontsize=14)
#    plt.ylabel('Longitude [$\degree$]', fontsize=14)
#    cbar = plt.colorbar()
#    plt.show()



