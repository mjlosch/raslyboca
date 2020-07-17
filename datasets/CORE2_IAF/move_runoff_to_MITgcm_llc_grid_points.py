#! /usr/bin/python
# coding: utf-8

# Move runoff points from lat lon grid to nearest point on MITgcm Lat-Lon-Cap grid 


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
     
    parser = argparse.ArgumentParser(description='interpolate (reanalysis) NetCDF data from a regular grid to MITgcm grid (by now, only llc tested), output in binary file (ready to be read by MITgcm)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input data file (NetCDF file)')                    
    parser.add_argument('outfile', help='Enter name for binary output file')  
    parser.add_argument('gridpath', help='Enter directory containing MITgcm grid (XC,YC)') 
    parser.add_argument('--parameter',default='runoff',help='parameter to interpolate and convert (must correspond to the name in the NetCDF file)')

    args = parser.parse_args()
    
    print 'Data file :',args.infile 
    print 'Grid directory:',args.gridpath 
    print 'parameter to interpolate:',args.parameter #}}}1

# Functions #{{{1

def closest_point(point,coordinates): #{{{2
    '''Find closest point'''
    
    index = np.argmin(cdist(point,coordinates),axis=1)#.argmin()

    return index

#}}}2

def move_points(parameter,infile,grid_path): #{{{2
    '''move variable points to closest points
       on MITgcm llc grid
    '''
   
    print '############### Move points to MITgcm grid ##############'
    data =  Dataset(infile,'r')
    variable = data.variables[parameter][:]
    lat = data.variables['yc'][:]
    lon = data.variables['xc'][:]
    area = data.variables['area'][:]

    volume_in = variable*area[None,:,:]/1000.  # kg/s/m^2 to m^3/s (divided by density)

    print 'VOLUME IN (m^3/s)',np.sum(volume_in)


    x_mit = rdmds(grid_path+'XC',astype=None)
    y_mit = rdmds(grid_path+'YC',astype=None)
    area2 = rdmds(grid_path+'RAC',astype=None)
    coordinates = np.vstack((y_mit.flat,x_mit.flat))


    var_new = np.zeros((len(variable[:,0,0]),len(x_mit[:,0]),len(x_mit[0,:])))
    
    for t in np.arange(0,len(var_new[:,0,0])):
        index_input = np.transpose(np.nonzero(volume_in[t,:,:]))
#        print 'timestep,number of non-zero points:',t,len(index_input)
      
#        print np.vstack((lat[index_input[:,0],index_input[:,1]],lon[index_input[:,0],index_input[:,1]])).shape,np.transpose(coordinates).shape
        index = closest_point(np.vstack((lat[index_input[:,0],index_input[:,1]],lon[index_input[:,0],index_input[:,1]])).T,coordinates.T)
#        print index.shape
#        print index[0]

        indices = np.unravel_index(index,(len(x_mit[:,0]),len(x_mit[0,:])))
        x_ind = indices[0]
        y_ind = indices[1]
#        print 'x',x_ind.shape,x_ind[0]

        for ind in np.arange(len(x_ind)):
            var_new[t,x_ind[ind],y_ind[ind]] = var_new[t,x_ind[ind],y_ind[ind]]+volume_in[t,index_input[ind,0],index_input[ind,1]]
#            print 'VAR2',var_new[t,x_ind[ind],y_ind[ind]],volume_in[t,index_input[ind,0],index_input[ind,1]]

            

#        for point in index_input:
#            i = point[0]
#            j = point[1]
##            print point, lat[i,j], lon[i,j]
#            
#            index = closest_point([[lat[i,j], lon[i,j]]],coordinates.T)
#            print index
#
#            indices = np.unravel_index(index,(len(x_mit[:,0]),len(x_mit[0,:])))
#            print indices
#
#            var_new[t,indices[0],indices[1]] = var_new[t,indices[0],indices[1]]+volume_in[t,i,j]
#            print 'VAR2',var_new[t,indices[0],indices[1]],volume_in[t,i,j]
        
#    print 'Old and new sum',np.sum(variable),np.sum(var_new)

    # convert from volume (m^3/s) to m/s
    area2[area2 == 0.0] = 1e-33
    var_new=var_new/area2[None,:,:]
   
    data.close()

    print 'VOLUME OUT',np.sum(var_new*area2[None,:,:])
    

    return var_new

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

def plot_interpolated_data(outfile,data): #{{{2
    '''Plot interpolated data before and after writing to file
       outfile: name of the output file
       data: interpolated array
    '''
    fig = plt.figure(figsize=(15,6))
    fig.clf()
    pcol(sq(llc.flat(data[6,:,:])),cmap='viridis')
#    pcol(llc.flat(t_new2[6,:,:]),cmap='viridis')
    plt.xlabel('Latitude [$\degree$]', fontsize=14)
    plt.ylabel('Longitude [$\degree$]', fontsize=14)
    cbar = plt.colorbar()

    fid = open(outfile,'r')
    v_new2 = np.fromfile(fid,np.float32)
    fid.close()

    v_new2 = np.reshape(v_new2,(12,1170,90))
    v_new2.byteswap(True)

    fig = plt.figure(figsize=(15,6))
    fig.clf()
    pcol(sq(llc.flat(v_new2[6,:,:])),cmap='viridis')
    plt.xlabel('Latitude [$\degree$]', fontsize=14)
    plt.ylabel('Longitude [$\degree$]', fontsize=14)
    cbar = plt.colorbar()
    plt.show()

#}}}1



if __name__=='__main__':
 
    v_new = move_points(args.parameter,args.infile,args.gridpath)

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


# Plot output (before and after writing) #{{{1
#    plot_interpolated_data(outname,v_new)

    fid = open(outname,'r')
    test = np.fromfile(fid,np.float32)
    fid.close()

    test.byteswap(True)
    grid = rdmds(args.gridpath+'XC',astype=None)

    test=np.reshape(test,(12,len(grid[:,0]),len(grid[0,:])))
    area = rdmds(args.gridpath+'RAC',astype=None)

    print 'VOLUME AFTER WRITING:',np.sum(test*area[None,:,:])
#    print 'READ:',np.sum(test[7,:,:]*area)
#    print 'READ:',np.sum(test[8,:,:]*area)


 #}}}1


