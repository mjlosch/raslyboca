#! /usr/bin/python
# coding: utf-8

# Move runoff points lying on land to the nearest point in the ocean (will only work for LLC grids)


# Elena Gerwing, Mar 2017                                                                                                
# Import modules #{{{1
import sys, os, glob

import numpy as np
from scipy import interpolate
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
     
    parser = argparse.ArgumentParser(description='Find points lying on land and move them to the nearest point in the ocean', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Enter input data file (on MITgcm grid, probably output from interpolate_to_MITgcm_llc_grid.py)')                    
    parser.add_argument('outfile', help='Enter name for binary output file')  
    parser.add_argument('gridpath', help='Enter directory containing MITgcm grid and mask (XC,YC,hFacC)') 
    parser.add_argument('--iterations',default=19,type=int,help='number of surrounding points in which to search for ocean point')
#    parser.add_argument('--extension',default=19,type=int,help='number of surrounding points in which to search for ocean point')
    parser.add_argument('--bathyfile',default='bathy_eccollc_90x50_min2pts_mod3.bin',help='bathymetry file (needs to have same dimension as grid)')
    args = parser.parse_args()
    
    print 'Data file :',args.infile 
    print 'Grid directory:',args.gridpath #}}}1

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

def list_of_movements(runoff,mask,iterations,bathy,extension=10): #{{{2
    '''Make list of points on land to move to ocean, 
       list contains source (on land) and destination (in ocean) coordinate
       runoff: runoff field (on MITgcm grid)
       mask: land-ocean-mask
       iterations: number of surrounding points in which to search for ocean-point
       bathy: bathymetry field
       extension: number of cells by which to extend current face
                in order to search for ocean points on other faces
    '''

    # array where only points on land are non-zero
    mask_inv = np.logical_not(mask)
    on_land = np.sum(runoff,axis=0)*mask_inv

    # reshape mds arrays to faces
    land_f = llc.faces(on_land)
    mask_f = llc.faces(mask)
    bathy_f = llc.faces(bathy)

    # rotate and relocate faces: move (differently shaped) Arctic face to the end of the list (4)
    # and rotate faces 2 and 3 in order to match orientation of faces 0 and 1
    mask_f = resort_faces(mask_f)
    land_f = resort_faces(land_f)
    bathy_f = resort_faces(bathy_f)

    moves = []
    # vector of order of faces, needed to find faces left and right of current face
    order = np.concatenate([[4],np.tile(np.array([0,1,2,3]),3),[4]])
    
    for face in np.arange(0,5):      # loop over faces
        k = 5+face                   # place in order vector
        if face == 4:
            k = -1

        # extend current face, add (10) surrounding cells at all boundaries
        mask_ext = embed_face(mask_f,face,order,k,extension)
        bathy_ext = embed_face(bathy_f,face,order,k,extension)
        land_ext = embed_face(land_f,face,order,k,extension)

#        plt.figure(3)
#        pcol(sq(mask_ext))
#        plt.show()
       
        # index of points lying on land
        index_land = np.transpose(np.nonzero(land_f[face]))
        index_land = index_land+extension       # account for extention of face

        for point in index_land:    # loop over points on land

            # search in surrounding cells for ocean-points with an increasing radius in each iteration
            for m in np.arange(1,iterations+1):
                if np.any(mask_ext[point[0]-m:point[0]+(m+1),point[1]-m:point[1]+(m+1)] == 1.):   # any point in ocean (mask = 1)

                    ocean = np.transpose(np.nonzero(mask_ext[point[0]-m:point[0]+(m+1),point[1]-m:point[1]+(m+1)])) # relative indices of all ocean-points to current point

                    out = point+ocean-m  # absolute position of ocean points
                    ocean_min = ocean[np.argmin(bathy_ext[out[:][:,0],out[:][:,1]])]  # find ocean-point with lowest elevation

                    point = point-extension                # account for extension of face (get coordinate without extension)
                    source = np.insert(point,0,face) # coordinate of point on land [face,y,x] 
                    outlet=point+ocean_min-m         # coordinate of ocean-point with lowest elevation [y,x]
#                    print 'out',outlet
            
                    # find coordinate if ocean-point is outside of current face #{{{3
                    if outlet[1] < 0:                                 # left of face
                        if face == 4:
                            outlet[1] = mask_f[order[k-1]].shape[0] + outlet[1]       # find correct x-coordinate
                            outlet = rotate_coordinate(outlet,1,np.rot90(mask_f[order[k-1]]))  # find coordinate on back-rotated left face
                        else:
                            outlet[1] = mask_f[order[k-1]].shape[1] + outlet[1]       # find correct x-coordinate
                        dest = np.insert(outlet,0,order[k-1])                         # coordinate of ocean-point [face,y,x]
                    elif outlet[1] >= mask_f[face].shape[1]:          # right of face
                        if face == 4:
                            outlet[1] = outlet[1] - mask_f[face].shape[1]
                            outlet = rotate_coordinate(outlet,3,np.rot90(mask_f[order[k-3]]))
                            dest = np.insert(outlet,0,order[k-3])
                        else:
                            outlet[1] = outlet[1] - mask_f[face].shape[1]
                            dest = np.insert(outlet,0,order[k+1])
                    elif outlet[0] < 0:                               # below face
                        outlet[0] = mask_f[order[k-4]].shape[0] + outlet[0]           # find correct y-coordinate
                        dest = np.insert(outlet,0,order[k-4])                         # coordinate of ocean-point [face,y,x]
                    elif outlet[0] >= mask_f[face].shape[0]:          # above face
                        outlet[0] = outlet[0] - mask_f[face].shape[0]
                        if face == 0:
                            dest = np.insert(outlet,0,order[0]) 
                        elif face == 4: 
                            outlet = rotate_coordinate(outlet,order[k-2],mask_f[order[k-2]])
                            dest = np.insert(outlet,0,order[k-2])  
                        else:
                            outlet = rotate_coordinate(outlet,order[k],mask_f[order[0]])   # back-rotate coordinate on Arctic face
                            dest = np.insert(outlet,0,order[0]) #}}}3
                    else:
                        dest = np.insert(outlet,0,order[k])   # coordinate of ocean-point [face,y,x] 
                    
                    move = list([source,dest])  # write source and destionation point as list elements
       
#                    print move

                    break

                if m == iterations:
                    point = point - extension
                    raise ValueError('ERROR: no wet point found for point {} on face {}. To solve problem the iteration could be increased'.format(point, face))
             
            moves.append(move)                  # append source and destionation point to list

    return moves #}}}2

def move(runoff,moves,mask):#{{{2
    '''move runoff from land- to ocean-point
       runoff: runoff field (on MITgcm grid)
       moves: list of source and destination points
       mask: land-ocean-mask
    '''
    
    area = rdmds(args.gridpath+'RAC')

    new_runoff = np.zeros_like(runoff)
    
    for month in np.arange(0,12):     # loop over months
        
        # move mds to list of faces and resort faces
        runoff_copy = np.copy(runoff[month,:,:])*area
        runoff_f = llc.faces(runoff_copy)
        runoff_f = resort_faces(runoff_f)

        for i in np.arange(len(moves)):  # loop over points to move

            # source and destination coordinates for face, y and x
            source_face=moves[i][0][0]
            source_y =  moves[i][0][1]
            source_x =  moves[i][0][2]
            dest_face = moves[i][1][0]
            dest_y =    moves[i][1][1]
            dest_x =    moves[i][1][2]
            
#            print '###########################################'
#            print 'source_coord',source_face,source_y,source_x
#            print 'dest_coord', dest_face,dest_y,dest_x 
#            print 'dest',runoff_f[dest_face][dest_y,dest_x]

            # add runoff of source cell to destination cell
            runoff_f[dest_face][dest_y,dest_x] = runoff_f[source_face][source_y,source_x] + runoff_f[dest_face][dest_y,dest_x]
#            print 'source', runoff_f[source_face][source_y,source_x] 
#            print 'new',runoff_f[dest_face][dest_y,dest_x] 


        # resort faces back and move list to mds array
        runoff_f = resort_back_faces(runoff_f)
        runoff_m = llc.faces2mds(runoff_f)
        
        runoff_m = runoff_m*mask    # delete runoff on land

        area[area == 0.0] = 1e-33
        runoff_m = runoff_m/area

        new_runoff[month,:,:] = runoff_m  # write runoff of month into array
        
    return new_runoff
#}}}2

def rotate_coordinate(coord,rotations,array): #{{{2
    '''calculate new coordinate after (back) rotation of face 
       coord: coordinate
       rotations: number how often face was rotated (in order to be extended to current face)
       array: rotated face (only needed for shape)
    '''

    new_coord = np.empty_like(coord)
    if rotations == 1:
        new_coord[0] = coord[1]
        new_coord[1] = array.shape[0]-1 - coord[0]
    elif rotations == 2:
        new_coord[0] = array.shape[1]-1 - coord[0]
        new_coord[1] = array.shape[0]-1 - coord[1]
    elif rotations == 3:
        new_coord[0] = array.shape[1]-1 - coord[1]
        new_coord[1] = coord[0]
    else:
        new_coord = coord

    return new_coord #}}}2

def resort_faces(field_faces): #{{{2

    field_faces[3] = np.rot90(field_faces[3])
    field_faces[4] = np.rot90(field_faces[4])
    field_faces[2] = np.rot90(field_faces[2],3)
    field_faces[2], field_faces[3], field_faces[4] = field_faces[3], field_faces[4], field_faces[2]

    return field_faces #}}}2

def resort_back_faces(field_faces): #{{{2
    
    field_faces[3], field_faces[4], field_faces[2] = field_faces[2], field_faces[3], field_faces[4] 
    field_faces[3] = np.rot90(field_faces[3],3)
    field_faces[4] = np.rot90(field_faces[4],3)
    field_faces[2] = np.rot90(field_faces[2])

    return field_faces #}}}2

def embed_face(field,face,order,k,n=10): #{{{2

    if face == 4:
        face1_rot = np.rot90(field[order[k-3]],3)
        face2_rot = np.rot90(field[order[k-2]],2)
        face3_rot = np.rot90(field[order[k-1]],1)
        field_ext = np.hstack((face3_rot[:,-n:],field[face],face1_rot[:,0:n]))
        field_ext = np.vstack((np.hstack((np.zeros((n,n)),field[order[k-4]][-n:,:],np.zeros((n,n)))),field_ext,np.hstack((np.zeros((n,n)),face2_rot[0:n,:],np.zeros((n,n))))))
    else:
        field_ext = np.hstack((field[order[k-1]][:,-n:],field[face],field[order[k+1]][:,:n]))
        face4_rot = np.rot90(field[order[0]],order[k])
        field_ext = np.vstack((np.hstack((np.zeros((n,n)),np.zeros((n,field[face].shape[1])),np.zeros((n,n)))),field_ext,np.hstack((np.zeros((n,n)),face4_rot[:n,:],np.zeros((n,n))))))
    
    return field_ext  #}}}2

#}}}1


if __name__=='__main__':
    
    mask = rdmds(args.gridpath+'hFacC')[0,:,:]
    runoff = readfield(args.infile,(12,len(mask[:,0]),len(mask[0,:])),np.float32)
    bathy = readfield(args.bathyfile,(len(mask[:,0]),len(mask[0,:])),np.float32)

    x_mit = rdmds(args.gridpath+'XC',astype=None)
    y_mit = rdmds(args.gridpath+'YC',astype=None)
    area = rdmds(args.gridpath+'RAC',astype=None)
    ind = np.unravel_index(np.argmax(runoff*area[None,:,:]),runoff.shape)
    print ind
    print 'Max', y_mit[ind[1],ind[2]],x_mit[ind[1],ind[2]] 

    extension = args.iterations

    moves = list_of_movements(runoff,mask,args.iterations,bathy,extension)

    new_runoff = move(runoff,moves,mask)

    ind = np.unravel_index(np.argmax(new_runoff*area[None,:,:]),new_runoff.shape)
    print ind
    print 'Max', y_mit[ind[1],ind[2]],x_mit[ind[1],ind[2]] 

    writefield(args.outfile,new_runoff)

    area = rdmds(args.gridpath+'RAC')

    print 'Old and new volume of runoff',np.sum(runoff*area[None,:,:]),np.sum(new_runoff*area[None,:,:])

