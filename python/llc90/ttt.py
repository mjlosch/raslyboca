#!/usr/bin/env python
import sys,os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq
import datetime

#m = Basemap(projection='moll',lon_0 = 0., resolution = 'c')
#m = Basemap(projection='hammer',lon_0 = 0., resolution = 'c')
#m = Basemap(projection='merc', resolution = 'c',llcrnrlon=-170,llcrnrlat=-85,urcrnrlon=190,urcrnrlat=85)
#m = Basemap(projection='npstere',lon_0 = 0., boundinglat=60, resolution = 'l') 
#m = Basemap(projection='npaeqd',lon_0 = 0., boundinglat=60, resolution = 'l') 
#m = Basemap(projection='spaeqd',lon_0 = 0., boundinglat=-50, resolution = 'l') 
#m = Basemap(projection='spstere',lon_0 = 0., boundinglat=-60, resolution = 'l') 

it2=[35795, 35857, 35913, 35975, 36035, 36097, 36157, 36219, 36281, 36341, 36403, 36463]
it2=[np.Inf]
#it2=[72743]
#it2=[72437]
#it2=[61, 119, 181, 241, 303, 363, 425, 487, 547, 609, 669, 731]

#it2=[108497,108581,108674,108764,108857,108947,109040,109133,109223,109316,109406,109499]
#it2=[108497]
#it2=[np.nan]

if len(it2)==1: 
    d,iter,meta=rdmds('diag2Dm',it2[0],returnmeta=True)
else: 
    d,iter,meta=rdmds('diag2Dm',it2,returnmeta=True)

mydir = os.getcwd().split(os.sep)[-1]

xg,yg=rdmds('XG'),rdmds('YG')
xc,yc=rdmds('XC'),rdmds('YC')

hf=rdmds('hFacC');

d=np.where(np.abs(d+999.)<1e-5,0.,d)

# #m = Basemap(projection='moll',lon_0 = 0., resolution = 'c')
# m = Basemap(projection='cyl',lon_0 = 10., resolution = 'c')
# #m = Basemap(projection='merc',lon_0 = 0., resolution = 'c')
# #m = Basemap(projection='ortho',lon_0 = 10., lat_0=90.,resolution = 'c')
# #m = Basemap(projection='cyl',
# #            llcrnrlon=-170,llcrnrlat=80, 
# #            urcrnrlon=190, urcrnrlat=90, resolution = 'c')
# #m = Basemap(projection='npstere',lon_0 = 0., boundinglat=60, resolution = 'l') 
# #m = Basemap(projection='splaea',lon_0 = 0., boundinglat=-50, resolution = 'l') 
# #m = Basemap(projection='stere',lon_0 = 0., lat_0=90., width = 1e7, height=1e7, boundinglat=-50, resolution = 'l') 
# plt.figure();
# plt.clf();
# #llc.pcol(xc,yc,(xc),m) #d[0,:,:])
# #llc.pcol(xg,yg,sq(d[1,:,:]*hf[0,:,:]),m) #,vmax=4.) #d[0,:,:])
# #llc.pcol(xg,yg,xc,m,edgecolors='k',facecolor='none') #d[0,:,:])
# #llc.pcol(xg,yg,sq(d[3,:,:]*hf[0,:,:]),m,vmax=4.,edgecolors='k',facecolor='none') #d[0,:,:])
# llc.pcol(xg,yg,sq(d[3,:,:]*hf[0,:,:]),m,vmax=4) #,edgecolors='k',facecolor='none')
# #llc.pcol(xg,yg,sq(d[1,:,:]*hf[0,:,:]),m) #d[0,:,:])
# m.drawmapboundary()
# #m.drawcoastlines(color='k')
# #m.fillcontinents(color='w')
# plt.colorbar(orientation='horizontal')
# #plt.colorbar()
# plt.show()
# sys.exit('bye bye')

vind,vmax =3, 6
#vind,vmax =4, 1
#vind,vmax =0, 1000
tmp = sq(xc)
fig=plt.figure();
m1 = Basemap(projection='cyl',lon_0 = 10., resolution = 'l')
m2 = Basemap(projection='npaeqd',lon_0 = 0., boundinglat= 60, resolution = 'l') 
m3 = Basemap(projection='spaeqd',lon_0 = 180., boundinglat=-50, resolution = 'l') 
refdate = 1
for i in range(len(iter)):

    if len(iter)==1: tmp = sq(d[vind,:,:])
    else:           tmp = sq(d[i,vind,:,:])

    daysSinceRefDate = iter[i]/3
    yearsSinceRefDate = daysSinceRefDate/365
    monthsSinceRefDate = round(np.mod(daysSinceRefDate,365)/(365./12.))
    mydate = "%04i/%02i"%(refdate + yearsSinceRefDate, monthsSinceRefDate)
    print "iteration[%i] = %i, %s"%(i,iter[i],mydate)
    fig.clf()
    fig.add_subplot(2,1,1)
    lh = llc.pcol(xg,yg,tmp,m1,vmax = vmax)
    plt.colorbar() #orientation='horizontal')
    m1.drawcoastlines(color='k')
    plt.title("%s: time step %06i, %s, min/max= %5.2f/%5.2f"%(
            mydir, iter[i], mydate, round(tmp.min(),2),round(tmp.max(),2)) )
    
    fig.add_subplot(2,2,3)
    lh = llc.pcol(xg,yg,tmp,m3,vmax = vmax)
#    m3.drawmapboundary(fill_color='w')
#    m3.fillcontinents(color='k',lake_color='g',alpha=0.1)
    m3.drawcoastlines(color='k')

    fig.add_subplot(2,2,4)
    lh = llc.pcol(xg,yg,tmp,m2,vmax = vmax)
#    m2.drawmapboundary(fill_color='w')
#    m2.fillcontinents(color='k',lake_color='g',alpha=0.1)
    m2.drawcoastlines(color='k')

    fname = fname = "diag2Dm_%02i_%010i"%(vind,iter[i])
#    plt.savefig(fname)
    plt.show()
#    plt.draw()
#    plt.pause(.001)

