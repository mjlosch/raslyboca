#!/usr/bin/env python
import sys,os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq


it2=[92, 176, 269, 359, 452, 542, 635, 728, 818, 911, 1001, 1094]
fldList = ['EXFuwind','EXFvwind','EXFpreci','EXFatemp','EXFaqh  ','EXFlwdn ','EXFswdn ']

fname = os.path.join('/uv/user/mlosch/judith/MITgcm/llc90/run01',"diagForc")

if len(it2)==1: f,iter,meta=rdmds(fname,it2[0],returnmeta=True)
else: f,iter,meta=rdmds(fname,it2,returnmeta=True)

gpath = "/uv/user/mlosch/judith/MITgcm/llc90/grid"
xg,yg=rdmds(os.path.join(gpath,'XG')),rdmds(os.path.join(gpath,'YG'))
xc,yc=rdmds(os.path.join(gpath,'XC')),rdmds(os.path.join(gpath,'YC'))
cs,sn=rdmds(os.path.join(gpath,'AngleCS')),rdmds(os.path.join(gpath,'AngleSN'))

hf=rdmds(os.path.join(gpath,'hFacC'))
msk=np.copy(hf[0,:,:])
msk[msk==0] = np.nan

m = Basemap(projection='nplaea',lon_0 = 0., boundinglat=60, resolution = 'l') 
#m = Basemap(projection='splaea',lon_0 = 0., boundinglat=-60, resolution = 'l') 
#m = Basemap(projection='cyl',lon_0 = 10., resolution = 'c')

wspd = np.sqrt(f[:,0,:,:]**2 + f[:,1,:,:]**2)

fig=plt.figure()

#fig.clf(); lh = llc.pcol(xg,yg,yg); plt.show()
#sys.exit()


#d3=rdmds('diag3Dm',iter)

i=0
for i in range(1):
    uu = np.mean(f[:,0,:,:],axis=0)
    vv = np.mean(f[:,1,:,:],axis=0) 
#    uu = np.mean(d3[:,2,1,:,:],axis=0)
#    vv = np.mean(d3[:,3,1,:,:],axis=0)
#    uu = np.ma.masked_array(uu,hf[0,:,:]<=0.0)
#    vv = np.ma.masked_array(vv,hf[0,:,:]<=0.0)
    spd=np.sqrt(uu*uu+vv*vv)
    u,v = Basemap.rotate_vector(m,uu*cs-vv*sn,uu*sn+vv*cs,xc,yc)
    
    fig.clf()
    lh = llc.pcol(xg,yg,sq(wspd[i,:,:]*hf[0,:,:]),m); 
    m.drawcoastlines(color='k')
    di = 1
    q =m.quiver(xc[::di,::di],yc[::di,::di],u[::di,::di],v[::di,::di],latlon=True,scale=200,
                  pivot = 'mid', edgecolor = 'none',headaxislength = 5) 
#                  pivot = 'mid',units = 'x', edgecolor = 'none',headaxislength = 5) 

    plt.show()

sys.exit()


i=0
fig.clf()
for vind in range(len(fldList)+1):
    fig.add_subplot(4,2,vind+1)
    if vind<len(fldList):
        tmp=sq(f[i,vind,:,:])
        fldName = fldList[vind]
    else:
        tmp=sq(wspd[i,:,:])
        fldName = "Wind speed"

    tmp = np.ma.masked_where(np.isnan(msk), tmp)

    lh = llc.pcol(xg,yg,sq(tmp),m); m.drawcoastlines(color='k')
#    lh = llc.pcol(xg,yg,sq(tmp))
    plt.colorbar() #orientation='horizontal')
    plt.title(fldName+', time step '+str(iter[i])+', min/max= '+str(round(tmp.min(),2))+"/"+str(round(tmp.max(),2)))

plt.show()
