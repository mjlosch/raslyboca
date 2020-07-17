#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
######################## -*- coding: utf-8 -*-
"""Usage: calc_transports.py INPUTFILE(S)
"""
import sys, os
from getopt import gnu_getopt as getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol
import datetime

# parse command-line arguments
try:
    optlist,args = getopt(sys.argv[1:], ':', ['verbose'])
    assert len(args) > 0
except (AssertionError):
    sys.exit(__doc__)

files=[]
if len(args)<2:
    from glob import glob
    for infile in glob(args[0]):
        files.append(infile)
else:
    files=args

mydir = os.getcwd().split(os.sep)[-1]

files.sort()
its = list(files)
for k, file in enumerate(files):
    its[k] = int(file.split('.')[1])

print("reading "+str(len(its))+" iterations: "+str(its))

bdir = '/work/ollie/mlosch/egerwing/llc270'
gdir = '/home/ollie/mlosch/MITgcm/MITgcm/llc270/grid'
xg,yg=rdmds(os.path.join(gdir,'XG')),rdmds(os.path.join(gdir,'YG'))
xc,yc=rdmds(os.path.join(gdir,'XC')),rdmds(os.path.join(gdir,'YC'))

drf = rdmds(os.path.join(gdir,'DRF'))
dxg,dyg=rdmds(os.path.join(gdir,'DXG')),rdmds(os.path.join(gdir,'DYG'))

rac = rdmds(os.path.join(gdir,'RAC'))

zg=-np.hstack([0.,np.cumsum(drf[:,0,0])])
zc=0.5*(zg[:-1]+zg[1:])

hf=rdmds(os.path.join(gdir,'hFacC'))
mskc=np.copy(hf)
mskc[mskc>0.]=1.
mskv=rdmds(os.path.join(gdir,'hFacS'))
mskv[mskv>0.]=1.
msku=rdmds(os.path.join(gdir,'hFacW'))
msku[msku>0.]=1.


refdate = datetime.datetime(1979,2,1,0,0)
def iter2date( myiter ):
    mydate  = refdate + datetime.timedelta(myiter/48.-1)
    return mydate

def date2iter( mydate ):
    myiter = ((mydate-refdate).days+1)*48
    return myiter

# atlantic mask
def atlantic_mask(mskc):
    # this is pure convenience
    fmsk=llc.faces(mskc)
    # mask southern ocean first
    js = 374
    nz,ny,nx = fmsk[0].shape
    fmsk[0][:,:374,172:]=0.
    fmsk[4][:,:188,527:]=0.
    # atlantic mask
    # this rough and removes part of the norwegian sea
    fmsk[0][:,:,180:]=0;
    fmsk[0][:,:380,175:]=0;
    # mediterranean sea
    fmsk[0][:,599:651,96:]=0;
    fmsk[0][:,651:662,125:]=0;
    # baltic and north sea
    fmsk[0][:,703:756,140:]=0;
    fmsk[0][:,756:768,172:]=0;
    fmsk[1][:]=0.
    fmsk[2][:]=0.
    fmsk[3][:]=0.
    # north pacific
    fmsk[4][:,:90,:]=0;
    # south pacific
    fmsk[4][:,:174,287:]=0;
    fmsk[4][:,:133,263:]=0;
    fmsk[4][:,:102,262]=0;
    # hudson bay etc
    fmsk[4][:,:160,:121]=0;
    fmsk[4][:,:189,26:121]=0;

    return llc.faces2mds(fmsk)

def atlantic_maskuv(msku,mskv):
    # this is pure convenience
    fmsku=llc.faces(msku)
    fmskv=llc.faces(mskv)
    # mask southern ocean first
    js = 374
    nz,ny,nx = fmsku[0].shape
    fmsku[0][:,:js,:]=0.
    fmsku[4][:,:,ny-js:]=0.
    fmskv[0][:,:js,:]=0.
    fmskv[4][:,:,ny-js:]=0.
    # atlantic mask
    iw = 176
    # this rough and removes part of the norwegian sea
    fmsku[0][:,:,iw:]=0;
    # mediterranean sea
    fmsku[0][:,599:651,96:]=0;
    fmsku[0][:,651:662,125:]=0;
    # baltic and north sea
    fmsku[0][:,703:756,140:]=0;
    fmsku[0][:,756:758,173:]=0;
    fmsku[1][:]=0.
    fmsku[2][:]=0.
    fmsku[3][:]=0.
    # north pacific
    fmsku[4][:,:90,:]=0;
    # south pacific
    fmsku[4][:,:174,287:]=0;
    fmsku[4][:,:133,263:]=0;
    # hudson bay etc
    fmsku[4][:,:160,:121]=0;
    fmsku[4][:,:189,26:121]=0;
    #
    # this rough and removes part of the norwegian sea
    fmskv[0][:,:,iw-1:]=0;
    # mediterranean sea
    fmskv[0][:,599:651,96:]=0;
    fmskv[0][:,651:662,125:]=0;
    # baltic and north sea
    fmskv[0][:,703:756,140:]=0;
    fmskv[0][:,756:758,172:]=0;
    fmskv[1][:]=0.
    fmskv[2][:]=0.
    fmskv[3][:]=0.
    # north pacific
    fmskv[4][:,:90,:]=0;
    # south pacific
    fmskv[4][:,:174,287:]=0;
    fmskv[4][:,:133,263:]=0;
    fmskv[4][:,90,262] = 0.
    fmskv[4][:,99:102,262] = 0.
    # hudson bay etc
    fmskv[4][:,:160,:121]=0;
    fmskv[4][:,:189,25:121]=0;

    return llc.faces2mds(fmsku), llc.faces2mds(fmskv)

def overturning(u,v,dxg,dyg,drf):

    ny,nx=dxg.shape
    nz   =drf.shape[0]
    dydr=np.tile(drf,(1,ny,nx))*np.tile(dyg.reshape((1,ny,nx)),(nz,1,1))
    dxdr=np.tile(drf,(1,ny,nx))*np.tile(dxg.reshape((1,ny,nx)),(nz,1,1))

    uf=llc.faces(u*dydr)
    vf=llc.faces(v*dxdr)

    vf0=np.sum(vf[0],axis=-1)
    vf1=np.sum(vf[1],axis=-1)
    uf3=np.roll(np.sum(-(uf[3])[:,:,::-1],axis=-2),1,axis=1)
    uf4=np.roll(np.sum(-(uf[4])[:,:,::-1],axis=-2),1,axis=1)
    psi = -np.cumsum((vf0+vf1+uf3+uf4)[::-1,:],axis=0)[::-1,:]

    return psi

def overturningw(w,rac):

    nz,ny,nx=w.shape
    dxdy=np.tile(rac.reshape((1,ny,nx)),(nz,1,1))
    wfx=llc.faces(w*dxdy)
    # integration over longitude
    wf = ( np.sum(wfx[0]+wfx[1],axis=-1)
         + np.sum(wfx[3]+wfx[4],axis=-2)[:,::-1] )
    # integration over latitudes
    psiw = np.cumsum(wf,axis=-1)

    return psiw

def barostream(u,v,dxg,dyg,drf):

    ny,nx=dxg.shape
    nz   =drf.shape[0]
    dydr=np.tile(drf,(1,ny,nx))*np.tile(dyg.reshape((1,ny,nx)),(nz,1,1))
    dxdr=np.tile(drf,(1,ny,nx))*np.tile(dxg.reshape((1,ny,nx)),(nz,1,1))
    ubf=llc.faces(np.sum(u*dydr,axis=-3))
    vbf=llc.faces(np.sum(v*dxdr,axis=-3))
    ub0=np.concatenate((ubf[0],ubf[1],
                        np.rot90(vbf[3],k=1),
                        np.rot90(vbf[4],k=1)),axis = -1)
    # shift arctic
    va=np.roll(vbf[2],-1,axis=0)
    va[-1,:]=0
    vb=np.hstack((np.rot90(-va,k=-1),np.zeros(va.shape),
                  np.rot90(vbf[2],k=1),np.zeros(va.shape)))[:nx//2,:]

    ub = np.concatenate((ub0,vb),axis=-2)
    psib=np.cumsum(np.ma.masked_array(ub,ub==0.),axis=0)
    return psib

#msku,mskv=atlantic_maskuv(msku,mskv)
#mskc=atlantic_mask(mskc)

ott=list(its)
dpt=list(its)
datex=[]
for k,i in enumerate(its):
    mydate = iter2date(i)
    d=rdmds(files[k].split('.')[0],i,rec=[2,3,4])
#    psi = overturning(d[0,:,:,:]*msku,d[1,:,:,:]*mskv,dxg,dyg,drf)
    psiw= overturningw(d[2,:,:,:]*mskc,rac)
    psib= barostream(d[0,:,:,:],d[1,:,:,:],dxg,dyg,drf)

    ott[k] = psiw[zg[:-1]<-1000,:][:,yg[:3*270,0]>40.].max()
    dpt[k] = psib[274,992]
    datex.append(mydate)
    print("iteration[%i] = %i, %s, ot=%f, dp=%f"%(k,i,mydate,
                                                  ott[k]*1e-6,dpt[k]*1e-6))


# clf(); contourf(yg[:3*270,0],zg[:-1],sq(psiw)*1e-6,np.linspace(-20,20,21)); colorbar(orientation='horizontal')
# clf(); pcol(yg[:3*270,0],zg,sq(psiw)*1e-6); colorbar(orientation='horizontal')
# clf(); pcol(sq(psib)*1e-6); colorbar(orientation='horizontal')

fig, ax = plt.subplots(2,sharex=True,squeeze=True) #,figsize=(8,10))

ax[0].plot(datex,np.array(ott)*1e-6)
ax[0].set_title('< 1000m-maximum of AMOC (Sv)')
ax[1].plot(datex,np.array(dpt)*1e-6)
ax[1].set_title('Drake Passage transport (Sv)')
locator = mdates.AutoDateLocator()
formatter = mdates.AutoDateFormatter(locator)
ax[1].xaxis.set_major_locator(locator)
ax[1].xaxis.set_major_formatter(formatter)

for bx in ax: bx.grid()

fig.show()
