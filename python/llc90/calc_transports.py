#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
######################## -*- coding: utf-8 -*-
"""Usage: calc_transports.py INPUTFILE(S)
"""
import sys, os
from getopt import gnu_getopt as getopt
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol

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

print "reading "+str(len(its))+" iterations: "+str(its)

xg,yg=rdmds('XG'),rdmds('YG')
xc,yc=rdmds('XC'),rdmds('YC')

drf = rdmds('DRF')
dxg,dyg=rdmds('DXG'),rdmds('DYG')

zg=-np.hstack([0.,np.cumsum(drf[:,0,0])])
zc=0.5*(zg[:-1]+zg[1:])


hf=rdmds('hFacC');
mskv=rdmds('hFacS');
mskv[mskv>0.]=1.
msku=rdmds('hFacW');
msku[msku>0.]=1.
# atlantic mask

def atlantic_mask(msku,mskv):
    # this is pure convenience
    fmsku=llc.faces(msku)
    fmskv=llc.faces(mskv)
    # mask southern ocean first
    fmsku[0][:,:123,:]=0.
    fmsku[4][:,:,148:]=0.
    fmskv[0][:,:123,:]=0.
    fmskv[4][:,:,147:]=0.
    # atlantic mask
    fmsku[0][:,:,57:]=0;
    fmsku[0][:,203:221,37:]=0;
    fmsku[0][:,208:210,33:37]=0;
    fmsku[1][:]=0.
    fmsku[2][:]=0.
    fmsku[3][:]=0.
    fmsku[4][:,:30,:]=0;
    fmsku[4][:,:51,91:]=0;
    fmsku[4][:,:59,123:]=0;
    fmsku[4][:,:42,85:]=0;
    fmsku[4][:,:54,:40]=0;
    #
    fmskv[0][:,:,57:]=0;
    fmskv[0][:,203:221,37:]=0;
    fmskv[0][:,209,33:37]=0;
    fmskv[1][:]=0.
    fmskv[2][:]=0.
    fmskv[3][:]=0.
    fmskv[4][:,:55,:40]=0;
    fmskv[4][:,:31,:]=0;
    fmskv[4][:,:36,84:]=0;
    fmskv[4][:,:51,92:]=0;
    fmskv[4][:,:60,123:]=0;
    fmskv[4][:,:45,86:]=0;

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

def barostream(u,v,dxg,dyg,drf):
    
    ny,nx=dxg.shape
    nz   =drf.shape[0]
    dydr=np.tile(drf,(1,ny,nx))*np.tile(dyg.reshape((1,ny,nx)),(nz,1,1))
    dxdr=np.tile(drf,(1,ny,nx))*np.tile(dxg.reshape((1,ny,nx)),(nz,1,1))
    ubf=llc.faces(np.sum(u*dydr,axis=-3))
    vbf=llc.faces(np.sum(v*dxdr,axis=-3))
    ub=np.concatenate((ubf[0],ubf[1],vbf[3][:,::-1].transpose(),
                       vbf[4][:,::-1].transpose()
                       ),axis = -1)
    psib=np.cumsum(np.ma.masked_array(ub,ub==0.),axis=0)

    return psib

msku,mskv=atlantic_mask(msku,mskv)
refdate = 1901

ott=list(its)
dpt=list(its)
for k,i in enumerate(its):
    daysSinceRefDate = i/3
    yearsSinceRefDate = daysSinceRefDate/365
    monthsSinceRefDate = round(np.mod(daysSinceRefDate,365)/(365./12.))
    mydate = "%04i/%02i"%(refdate + yearsSinceRefDate, monthsSinceRefDate)
    d=rdmds(files[k].split('.')[0],i,rec=[4,5])
    psi = overturning(d[0,:,:,:]*msku,d[1,:,:,:]*mskv,dxg,dyg,drf)
    psib= barostream(d[0,:,:,:],d[1,:,:,:],dxg,dyg,drf)

    ott[k] = psi[zg[:-1]<-1000,:][:,yg[:270,0]>40.].max()
    dpt[k] = psib[93,333]
    print "iteration[%i] = %i, %s, ot=%f, dp=%f"%(k,i,mydate,
                                                  ott[k]*1e-6,dpt[k]*1e-6)


fig=plt.figure()
fig.clf()
ax = fig.add_subplot(211)
ax.plot(its,np.array(ott)*1e-6)
ax.grid()
plt.title('< 1000m-maximum of AMOC (Sv)')
bx = fig.add_subplot(212)
bx.plot(its,np.array(dpt)*1e-6)
bx.grid()
plt.title('Drake Passage transport (Sv)')
plt.xlabel('timestep')

plt.show()
