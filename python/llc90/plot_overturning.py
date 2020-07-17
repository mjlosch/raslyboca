#!/usr/bin/env python
import sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol

files = sorted(glob.glob("diag3Dm.*.data"))

usedfiles = files[-12:]
#usedfiles = files[:12]
its = list(usedfiles)
for k, file in enumerate(usedfiles):
    its[k] = int(file.split('.')[1])

#its=[35795, 35857, 35913, 35975, 36035, 36097, 36157, 36219, 36281, 36341, 36403, 36463]
#its=[np.Inf]
#its=[72743]
#its=[72437]
#its=[61, 119, 181, 241, 303, 363, 425, 487, 547, 609, 669, 731]

#its=[108497,108581,108674,108764,108857,108947,109040,109133,109223,109316,109406,109499]
#its=[108497]
#its=[np.nan]
#its=[911]

print "reading iterations "+str(its)

if len(its)==1: 
    d,iter,meta=rdmds('diag3Dm',its[0],rec=[4,5],returnmeta=True)
else:
    d,iter,meta=rdmds('diag3Dm',its,rec=[4,5],returnmeta=True)

if len(d.shape)==5: d = np.mean(d,axis=0)

myyear = np.round(np.mean(iter)*8*3600./(365*86400.)) + 1901
mydir = os.getcwd().split(os.sep)[-1]

xg,yg=rdmds('XG'),rdmds('YG')
xc,yc=rdmds('XC'),rdmds('YC')

drf = rdmds('DRF')
dxg,dyg=rdmds('DXG'),rdmds('DYG')

zg=-np.hstack([0.,np.cumsum(drf[:,0,0])])
zc=0.5*(zg[:-1]+zg[1:])


hf=rdmds('hFacC');
mskv=rdmds('hFacS');
msku=rdmds('hFacW');

mskv[mskv>0.]=1.
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
    # psi = -(np.cumsum((np.sum(vf[0]*mskv[0],axis=-1)
    #                    + np.sum(vf[1]*mskv[1],axis=-1)
    #                    + np.sum(-(uf[3]*msku[3])[:,:,::-1],axis=-2)
    #                    + np.sum(-(uf[4]*msku[4])[:,:,::-1],axis=-2)
    #                    )[::-1,:],axis=0)
    #         )[::-1,:]

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

psi = overturning(d[-2,:,:,:]*msku,d[-1,:,:,:]*mskv,dxg,dyg,drf)

psib= barostream(d[-2,:,:,:],d[-1,:,:,:],dxg,dyg,drf)

fig=plt.figure()
fig.clf()
plt.pcolormesh(yg[:270,0],zc,sq(psi)*1e-6); 
#plt.clim([-10,10])
plt.colorbar(orientation='horizontal')
plt.title("%s, %04i: overturning stream function (Sv)"%(mydir,myyear))
    
fig=plt.figure()
fig.clf()
plt.pcolormesh(sq(psib)*1e-6); 
plt.colorbar(orientation='horizontal')
plt.title("%s, %04i: overturning stream function (Sv)"%(mydir,myyear))

plt.show()
# sys.exit('bye bye')
