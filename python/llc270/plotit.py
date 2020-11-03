# for/scratch/users/mlosch/MITgcm/MITgcm/verification/global_oce_in_p2/
import matplotlib.pyplot as plt
import numpy as np
from MITgcmutils import rdmds, llc
from myutils import *
import cmocean.cm as cm
import os

printit=False
bdir = '/work/ollie/mlosch/raslyboca/llc270'
# myiter = np.Inf
# myiter2 = [53946, 55434, 56874, 58362, 59802, 61290,
#            62778, 64218, 65706, 67146, 68634, 70122]
# myiter2 = [115386, 116826, 118314, 119754, 121242, 122730,
#            124074, 125562, 127002, 128490, 129930, 131418]
# myiter2 = [144528, 146016, 147456, 148944, 150432, 151872,
#            153360, 154800, 156288, 157776, 159168, 160656]
# myiter2 = [305280, 306768, 308256, 309696, 311184, 312624,
#            314112, 315600, 316944, 318432, 319872, 321360]
# myiter = myiter2[-1]

# myiter = np.Inf; myiter2 = myiter
# # 1984-12-31
# myiter = 103728; myiter2 = myiter
# # 1989-12-13
# myiter = 191376; myiter2 = myiter
# # 1994-05-31
# myiter = 268752; myiter2 = myiter
# # 1997-05-31
# myiter = 191376; myiter2 = myiter
# # 1997-05-31
# myiter = 321360; myiter2 = myiter
# 1998-05-31
myiter = 338880; myiter2 = myiter
# # 1998-12-31
# myiter = 349152; myiter2 = myiter

zdir = os.path.join(bdir,'runz00')
rdir = os.path.join(bdir,'runp00')
gdir = '/home/ollie/mlosch/MITgcm/MITgcm/llc270/grid'
xg = rdmds(os.path.join(gdir,'XG'))
yg = rdmds(os.path.join(gdir,'YG'))
rac= rdmds(os.path.join(rdir,'RAC'))
rf = rdmds(os.path.join(rdir,'RF'))[:,0,0]
rc = rdmds(os.path.join(rdir,'RC'))[:,0,0]
drf= rdmds(os.path.join(rdir,'DRF'))
hf = rdmds(os.path.join(rdir,'hFacC'))
zf = rdmds(os.path.join(zdir,'RF'))[:,0,0]
zc = rdmds(os.path.join(zdir,'RC'))[:,0,0]
dzf= rdmds(os.path.join(zdir,'DRF'))
hfz= rdmds(os.path.join(rdir,'hFacC'))

dpz= rdmds(os.path.join(zdir,'Depth'))
dpp= rdmds(os.path.join(rdir,'Depth'))

zdiry = os.path.join(zdir,'20years')
rdiry = os.path.join(rdir,'20years')

rhoConst=1035.
gravity = 9.81
grho = rhoConst*gravity

nz,ny,nx=hf.shape

bp0  = np.sum(hf*np.tile(drf,[1,ny,nx]),axis=0)

mskp = np.copy(hf)
mskp[mskp>0]=1.
# mskz = np.copy(hfz)
# mskz[mskz>0]=1.

#mskp = np.ones(mskp.shape);
#mskz = np.ones(mskz.shape)

lmsk = mskp[-1,:,:]

#dp = rdmds(os.path.join(rdiry,'diag2D'),myiter)
#dz = rdmds(os.path.join(zdiry,'diag2D'),myiter)

sip,myiter,msidata=rdmds(os.path.join(rdiry,'SIdiags'),myiter,returnmeta=True)
d2p,myiter2,m2data=rdmds(os.path.join(rdiry,'diags2D'),myiter2,returnmeta=True)
d3p,myiter,m3data=rdmds(os.path.join(rdiry,'diags3D'),myiter,returnmeta=True)


siz=rdmds(os.path.join(zdiry,'SIdiags'),myiter)
d2z=rdmds(os.path.join(zdiry,'diags2D'),myiter2)
d3z=rdmds(os.path.join(zdiry,'diags3D'),myiter)

if len(myiter2)>1:
    d2p=d2p.mean(axis=0)
    d2z=d2z.mean(axis=0)

def calc_std(fld2,fld):
    var = fld2-fld**2
    return np.where(var<0.,0,np.sqrt(var))

fig, ax = plt.subplots(5,2,sharex=True,sharey=True,squeeze=True,figsize=(8,10))

csf = np.copy(ax)

ivar = 1 # sea ice thickness
plt.sca(ax[0,0])
csf[0,0]=llc.pcol(xg,yg,sq(sip[ivar,:,:]))
plt.sca(ax[0,1])
csf[0,1]=llc.pcol(xg,yg,sq(siz[ivar,:,:]))
ivar = 4 # MXLDEPTH
plt.sca(ax[1,0])
csf[1,0]=llc.pcol(xg,yg,sq(d2p[ivar,:,:])/grho)
plt.sca(ax[1,1])
csf[1,1]=llc.pcol(xg,yg,sq(d2z[ivar,:,:]))
# sea level
# substract mean sea level
etam = (d2p[2,:,:]*lmsk*rac).sum()/(rac*lmsk).sum()
etap = (sq(d2p[2,:,:])-etam)/gravity
#etam = 0.
plt.sca(ax[2,0])
csf[2,0]=llc.pcol(xg,yg,etap)
plt.sca(ax[2,1])
etazm = (d2z[0,:,:]*lmsk*rac).sum()/(rac*lmsk).sum()
etaz =sq(d2z[0,:,:])-etazm
csf[2,1]=llc.pcol(xg,yg,etaz)
# sea level variability
sl2p = calc_std(d2p[3,:,:],d2p[2,:,:])*lmsk/gravity
sl2z = calc_std(d2z[1,:,:],d2z[0,:,:])*lmsk
# surface temperature
# ivar=0
# sl2p = d3p[ivar,-1,:,:]
# sl2z = d3z[ivar,0,:,:]
plt.sca(ax[3,0])
csf[3,0]=llc.pcol(xg,yg,sq(sl2p))
plt.sca(ax[3,1])
csf[3,1]=llc.pcol(xg,yg,sq(sl2z))
# bottom pressure anomaly
bp2p = calc_std(d2p[1,:,:],d2p[0,:,:])*lmsk
bp2z = calc_std(d2z[3,:,:],d2z[2,:,:])*lmsk*rhoConst
plt.sca(ax[4,0])
csf[4,0]=llc.pcol(xg,yg,sq(bp2p)/grho)
plt.sca(ax[4,1])
csf[4,1]=llc.pcol(xg,yg,sq(bp2z)/grho)

ax[0,0].title.set_text('p-coord: timestep %7i'%myiter[0])
ax[0,1].title.set_text('z-coord: timestep %7i'%myiter[0])

varname = ['SIheff (m)','MXLDEPTH (m)','sea level (m)','std(sla)','std(bp) (m)']
clims = [(0,2),(0,1000),(-2.5,1),(0,.1),(0,.06)]
#varname = ['SIheff (m)','MXLDEPTH (m)','sea level (m)','SST(degC),'std(bp) (Pa)'']
#clims = [(0,4),(0,1000),(-2.5,1),(-1.8,30.),(0,200)]
for k, bx in enumerate(ax):
    for kk in range(2):
        # ax[k,kk].set_xlim((-180.,180))
        # ax[k,kk].set_ylim((-90.,90.))
        # ax[k,kk].axis('image')
        plt.colorbar(csf[k,kk][0],ax=ax[k,kk],orientation='vertical')
        for ccsf in csf[k,kk]:
            ccsf.set_clim(clims[k])
    # for ccsf in bcsf:
    #     ccsf.set_clim(clims[k])
    # acsf = csf[k,0]
    # bcsf = csf[k,1]
    ax[k,0].set_ylabel(varname[k])

fig.show()

def section2d(yg, fld3d):
    nr,ny,nx=fld3d.shape
    iy0=140
    iy1=3*nx
    ix = 43
    y, f = yg[iy0:iy1,ix], fld3d[:,iy0:iy1,ix]
    return y, f

f3d, a3d = plt.subplots(5,2,sharex=True,sharey=True,squeeze=True,figsize=(8,10))
csf3d = np.copy(a3d)
ivar = 0 # theta
ysec,fsec = section2d(yg,d3p[ivar,:,:,:])
csf3d[0,0]=a3d[0,0].pcolormesh(ysec,-rf/grho,sq(fsec),cmap=cm.thermal)
ysec,fsec = section2d(yg,d3z[ivar,:,:,:])
csf3d[0,1]=a3d[0,1].pcolormesh(ysec,zf,sq(fsec),cmap=cm.thermal)
ivar = 1 # salinity
ysec,fsec = section2d(yg,d3p[ivar,:,:,:])
csf3d[1,0]=a3d[1,0].pcolormesh(ysec,-rf/grho,sq(fsec),cmap=cm.haline)
ysec,fsec = section2d(yg,d3z[ivar,:,:,:])
csf3d[1,1]=a3d[1,1].pcolormesh(ysec,zf,sq(fsec),cmap=cm.haline)
ivar = 2 # uvel
ysec,fsec = section2d(yg,d3p[ivar,:,:,:])
csf3d[2,0]=a3d[2,0].pcolormesh(ysec,-rf/grho,sq(fsec),cmap=cm.delta)
ysec,fsec = section2d(yg,d3z[ivar,:,:,:])
csf3d[2,1]=a3d[2,1].pcolormesh(ysec,zf,sq(fsec),cmap=cm.delta)
wp = d3p[4,...]; mwp = np.copy(wp); mwp[mwp!=0]=1.
wz = d3z[4,...]; mwz = np.copy(wz); mwz[mwz!=0]=1.
# # this is complicated to account for the shift in numbering of faces
mwp[1:,:] = mwp[1:,:]*mwp[:-1,:]; mwp[0,:]=0
mwz[0,:] = 0
# reset masks for test purposes
#mwp = np.ones(mwp.shape); mwz = np.ones(mwz.shape);
efac = 1e3
#efac = 1.
ivar = 5 # GGL90TKE
# ivar = 7 # IDEMIX_E
# tkep = np.log10(sq(ocp[ivar,...]*mwp))
# tkez = np.log10(sq(ocz[ivar,...]*mwz))
tkep = sq(d3p[ivar,...]*mwp)
tkez = sq(d3z[ivar,...]*mwz)
ysec,fsecp = section2d(yg,tkep)
csf3d[3,0]=a3d[3,0].pcolormesh(ysec,-rc[::-1]/grho,efac*fsecp[::-1,:],
                            cmap=cm.haline)
ysec,fsecz = section2d(yg,tkez)
csf3d[3,1]=a3d[3,1].pcolormesh(ysec, zc,           efac*fsecz[1:,:],
                            cmap=cm.haline)
# ivar = 7 # IDEMIX_E
# idmp = sq(d3p[ivar,...]*mwp)
# idmz = sq(d3z[ivar,...]*mwz)
# ysec,fsecp = section2d(yg,idmp)
# csf3d[4,0]=a3d[4,0].pcolormesh(ysec,-rc[::-1]/grho,efac*fsecp[::-1,:],
#                             cmap=cm.haline)
# ysec,fsecz = section2d(yg,idmz)
# csf3d[4,1]=a3d[4,1].pcolormesh(ysec, zc,           efac*fsecz[1:,:],
#                             cmap=cm.haline)
ivar = 6 # GGL90Kr
kp = sq(d3p[ivar,...]*mwp)/grho**2
kz = sq(d3z[ivar,...]*mwz)
ysec,fsecp = section2d(yg,kp)
csf3d[4,0]=a3d[4,0].pcolormesh(ysec,-rc[::-1]/grho,np.log10(fsecp[::-1,:]),
                             cmap=cm.haline)
ysec,fsecz = section2d(yg,kz)
csf3d[4,1]=a3d[4,1].pcolormesh(ysec, zc,           np.log10(fsecz[1:,:]),
                            cmap=cm.haline)

a3d[0,0].title.set_text('p-coord: timestep %7i'%myiter[0])
a3d[0,1].title.set_text('z-coord: timestep %7i'%myiter[0])

varname = ['THETA (degC)','salinity','UVEL (m/s)','TKE (10$^{-3}$ m$^2$/s$^2$)',
            'log10(Kr) (log10(m$^2$/s))']
clims = [(-1.8,30),(30.34,37.6),(-.07,.07),(0,2e-4*efac),(-4.3,0)]
# varname = ['THETA (degC)','salinity','UVEL (m/s)','TKE (10$^{-3}$ m$^2$/s$^2$)',
#             'IDEMIX_E (10$^{-3}$ m$^2$/s$^2$)']
# clims = [(-1.8,30),(30.34,37.6),(-.07,.07),(0,2e-4*efac),(0,2e-4*efac)]
#clims = [(-1.8,30),(30.34,37.6),(-.07,.07),(-6,-3),(-4.3,0)]
for k, bx in enumerate(a3d[:,:]):
    plt.colorbar(csf3d[k,0],ax=a3d[k,0],orientation='vertical')
    plt.colorbar(csf3d[k,1],ax=a3d[k,1],orientation='vertical')
    a3d[k,0].set_ylabel(varname[k])
#    if k!=3: csf3d[k,0].set_clim(clims[k])
#    if k!=3: csf3d[k,1].set_clim(clims[k])
    csf3d[k,0].set_clim(clims[k]); csf3d[k,1].set_clim(clims[k])
    a3d[k,0].set_ylim((-5200,0));
    a3d[k,1].set_ylim((-5200,0))

# # # csf3d0=a3d[0].imshow(sq(dp[ivar,:,:,ii]*mskp[:,:,ii])[::-1,:])
# # # csf3d1=a3d[1].imshow(sq(dz[ivar,:,:,ii]*mskz[:,:,ii]))
# # # for b in a3d:
# # #     b.grid(b=True)

f3d.show()

fcmp,ax=plt.subplots(2,2,sharex=True,sharey=True,squeeze=True,figsize=(8,6))
csf = np.copy(ax)
# mean sea level
plt.sca(ax[0,0])
csf[0,0]=llc.pcol(xg,yg,etap-etaz,cmap=cm.delta)
# mixed layer depth
ivar = 4 # MXLDEPTH
dmld = sq(d2p[ivar,:,:])/grho-sq(d2z[ivar,:,:])
plt.sca(ax[0,1])
csf[0,1]=llc.pcol(xg,yg,dmld,cmap=cm.delta)
# sea level anomaly
plt.sca(ax[1,0])
csf[1,0]=llc.pcol(xg,yg,sl2p-sl2z,cmap=cm.delta)
# bottom pressure anomaly
plt.sca(ax[1,1])
csf[1,1]=llc.pcol(xg,yg,(bp2p-bp2z)/grho,cmap=cm.delta)

varname = ['mean sea level (m)','mixed layer depth (m)'
           ,'sea level variability (m)','bottom pressure variability (m)']
clims = [(-.1,.1),(-100,100),(-.05,.05),(-.01,.01)]
# min/max values
minmax = [[(etap-etaz).min(), (etap-etaz).max()],
          [dmld.min(), dmld.max()],
          [(sl2p-sl2z).min(), (sl2p-sl2z).max()],
          [((bp2p-bp2z)/grho).min(), ((bp2p-bp2z)/grho).max()]]

for k, bx in enumerate(ax):
    for kk in range(2):
        if k == 0:
            cax =plt.colorbar(csf[k,kk][0],ax=ax[k,kk],orientation='horizontal',
                              extend='both',pad = 0.05)
        else:
            cax =plt.colorbar(csf[k,kk][0],ax=ax[k,kk],orientation='horizontal',
                              extend='both')
        for ccsf in csf[k,kk]:
            ccsf.set_clim(clims[k*ax.shape[0]+kk])

for k, bx in enumerate(ax.flatten()):
    bx.set_title("%s \n range = (%4.2f,%4.2f)"%(varname[k],
                                              minmax[k][0],minmax[k][1]))

fcmp.show()


if printit:
    mydpi = 300
    fig.savefig('cmp2d.%10.10i.png'%myiter[0],dpi=mydpi)
    f3d.savefig('cmp3d.%10.10i.png'%myiter[0],dpi=mydpi)
    fcmp.savefig('diff2d.%10.10i.png'%myiter[0],dpi=mydpi)
