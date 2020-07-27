#!/usr/bin/env python
# coding: utf-8

# # Load modules

# In[1]:


import numpy as np
import matplotlib.pylab as plt
import os
import sys
sys.path.append('../../utils/python/MITgcmutils/')
from MITgcmutils import rdmds,llc
#import xmitgcm
from mpl_toolkits.basemap import Basemap
import datetime

get_ipython().run_line_magic('matplotlib', 'notebook')


# # Read model results for both configurations

# In[42]:


# Read data
dir_th           = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test/'
dir_thad         = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test_thad/output_oldcode'
dir_thad_nc      = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test_thad/'
dir_thad_nc_f    = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test_thad_flood/'
dir_thad_nc_pl   = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test_thad_plume/'
dir_thad_nc_pl_f = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/test_thad_flood_plume/'

fname = 'diag2Dy'

t = np.zeros((6,60,1170,90))
a = np.zeros((6,60,1170,90))

dates = []

for idic,dici in enumerate([dir_th,dir_thad,dir_thad_nc,dir_thad_nc_f,dir_thad_nc_pl,dir_thad_nc_pl_f]):
    flist = [ifile for ifile in os.listdir(os.path.join(dici)) if ifile.startswith(fname) and ifile.endswith('.data')]
    flist.sort()
    steps = [int(ifile.split('.')[-2]) for ifile in flist]
    for istep,stepi in enumerate(steps):
        [a[idic,istep,:,:],t[idic,istep,:,:]] = rdmds(os.path.join(dici,fname),stepi)[2:4,:,:]
        if idic==0:
            dates.append(datetime.datetime(1948,1,1,0,0,0)+datetime.timedelta(0,stepi*28800))


# In[22]:


# Read grid
dir_grid = '/home/ollie/nhutter/raslyboca/MITgcm/llc90/grid'
xg  = rdmds(os.path.join(dir_grid,'XG'))
yg  = rdmds(os.path.join(dir_grid,'YG'))
rac = rdmds(os.path.join(dir_grid,'RAC'))


# In[23]:


# Mask northern and southern hemisphere
mask_sh = (yg<0)
mask_nh = (yg>0)


# In[33]:


# Map projections
mnh = Basemap(projection='npaeqd',boundinglat=60,lon_0=0,resolution='l')
msh = Basemap(projection='spaeqd',boundinglat=-50,lon_0=0,resolution='l')


# In[44]:


# Difference plots
fig = plt.figure(figsize=(10,11))

labels = ['old si code', 'old adx code', 'new adx code', 'n a c + flooding', 'n a c + plumes', 'n a c + flood + plumes']

for i in range(5):
    # Differences in Concentration NH
    fig.add_subplot(5,4,1+i*4)
    llc.pcol(xg,yg,np.nanmean(a[0,:,:,:]-a[i+1,:,:,:],axis=0)*100,mnh,cmap=plt.get_cmap('coolwarm'),vmax=10,vmin=-10)
    mnh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff conc [%] NH')
    plt.gca().set_ylabel(labels[i+1])

    # Differences in Concentration SH
    fig.add_subplot(5,4,2+i*4)
    llc.pcol(xg,yg,np.nanmean(a[0,:,:,:]-a[i+1,:,:,:],axis=0)*100,msh,cmap=plt.get_cmap('coolwarm'),vmax=50,vmin=-50)
    msh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff conc [%] SH')

    # Differences in Thickness NH
    fig.add_subplot(5,4,3+i*4)
    llc.pcol(xg,yg,np.nanmean(t[0,:,:,:]-t[i+1,:,:,:],axis=0)*100,mnh,cmap=plt.get_cmap('coolwarm'),vmax=50,vmin=-50)
    mnh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff thic [cm] NH')

    # Differences in Concentration SH
    fig.add_subplot(5,4,4+i*4)
    llc.pcol(xg,yg,np.nanmean(t[0,:,:,:]-t[i+1,:,:,:],axis=0)*100,msh,cmap=plt.get_cmap('coolwarm'),vmax=50,vmin=-50)
    msh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff thic [cm] SH')


# In[53]:


# Difference plots
fig = plt.figure(figsize=(10,12))
it=9+4*12
fi = 1

for fi in range(6):
    # Differences in Concentration NH
    fig.add_subplot(6,4,1+fi*4)
    llc.pcol(xg,yg,a[fi,it,:,:]*100,mnh,cmap=plt.get_cmap('viridis'),vmax=100,vmin=0)
    mnh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff conc [%] NH')
    plt.gca().set_ylabel(labels[fi])

    # Differences in Concentration SH
    fig.add_subplot(6,4,2+fi*4)
    llc.pcol(xg,yg,a[fi,it,:,:]*100,msh,cmap=plt.get_cmap('viridis'),vmax=100,vmin=0)
    msh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff conc [%] SH')

    # Differences in Thickness NH
    fig.add_subplot(6,4,3+fi*4)
    llc.pcol(xg,yg,t[fi,it,:,:]*100,mnh,cmap=plt.get_cmap('viridis'),vmax=200,vmin=0)
    mnh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff thic [cm] NH')

    # Differences in Concentration SH
    fig.add_subplot(6,4,4+fi*4)
    llc.pcol(xg,yg,t[fi,it,:,:]*100,msh,cmap=plt.get_cmap('viridis'),vmax=200,vmin=0)
    msh.drawcoastlines()
    plt.colorbar()
    plt.gca().set_title('diff thic [cm] SH')


# In[46]:


# Time series plots
fig,ax = plt.subplots(4,1,sharex=True,figsize=(10,7))

mask_nh.shape
np.rollaxis(np.rollaxis((a[0,:,:,:]*rac),1,0),2,1)[mask_nh].shape

for i in range(6):
    # Extent NH
    ax[0].plot(dates,np.rollaxis(np.rollaxis((a[i,:,:,:]*rac),1,0),2,1)[mask_nh].sum(axis=0),label=labels[i],linestyle='-')

    # Extent SH
    ax[1].plot(dates,np.rollaxis(np.rollaxis((a[i,:,:,:]*rac),1,0),2,1)[mask_sh].sum(axis=0),label=labels[i],linestyle='-')

    # Volume NH
    ax[2].plot(dates,np.rollaxis(np.rollaxis((t[i,:,:,:]*rac),1,0),2,1)[mask_nh].sum(axis=0),label=labels[i],linestyle='-')

    # Volume SH
    ax[3].plot(dates,np.rollaxis(np.rollaxis((t[i,:,:,:]*rac),1,0),2,1)[mask_sh].sum(axis=0),label=labels[i],linestyle='-')

# Extent NH
ax[0].set_ylabel('extent NH')

# Extent SH
ax[1].set_ylabel('extent SH')

# Volume NH
ax[2].set_ylabel('volume NH')

# Volume SH
ax[3].set_ylabel('volume SH')
ax[3].legend()


# In[17]:


#
ds = xmitgcm.open_mdsdataset(dir_th, geometry='llc')


# In[55]:


plt.figure()
llc.pcol(xg,yg,np.nanmean(t[1,:,:,:]-t[2,:,:,:],axis=0)*100,msh,cmap=plt.get_cmap('coolwarm'),vmax=50,vmin=-50)
msh.drawcoastlines()
plt.colorbar()

np.all(t[1,:,:,:]==t[2,:,:,:])
