import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.animation import FuncAnimation
import os
from mpl_toolkits.basemap import Basemap

year = 2018

fields = ['u10','v10','q','ssrd','strd','t2m']

def readfield(fname,dims):
    import sys
    """Call signatures::

    readfield(filename, dims, numpy.datatype)

    Read unblocked binary data with dimentions "dims".
    """

    try:
        fid = open(fname,"rb")
    except:
        sys.exit( fname+": no such file or directory")
    else:
        v   = np.fromfile(fid,'>f4')
        fid.close()

        if   len(v) == np.prod(dims):     v = v.reshape(dims)
        elif len(v) == np.prod(dims[1:]): v = v.reshape(dims[1:])
        else:
            errstr = (  "dimensions do not match: \n len(data) = " + str(len(v))
                        + ", but prod(dims) = " + str(np.prod(dims)) )
            raise RuntimeError(errstr)

    return v

# erai grid
rEarth = 6371e3
lon = np.linspace(0.,360.,480)
lat = np.linspace(-90.,90.,241)
dlat = 0.75*np.pi/180.*np.cos(lat*np.pi/180.)*rEarth
dlat[0],dlat[-1]=dlat[1],dlat[-2]
dlon = 0.75*np.pi/180.*np.ones(480,)*rEarth
dx,dy=np.meshgrid(dlon,dlat)
lon2,lat2=np.meshgrid(lon,lat)
rac=dx*dy

nt = 365*4
if year==2019: nt = 244*4
u=readfield('u10_ERAi_6hourly_'+str(year),[nt,241,480])
v=readfield('v10_ERAi_6hourly_'+str(year),[nt,241,480])

# estimate divergence and vorticity:
dudi = (u-np.roll(u,-1,axis=-1))/dx
dvdj = (v-np.roll(v,-1,axis=-2))/dy
dudj = (u-np.roll(u,-1,axis=-2))/dy
dvdi = (v-np.roll(v,-1,axis=-1))/dx

div = dudi+dvdj
vor = dvdi-dudj
# div=readfield('q_ERAi_6hourly_'+str(year),[nt,241,480])
# vor=readfield('t2m_ERAi_6hourly_'+str(year),[nt,241,480])

# we need an arctic map

# the first plot is only created to set up a map, which is rotated in the
# wrong way and grab the corner point coordinates
lon0 = -160
lat0 = 83
m =Basemap(projection='stere',lat_0 = lat0, lon_0 = lon0,
           width = 4e6, height =4e6,
           resolution = 'l')
f0, a0 =plt.subplots(1,1)
m.fillcontinents()
# inverse coordinate transformation to arrive at units of meters
urc = m(a0.get_xlim()[0],a0.get_ylim()[0],inverse=True)
llc = m(a0.get_xlim()[1],a0.get_ylim()[1],inverse=True)
print('coordinates of upper right corner ({}, {})'.format(urc[0],urc[1]))
print('coordinates of lower left  corner ({}, {})'.format(llc[0],llc[1]))
plt.close(f0)

# now use the flipped corner coordinates to set up the proper projection
m =Basemap(projection='stere',lat_0 = lat0, lon_0 = lon0,
           llcrnrlon = llc[0], llcrnrlat = llc[1],
           urcrnrlon = urc[0], urcrnrlat = urc[1],
           resolution = 'l')

x,y=m(lon2,lat2)

k = 0

f1, ax =plt.subplots(1,2)
csf0=ax[0].pcolormesh(x,y,div[k,:,:],vmin=-1e-3,vmax=1e-3)
csf1=ax[1].pcolormesh(x,y,vor[k,:,:],vmin=-1e-3,vmax=1e-3)
plt.colorbar(csf0,ax=ax[0],orientation='horizontal')
plt.colorbar(csf1,ax=ax[1],orientation='horizontal')
th0 = ax[0].set_title('wind divergence')
th1 = ax[1].set_title('wind vorticity')

for a1 in ax:
    m.fillcontinents(ax=a1)
    m.fillcontinents(ax=a1)
    m.drawparallels(np.arange(50.,80,10.),ax=a1)
    m.drawmeridians(np.arange(0.,360.,30.),ax=a1)


refdate = datetime.datetime(year,1,1)

def animate(t):
    csf0.set_array(div[t,:,:].ravel())
    csf1.set_array(vor[t,:,:].ravel())
    mydate = refdate + datetime.timedelta(days=t/4)
    th0.set_text('wind divergence at '+mydate.strftime("%b %d %Y %H:%M:%S"))
    th1.set_text('wind vorticity at '+mydate.strftime("%b %d %Y %H:%M:%S"))

anim = FuncAnimation(f1, animate, interval=10, frames=range(250*4,290*4))

plt.show()

#f1.show()
