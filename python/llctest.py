import sys
import numpy as np
import matplotlib.pyplot as plt

from MITgcmutils import rdmds
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='npstere',lon_0 = 0., boundinglat=60, resolution = 'l')
m = Basemap(projection='moll',lon_0 = 0., resolution = 'c')

d = rdmds('diag2Dm',[72319,72381,72437,72499,72559,72621,72681,72743,
                     72805,72865,72927,72987])
d3 = rdmds('diag3Dm',np.Inf)
xc,yc=m(rdmds('XC'),rdmds('YC'))
hf = rdmds('hFacC')

def llc_contourf(*arguments, **kwargs):

    arglen = len(arguments)
    h = []
    if arglen >= 3:
        data = arguments[2].flatten()
        x = arguments[0].flatten()
        y = arguments[1].flatten()
        data = np.ma.masked_array(data)
        data = np.ma.masked_where(data == 0., data)
        if data.min() != 0.:
            i = data != 0
            x = arguments[0].flatten()[i]
            y = arguments[1].flatten()[i]
            data = data[i]
            
        if arglen == 3:
            h = plt.tricontourf(x, y, data, **kwargs)
        elif arglen == 4:
            h = plt.tricontourf(x, y, data, arguments[3], **kwargs)
        
    else:
        print "wrong number of arguments"
        print "need at least x,y,fld"
        sys.exit(__doc__)

    return h

plt.clf()
#h  = llc_contourf(xc,yc,d[0,2,:,:],levels=np.arange(20)*6./20+1e-4,extend='max')
#h  = llc_contourf(xc,yc,d[0,1,:,:],levels=np.arange(20)/20.+0.05)
h  = llc_contourf(xc,yc,d3[0,0,:,:],20,extend='both')
m.drawmapboundary(fill_color='w'); 
#m.drawcoastlines(color='grey');
m.fillcontinents(color = 'grey')
m.drawmeridians(np.arange(0, 360, 30));
plt.colorbar(orientation='horizontal')
m.drawparallels(np.arange(-90, 90, 30))

plt.show()

