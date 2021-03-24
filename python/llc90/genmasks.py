import numpy as np
from MITgcmutils import rdmds, wrmds, llc
import os

# need to adjust this
# gdir = '../../grid'
# depth = rdmds(os.path.join(gdir,'Depth'))

nx = 90
ny = 1170
nl = 4

# four levels clockwise around the arctic
# level 1: Bering Strait
# level 2: Barents Sea
# level 3: Fram Strait
# level 4: Canadian Arctic Archipelago

masku = np.zeros((nl,ny,nx),dtype='float32')
maskv = np.zeros((nl,ny,nx),dtype='float32')

maskuf=llc.faces(masku)
maskvf=llc.faces(maskv)

# level 1: Bering Strait
maskuf[3][0,48:50,13] =  1.

# level 2: Barents Sea
maskuf[2][1,26:36,13] = -1.
maskvf[2][1,26,:13]   =  1.

# level 3: Fram Strait
maskuf[2][2,43:59,18] = -1.

# level 4: Canadian Arctic Archipelago
maskvf[2][3,63,35:38] =  1.
maskvf[2][3,71,51:59] =  1.
maskvf[2][3,72,60:64] =  1.
maskvf[2][3,76,69:72] =  1.
maskuf[2][3,74:76,69] = -1.

masku=llc.faces2mds(maskuf)
maskv=llc.faces2mds(maskvf)

wrmds('maskU_arctic_sections',masku)
wrmds('maskV_arctic_sections',masku)
