# ERAiterim data ends on Aug 31, 2019. In order to simulate to the end of this
# we need an extra forcing record. We just add Aug 31, 2019 a second time
# and pretend it is Sep 01, 2019.
# this script can only be run one time, because it assumes a fixed
# time record length that is then changed in this script.

import numpy as np
import os



ncdir = '/work/ollie/clidyn/forcing/erai/raw'
input_variables = ['q','strd','ssrd','tp','t2m', 'u10', 'v10']
#input_variables = ['tp','t2m', 'u10', 'v10']
#input_variables = ['q']

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

def writefield(fname,data):
    import sys
    print('writing '+fname)
    if True:
        print('not writing anything for safety')
        pass
    else:
        fid = open(fname,"wb")
        data.astype('>f4').tofile(fid)
        fid.close()

for v in input_variables:
    fname = v+'_ERAi_6hourly_2019'
    fname_old = fname+'_orig'

    print(v)
    fld = readfield(fname,[243*4,241,480])
    writefield(fname_old,fld)

    print(fld.shape)
    fld = np.append(fld,fld[-4:,:,:],axis=0)
    print(fld.shape)

    writefield(fname,fld)
