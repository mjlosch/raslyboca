# compute net precipitation into llc90 and into llc270 ocean

import numpy as np
from MITgcmutils import rdmds
from netCDF4 import Dataset
import os
from myutils import *

mdir = '/home/ollie/mlosch/MITgcm/MITgcm'
rac_llc90 = rdmds(os.path.join(mdir,'llc90/grid','RAC'))
hf_llc90 = rdmds(os.path.join(mdir,'llc90/grid','hFacC'),lev=[0])
rac_llc270 = rdmds(os.path.join(mdir,'llc270/grid','RAC'))
hf_llc270 = rdmds(os.path.join(mdir,'llc270/grid','hFacC'),lev=[0])

a90  = (rac_llc90*hf_llc90).sum()
a270 = (rac_llc270*hf_llc270).sum()

# precip from an 30year llc90 spinup run
bdir = '/work/ollie/mlosch/raslyboca'
m90 = 'MITgcm_elena/llc90/spinup_era_interim/monitor_exf.0000000000.t001.nc'
ds90 = Dataset(os.path.join(bdir,m90))

# precip from a 20year llc270 spinup run
m270 = 'llc270/runz00/20years/stdout.*'
globfiles = os.path.join(bdir,m270)
files=[]
if len(args)<2:
    from glob import glob
    for infile in glob(globfiles):
        files.append(infile)


def get_output (fnames, mystring):
    """parse fname and get some numbers out"""
    timev = []
    myvar = []
    if mystring[:3]=='exf':
        timename = 'exf_time_sec'
    elif mystring[:3]=='sea':
        timename = 'seaice_time_sec'
    else:
        timename = 'time_secondsf'
    for fname in fnames:
        try:
            f = open(fname)
        except:
            print(fname + " does not exist, continuing")
        else:
            for line in f:
                if timename in line:
                    ll = line.split()
                    timev.append(float(ll[-1].replace('D','e')))
                elif mystring in line:
                    ll = line.split()
                    myvar.append(float(ll[-1].replace('D','e')))

            f.close()

    # reverse order
    timevs=np.asarray(timev[::-1])
    myvars=np.asarray(myvar[::-1])
    # This sorts again in ascending order and returns the index of
    # the first occurrence of duplicates. Because we have reverted the order
    # before, in this way we use the values at the beginning of a pickup run
    # rather than the overlapping values of the previous (potentially crashed)
    # run
    timevs, isort = np.unique(timevs,return_index=True)
    myvars=myvars[isort]

    return timevs, myvars
# done

etimsec, prc = get_output(files,'exf_precip_mean')

plt.clf()
plt.plot(ds90['exf_time_sec'][:],ds90['exf_precip_mean'][:])
plt.plot(etimsec,prc)
