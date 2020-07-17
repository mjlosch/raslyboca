#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol

d=rdmds('diag3Dm',[120550, 120554])
print d.shape

df=llc.faces(d[:,:,:,:,:])

print len(df)
#print len(df[0])
for k in range(len(df)): print k, df[k].shape
