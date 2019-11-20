"""
print u and v vectors for mslice/horace. Use in conjunction with make_d_spacing.py
"""
from __future__ import print_function
import numpy as np
wksname='SingleCrystalPeakTable'
a=3.78
b=3.78
c=13.2
alpha=90
beta=90
gamma=90
CalculateUMatrix(PeaksWorkspace=wksname,a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
h=mtd[wksname]
uV=h.sample().getOrientedLattice().getuVector()
vV=h.sample().getOrientedLattice().getvVector()
print('u=%s' %(np.array([uV.X(),uV.Y(),uV.Z()])/uV.norm()))
print('v=%s'  %(np.array([vV.X(),vV.Y(),vV.Z()])/vV.norm()))
