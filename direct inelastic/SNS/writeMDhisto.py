"""
Saves an md histo workspace with wname (2D slice from SliceViewer) into an ascii file
"""
from numpy import *
from string import *

wname='nxspe_md_rebinned'
filename='/home/3y9/Desktop/file.csv'

w=mtd[wname]
errsq=w.getErrorSquaredArray()
signal=w.getSignalArray()
nev=w.getNumEventsArray()
data=signal/nev
err=sqrt(errsq)/nev

d0=w.getDimension(0)
d1=w.getDimension(1)

f = open(filename, 'w')
f.write(d0.getName()+', '+str(d0.getMinimum())+', '+str(d0.getMaximum())+', '+str(d0.getNBins())+'\n')
f.write(d1.getName()+', '+str(d1.getMinimum())+', '+str(d1.getMaximum())+', '+str(d1.getNBins())+'\n')
for i in range(d1.getNBins()):
	s=array_str(data[:,i],1000)
	f.write( ','.join(s.strip('[').strip(']').split())+'\n')
f.close()
