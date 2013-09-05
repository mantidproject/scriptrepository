"""
Generate a list of runs where some log values obey certain criteria. Input can be from log.py 
"""
from string import *
from numpy import *
f=open('/SNS/SEQ/IPTS-6179/shared/2K_50meV_complete/complete_log_2K_50meV.txt','r')
info=f.read().split()
f.close()
info=info[4:]
listsize=array(info).size/4
runs=arange(listsize)
temps=arange(listsize)
angs=arange(listsize)
for i in arange(listsize):
	runs[i]=info[4*i]
	temps[i]=float(info[4*i+1])
	angs[i]=int(round(float(info[4*i+2])))
runlist=[]
for i in arange(-90,39,1):
	runlist.append(list(runs[where((angs==i))]))#(temps>70) &
print runlist
