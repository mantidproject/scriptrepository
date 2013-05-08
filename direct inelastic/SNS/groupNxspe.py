"""
Reduce file size of NXSPE files by grouping pixels togeter. 
"""
from string import *
from numpy import *

# 40 meV
r_0 = [33008,33009]  
r_180 = [33383,33384]
# -0 < Omega CCR < -90
r1 = range(33012,33193)
r2 = range(33385,33566)
# -30 < Omega CCR < -60
r3 = range(33573,33635)
# -90 < Omega CCR < -180
r4 = range(33196,33377)
All_runs = r1+r2+r3+r4 + r_0 + r_180

# 50 meV
#All_runs = [33194,33195,33381,33382,33572,33568,33571,33569]


for run in All_runs:
	filename='/SNS/SEQ/IPTS-7672/shared/autoreduce/SEQ_'+str(run)+'_autoreduced.nxspe'
        outname='/SNS/SEQ/IPTS-7672/shared/40meV/Empty_CanSubtracted/SEQ_'+str(run)+'.nxspe'
	w=LoadNXSPE(filename)
	LoadInstrument(w,Filename='/opt/Mantid/instrument/SEQUOIA_Definition.xml')
	w1=GroupDetectors(w,'/SNS/SEQ/IPTS-7672/shared/4x2_groupDet.xml',Behaviour='Average')
	efixed=w.getRun()['Ei'].value
	psi=w.getRun()['psi'].value
	SaveNXSPE(w1,outname,efixed,psi,1)

