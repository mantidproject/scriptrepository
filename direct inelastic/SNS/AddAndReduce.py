"""
Simple script to add runs together, reduce using DGSReduction, and save to NXSPE format.
It assumes that the preprocessed vanadium for normalization was saved as van.nx5
"""

from numpy import *
from string import *

path='/SNS/SEQ/IPTS-7446/data/'
#runs=range(31225 , 31228) #first to last+1
runs=[31229,31230,31231,31232,31233,31234,31235,31236,31237,31238,31239,31242,31243]
filenames=''
for r in runs:
	filenames+=path+'SEQ_'+str(r)+'_event.nxs+'
filenames=filenames[:-1]

Load(Filename=filenames,OutputWorkspace='input')

config['default.facility']="SNS"
DgsReduction(
             SampleInputWorkspace='input',
             OutputWorkspace="output",
             IncidentBeamNormalisation="ByCurrent",
             DetectorVanadiumInputFile="/SNS/SEQ/IPTS-7446/shared/autoreduce/van.nx5",
            )
SaveNXSPE("output","/SNS/SEQ/IPTS-7446/shared/3y9/25K250meV.nxspe",mtd["output"].getRun()['Ei'].value,0,'1')
