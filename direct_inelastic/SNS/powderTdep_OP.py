"""
Powder line intensity - temperature dependence
"""
from numpy import *
w=Load('/SNS/SEQ/IPTS-4273/data/SEQ_10235_event.nxs+/SNS/SEQ/IPTS-4273/data/SEQ_10236_event.nxs')
tomask=[]

angmin=19
angmax=21

for i in range(w.getNumberHistograms()):
	ang=degrees(w.detectorTwoTheta(w.getDetector(i)))
	if (ang<angmin) or (ang>angmax):
		tomask.append(i)
		
MaskDetectors(w,tomask)

GenerateEventsFilter(InputWorkspace='w',OutputWorkspace='temp_split',InformationWorkspace='info',UnitOfTime='Nanoseconds',LogName='SampleTemp',MinimumLogValue='7',MaximumLogValue='300',LogValueInterval='5')
FilterEvents(InputWorkspace='w',OutputWorkspaceBaseName='split',InformationWorkspace='info',SplitterWorkspace='temp_split',FilterByPulseTime='1',GroupWorkspaces='1')
DeleteWorkspace('split_unfiltered')
sum=SumSpectra('split')
sum=NormaliseByCurrent(sum)

temperature=[]
intensity=[]
error=[]

for name in sum.getNames():
	si=mtd[name]
	if si.run().getProtonCharge()>1:
		intensity.append(si.extractY()[0][0])
		error.append(si.extractE()[0][0])
		temperature.append(si.run()['SampleTemp'].getStatistics().mean)

final=CreateWorkspace(temperature,intensity,error)

plotSpectrum(final,0,1)
SaveAscii(final,'/SNS/SEQ/IPTS-4273/shared/int_OP.dat')
