"""
Quick temperature scan
"""
Load(Filename='/SNS/SEQ/IPTS-6057/data/SEQ_23397_event.nxs',OutputWorkspace='SEQ_23397_event',LoadMonitors='1')
ConvertUnits(InputWorkspace='SEQ_23397_event',OutputWorkspace='d_event',Target='dSpacing')
DeleteWorkspace('SEQ_23397_event')
Thi=np.arange(250.0,0,-10.0)
Tlo=Thi-10.0
for idx in range(len(Thi)):
	FilterByLogValue(InputWorkspace='d_event',OutputWorkspace='tempT',LogName='SampleTemp',MinimumValue=str(Tlo[idx]),MaximumValue=str(Thi[idx]),TimeTolerance='0.0089999999999999993',LogBoundary='Left')
	Rebin(InputWorkspace='tempT',OutputWorkspace='tempT',Params='0.5,0.005,7',PreserveEvents='0')
	OWS_name=T_%d_%d %(Tlo[idx],Thi[idx])
	SumSpectra(InputWorkspace='tempT',OutputWorkspace=OWS_name,IncludeMonitors='0')
	fname=OWS_name+'.nxs'
	SaveNexus(OWS_name,fname)
	DeleteWorkspace(OWS_name)



