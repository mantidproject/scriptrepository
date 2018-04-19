# runlist=range(20764,20795)
runlist=[20721,20722,20730,20723,20731,20724,20732,20725,20733,20726,20727,20728,20729]
wksplist=''
for run in (runlist):
	Load(Filename='//emu/data\Cycle 10_1/emu000'+str(run)+'.nxs',OutputWorkspace='raw')
	AsymmetryCalc(InputWorkspace='raw',OutputWorkspace='asym',ForwardSpectra='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48',BackwardSpectra='49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96',Alpha='1.0')
	CropWorkspace(InputWorkspace='asym',OutputWorkspace='trim',XMin=0.1,XMax=33.0)
#	Rebunch(InputWorkspace='trim',OutputWorkspace='bunch',NBunch='5')
	# c1=0.2276 usually
	ScaleX(InputWorkspace='trim',Factor=-0.113,Operation='Add',OutputWorkspace='corr')
	Rebin(InputWorkspace='corr',Params='0,0.016,128.0',OutputWorkspace='corr') # zero pad at end
	ExponentialCorrection(InputWorkspace='corr',OutputWorkspace='corr',C1='0.22',Operation='Multiply')
	RealFFT(InputWorkspace='corr',OutputWorkspace='freq'+str(run),IgnoreXBins='1')
	if(wksplist != ''):
		wksplist=wksplist+','
	wksplist=wksplist+('freq'+str(run))
# FFT puts Modulus in spectrum 2
ConjoinSpectraNumAx(InputWorkspaces=wksplist,OutputWorkspace='fftmap0',WorkspaceIndex=2,LabelUsing='Temp_Sample',LabelValue='Mean')

