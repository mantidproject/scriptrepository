# runlist=range(20764,20795)
runlist=list(range(130292,130308))
wksplist=''
for run in (runlist):
	Load(Filename='//hifi/data/hifi00'+str(run)+'.nxs',OutputWorkspace='raw')
	AsymmetryCalc(InputWorkspace='raw',OutputWorkspace='asym',ForwardSpectra=list(range(3,17))+list(range(49,65)),BackwardSpectra=list(range(17,49)),Alpha='1.0')
	CropWorkspace(InputWorkspace='asym',OutputWorkspace='trim',XMin=3.9,XMax=33.0)
#	Rebunch(InputWorkspace='trim',OutputWorkspace='bunch',NBunch='5')
	# c1=0.2276 usually
	ScaleX(InputWorkspace='trim',Factor=-3.7,Operation='Add',OutputWorkspace='corr')
	#Rebin(InputWorkspace='corr',Params='0,0.016,128.0',OutputWorkspace='corr') # zero pad at end
	ExponentialCorrection(InputWorkspace='corr',OutputWorkspace='corr',C1='0.22',Operation='Multiply')
	RealFFT(InputWorkspace='corr',OutputWorkspace='freq'+str(run),IgnoreXBins='1')
	if(wksplist != ''):
		wksplist=wksplist+','
	wksplist=wksplist+('freq'+str(run))
# FFT puts Modulus in spectrum 2
ConjoinSpectraNumAx(InputWorkspaces=wksplist,OutputWorkspace='fftmap0',WorkspaceIndex=2,LabelUsing='sample_magn_field',LabelValue='Mean')
