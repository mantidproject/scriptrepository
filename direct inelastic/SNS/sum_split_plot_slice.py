"""
Simple script for continuous rotations
"""
Load(Filename='HYS_13656-13674',OutputWorkspace='sum')
FilterByLogValue(InputWorkspace='sum',OutputWorkspace='sum1',LogName='s1',MinimumValue='0',MaximumValue='44.5',LogBoundary='Left')
GenerateEventsFilter(InputWorkspace='sum1',OutputWorkspace='splboth',InformationWorkspace='info',UnitOfTime='Nanoseconds',LogName='s1',MaximumLogValue='44.5',LogValueInterval='2')	
FilterEvents(InputWorkspace='sum1',OutputWorkspaceBaseName='split',InformationWorkspace='info',SplitterWorkspace='splboth',FilterByPulseTime='1',GroupWorkspaces='1')				
DeleteWorkspace('split_unfiltered')
CompressEvents('split',0.1,OutputWorkspace='splitc')
grouping_file='/SNS/HYS/IPTS-8361/shared/3y9/group4x1.xml'
DgsReduction(SampleInputWorkspace='splitc',IncidentBeamNormalisation='ByCurrent',OutputWorkspace='reduced',GroupingFile=grouping_file,TimeIndepBackgroundSub ='1',TibTofRangeStart =10400,TibTofRangeEnd =12400,IncidentEnergyGuess=50)
SetGoniometer('reduced',Axis0="s1,0,1,0,1")
SetUB('reduced',5.823,6.475,3.186,90,90,90,'0,1,0','0,0,1')
ConvertToMD(InputWorkspace='reduced',OutputWorkspace='md',QDimensions='Q3D',QConversionScales='HKL',MinValues='-0.5,-3,-5,-10',MaxValues='0.5,6,2,45')
MergeMD(InputWorkspaces='md',OutputWorkspace='merged')
BinMD(InputWorkspace='merged',AxisAligned='0',BasisVector0='[H,0,0],in 1.079 A^-1,1,0,0,0',BasisVector1='[0,K,0],in 0.97 A^-1,0,1,0,0',BasisVector2='[0,0,L],in 1.972 A^-1,0,0,1,0',BasisVector3='DeltaE,DeltaE,0,0,0,1',OutputExtents='-3,3,-2,6,-4,-1.5,-3,3',OutputBins='1,100,100,1',Parallel='1',OutputWorkspace='slice')
sv=plotSlice("slice", xydim=[1,2], slicepoint=[0.0,0,0] ) 
sv.setColorScaleMin(0)
sv.setColorScaleMax(10)
