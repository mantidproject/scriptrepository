# Some useful comments
ws = CreateMDWorkspace(Dimensions=2, Extents=[-10,10,-10,10], Names='A,B', Units='U,U')
FakeMDEventData(InputWorkspace=ws, PeakParams='100000,-5,0,1')
FakeMDEventData(InputWorkspace=ws, PeakParams='100000,5,0,1')
histogram = BinMD(InputWorkspace='a', AlignedDim0='A,-10,10,100', AlignedDim1='B,-10,10,100', OutputExtents='-10,10,-10,10,-10,10', OutputBins='10,10,10')
#plotSlice(histogram)
smoothed = SmoothMD(InputWorkspace=histogram, WidthVector=5, Function='Hat')
#plotSlice(smoothed)
