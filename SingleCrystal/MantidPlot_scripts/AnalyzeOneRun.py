# Initial analysis to find peaks and to obtain lattice constants and UB matrix.
#
LoadEventNexusDialog(Filename = '/SNS/MANDI/IPTS-8776/data/MANDI_4362_event.nxs',
    OutputWorkspace = 'Event_ws', FilterByTofMin=500, FilterByTofMax=16500, Precount=False, 
    Enable = 'Filename, OutputWorkspace, FilterByTofMin, FilterByTofMax')
ConvertToDiffractionMDWorkspaceDialog(InputWorkspace='event_ws', OutputWorkspace='MD_ws', 
    OneEventPerBin=False, LorentzCorrection=True, SplitThreshold=50, MaxRecursionDepth=11, 
    Enable ='InputWorkspace, OutputWorkspace, LorentzCorrection')
FindPeaksMDDialog(InputWorkspace='MD_ws', OutputWorkspace='Peaks', MaxPeaks = 500, 
    Enable = 'InputWorkspace, OutputWorkspace, MaxPeaks')
FindUBUsingFFTDialog(PeaksWorkspace='Peaks', MinD=5, MaxD=13, Tolerance=0.12, Enable = 'PeaksWorkspace, MinD, MaxD, Tolerance')
IndexPeaksDialog(PeaksWorkspace='Peaks', Tolerance = 0.12, RoundHKL=True, Enable = 'PeaksWorkspace, Tolerance, RoundHKL')
SaveIsawUBDialog(InputWorkspace = 'Peaks', Enable = 'InputWorkspace')
SaveIsawPeaksDialog(InputWorkspace = 'Peaks', Enable = 'InputWorkspace')

# Create MD workspace with Lorentz correction turned off for integration.
#
ConvertToDiffractionMDWorkspaceDialog(InputWorkspace='event_ws', OutputWorkspace='MD_ws', 
    OneEventPerBin=False, LorentzCorrection=False, SplitThreshold=50, MaxRecursionDepth=11, 
    Enable ='InputWorkspace, OutputWorkspace, LorentzCorrection')
    
IntegratePeaksMDDialog(InputWorkspace = 'MD_ws', PeaksWorkspace = 'Peaks',  OutputWorkspace = 'IntegratedPeaks_MD', PeakRadius = 0.10, 
    BackgroundInnerRadius = 0.10, BackgroundOuterRadius = 0.12, Cylinder = False,
    Enable = 'PeakRadius, BackgroundInnerRadius, BackgroundOuterRadius' )
    
SaveIsawPeaksDialog(InputWorkspace = 'IntegratedPeaks_MD', Enable = 'InputWorkspace')
    
IntegratePeaksHybridDialog(InputWorkspace = 'MD_ws', PeaksWorkspace = 'Peaks', OutputWorkspace = 'IntegratedPeaks_Hybrid',
    NumberOfBins = 20, BackgroundOuterRadius = 0.15, OutputWorkspaces = 'MDHisto',
    Enable = 'NumberOfBins, BackgroundOuterRadius' )
    
SaveIsawPeaksDialog(InputWorkspace = 'IntegratedPeaks_Hybrid', Enable = 'InputWorkspace')

# BinMDDialog(InputWorkspace = 'MD_ws', OutputWorkspace = 'BinMD_ws')
# IntegratePeaksUsingClustersDialog(InputWorkspace = 'BinMD_ws', PeaksWorkspace = 'Peaks', Threshold = 10, 
    # OutputWorkspace = 'IntegratedPeaks_Cluster', OutputWorkspaceMD = 'MD_out')

