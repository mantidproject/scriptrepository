# Load single crystl white beam data from MERLIN - Sr3 Fe2 O7 sample

LoadRaw(Filename='MER16000.raw', OutputWorkspace='MER16000')

# Lattice parameters input

a=3.85
b=3.85
c=20.15

u=UnitCell(a,b,c,90,90,90)

# Print d spacing or angle between reflections - examples:

print u.d(1,1,0)
print u.d(1,0,0)

print u.recAngle(1,1,0,1,0,0)

# If you want to convert to d spacing : 
#ConvertUnits(InputWorkspace='MER16000', OutputWorkspace='MER16000d', Target='dSpacing')

# To create a table of peask, go to instrument view of the data set in tof and create a table with peak positions picking peaks using the tab Pick-->add single crystal peak
# Zoom and pick some of them. The mini-plot will display the spectrum of this detector. 
#Click on the mini-plot at the position of the peak in time-of-flight. 
#If it is the first peak Mantid will create a PeaksWorkspace named "SingleCrystalPeakTable" and add a peak to it. 
#A peak marker will be displayed over the instrument. Repeating these steps will add new peaks to the same table and display them. Then find the UBmatrix:

FindUBUsingLatticeParameters('SingleCrystalPeakTable',3.85,3.85,20.15,90,90,90,NumInitial=4,Tolerance=0.3)

# Numinitial : selects only the first n peaks in the peak table (4)
# Tolerance of Q: 0.3 is generous

#Then run and use PredictPeaks:  leave HKLPeakworkspce empty; 

PredictPeaks('SingleCrystalPeakTable',WavelengthMin=1, WavelengthMax=5, MinDSpacing=0.8, MaxDSpacing=10, OutputWorkspace='Prediction')

# drug and drop the output workspace 'Prediction' 
# on the Instrument view and it will show the predicted peaks

# You can save the UB matrix 

SaveIsawUB('SingleCrystalPeakTable','UBmatrix')

# Or copy the UB matrix and attach it to the workspace using "

CopySample(InputWorkspace='SingleCrystalPeakTable', OutputWorkspace='MER16000', CopyEnvironment=False, CopyOrientationOnly=True)

# First input is the 'SingleCrystalPeakTable' 
# and the second is the sample run you want to copy the UB matrix on

#To convert UB to u,v (use the script interpreter window):

wk=mtd['MER16000']
uVector=wk.sample().getOrientedLattice().getuVector()
vVector=wk.sample().getOrientedLattice().getvVector()
print uVector,vVector

# You can set different axix of the goniometer using the command: SetGoniometer(Workspace,name, x,y,z, 1/-1 (1 for ccw, -1 for cw rotation)). A number of degrees can be used instead of name. 
#X, Y, Z components of the vector of the axis of rotation. Right-handed coordinates with +Z=beam direction; +Y=Vertically up (against gravity); +X to the left.
 #The sense of rotation as 1 or -1: 1 for counter-clockwise, -1 for clockwise rotation.

SetGoniometer(Workspace='MER16000', Axis0='10,0,1,0,1')

# Now predict the peaks for the new rotation

PredictPeaks(InputWorkspace='MER16000', CalculateStructureFactors=True, OutputWorkspace='rotated')

# To predict the position of magnetic peaks you can use PredictFractionalPeak 

PredictFractionalPeaks(Peaks='Prediction',FracPeaks='PredictionFrac',HOffset='0.5',KOffset='0',LOffset='0')

# To convert in HKL space use the workspace with the UB matrix and ConvertToDiffractionMDWorkspace output HKL and look with Slice Viewer if you are able to!!!
# To look with Slice Viewer, unzoom all, then enable rebin and auto rebin and you can overlay peaks usig the tab on Slice Viewer View-->Peak Overlay

ConvertToDiffractionMDWorkspace(InputWorkspace='MER16000', OutputWorkspace='MER16000_HKL', OutputDimensions='HKL')



