"""Model:  Resolution x (Delta + Lorentzian + Lorentzian) + LinearBackground
 - Sequential fitting
 - Plot several quantities
"""
from __future__ import print_function

import re
from os import path
import numpy
import scipy.stats as stc
import matplotlib.pyplot as plt

data_dir="/projects/development/QENSmodelling/benchmark"  #substitute with your own directory
#data_dir='/Users/far/Documents/ORNL/GROUP LEADER Duties/BASIS/2016-03-30 Data'

resolution_run=53554
signal_run=53425

# Make sure Peak functions are not set to zero for high energies
mantid.config["curvefitting.peakRadius"]="1000"

#Load resolution and signal files (reduced files)
#resolution=Load(Filename=path.join(data_dir,"BASIS_"+str(resolution_run)+"_sqw.nxs"))
#signal=Load(path.join(data_dir,"BASIS_"+str(signal_run)+"_sqw"))
resolution=LoadDaveGrp(Filename=path.join(data_dir,"BASIS_"+str(resolution_run)+"_1run_divided.dat"), isMicroEV=1)
signal=LoadDaveGrp(path.join(data_dir,"BASIS_"+str(signal_run)+"_1run_divided.dat"), isMicroEV=1)

#Define the model:  Resolution x (Delta + Lorentzian + Lorentzian) + LinearBackground
function_string=\
"""(composite=Convolution,FixResolution=true,NumDeriv=true;
    name=TabulatedFunction,Workspace=resolution,WorkspaceIndex=0,Scaling=1,Shift=0,XScaling=1,ties=(XScaling=1);
     (name=DeltaFunction,Height=1,Centre=0;
      name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0,constraints=(0<Amplitude,0<FWHM);
      name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0,constraints=(0<Amplitude,0<FWHM);
      ties=(f0.Centre=f1.PeakCentre=f2.PeakCentre)

     );
    name=LinearBackground,A0=0,A1=0
   )"""
   
function_string=re.sub(r"[\s\t\n]", "", function_string)  #remove spaces,tabs,newlines
print(function_string)
#Sequential fitting (more info in http://docs.mantidproject.org/nightly/algorithms/PlotPeakByLogValue-v1.html)
inputs=';'.join([signal.name() + ',i%d' % i for i in range(signal.getNumberHistograms())])
fitSeq=PlotPeakByLogValue(inputs,
                          OutputWorkspace="fitSeq",
                          Function=function_string,
                          startX=-0.12, endX=0.52,
                          Fittype="Sequential",
                          passWSIndexToFunction=1,
                          CreateOutput=1,
                          OutputCompositeMembers=1,
                          ConvolveMembers=1)

##################### Fitting has ended here. Some analysis now #####################

#Find the Area, FWHM, skew, and kurtosis for each spectrum of the signal, then plot
dE = signal.dataX(0)[1]-signal.dataX(0)[0]  #energy bin
s=list()
for i in range(signal.getNumberHistograms()):
    y = signal.dataY(i)  # intensities along the energy for spectrum with workspace index "i"
    s.append([dE*y.sum(), 2.355*stc.moment(y,2), stc.skew(y), stc.kurtosis(y)])
s=numpy.array(s).transpose()
qvalues=0.1+signal.getAxis(1).extractValues()[:-1]
#for i in range(4): plot(qvalues,s[i]) #plot mean, FWHM, swew, and kurtosis, in this order

"""Viewing the results:
The sequential fitting produces the following workspaces:
  fitSeq_Workspaces: a set of data workspaces containing the curves (data,model,residuals,components) for each spectrum
  fitSeq_Parameters: a set of tables containin the optimized fitting parameters for each spectrum
  fitSeq: a composite table of all tables of the previous spectrum
"""

#Plotting one of the fits:
i=8 #the last index (remember they start from zero, not one)
workspace_name="signal_"+str(i)+"_Workspace"
plotSpectrum(workspace_name,indices=[0,1,2]) #plot data,fit, and residuals

#Quick plot together the fits for all Q's (this is always the same, could be packaged into some function)
ws=ExtractSpectra("signal_0_Workspace",StartWorkspaceIndex=0,EndWorkspaceIndex=1)
for i in range(1, signal.getNumberHistograms()):
    ws2 = ExtractSpectra("signal_"+str(i)+"_Workspace",StartWorkspaceIndex=0,EndWorkspaceIndex=1)
    ws = AppendSpectra(ws,ws2)
plotSpectrum(ws, indices=list(range(ws.getNumberHistograms())),waterfall=1)

#Prepare output for plotting the parameters (this is always the same, could be packaged into some function)
column_names=fitSeq.getColumnNames()  #names of all fitting parameters
print(column_names)
n_parameters = (len(column_names)-2)/2
x=fitSeq.column(0)*n_parameters
y=list()
for i in range(1,len(column_names)-1,2): y += fitSeq.column(i)
yerr=list()
for i in range(2,len(column_names)-1,2): yerr += fitSeq.column(i)
print(n_parameters,len(x),len(y),len(yerr))
parms=CreateWorkspace(x,y,yerr,NSpec=n_parameters,UnitX="MomentumTransfer")

#Plot the FWHM's
for parm_name in ("f1.f1.FWHM", "f1.f2.FWHM"):
    plotSpectrum(parms, (column_names.index(parm_name)-1)/2)

#What is the EISF?
EISF = numpy.array(fitSeq.column("f1.f0.Height")) # intensity of the elastic line
intensity = numpy.copy(EISF) # this will be the total intensity
for parm_name in ("f1.f1.Amplitude", "f1.f2.Amplitude"):
    intensity += numpy.array(fitSeq.column(parm_name))
EISF = EISF/intensity #EISF as a fraction of the total intensity
plot(fitSeq.column("axis-1"),EISF)

"""If you want a single lorentzian instead of two, you could either:
  - Set the last of the Lorentzians to always have zero amplitude:
            name=Lorentzian,Amplitude=0,PeakCentre=0,FWHM=0,ties=(Amplitude=0)
  - or remove the last lorentzian line from the model:
            name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0;
            name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0
    becomes:
            name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0
"""
