'''
  For problems with this script, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

  Script to global fit of QENS data to the Fourier transform of a stretched exponential
  This script should be run in the "Script Window" of MantidPlot
  Global parameter: beta exponent in the stretched exponential

  Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFT ) + LinearBackground + BackgroundFile
    with 0<x<1 is the fraction of the elastic intensity

  Parameter Beta of fit function StretchedExFT is the same for all spectra. All other fitting
  parameters are different for each spectrum
'''
import re
from copy import copy
import numpy as np
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from seqFitStretchedExFT import sequentialFit

"""
   Below are the variables that can be changed by the user
"""
data_name="data"  # Name of the workspace containing the QENS signal
resolution_name="resolution"  # Name of the workspace containing the resolution
background_name="background"  # Name of the workspace containing the background

# Energy range over which we do the fitting.
minE = -0.1  # Units are in meV
maxE =  0.1

# Do the fit only on these workspace indexes (Note: the index begins at zero, not one)
#selected_wi = [ 1, 2, 3]
selected_wi = None   # uncomment this line if your want to select all spectra

# Initial guess for the HIGHEST Q of the sequential fit. A guess can be obtained by
# running MantidPlot interactively
initguess = { 'f0.f1.f0.Height' :    0.0,   # intensity due to elastic line
              'f0.f1.f1.Height' :    0.4,   # intensity of the stretched exponential
              'f0.f1.f1.Tau'    :   25.0,   # relaxation time
              'f0.f1.f1.Beta'   :    0.6,   # stretching exponent
}

# Option to select an initial value of stretching exponent beta for the global fit. 
# If "None" is selected, the initial value will be average over all the values
# that beta takes during the sequential fit of the spectra
#initial_beta=None
initial_beta=0.61

# Settings for the minimizer. See the "Fit" algorithm in the documentation
minimizer="FABADA"  # slow, but more reliable than "Levenberg-Marquardt"
maxIterations=1000

"""
   Beginning here, the user does not need to change anything
"""

# Get handle to workspaces
data=mtd[data_name]
resolution=mtd[resolution_name]
background=mtd[background_name]

# Extract Q-values
vertical_axis = data.getAxis(1)
qvalues = vertical_axis.extractValues()
if not selected_wi:
    selected_wi=range(len(qvalues))

""" Below is the model cast as a template string suitable for the
    Fit algoritm of Mantid. You can obtain similar string by setting up a model
    in the "fit wizard" of MantidPlot and then
    "Manage Setup" --> "Copy to Clipboard". This actions will save the model
    as a string which you can later paste onto this script.
"""
fitstring_template ="""
    (composite=Convolution,FixResolution=false,NumDeriv=true;
        name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=_IQ_,
            Scaling=1,Shift=0,XScaling=1,ties=(Scaling=1,XScaling=1);
        (
         name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height),
             ties=(Centre=0);
         name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
             ties=(Centre=0)
        );
    );
    name=LinearBackground,A0=0.0,A1=0.0;
    name=TabulatedFunction,Workspace=_BACKGROUND_,WorkspaceIndex=_IQ_,Scaling=1,Shift=0,XScaling=1,ties=(XScaling=1)"""
fitstring_template = re.sub('[\s+]', '', fitstring_template)  # remove whitespaces and such

print "\n#######################\nRunning a sequential fit to obtain a good initial guess\n#######################"
seqOutput = sequentialFit(resolution, data, background, fitstring_template, initguess, [minE, maxE], qvalues, selected_wi)

# Since we are going to tie parameter Beta, find the average and use this number as initial guess
betas = list()
for i in range(len(seqOutput["funcStrings"])):
    match = re.search("Beta=(\d+\.\d+)", seqOutput["funcStrings"][i])
    betas.append(float(match.groups()[0]))
average_beta = "Beta="+str(sum(betas)/len(betas))
if initial_beta:
    average_beta="Beta="+str(initial_beta)
for i in range(len(seqOutput["funcStrings"])):
    seqOutput["funcStrings"][i] = re.sub("Beta=\d+\.\d+",average_beta, seqOutput["funcStrings"][i])

print seqOutput["funcStrings"]

# Merge models for each spectra
global_model= 'composite=MultiDomainFunction,NumDeriv=true;'
for funcString in seqOutput["funcStrings"]:
    global_model += "(composite=CompositeFunction,NumDeriv=true,$domains=i;{0});\n".format(funcString)

# Insert the ties for beta parameter and relate spectra to domains
ties = list()
spectra_domain_relation = dict()
domain_index = 0
for iq in selected_wi:
    di = str(domain_index) 
    if domain_index == 0:
        spectra_domain_relation.update({"InputWorkspace": data.name(), "WorkspaceIndex": str(iq), "StartX": str(minE), "EndX": str(maxE)})
        #spectra_domain_relation.update({"InputWorkspace_"+di: data.name(), "WorkspaceIndex_"+di: str(iq), "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
    else:
        spectra_domain_relation.update({"InputWorkspace_"+di: data.name(), "WorkspaceIndex_"+di: str(iq), "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
        ties.append("f{0}.f0.f1.f1.Beta".format(domain_index))
    domain_index += 1   
global_model += "ties=({0}=f0.f0.f1.f1.Beta)".format('='.join(ties))

# Carry out the fit
print "#######################\nRunning the global fit\n#######################"
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True,
    Minimizer=minimizer, MaxIterations=maxIterations,
    **spectra_domain_relation)

# Save Q-dependencies of the optimized parameters
parameters_workspace = mtd[output_workspace+"_Parameters"]
other = dict()
other_error = dict()
other["qvalues"] = [qvalues[iq] for iq in selected_wi]
other_error["qvalues"] = [0.0,]*len(selected_wi)
for row in parameters_workspace:
    name = row["Name"]
    matches = re.search("^f(\d+)\.(.*)", row['Name']) # for instance, f3.f0.f1.f0.Height
    if matches:
        iq, name = matches.groups()  # for instance, 3 and f0.f1.f0.Height
        if name in other.keys():
            other[name].append(row["Value"])
            other_error[name].append(row["Error"])
        else:
            other[name]=[row["Value"],]
            other_error[name]=[row["Error"],]
    if name=="Cost function value":
        chi2 = row["Value"]

# Calculate the EISF and associated error with error propagation
a=np.array(other["f0.f1.f0.Height"]); b=np.array(other["f0.f1.f1.Height"])
ae=other_error["f0.f1.f0.Height"];  be=other_error["f0.f1.f1.Height"]
other["EISF"] = (a/(a+b)).tolist()
other_error["EISF"]=(np.sqrt(b*ae*ae+a*be*be)/(a+b)).tolist()  # Used error propagation formula

# Save dependence versus Q
dataY = other["EISF"] + other["f0.f1.f1.Tau"] + other["f0.f1.f1.Beta"]
dataX = other["qvalues"] * 3
glofit_Qdependencies = CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY,
                                       NSpec=3, WorkspaceTitle="Q-dependence of parameters",
                                       VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "tau", "beta"])
# Save dependence versus Q**2
dataX = np.array(dataX) ** 2
glofit_Q2dependencies = CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, NSpec=3,
                                        WorkspaceTitle="Q squared-dependence of parameters",
                                        VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "tau", "beta"])
print "Chi-square of global fit=",chi2


