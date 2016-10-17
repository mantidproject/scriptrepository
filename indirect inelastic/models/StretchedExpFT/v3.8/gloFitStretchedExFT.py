'''
  For problems with this script, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

  Script to global fit of QENS data to the Fourier transform of a stretched exponential
  This script should be run in the "Script Window" of MantidPlot
  Global parameter: beta exponent in the stretched exponential

  Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFT ) + LinearBackground
    with 0<x<1 is the fraction of the elastic intensity

  Parameter Beta of fit function StretchedExFT is the same for all spectra. All other fitting
  parameters are different for each spectrum
'''
import re
from copy import copy
from seqFitStretchedExFT import sequentialFit
import sys

"""
   Below are the variables that can be changed by the user
"""
data_name="data"  # Name of the workspace containing the QENS signal
resolution_name="resolution"  # Name of the workspace containing the resolution

# Energy range over which we do the fitting.
minE = -0.1  # Units are in meV
maxE =  0.1

# Do the fit only on these workspace indexes (Note: the index begins at zero, not one)
selected_wi = [ 0, 1, 2, 3, 4]
#selected_wi = None   # uncomment this line if your want to select all spectra

# Initial guess for the lowest Q. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f1.f0.Height' :    0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.Height' :    0.9,   # This has to be 1-f0.f1.f0.Height
              'f0.f1.f1.Tau'    : 1000.0,   # tau or relaxation time
              'f0.f1.f1.Beta'   :    0.9,   # exponent
              'f1.A0'           :    0.0,   # intercept background
              'f1.A1'           :    0.0,   # slope background
}


"""
   Beginning here, the user does not need to change anything
"""

# Get handle to workspaces
data=mtd[data_name]
resolution=mtd[resolution_name]

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
fitstring_template = """(composite=Convolution,FixResolution=true,NumDeriv=true;
name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=_IQ_,
Scaling=1,Shift=0,XScaling=1,ties=(Scaling=1,Shift=0,XScaling=1);
(name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height);
name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
constraints=(0<Tau,0<Beta);
ties=(f1.Centre=f0.Centre)));
name=LinearBackground,A0=0,A1=0"""

fitstring_template = re.sub('[\s+]', '', fitstring_template)  # remove whitespaces and such

print "\n#######################\nRunning a sequential fit to obtain a good initial guess\n#######################"
seqOutput = sequentialFit(resolution, data, fitstring_template, initguess, [minE, maxE], qvalues, selected_wi)
# Since we are going to tie parameter Beta, find the average and use this number as initial guess
betas = list()
for i in range(len(seqOutput["funcStrings"])):
    match = re.search("Beta=(\d+\.\d+)", seqOutput["funcStrings"][i])
    betas.append(float(match.groups()[0]))
average_beta = "Beta="+str(sum(betas)/len(betas))
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
    if domain_index == 0:
        spectra_domain_relation.update({"InputWorkspace": data.name(), "WorkspaceIndex": str(iq), "StartX": str(minE), "EndX": str(maxE)})
    else:
        di = str(domain_index)
        spectra_domain_relation.update({"InputWorkspace_"+di: data.name(), "WorkspaceIndex_"+di: str(iq), "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
        ties.append("f{0}.f0.f1.f1.Beta".format(domain_index))
    domain_index += 1
global_model += "ties=({0}=f0.f0.f1.f1.Beta)".format('='.join(ties))

# Carry out the fit
print "#######################\nRunning the global fit\n#######################"
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True,
    Minimizer="FABADA", MaxIterations=1000,
    **spectra_domain_relation)

# Save Q-dependencies of the optimized parameters
parameters_workspace = mtd[output_workspace+"_Parameters"]
shorties = {"f0.f1.f0.Height": "EISF", "f0.f1.f1.Tau": "tau", "f0.f1.f1.Beta": "beta"}
other = dict()
for shorty in shorties.values():
    other[shorty] = list()
other["qvalues"] = [qvalues[iq] for iq in selected_wi]
names = shorties.keys()
for row in parameters_workspace:
    # name of the fitting parameter for this particular row
    matches = re.search("^f(\d+)\.(.*)", row['Name']) # for instance, f0.f0.f1.f0.Height
    if matches:
        iq, name = matches.groups()
        if name in names:
            other[shorties[name]].append(row["Value"])
dataY = other["EISF"] + other["tau"] + other["beta"]
dataX = other["qvalues"] * 3
# Save dependence versus Q
glofit_Qdependencies = CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY,
                                       NSpec=3, WorkspaceTitle="Q-dependence of parameters",
                                       VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "tau", "beta"])
# Save dependence versus Q**2
dataX = np.array(dataX) ** 2
glofit_Q2dependencies = CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, NSpec=3,
                                        WorkspaceTitle="Q squared-dependence of parameters",
                                        VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "tau", "beta"])

