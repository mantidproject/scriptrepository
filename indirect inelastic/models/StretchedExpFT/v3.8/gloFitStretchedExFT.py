'''
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

# Data directory (update with your own)
datadir = "/home/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/models/StretchedExpFT"
data_davegroup_filename = "data.dat"
resolution_davegroup_filename = "resolution.dat"

""" Below is the model cast as a template string suitable for the
    Fit algoritm of Mantid. You can obtain similar string by setting up a model
    in the "fit wizard" of MantidPlot and then
    "Manage Setup" --> "Copy to Clipboard". This actions will save the model
    as a string which you can later paste onto this script.
"""
fitstring_template = """(composite=Convolution,FixResolution=false,NumDeriv=true;
  name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=_IQ_,
  Scaling=f0.f0.Scaling,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1);
  (name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height<1);
   name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
   constraints=(0<Tau,0<Beta);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=f1.A0,A1=f1.A1"""

fitstring_template = re.sub('[\s+]', '', fitstring_template)  # remove whitespaces and such

# Load the data. We assume the format is DAVE group file.
#  Use the "LoadDaveGrp" algorithm
data = msapi.LoadDaveGrp( Filename = '{0}/{1}'.format(datadir, data_davegroup_filename),
                          OutputWorkspace = 'data', XAxisUnits = 'DeltaE', IsMicroEV = 1 )

# Alternatively, we use the "LoadNexus" algorithm if we have the reduced data
# as a Nexus file
# LoadNexus( Filename = '{0}/BASIS_17706_1run_divided.nxs'.format( datadir ),
#             OutputWorkspace = 'data' )

# Load the resolution
resolution = LoadDaveGrp( Filename = '{0}/resolution.dat'.format( datadir ),
    OutputWorkspace = 'resolution',
    XAxisUnits = 'DeltaE',
    IsMicroEV = 1 )

# Extra information
# list of Q-values
qvalues = [ 0.275, 0.425, 0.575, 0.725, 0.875, 1.025,
            1.175, 1.325, 1.475, 1.625, 1.775, 1.925 ]

# Do the fit only on these workspace indexes
selected_wi = [ 1, 2, 3, 4] # select a few workspace indexes
#selected_wi = range(0,len(qvalues))  # select all spectra

# Energy range over which we do the fitting.
#  You can edit this to change these boundaries.
minE = -0.1  # Units are in meV
maxE =  0.1

# Initial guess for the lowest Q. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f0.Scaling'   :     0.5,   # Overall intensity
              'f0.f1.f0.Height' :    0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.Height' :    0.9,   # This has to be 1-f0.f1.f0.Height
              'f0.f1.f1.Tau'    :  500.0,   # tau or relaxation time
              'f0.f1.f1.Beta'   :    0.9,   # exponent
              'f1.A0'           :    0.0,   # intercept background
              'f1.A1'           :    0.0,   # slope background
}

# Run sequential fit to obtain good initial guess for each spectrum
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
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True, MaxIterations=500, **spectra_domain_relation)

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

