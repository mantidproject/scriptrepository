'''
  For problems with this script, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

   Script to global fit of QENS data to the Fourier transform of a stretched exponential
   This script should be run in the "Script Window" of MantidPlot

   Fitting model. In this case:
       Convolution( A*Resolution, EISF*Delta + (1-EISF)*StretchedExFTTauQ ) + LinearBackground
   with 0<EISF<1 is the fraction of the elastic intensity
   and the relaxation time follows a power-law: Tau = Taumax * (Qmin/Q)**Alpha

  Parameters Beta, TauMax, and Alpha of fit function StretchedExFTTauQ are the same for all spectra.
  All other fitting parameters are different for each spectrum.
'''
from __future__ import print_function

import re
from copy import copy
import sys

# Directory holding the scripts and the fit function
scriptdir = "/home/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/models/StretchedExpFTTauQ/v3.8"
sys.path.append(scriptdir)
from seqFitStretchedExFTTauQ import sequentialFit


# Data directory (update with your own)
datadir = "/home/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/BASIS/models/StretchedExpFT"
data_davegroup_filename = "data.dat"
resolution_davegroup_filename = "resolution.dat"

""" Below is the model cast as a template string suitable for the
    Fit algoritm of Mantid. You can obtain similar string by setting up a model
    in the "fit wizard" of MantidPlot and then
    "Manage Setup" --> "Copy to Clipboard". This actions will save the model
    as a string which you can later paste onto this script.
"""

fitstring_template = """
(composite=Convolution,FixResolution=false,NumDeriv=true;
  name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=_IQ_,
  Scaling=f0.f0.Scaling,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1);
  (name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height<1);
   name=StretchedExpFTTauQ,Height=f0.f1.f1.Height,TauMax=f0.f1.f1.TauMax,
        Alpha=f0.f1.f1.Alpha,Beta=f0.f1.f1.Beta,Centre=0,Qmin=_QMIN_,Q=_Q_,
   constraints=(0<TauMax,0<Beta,0<Alpha);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=f1.A0,A1=f1.A1"""

fitstring_template = re.sub('[\s+]', '', fitstring_template)  # remove whitespaces and such

# Load the data. We assume the format is DAVE group file.
#  Use the "LoadDaveGrp" algorithm
data = LoadDaveGrp( Filename = '{0}/{1}'.format(datadir, data_davegroup_filename),
                    OutputWorkspace = 'data', XAxisUnits = 'DeltaE', IsMicroEV = 1 )

# Alternatively, we use the "LoadNexus" algorithm if we have the reduced data
# as a Nexus file
# LoadNexus( Filename = '{0}/BASIS_17706_1run_divided.nxs'.format( datadir ),
#             OutputWorkspace = 'data' )

# Load the resolution
resolution = LoadDaveGrp( Filename = '{0}/{1}'.format(datadir, resolution_davegroup_filename),
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

qmin = qvalues[selected_wi[0]] # Pick qmin as the first of the selected workspace indexes
fitstring_template = fitstring_template.replace("_QMIN_", str(qmin))

# Initial guess for the lowest Q. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f0.Scaling'   :   1.0,   # Overall intensity
              'f0.f1.f0.Height' :   0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.Height' :   0.9,   # This has to be 1-f0.f1.f0.Height
              'f0.f1.f1.TauMax' : 500.0,   # maximum relaxation time
              'f0.f1.f1.Alpha'  :   2.0,   #
              'f0.f1.f1.Beta'   :   0.9,   # exponent
              'f1.A0': 0.0,  # intercept background
              'f1.A1': 0.0,  # slope background
}

# Run sequential fit to obtain good initial guess for each spectrum
seqOutput = sequentialFit(resolution, data, fitstring_template, initguess, [minE, maxE], qvalues, selected_wi)

# Since we are going to tie parameters Beta, TauMax, and Alpha, find the averages
# and use these number as initial guesses
individuals = {"Beta": list(), "TauMax": list(), "Alpha": list()}
# collect the individual values at each spectrum
for i in range(len(seqOutput["funcStrings"])):
    for name in list(individuals.keys()):
        match = re.search("{0}=(\d+\.\d+)".format(name), seqOutput["funcStrings"][i])
        individuals[name].append(float(match.groups()[0]))
# calculate the averages
averages = dict()
for name in list(individuals.keys()):
    averages[name] = "{0}={1}".format(name, str(sum(individuals[name])/len(individuals[name])))
# substitute the individual values with the averages
for i in range(len(seqOutput["funcStrings"])):
    for name in list(individuals.keys()):
        seqOutput["funcStrings"][i] = re.sub("{0}=\d+\.\d+".format(name),averages[name],
                                             seqOutput["funcStrings"][i])
print(seqOutput["funcStrings"])

# Merge models for each spectra
global_model= 'composite=MultiDomainFunction,NumDeriv=true;'
for funcString in seqOutput["funcStrings"]:
    global_model += "(composite=CompositeFunction,NumDeriv=true,$domains=i;{0});\n".format(funcString)

# Insert the ties for Beta, TauMax, and Alpha parameters
ties = {"Beta": "", "TauMax": "", "Alpha": ""}
for name in list(ties.keys()):
    for i in reversed(list(range(len(selected_wi)))):  
        ties[name] += "f{0}.f0.f1.f1.{1}=".format(i, name)
all_ties = "ties=(" + ','.join([ties[name].strip("=") for name in list(ties.keys())]) + ')'
global_model += all_ties
print(global_model)
#Relate spectra to domains
spectra_domain_relation = dict()
domain_index = 0
for iq in selected_wi:
    if domain_index == 0:
        spectra_domain_relation.update({"InputWorkspace": data.name(), "WorkspaceIndex": str(iq),
                                        "StartX": str(minE), "EndX": str(maxE)})
    else:
        di = str(domain_index)
        spectra_domain_relation.update({"InputWorkspace_"+di: data.name(), "WorkspaceIndex_"+di: str(iq),
                                        "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
    domain_index += 1

# Carry out the fit
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True, MaxIterations=500, **spectra_domain_relation)

# Retrieve the optimized parameters ordered by momentum transfer
parameters_workspace = mtd[output_workspace+"_Parameters"]
shorties = {"f0.f1.f0.Height": "EISF", "f0.f1.f1.TauMax": "taumax",
            "f0.f1.f1.Alpha": "alpha", "f0.f1.f1.Beta": "beta", }
other = dict()
for shorty in list(shorties.values()):
    other[shorty] = list()
other["qvalues"] = [qvalues[iq] for iq in selected_wi]
names = list(shorties.keys())
for row in parameters_workspace:
    # name of the fitting parameter for this particular row
    matches = re.search("^f(\d+)\.(.*)", row['Name']) # for instance, f0.f0.f1.f0.Height
    if matches:
        iq, name = matches.groups()
        if name in names:
            other[shorties[name]].append(row["Value"])
# calculate tau
other["tau"] = (other['taumax'][0]*(qmin/np.array(other["qvalues"]))**other["alpha"][0]).tolist()
print(other)
# Save Q-dependencies of the optimized parameters, and also tau
nspectra = 5
dataY = other["EISF"] + other['taumax'] + other['alpha'] + other['tau'] + other['beta']
dataX = other["qvalues"] * nspectra
# Save dependence versus Q
glofit_Qdependencies = CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY, NSpec=nspectra,
                                       WorkspaceTitle="Q-dependence of parameters", VerticalAxisUnit="Text",
                                       VerticalAxisValues=["EISF", "taumax", "alpha", "tau", "beta"])
dataX = np.array(dataX) ** 2
# Save dependence versus Q**2
glofit_Q2dependencies = CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, NSpec=nspectra,
                                        WorkspaceTitle="Q squared-dependence of parameters",
                                        VerticalAxisUnit="Text",
                                        VerticalAxisValues=["EISF", "tauMax", "alpha", "tau", "beta"])



