'''
  For problems with this script, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

  Script to global fit of QENS data to the jump-diffusion model by Teixeira
  This script should be run in the "Script Window" of MantidPlot

  Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*TeixeiraWaterSQE ) + LinearBackground
    with 0<x<1 is the fraction of the elastic intensity

  Parameters DiffCoeff and Tau (residence time) of fit function TeixeiraWaterSQE are
  the same for all spectra. All other fitting  parameters are different for each spectrum
'''
import re
from copy import copy
import sys
import numpy as np
from os.path import join as pjoin

"""
   Below are the variables that can be changed by the user
"""
data_name="irs26176_graphite002_red"  # Name of the workspace containing the QENS signal
resolution_name="irs26173_graphite_res"  # Name of the workspace containing the resolution
# Energy range over which we do the fitting.
minE = -0.4  # Units are in meV
maxE =  0.4
# Do the fit using only these workspaces
#selected_wi = [ 0, 3, 4, 7] # select a few workspace indexes
selected_wi = None  # use all spectra for the fit

# Get handle to workspaces
data=mtd[data_name]
resolution=mtd[resolution_name]

# Find out the Q-values from the loaded data
Qs=GetQsInQENSData(data)
nQ=len(Qs)
if not selected_wi:
    selected_wi = range(0,len(Qs))

""" Below is the model cast as a template string suitable for the
    Fit algoritm of Mantid. You can obtain similar string by setting up a model
    in the "fit wizard" of MantidPlot and then
    "Manage Setup" --> "Copy to Clipboard". This actions will save the model
    as a string which you can later paste onto this script.

    Our initial guesses are DiffCoeff=1.0 and Tau=1.0
"""
single_model_template="""(composite=Convolution,NumDeriv=true;
name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=0,Scaling=1,Shift=0,XScaling=1;
(name=DeltaFunction,Height=1.5,Centre=0;
name=TeixeiraWaterSQE,Q=_Q_,Height=1.0,Tau=1.0,DiffCoeff=1.0,Centre=0;
ties=(f1.Centre=f0.Centre)));
name=LinearBackground,A0=0,A1=0"""
single_model_template = re.sub('[\s+]', '', single_model_template)  # remove whitespaces and such

nWi=len(selected_wi)  # number of selected spectra for fitting

# Create the string representation of the global model (for selected spectra):
global_model="composite=MultiDomainFunction,NumDeriv=true;"
for wi in selected_wi:
    Q=Qs[wi]
    single_model = single_model_template.replace("_Q_", str(Q))  # insert Q-value
    single_model = single_model.replace("_RESOLUTION_", resolution_name)
    global_model += "(composite=CompositeFunction,NumDeriv=true,$domains=i;{0});\n".format(single_model)

# Tie DiffCoeff and Tau since they are global parameters
ties=['='.join(["f{0}.f0.f1.f1.DiffCoeff".format(di) for di in reversed(range(nWi))]),
    '='.join(["f{0}.f0.f1.f1.Tau".format(wi) for wi in reversed(range(nWi))]) ]
global_model += "ties=("+','.join(ties)+')'  # tie Radius

# Now relate each domain(i.e. spectrum) to each single-spectrum model
domain_model = dict()
domain_index = 0
for wi in selected_wi:
    if domain_index == 0:
        domain_model.update({"InputWorkspace": data.name(), "WorkspaceIndex": str(wi),
                             "StartX": str(minE), "EndX": str(maxE)})
    else:
        di = str(domain_index)
        domain_model.update({"InputWorkspace_"+di: data.name(), "WorkspaceIndex_"+di: str(wi),
                             "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
    domain_index += 1

# Invoke the Fit algorithm using global_model and domain_model:
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True, MaxIterations=500, **domain_model)

# Save Q-dependencies of the optimized parameters
parameters_workspace = mtd[output_workspace+"_Parameters"]
aliases = {"f0.f1.f0.Height": "Ha", "f0.f1.f1.Height": "Hb", "f0.f1.f1.Tau": "tau", "f0.f1.f1.DiffCoeff": "diffCoeff"}
optimals = dict()
for alias in aliases.values():
    optimals[alias] = list()
names = aliases.keys()
for row in parameters_workspace:
    # name of the fitting parameter for this particular row
    matches = re.search("^f(\d+)\.(.*)", row['Name']) # for instance, f0.f0.f1.f0.Height
    if matches:
        iq, name = matches.groups()
        if name in names:
            optimals[aliases[name]].append(row["Value"])
# Calculate EISF as the fraction of the total intensity multiplying the Delta-Dirac
# (The data is water at ambient conditions, thus EISF turns out very small)
optimals["EISF"]=  [optimals["Ha"][i]/(optimals["Ha"][i]+optimals["Hb"][i]) for i in range(nWi)]
optimals["qvalues"]=[Qs[wi] for wi in selected_wi]  # list the Q-values for the selected spectra
# Save parameters Q-dependence to workspace "glofit_Qdependencies"
dataY = optimals["EISF"] + optimals["diffCoeff"] + optimals["tau"]
dataX = optimals["qvalues"] * 3
glofit_Qdependencies = CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY,
                                       NSpec=3, WorkspaceTitle="Q-dependence of parameters",
                                       VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "diffCoeff", "tau"])
# Save parameters Q-dependence to workspace "glofit_Qdependencies"
dataX = np.array(dataX) ** 2
glofit_Q2dependencies = CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, NSpec=3,
                                        WorkspaceTitle="Q squared-dependence of parameters",
                                        VerticalAxisUnit="Text", VerticalAxisValues=["EISF", "diffCoeff", "tau"])


