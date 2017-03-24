import numpy as np

data=Load("irs26176_graphite002_red")
resolution=Load("irs26173_graphite002_res.nxs")
Emin=-0.4;  Emax=0.4

# Find out the Q-values
numHist = data.getNumberHistograms()
Qs=list()
for idx in range(numHist):
    detector = data.getDetector(idx) 
    efixed = data.getEFixed(detector.getID())  # in meV
    wavelength=9.044567/np.sqrt(efixed)  # in Angstroms
    usignTheta = 0.5 * data.detectorTwoTheta(detector)
    Q = (4*np.pi/wavelength)*np.sin(usignTheta)  # in inverse Angstroms
    Qs.append(Q)
nQ=len(Qs)

# This is the template fitting model for each spectrum (each Q-value):
# Our initial guesses are length=1.0 and  tau=10
single_model_template="""(composite=Convolution,NumDeriv=true;
name=TabulatedFunction,Workspace=resolution,WorkspaceIndex=0,Scaling=1,Shift=0,XScaling=1;
(name=DeltaFunction,Height=1.0,Centre=0;
name=TeixeiraWaterSQE,Q=_Q_,Height=1.0,Tau=1.25,DiffCoeff=2.25,Centre=0;
ties=(f1.Centre=f0.Centre)));
name=LinearBackground,A0=0,A1=0"""
# Now create the string representation of the global model (all spectra, all Q-values):
global_model="composite=MultiDomainFunction,NumDeriv=true;"
wi=0
for Q in Qs:
    single_model = single_model_template.replace("_Q_", str(Q))  # insert Q-value
    print single_model
    global_model += "(composite=CompositeFunction,NumDeriv=true,$domains=i;{0});\n".format(single_model)
    wi+=1
# The Length and Tau are the same for all spectra, thus tie them:
ties=['='.join(["f{0}.f0.f1.f1.DiffCoeff".format(wi) for wi in reversed(range(nQ))]),
    '='.join(["f{0}.f0.f1.f1.Tau".format(wi) for wi in reversed(range(nQ))]) ]
global_model += "ties=("+','.join(ties)+')'  # tie Radius

# Now relate each domain(i.e. spectrum) to each single model
domain_model=dict()
for wi in range(nQ):
    if wi == 0:
        domain_model.update({"InputWorkspace": data.name(), "WorkspaceIndex": str(wi),
            "StartX": str(Emin), "EndX": str(Emax)})
    else:
        domain_model.update({"InputWorkspace_"+str(wi): data.name(), "WorkspaceIndex_"+str(wi): str(wi),
            "StartX_"+str(wi): str(Emin), "EndX_"+str(wi): str(Emax)})
# Invoke the Fit algorithm using global_model and domain_model:
output_workspace = "glofit_"+data.name()
Fit(Function=global_model, Output=output_workspace, CreateOutput=True, MaxIterations=500, **domain_model)

# Extract Height, Radius, and Tau from workspace glofit_data_Parameters, the output of Fit:
nparms=0
parameter_ws = mtd[output_workspace+"_Parameters"]
for irow in range(parameter_ws.rowCount()):
    row = parameter_ws.row(irow)
    if row["Name"]=="f0.f0.f1.f1.DiffCoeff":
        DiffCoeff=row["Value"]
        nparms+=1
    elif row["Name"]=="f0.f0.f1.f1.Tau":
        Tau=row["Value"]
        nparms+=1
    elif row["Name"]=="Cost function value":
        chi2=row["Value"]
        nparms+=1
    if nparms==3:
        break  # We got the three parameters we are interested in
print "DiffCoeff=",DiffCoeff," Tau=",Tau, "Chi square=",chi2

