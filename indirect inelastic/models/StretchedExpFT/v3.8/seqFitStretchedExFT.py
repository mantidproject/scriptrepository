"""
Script to sequential fit of QENS data to a the Fourier trasnform of a stretched exponential
This script should be run in the "Script Window" of MantidPlot
"""
import re
from copy import copy
import numpy as np
import mantid.simpleapi as msapi

# Data directory (User, update with your own)
datadir = "/home/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/models/StretchedExpFT"     #update this with your own directory
data_davegroup_filename = "data.dat"
resolution_davegroup_filename = "resolution.dat"

""" Fitting model. In this case:
        Convolution( A*Resolution, EISF*Delta + (1-EISF)*StretchedExFT ) + LinearBackground
    with 0<EISF<1 is the fraction of the elastic intensity
"""

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
   name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
   constraints=(0<Tau,0<Beta);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=f1.A0,A1=f1.A1"""

'''
############
# If you want the same model with beta fixed to one, use this template string instead
############
fitstring_template = """
(composite=Convolution,FixResolution=false,NumDeriv=true;
  name=TabulatedFunction,Workspace=_RESOLUTION_,WorkspaceIndex=_IQ_,
  Scaling=f0.f0.Scaling,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1);
  (name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height<1);
   name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
   constraints=(0<Tau),ties=(Beta=1);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=f1.A0,A1=f1.A1"""
'''
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
resolution = msapi.LoadDaveGrp( Filename = '{0}/{1}'.format(datadir, resolution_davegroup_filename),
                                OutputWorkspace = 'resolution',
                                XAxisUnits = 'DeltaE',
                                IsMicroEV = 1 )

# Extra information
# list of Q-values
qvalues = [ 0.275, 0.425, 0.575, 0.725, 0.875, 1.025,
            1.175, 1.325, 1.475, 1.625, 1.775, 1.925 ]

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

# Do the fit only on these workspace indexes
selected_wi = [ 1, 2, 3, 4] # select a few workspace indexes
#selected_wi = range(0,len(qvalues))  # select all spectra

def sequentialFit(resolution, data, fitstring_template, initguess, erange, qvalues, selectedwi):
    """
    Carry out the sequential fitting
    :param resolution: workspace containing the resolution function
    :param data: workspace containing the QENS signal
    :param fitstring_template: string defining the model
    :param initguess: dictionary with guess for parameters of the first workspace selected
    :param erange: list with minimum and maximum energy values
    :param qvalues: list of momentum transfer values
    :param selectedwi: list of indexes defining which values of qvalues we keep
    :return: list of function strings with the optimized parameters
    """
    dataName = data.name()
    fitstring_template = fitstring_template.replace("_RESOLUTION_", resolution.name())
    minE, maxE = erange  # Units are in meV
    nq = len(qvalues)
    names = initguess.keys()  # Store the names of the parameters in a list
    chi2 = []     # store the Chi-square values of the fits for every Q
    results = []  # we will store in this list the fitted parameters for every Q
    funcStrings = []  # store the fitstring with the optimized values

    # "iq" is the workspace index variable
    # (Note: the first value of iq is zero, not one)
    for iq in range(nq):
            # Insert blank data if the workspace is not among the selected ones
            if iq not in selected_wi:
                    chi2.append(None)
                    results.append(None)
                    continue  # skip to next workspace index
            # Update the model template string with the particular workspace index
            fitstring = fitstring_template.replace( '_IQ_', str(iq) )
            funcString = fitstring_template.replace( '_IQ_', str(iq) )

            # Update the model template string with the initial guess
            for key, value in initguess.items():
                fitstring = fitstring.replace( key, str( value ) )

            # Call the Fit algorithm using the updated "fitstring".
            print fitstring+'\n'
            msapi.Fit( fitstring, InputWorkspace=dataName, WorkspaceIndex=iq,
                       CreateOutput=1, startX=minE, endX=maxE, MaxIterations=500 )

            # As a result of the fit, three workspaces are created:
            # dataName+"_Parameters" : optimized parameters and Chi-square
            # dataName+"_NormalisedCovarianceMatrix" : correlations between parameters
            # dataName+"_Workspace"       : data, fit, residuals, and model

            # Update the initial guess for the next fit using the results of the
            # current fit. We do this by scanning every row in dataName+"_Parameters"
            params_workspace = msapi.mtd[dataName+'_Parameters']
            for row in params_workspace:
                    # name of the fitting parameter for this particular row
                    name = row['Name']
                    if name in names:
                            # update the value of the fitting parameter
                            initguess[name] = row['Value']
                    if name == 'Cost function value':  # store Chi-square separately
                            chi2.append(row['Value'])
            # Save the string representation with optimal values
            for key, value in initguess.items():
                funcString = funcString.replace(key, str(value))
            funcStrings.append(funcString)
            
            # z is a python dictionary containing the fitted parameters for this Q
            z = copy( initguess )
            # add also the workspace index and the particular Q value
            z.update({'iq': iq, 'qvalue': qvalues[iq]})

            # Append the dictionary to the list of results
            results.append( z )

            # Strangely, we need a non-vanishing value of f0.f1.f0.Height as initial guess
            initguess['f0.f1.f0.Height'] = max(0.1, initguess['f0.f1.f0.Height'])
            initguess['f0.f1.f1.Height'] = 1 - initguess['f0.f1.f0.Height']

            # Special treatment for the tau parameter!
            # tau decreases very fast with increasing Q,
            # thus a better initial guess is to divide the
            # tau from the previous Q by some number (here 2.0)
            # as the initial guess for the next Q
            initguess['f0.f1.f1.Tau' ] = initguess['f0.f1.f1.Tau' ]  / 2.0

            # The workspaces resulting from the current fit
            # will be overwritten in the next iteration. This is OK,
            # but if we want to keep them we should rename then by
            # tagging them with the workspace index.
            msapi.RenameWorkspace( InputWorkspace = dataName+'_Parameters',
                                   OutputWorkspace = "seqfit_"+dataName+'_{0}_Parameters'.format(iq) )
            msapi.RenameWorkspace( InputWorkspace = dataName+'_NormalisedCovarianceMatrix',
                                   OutputWorkspace= "seqfit_"+dataName+'_{0}_NormalisedCovarianceMatrix'.format(iq))
            msapi.RenameWorkspace( InputWorkspace = dataName+'_Workspace',
                                   OutputWorkspace = "seqfit_"+dataName+'_{0}_Workspace'.format(iq) )


    # Save some selected parameters to a string.
    # Remember that only the selected spectra contain data.
    buffer = '#Sequential fit summary\n#  Q    Chi2   x   Tau(ps) Beta\n'
    print buffer,
    for iq in range( nq ):
            if iq not in selected_wi:
                continue # skip to next workspace index
            # python dictionary holding  results for the fit of this workspace index
            result = results[ iq ]
            line = '{0} {1:6.2f} {2:5.3f} {3:6.2f} {4:5.3f}\n'.format( result['qvalue'],
                chi2[ iq ],
                result['f0.f1.f0.Height'],
                result['f0.f1.f1.Tau'],
                result['f0.f1.f1.Beta'])
            print line,
            buffer += line

    # Save Q-dependence of parameters to a Workspace
    other = {}
    for iq in range( nq ):
        if iq not in selected_wi:
            continue # skip to next workspace index
        # python dictionary holding  results for the fit of this workspace index
        result = results[ iq ]
        for key,value in result.items():
            if key not in other.keys():
                other[key] = [value,]
            else:
                other[key].append(value)
        if "Chi2" not in other.keys():
            other["Chi2"] = [chi2[ iq ],]
        else:
            other["Chi2"].append(chi2[ iq ])
    nspectra = 4
    dataY = other["f0.f1.f0.Height"] + other['f0.f1.f1.Tau'] + other['f0.f1.f1.Beta'] + other["Chi2"]
    dataX = other["qvalue"] * nspectra
    # Save dependence versus Q
    seqfit_Qdependencies = msapi.CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY,
                                                 NSpec=nspectra, WorkspaceTitle="Q-dependence of parameters",
                                                 VerticalAxisUnit="Text",
                                                 VerticalAxisValues=["EISF", "Tau", "Beta", "Chi2"])
    # Save dependence versus Q**2
    dataX = np.array(dataX)**2
    seqfit_Q2dependencies = msapi.CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, NSpec=nspectra,
                                                  WorkspaceTitle="Q squared-dependence of parameters",
                                                  VerticalAxisUnit="Text",
                                                  VerticalAxisValues=["EISF", "Tau", "Beta", "Chi2"])

    return {"funcStrings": funcStrings,}

if __name__ == "__main__":
    sequentialFit(resolution, data, fitstring_template, initguess, [minE, maxE], qvalues, selected_wi)
