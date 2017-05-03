"""
  For problems with this script, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

Script to sequential fit of QENS data to a the Fourier trasnform of a stretched exponential
This script should be run in the "Script Window" of MantidPlot
"""
import re
from copy import copy
import numpy as np
import mantid.simpleapi as msapi


"""
   Below are the variables that can be changed by the user
"""
data_name="data"  # Name of the workspace containing the QENS signal
resolution_name="resolution"  # Name of the workspace containing the resolution

# Energy range over which we do the fitting.
minE = -0.1  # Units are in meV
maxE =  0.1

# Do the fit only on these workspace indexes (Note: the index begins at zero, not one)
selected_wi = [0, 1, 2, 3, 4]
#selected_wi = None   # uncomment this line if your want to select all spectra

# Initial guess for the lowest selected index. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f1.f0.Height' :    0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.Height' :    0.9,   # This has to be 1-f0.f1.f0.Height
              'f0.f1.f1.Tau'    : 1000.0,   # tau or relaxation time
              'f0.f1.f1.Beta'   :    0.9,   # exponent
              'f1.A0'           :    0.0,   # intercept background
              'f1.A1'           :    0.0,   # slope background
}

# Settings for the minimizer. See the "Fit" algorithm in the documentation
minimizer="FABADA"  # slow, but more reliable than "Levenberg-Marquardt"
maxIterations=5000


"""
   Beginning here, the user does not need to change anything
"""

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
    errors = []  # we will store in this list the errors of fitted parameters for every Q
    funcStrings = []  # store the fitstring with the optimized values

    # "iq" is the workspace index variable
    # (Note: the first value of iq is zero, not one)
    for iq in range(nq):
            # Insert blank data if the workspace is not among the selected ones
            if iq not in selectedwi:
                    chi2.append(None)
                    results.append(None)
                    errors.append(None)
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
                       CreateOutput=1, startX=minE, endX=maxE,
                       Minimizer=minimizer, MaxIterations=maxIterations )

            # As a result of the fit, three workspaces are created:
            # dataName+"_Parameters" : optimized parameters and Chi-square
            # dataName+"_NormalisedCovarianceMatrix" : correlations between parameters
            # dataName+"_Workspace"       : data, fit, residuals, and model

            # Update the initial guess for the next fit using the results of the
            # current fit. We do this by scanning every row in dataName+"_Parameters"
            params_workspace = msapi.mtd[dataName+'_Parameters']
            error = dict()
            for row in params_workspace:
                    # name of the fitting parameter for this particular row
                    name = row['Name']
                    if name in names:
                            # update the value of the fitting parameter
                            initguess[name] = row['Value']
                            error[name] = row['Error']
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
            error.update({'iq': 0, 'qvalue': 0})
            # Append dictionaries to the list of results and associated errors
            results.append( z )
            errors.append(error)

            # Strangely, we need a non-vanishing value of f0.f1.f0.Height as initial guess
            #initguess['f0.f1.f0.Height'] = max(0.1, initguess['f0.f1.f0.Height'])
            #initguess['f0.f1.f1.Height'] = 1 - initguess['f0.f1.f0.Height']

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
                                   OutputWorkspace = "seqfit_"+dataName+'_Q{0}_Parameters'.format(qvalues[iq]) )
            msapi.RenameWorkspace( InputWorkspace = dataName+'_NormalisedCovarianceMatrix',
                                   OutputWorkspace= "seqfit_"+dataName+'_Q{0}_NormalisedCovarianceMatrix'.format(qvalues[iq]))
            msapi.RenameWorkspace( InputWorkspace = dataName+'_Workspace',
                                   OutputWorkspace = "seqfit_"+dataName+'_Q{0}_Workspace'.format(qvalues[iq]) )
            msapi.RenameWorkspace( InputWorkspace = "PDF",
                                   OutputWorkspace = "seqfit_"+dataName+'_Q{0}_PDF'.format(qvalues[iq]) )


    # Save some selected parameters to a string.
    # Remember that only the selected spectra contain data.
    buffer = '#Sequential fit summary\n#  Q    Chi2 Tau(ps)+-error  Beta+-error\n'
    print buffer,
    for iq in range( nq ):
            if iq not in selectedwi:
                continue # skip to next workspace index
            # python dictionary holding  results for the fit of this workspace index
            result = results[ iq ]
            error = errors[iq]
            line = '{0} {1:6.2f} {2:7.2f} {3:6.2f} {4:5.3f} {5:4.3f}\n'.format( result['qvalue'],
                chi2[ iq ],
                result['f0.f1.f1.Tau'], error['f0.f1.f1.Tau'],
                result['f0.f1.f1.Beta'], error['f0.f1.f1.Beta'])
            print line,
            buffer += line

    # Save Q-dependence of parameters to a Workspace
    other = {}
    other_error = {}
    for iq in range( nq ):
        if iq not in selectedwi:
            continue # skip to next workspace index
        # python dictionary holding  results for the fit of this workspace index
        result = results[ iq ]
        error = errors[iq]
        for key,value in result.items():
            if key not in other.keys():
                other[key] = np.array([value,])
                other_error[key] = np.array([error[key],])
            else:
                other[key] = np.append(other[key], value)
                other_error[key] = np.append(other_error[key], error[key])
        if "Chi2" not in other.keys():
            other["Chi2"] = np.array([chi2[iq],])
            other_error["Chi2"] = np.array([0.0,])  # no error in the optimal Chi2
        else:
            other["Chi2"] = np.append(other["Chi2"], chi2[ iq ])
            other_error["Chi2"] = np.append(other_error["Chi2"],0.0)
    # Get the EISF and associated error with error propagation
    a=other["f0.f1.f0.Height"]; b=other["f0.f1.f1.Height"]
    ae=other_error["f0.f1.f0.Height"];  be=other_error["f0.f1.f1.Height"]
    other["EISF"] = a/(a+b)
    other_error["EISF"]=np.sqrt(b*ae*ae+a*be*be)/(a+b)  # Used error propagation formula
    # Plot the Q-dependenced of the optimal parameters
    nspectra = 4
    dataY = np.concatenate((other["EISF"],other['f0.f1.f1.Tau'],other['f0.f1.f1.Beta'],other["Chi2"]))
    dataE = np.concatenate((other_error["EISF"],other_error['f0.f1.f1.Tau'],
        other_error['f0.f1.f1.Beta'],other_error["Chi2"]))
    dataX = np.tile(other["qvalue"], nspectra)
    # Save dependence versus Q
    seqfit_Qdependencies = msapi.CreateWorkspace(DataX=dataX, UnitX="MomentumTransfer", DataY=dataY, DataE=dataE,
                                                 NSpec=nspectra, WorkspaceTitle="Q-dependence of parameters",
                                                 VerticalAxisUnit="Text",
                                                 VerticalAxisValues=["EISF", "Tau", "Beta", "Chi2"])
    # Save dependence versus Q**2
    dataX = np.array(dataX)**2
    seqfit_Q2dependencies = msapi.CreateWorkspace(DataX=dataX, UnitX="QSquared", DataY=dataY, DataE=dataE,
                                                  NSpec=nspectra, WorkspaceTitle="Q squared-dependence of parameters",
                                                  VerticalAxisUnit="Text",
                                                  VerticalAxisValues=["EISF", "Tau", "Beta", "Chi2"])

    return {"funcStrings": funcStrings,}

if __name__ == "__main__":
    
    # Get handle to workspaces
    data=mtd[data_name]
    resolution=mtd[resolution_name]

    # Extract Q-values
    vertical_axis = data.getAxis(1)
    qvalues = vertical_axis.extractValues()
    if not selected_wi:
        selected_wi=range(len(qvalues))

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
            Scaling=1,Shift=0,XScaling=1,ties=(Scaling=1,XScaling=1);
        (
         name=DeltaFunction,Height=f0.f1.f0.Height,Centre=0,constraints=(0<Height),ties=(Centre=0);
         name=StretchedExpFT,Height=f0.f1.f1.Height,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
             constraints=(0<Tau,0<Beta),ties=(Centre=0)
        );
    );
    name=LinearBackground,A0=0.0,A1=0.0"""
    fitstring_template = re.sub('[\s+]', '', fitstring_template)  # remove whitespaces and such
    sequentialFit(resolution, data, fitstring_template, initguess, [minE, maxE], qvalues, selected_wi)
