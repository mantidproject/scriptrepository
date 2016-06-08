''' 
                  Script to fit QENS data to a the Fourier trasnform of a stretched exponential
                  This script should be run in the "Script Window" of MantidPlot
'''
from copy import copy

# Data directory
datadir = "/SNS/users/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/BASIS/models/StretchedExpFT"     #update this with your own directory
#datadir = "/SNS/users/jbq/test"  #update this with your own directory

""" Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFT ) + LinearBackground
    with 0<x<1 is the fraction of the elastic intensity
"""

""" Below is the model cast as a template string suitable for the
    Fit algoritm of Mantid. You can obtain similar string by setting up a model
    in the "fit wizard" of MantidPlot and then
    "Manage Setup" --> "Copy to Clipboard". This actions will save the model
    as a string which you can later paste onto this script.
"""
fitstring_template = """
(composite=Convolution,FixResolution=false,NumDeriv=true;
  name=TabulatedFunction,Workspace=resolution,WorkspaceIndex=1,
  Scaling=f0.f0.Scaling,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1);
  (name=DeltaFunction,Height=0.1,Centre=0,constraints=(0<Height<1);
   name=StretchedExpFT,Height=0.9,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
   constraints=(0<Tau,0<Beta);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=0,A1=0"""

'''
############
# If you want the same model with beta fixed to one, use this template string instead
############
fitstring_template = """
(composite=Convolution,FixResolution=false,NumDeriv=true;
  name=TabulatedFunction,Workspace=resolution,WorkspaceIndex=1,
  Scaling=f0.f0.Scaling,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1);
  (name=DeltaFunction,Height=0.1,Centre=0,constraints=(0<Height<1);
   name=StretchedExpFT,Height=0.9,Tau=f0.f1.f1.Tau,Beta=f0.f1.f1.Beta,Centre=0,
   constraints=(0<Tau),ties=(Beta=1);
   ties=(f1.Height=1-f0.Height,f1.Centre=f0.Centre)
  )
);name=LinearBackground,A0=0,A1=0"""
'''
fitstring_template = fitstring_template.strip(' \t\n\r')  # remove whitespaces and such
print fitstring_template

# Load the data. We assume the format is DAVE group file.
#  Use the "LoadDaveGrp" algorithm
LoadDaveGrp( Filename = '{0}/data.dat'.format( datadir ),
            OutputWorkspace = 'data', XAxisUnits = 'DeltaE', IsMicroEV = 1 )

# Alternatively, we use the "LoadNexus" algorithm if we have the reduced data
# as a Nexus file
# LoadNexus( Filename = '{0}/BASIS_17706_1run_divided.nxs'.format( datadir ),
#             OutputWorkspace = 'data' )

# Load the resolution
LoadDaveGrp( Filename = '{0}/resolution.dat'.format( datadir ),
             OutputWorkspace = 'resolution',
             XAxisUnits = 'DeltaE',
             IsMicroEV = 1 )

# Extra information
# list of Q-values
qvalues = [ 0.275, 0.425, 0.575, 0.725, 0.875, 1.025,
            1.175, 1.325, 1.475, 1.625, 1.775, 1.925 ]
nq = len( qvalues )

# Energy range over which we do the fitting.
#  You can edit this to change these boundaries.
minE = -0.1  # Units are in meV
maxE =  0.1

# Initial guess for the lowest Q. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f0.Scaling'   :   1.0,   # Overall intensity
              'f0.f1.f0.Height' :   0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.Tau'    : 500.0,   # tau or relaxation time
              'f0.f1.f1.Beta'   :   1.0,   # exponent
}

names = initguess.keys()  # Store the names of the parameters in a list

chi2 = []     # store the Chi-square values of the fits for every Q
results = []  # we will store in this list the fitted parameters for every Q

# Do the fit only on these workspace indexes
selected_wi = [ 1, 2, 3, 4]

# "iq" is the workspace index variable
# (Note: the first value of iq is zero, not one)
for iq in range(nq):
        # Insert blank data if the workspace is not among the selected ones
        if iq not in selected_wi:
                chi2.append(None)
                results.append(None)
                continue  # skip to next workspace index 
        # Update the model template string with the particular workspace index
        fitstring = fitstring_template.replace( 'iQ', str(iq) )

        # Update the model template string with the initial guess
        for key, value in initguess.items():
                fitstring = fitstring.replace( key, str( value ) )

        # Call the Fit algorithm using the updated "fitstring".
        print fitstring+'\n'
        Fit( fitstring, InputWorkspace = 'data', WorkspaceIndex = iq,
        CreateOutput = 1, startX = minE, endX = maxE,
        MaxIterations=200 )

        # As a result of the fit, three workspaces are created:
        # "data_Parameters" : optimized parameters and Chi-square
        # "data_NormalisedCovarianceMatrix" : correlations between parameters
        # "data_Workspace"       : data, fit, residuals, and model

        # Update the initial guess for the next fit using the results of the
        # current fit. We do this by scanning every row in "data_Parameters"
        params_workspace = mtd[ 'data_Parameters' ]
        for row in params_workspace:
                # name of the fitting parameter for this particular row
                name = row[ 'Name' ]
                if name in names:
                        # update the value of the fitting parameter
                        initguess[ name ] = row[ 'Value' ]
                if name == 'Cost function value':  # store Chi-square separately
                        chi2.append( row[ 'Value' ] )

        # z is a python dictionary containing the fitted parameters for this Q
        z = copy( initguess )
        # add also the workspace index and the particular Q value
        z.update( {'iq' : iq, 'qvalue' : qvalues[ iq ] } )

        # Append the dictionary to the list of results
        results.append( z )

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
        RenameWorkspace( InputWorkspace = 'data_Parameters',
                         OutputWorkspace = 'data_{0}_Parameters'.format(iq) )
        RenameWorkspace( InputWorkspace = 'data_NormalisedCovarianceMatrix',
                         OutputWorkspace='data_{0}_NormalisedCovarianceMatrix'.format(iq))
        RenameWorkspace( InputWorkspace = 'data_Workspace',
                         OutputWorkspace = 'data_{0}_Workspace'.format(iq) )


# Save some selected parameters to a string.
# Remember that only the selected spectra contain data.
buffer = '#  Q    Chi2   x   Tau(ps) Beta\n'
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

# Save the string to file "results_StretchedExpFT.dat" whithin the data dir
open( '{0}/results_StretchedExpFT.dat'.format(datadir), 'w' ).write( buffer )
