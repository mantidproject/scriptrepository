''' Script to fit QENS data to a the Fourier trasnform of a stretched exponential
    Relaxation time follows a power-law: tau = taumax * (Qmin/Q)**alpha
    This script should be run in the "Script Window" of MantidPlot
'''

from __future__ import (absolute_import, division, print_function)
from copy import copy
import numpy as np
import scipy.constants
from scipy.fftpack import fft, fftfreq
from scipy.special import gamma
from mantid.api import IFunction1D, FunctionFactory

# Data directory (update with your own)
datadir = "/home/jbq/repositories/mantidproject/scriptrepository/indirect inelastic/BASIS/models/StretchedExpFTTauQ"

""" Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFTTauQ ) + LinearBackground
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
   name=StretchedExpFTTauQ,Height=0.9,TauMax=f0.f1.f1.TauMax,Alpha=f0.f1.f1.Alpha,Beta=f0.f1.f1.Beta,Centre=0,Qmin=_Qmin_,Q=_Q_,
   constraints=(0<TauMax,0<Beta,0<Alpha);
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
   name=StretchedExpFTTauQ,Height=0.9,TauMax=f0.f1.f1.TauMax,Alpha=f0.f1.f1.Alpha,Beta=f0.f1.f1.Beta,Centre=0,Qmin=_Qmin_,Q=_Q_,
   constraints=(0<TauMax,0<Alpha),ties=(Beta=1);
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

# Do the fit only on these workspace indexes, and select Qmin
selected_wi = [ 1, 2, 3, 4]
qmin = qvalues[selected_wi[0]] # Pick qmin as the first of the selected workspace indexes
fitstring_template = fitstring_template.replace("_QMIN_", str(qmin))

# Energy range over which we do the fitting.
#  You can edit this to change these boundaries.
minE = -0.1  # Units are in meV
maxE =  0.1

# Initial guess for the lowest Q. A guess can be obtained by
# running MantidPlot interactively just for the first Q
initguess = { 'f0.f0.Scaling'   :   1.0,   # Overall intensity
              'f0.f1.f0.Height' :   0.1,   # intensity fraction due to elastic line
              'f0.f1.f1.TauMax' : 500.0,   # maximum relaxation time
              'f0.f1.f1.Alpha'  :   3.0,   #
              'f0.f1.f1.Beta'   :   1.0,   # exponent
}

"""Here is the definition of the stretched exponential fit funciton with a power-law
    dependence for the relaxation time
"""
class StretchedExpFTTaoQ(IFunction1D):
    # Class variables
    _planck_constant = scipy.constants.Planck/scipy.constants.e*1E15  # meV*psec

    # pylint: disable=super-on-old-class
    def __init__(self):
        """declare some constants"""
        super(StretchedExpFTTaoQ, self).__init__()
        self._parmList = list()

    def category(self):
        return 'QuasiElastic'

    def init(self):
        """Declare parameters that participate in the fitting"""
        # Active fitting parameters
        self.declareParameter('Height', 0.1, 'Intensity at the origin')
        self.declareParameter('TauMax', 100.0, 'Maximum relaxation time')
        self.declareParameter('Alpha', 2.0, 'Exponent for power-law of Relaxation time dependence')
        self.declareParameter('Beta', 1.0, 'Stretching exponent')
        self.declareParameter('Centre', 0.0, 'Centre of the peak')
        # Keep order in which parameters are declared. Should be a class
        # variable but we initialize it just below parameter declaration
        # to make sure we follow the order.
        self._parmList = ['Height', 'TauMax', 'Alpha', 'Beta', 'Centre']
        # Not fitting parameters, but fixed attributes
        self.declareAttribute("Qmin", 0.3, 'Momentum transfer for maximum observed relaxation')
        self.declareAttribute("Q", 0.3, 'Momentum transfer')

    def setAttributeValue(self, name, value):
        """
        This is called by the framework when an attribute is passed to Fit and its value set.
        It's main use is to store the attribute value on the object once to avoid
        repeated calls during the fitting process
        """
        if name == "Qmin":
            self._qmin = value
        elif name == "Q":
            self._q = value

    def validateParams(self):
        """Check parameters are positive"""
        height = self.getParameterValue('Height')
        taumax = self.getParameterValue('TauMax')
        alpha = self.getParameterValue('Alpha')
        beta = self.getParameterValue('Beta')
        Centre = self.getParameterValue('Centre')
        qmin = self.getAttributeValue("Qmin")
        q = self.getAttributeValue("Q")

        for _, value in {'Height': height, 'TauMax': taumax, 'Alpha': alpha,
                         'Beta': beta, 'Qmin': qmin, 'Q': q}.items():
            if value <= 0:
                return None
        return {'Height': height, 'TauMax': taumax, 'Alpha': alpha, 'Beta': beta,
                'Centre': Centre, 'Qmin': qmin, 'Q': q}

    def function1D(self, xvals, **optparms):
        """ Fourier transform of the Symmetrized Stretched Exponential

        The Symmetrized Stretched Exponential with a power-law Q-dependent relaxation time:
                height * exp( - |t/tau|**beta )
                tau = taumax*(Qmin/Q)**alpha

        The Fourier Transform:
            F(E) \int_{-infty}^{infty} (dt/h) e^{-i2\pi Et/h} f(t)
            F(E) is normalized:
                \int_{-infty}^{infty} dE F(E) = 1
        """
        p = self.validateParams()
        if not p:
            # return zeros if parameters not valid

        tau = p['TauMax'] * (p['Qmin'] / p['Q']) ** p['alpha']
        ne = len(xvals)
        # energy spacing. Assumed xvals is a single-segment grid
        # of increasing energy values
        de = (xvals[-1] - xvals[0]) / (ne-1)
        erange = 2*max(abs(xvals))
        dt = 0.5*StretchedExpFT._planck_constant/erange  # spacing in time
        tmax = StretchedExpFT._planck_constant/de  # maximum reciprocal time
        # round to an upper power of two
        nt = 2**(1+int(np.log(tmax/dt)/np.log(2)))
        sampled_times = dt * np.arange(-nt, nt)
        decay = np.exp(-(np.abs(sampled_times)/tau)**p['Beta'])
        # The Fourier transform introduces an extra factor exp(i*pi*E/de),
        # which amounts to alternating sign every time E increases by de,
        # the energy bin width. Thus, we take the absolute value
        fourier = np.abs(fft(decay).real)  # notice the reverse of decay array
        fourier /= fourier[0]  # set maximum to unity
        # Normalize the integral in energies to unity
        fourier *= 2*tau*gamma(1./p['Beta'])/(p['Beta']*StretchedExpFT._planck_constant)
        # symmetrize to negative energies
        fourier = np.concatenate([fourier[nt:], fourier[:nt]])  # increasing ordering
        # Find corresponding energies
        energies = StretchedExpFT._planck_constant * fftfreq(2*nt, d=dt)  # standard ordering
        energies = np.concatenate([energies[nt:], energies[:nt]])  # increasing ordering
        transform = p['Height'] * np.interp(xvals-p['Centre'], energies, fourier)
        return transform

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(StretchedExpFTTauQ)

names = initguess.keys()  # Store the names of the parameters in a list

chi2 = []     # store the Chi-square values of the fits for every Q
results = []  # we will store in this list the fitted parameters for every Q

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
        fitstring = fitstring.replace("_Q_", str(qvalues[iq])) # momentum transfer

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
        line = '{0} {1:6.2f} {2:5.3f} {3:6.2f} {4:5.3f} {5:5.3f}\n'.format( result['qvalue'],
            chi2[ iq ],
            result['f0.f1.f0.Height'],
            result['f0.f1.f1.TauMax'],
            result['f0.f1.f1.Alpha'],
            result['f0.f1.f1.Beta'])
        print line,
        buffer += line

# Save the string to file "results_StretchedExpFT.dat" whithin the data dir
open( '{0}/results_StretchedExpFTTauQ.dat'.format(datadir), 'w' ).write( buffer )
