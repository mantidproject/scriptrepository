from __future__ import (absolute_import, division, print_function)
from mantid.api import IFunction1D
import numpy as np
import scipy.constants
from scipy.fftpack import fft, fftfreq
from scipy.special import gamma

class StretchedExpFTTauQ(IFunction1D):
    # Class variables
    _planck_constant = scipy.constants.Planck/scipy.constants.e*1E15  # meV*psec

    # pylint: disable=super-on-old-class
    def __init__(self):
        """declare some constants"""
        super(StretchedExpFTTauQ, self).__init__()
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
        self.declareAttribute("Qmin", 0.3)  #Momentum transfer for maximum observed relaxation
        self.declareAttribute("Q", 0.3)  #Momentum transfer

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
            return np.zeros(len(xvals), dtype=float)

        tau = p['TauMax'] * (p['Qmin'] / p['Q']) ** p['Alpha']
        ne = len(xvals)
        # energy spacing. Assumed xvals is a single-segment grid
        # of increasing energy values
        de = (xvals[-1] - xvals[0]) / (ne-1)
        erange = 2*max(abs(xvals))
        dt = 0.5*StretchedExpFTTauQ._planck_constant/erange  # spacing in time
        tmax = StretchedExpFTTauQ._planck_constant/de  # maximum reciprocal time
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
        fourier *= 2*tau*gamma(1./p['Beta'])/(p['Beta']*StretchedExpFTTauQ._planck_constant)
        # symmetrize to negative energies
        fourier = np.concatenate([fourier[nt:], fourier[:nt]])  # increasing ordering
        # Find corresponding energies
        energies = StretchedExpFTTauQ._planck_constant * fftfreq(2*nt, d=dt)  # standard ordering
        energies = np.concatenate([energies[nt:], energies[:nt]])  # increasing ordering
        transform = p['Height'] * np.interp(xvals-p['Centre'], energies, fourier)
        return transform

