from __future__ import (absolute_import, division, print_function)

import os
import pytest
from distutils.version import LooseVersion as version
from scipy.fftpack import fft, fftfreq
from scipy.special import gamma
from scipy import constants
import numpy as np
from numpy.testing import assert_allclose
import lmfit
from lmfit.models import (Model, index_of)
from lmfit.lineshapes import gaussian, lorentzian

planck_constant = constants.Planck / constants.e * 1E15  # meV*psec

def strexpft(x, amplitude=1.0, center=0.0, tau=10.0, beta=1.0):
    r"""Fourier transform of the symmetrized stretched exponential

    .. math::

        S(E) = A \int_{-\infty}^{\infty} dt/h e^{-i2\pi(E-E_0)t/h} e^{|\frac{x}{\tau}|^\beta}

    Normalization and maximum at :math:`E=E_0`:

    .. math::

        \int_{-\infty}^{\infty} dE S(E) = A
        max(S) = A \frac{\tau}{\beta} \Gamma(\beta^{-1})

    Uses scipy.fftpack.fft for the Fourier transform

    Parameters
    ----------
    x: sequence of floats
        domain of the function, energy
    amplitude : float
        Integrated intensity of the curve
    center : float
        position of the peak
    tau: float
        relaxation time.
    beta: float
        stretching exponent

    If the time units are picoseconds, then the energy units are mili-eV.

    Returns
    -------
    values: sequence of floats
        function over the domain
    """
    ne = len(x)
    # energy spacing. Assumed x is a grid of increasing energy values
    refine_factor = 16  # for better calculation of the fourier transform
    de = (x[-1] - x[0]) / (refine_factor * (ne - 1))
    erange = 2 * max(abs(x))
    dt = 0.5 * planck_constant / erange  # spacing in time
    tmax = planck_constant / de  # maximum reciprocal time
    # round to an upper power of two
    nt = 2 ** (1 + int(np.log(tmax / dt) / np.log(2)))
    sampled_times = dt * np.arange(-nt, nt)
    decay = np.exp(-(np.abs(sampled_times) / tau) ** beta)
    # The Fourier transform introduces an extra factor exp(i*pi*E/de),
    # which amounts to alternating sign every time E increases by de,
    # the energy bin width. Thus, we take the absolute value
    fourier = np.abs(fft(decay).real)  # notice the reverse of decay array
    fourier /= fourier[0]  # set maximum to unity
    # Normalize the integral in energies to unity
    fourier *= 2 * tau * gamma(1. / beta) / (beta * planck_constant)
    # symmetrize to negative energies
    fourier = np.concatenate(
        [fourier[nt:], fourier[:nt]])  # increasing ordering
    # Find energy values corresponding to the fourier values
    energies = planck_constant * fftfreq(2 * nt, d=dt)  # standard ordering
    energies = np.concatenate(
        [energies[nt:], energies[:nt]])  # increasing ordering
    # Interpolate at the requested energy values x, shifted by center.
    return amplitude * np.interp(x - center, energies, fourier)


class StretchedExponentialFTModel(Model):
    r"""Fourier transform of the symmetrized stretched exponential

    .. math::

        S(E) = A \int_{-\infty}^{\infty} dt/h e^{-i2\pi(E-E_0)t/h} e^{|\frac{x}{\tau}|^\beta}

    Normalization and maximum at :math:`E=E_0`:

    .. math::

        \int_{-\infty}^{\infty} dE S(E) = A
        max(S) = A \frac{\tau}{\beta} \Gamma(\beta^{-1})

    Uses scipy.fftpack.fft for the Fourier transform

    Fitting parameters:
        - integrated intensity ``amplitude`` :math:`A`
        - position of the peak ``center`` :math:`E_0`
        - nominal relaxation time ``tau``` :math:`\tau` in picoseconds
        - stretching exponent ``beta`` :math:`\beta`

    If the time unit is picoseconds, then the reciprocal energy unit is mili-eV
    """

    def __init__(self, independent_vars=['x'], prefix='', missing=None,
                 name=None,  **kwargs):
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars})
        super(StretchedExponentialFTModel, self).__init__(strexpft, **kwargs)
        self.set_param_hint('amplitude', min=0)
        self.set_param_hint('tau', min=0)
        self.set_param_hint('beta', min=0)

    if version(lmfit.__version__) > version('0.9.5'):
        __init__.__doc__ = lmfit.models.COMMON_INIT_DOC

    def guess(self, y, x=None, **kwargs):
        r"""Guess starting values for the parameters of a model.

        Parameters
        ----------
        y : sequence of floats
            Intensities
        x : sequence of floats
            energy values
        kwargs : dict
            additional optional arguments, passed to model function.

        Returns
        -------
        Parameters
            parameters with guessed values
        """
        amplitude = 1.0
        center = 0.0
        tau = 10.0
        beta = 1.0
        if x is not None:
            center = x[index_of(y, max(y))]  # assumed peak within domain x
            # Assumptions:
            #   1. large dynamic range, function decays fast in domain x
            #   2. x-values are equispaced
            amplitude = sum(y)* (max(x)-min(x))/len(x)
            tau = max(y) / amplitude  # assumed beta==1.0
        return self.make_params(amplitude=amplitude,
                                center=center,
                                tau=tau,
                                beta=beta)


# Test that FT{Gaussian} = Gaussian, and FT{exponential} = Lorentzian
x = np.arange(-0.1, 0.5, 0.0004)  # energy domain, in meV
# items are (tau, beta, intensities). Assumed that tau unit is picoseconds
trios = [(20.0, 2.0,
          gaussian(x, amplitude=1.0, center=0.0,
                   sigma=planck_constant/(np.sqrt(2.0)*20.0*np.pi))),
         (100.0, 1.0,
          # sigma below is actually the HWHM
          lorentzian(x, amplitude=1.0, center=0.0,
                     sigma=planck_constant / (2 * np.pi * 100.0)))
         ]

@pytest.mark.parametrize('tau, beta, y', trios)
def test_lineshapes(tau, beta, y):
    model = StretchedExponentialFTModel()
    params = model.guess(y, x=None)
    # Stray away from optimal parameters
    [params[name].set(value=val) for name, val in
     dict(amplitude=3.0, center=0.0002, tau=50.0, beta=1.5).items()]
    r = model.fit(y, params, x=x)
    all_params = 'tau beta amplitude center'.split()
    assert_allclose([tau, beta, 1.0, 0.0],
                    [r.params[p].value for p in all_params],
                    rtol=0.1, atol=0.000001)


if __name__ == '__main__':
    pytest.main([os.path.abspath(__file__)])
