Fitting Models
================

Python scripts for fitting of with different QENS models

* StretchedExpFT  
Fourier transform of the stretched exponential. The model is:  
S(Q,E) = Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFT ) + LinearBackground  
with StretchedExFT = height * exp( - |t/tau|**beta )

* StretchedExpFTTauQ - as StretchedExpFT, but with a power-law dependence for the relaxation time.  
The model is:  
Convolution( A*Resolution, x*Delta + (1-x)*StretchedExFTTauQ ) + LinearBackground  
with StretchedExFTTauQ = height * exp( - |t/tau|**beta )  
and tau = taumax*(Qmin/Q)**alpha


 
