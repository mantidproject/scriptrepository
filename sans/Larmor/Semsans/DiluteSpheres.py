import numpy as np
from scipy.special import kn
from mantid.api import IFunction1D, FunctionFactory

import numpy as np

class DiluteSpheres(IFunction1D):

    def category(self):
        return "SESANS"
        
    def init(self):
        self.declareParameter("Radius", 1.0)
        self.declareParameter("Const", 0.003)
        
    def function1D(self, x):
        R = self.getParameterValue("Radius")
        C = self.getParameterValue("Const")
        
        z = np.array(x/R, dtype=np.complex128)
        
        first = np.sqrt(1-(z/2)**2)*(1+1/8.0*z**2)
        second = 0.5*z**2*(1-(z/4)**2)*np.log(z/(2+np.sqrt(4-z**2)))
        
        return C*np.array(np.real(first+second-1), dtype=np.float64)
        
FunctionFactory.subscribe(DiluteSpheres)