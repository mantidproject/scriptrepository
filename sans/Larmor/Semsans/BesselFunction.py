import numpy as np
from scipy.special import kn
from mantid.api import IFunction1D, FunctionFactory

import numpy as np

class BesselFunction(IFunction1D):

    def category(self):
        return "SESANS"
        
    def init(self):
        self.declareParameter("length_scale", 0.2)
        self.declareParameter("Const", 1)
        
    def function1D(self, x):
        c = self.getParameterValue("length_scale")
        C = self.getParameterValue("Const")
        kappa = c*x**2
        
        return C*kappa*kn(1, kappa)
        
FunctionFactory.subscribe(BesselFunction)