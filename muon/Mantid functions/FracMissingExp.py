from mantid.api import *
import math
import numpy 

class FracMissingExp(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.3)
        self.declareParameter("alpha", 0.666)
        self.declareParameter("lam", 0.5)
        self.declareParameter("Abg", 0.0)     
 
  def function1D(self, x):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        alpha = self.getParameterValue("alpha")
        lam = self.getParameterValue("lam")
        Abg= self.getParameterValue("Abg")

        return A0*alpha+A0*(1-alpha)*numpy.exp(-lam*x)+Abg

FunctionFactory.subscribe(FracMissingExp)