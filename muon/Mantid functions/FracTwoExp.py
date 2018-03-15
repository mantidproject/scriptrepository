from mantid.api import *
import math
import numpy 

class FracTwoExp(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.3)
        self.declareParameter("alpha", 0.666)
        self.declareParameter("lambda1", 5)
        self.declareParameter("lambda2", 0.5)
        self.declareParameter("Abg", 0.0)     
 
  def function1D(self, x):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        alpha = self.getParameterValue("alpha")
        lambda1 = self.getParameterValue("lambda1")
        lambda2 = self.getParameterValue("lambda2")
        Abg= self.getParameterValue("Abg")

        return A0*alpha*numpy.exp(-lambda1*x)+A0*(1-alpha)*numpy.exp(-lambda2*x)+Abg

FunctionFactory.subscribe(FracTwoExp)