from mantid.api import *
import math
import numpy 

class FracTwoExp(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.3)
        self.declareParameter("Af", 0.666)
        self.declareParameter("lambda1", 5)
        self.declareParameter("lambda2", 0.5)
        self.declareParameter("Abg", 0.0)     
 
  def function1D(self, x):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        A = self.getParameterValue("A")
        lambda1 = self.getParameterValue("lambda1")
        lambda2 = self.getParameterValue("lambda2")
        Abg= self.getParameterValue("Abg")

        return Af*alpha*numpy.exp(-lambda1*x)+A0*(A0-Af)*numpy.exp(-lambda1*x)+Abg

FunctionFactory.subscribe(FracTwoExp)