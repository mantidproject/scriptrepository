from mantid.api import *
import math
import numpy 

class GausOscSC(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.3)
        self.declareParameter("sigma", 0.2)
        self.declareParameter("field", 0.0)
        self.declareParameter("phi", 0.0)
 
  def getOmega(self,field):
       return field*0.1355
 
  def function1D(self, x):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        sigma = self.getParameterValue("sigma")
        field = self.getParameterValue("field")
        phi= self.getParameterValue("phi")
        omega = self.getOmega(field)

        return A0*numpy.exp(-sigma*sigma*x*x/2.)*numpy.cos(omega*x+phi)

FunctionFactory.subscribe(GausOscSC)