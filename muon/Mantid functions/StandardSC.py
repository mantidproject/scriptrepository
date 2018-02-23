from mantid.api import *
import math
import numpy 

class StandardSC(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.16)
        self.declareParameter("sigma", 0.2)
        self.declareParameter("field", 300.0)
        self.declareParameter("phi", 0.0)
        self.declareParameter("Abg", 0.1)
 
  def getOmega(self,field):
       return field*0.1355
 
  def function1D(self, xvals):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        sigma = self.getParameterValue("sigma")
        field = self.getParameterValue("field")
        phi= self.getParameterValue("phi")
        omega = self.getOmega(field)
        Abg = self.getParameterValue("Abg")

        return A0*numpy.exp(-sigma*sigma*xvals*xvals/2.)*numpy.cos(omega*xvals+phi)+Abg*numpy.cos(omega*xvals+phi)

FunctionFactory.subscribe(StandardSC)