from mantid.api import *
import math
import numpy 

class StandardSC(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
        self.declareParameter("A0", 0.16)
        self.declareParameter("sigma", 0.2)
        self.declareParameter("fieldSC", 295.0)
        self.declareParameter("fieldBG", 300.0)
        self.declareParameter("phi", 0.0)
        self.declareParameter("Abg", 0.1)
 
  def getFrequency(self,field):
       return (field*0.0851615)/(2.0*3.1415927)
 
  def function1D(self, xvals):
        # Access current values during the fit
        A0 = self.getParameterValue("A0")
        sigma = self.getParameterValue("sigma")
        fieldSC = self.getParameterValue("fieldSC")
        fieldBG = self.getParameterValue("fieldBG")
        phi= self.getParameterValue("phi")
        freqSC = self.getFrequency(fieldSC)
        freqBG = self.getFrequency(fieldBG)
        Abg = self.getParameterValue("Abg")

        return A0*numpy.exp(-sigma*sigma*xvals*xvals/2.)*numpy.cos(2.0*3.1415927*freqSC*xvals+phi)+Abg*numpy.cos(2.0*3.1415927*freqBG*xvals+phi)

FunctionFactory.subscribe(StandardSC)