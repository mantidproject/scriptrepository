from mantid.api import *
import math
import numpy as np

class LFradical(IFunction1D):
   # Function to fit the reposarisation of anisotropic radical species in applied longitudinal fields (in gauss)
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("AmpR",1.0)
    self.declareParameter("HF",1585.0)
    self.declareParameter("Aniso",0.0)
       
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    AmpR = self.getParameterValue("AmpR")
    HF = self.getParameterValue("HF")
    Aniso = self.getParameterValue("Aniso")
    return A0+AmpR*(1.0+2.0*((x/HF)**2.0))/(Aniso+2.0*((x/HF)**2.0))
   
  def der0(self,A0,AmpR,HF,Aniso,x):	
	return 1.0
  def der1(self,A0,AmpR,HF,Aniso,x):
      return (1.0+2.0*((x/HF)**2.0))/(Aniso+2.0*((x/HF)**2.0))
  def der2(self,A0,AmpR,HF,Aniso,x):
      return AmpR*(4.0*(x**2.0)*HF*(1.0-Aniso))/((Aniso*(HF**2.0)+2.0*(x**2.0))**2.0)
  def der3(self,A0,AmpR,HF,Aniso,x):
      return -1.0*AmpR*((HF**2.0)+2.0*(x**2.0))*(HF**2.0)/((Aniso*(HF**2.0)+2.0*(x**2.0))**2.0)

   
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    AmpR = self.getParameterValue("AmpR");
    HF = self.getParameterValue("HF")
    Aniso = self.getParameterValue("Aniso")
    
    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,AmpR,HF,Aniso,x))
      jacobian.set(i,1, self.der0(A0,AmpR,HF,Aniso,x))
      jacobian.set(i,2, self.der0(A0,AmpR,HF,Aniso,x))
      jacobian.set(i,3, self.der0(A0,AmpR,HF,Aniso,x))


# jacobian.set(index =i, parameter number,  d/d parameter [function]) 
  #def activeParameter(self, index):
  #  param_value = self.getParameterValue(index)
  #  if index == 2: #Sigma. Actually fit to 1/(sigma^2) for stability
  #    return 1./math.pow(param_value,2)
  #  else:
  #    return param_value

  #def setActiveParameter(self, index, value):
  #  param_value = value
  #  explicit = False
  #  if index == 2:
  #    param_value = math.sqrt(math.fabs(1.0/value))
  #  else:
  #    param_value = value
  #    # Final explicit argument is required to be false here
  #    self.setParameter(index, param_value, False) 
  
  def A0(self):
    return self.getParameterValue("A0")
  def AmpR(self):
    return self.getParameterValue("AmpR")
  def HF(self):
    return self.getParameterValue("HF")
  def Aniso(self):
    return self.getParameterValue("Aniso")


  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)
  def setAmpR(self, new_AmpR):
    self.setParameter("AmpR", new_AmpR)
  def setHF(self, new_HF):
    self.setParameter("HF", new_HF)
  def setAniso(self, new_Aniso):
    self.setParameter("Aniso", new_Aniso)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(LFradical)