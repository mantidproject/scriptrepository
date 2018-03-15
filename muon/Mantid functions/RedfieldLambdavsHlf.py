from mantid.api import *
import math
import numpy as np

class RedfieldCutoff(IFunction1D):
   # fits Lambda(H_LF), applies when the distribution of the magnitude of H_loc and of spin fluctuation frequencies is homogeneous in time and space
   # from Suzuki et al., Journal of Physics: Conference Series, 502 (2014) 012041
  
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("Hloc",0.1)
    self.declareParameter("Hlf",0.1)
    self.declareParameter("tau",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    Hloc = self.getParameterValue("Hloc")  # fields in G 
    Hlf = self.getParameterValue("Hlf")
    tau = self.getParameterValue("tau") 
    gamu = 13.55; # kHz/G
    return A0*(2.0*(gamu**2.0)*(Hloc**2.0)*tau)/(1.0+(gamu**2.0)*(Hlf**2.0)*(tau**2.0))

	
  def der0(self,A0,Hloc,Hlf,tau,x)
     return (2.0*(gamu**2.0)*(Hloc**2.0)*tau)/(1.0+(gamu**2.0)*(Hlf**2.0)*(tau**2.0))
  def der1(self,A0,Hloc,Hlf,tau,x)	 
      return (4.0*Hloc*tau*(gamu**2.0))/(1.0+(gamu**2.0)*(Hlf**2.0)*(tau**2.0))
  def der2(self,A0,Hloc,Hlf,tau,x)
      return -(4.0*(Hloc**2.0)*Hlf*(tau**3.0)*(gamu**4.0))/np.power((1.0+(gamu**2.0)*(Hlf**2.0)*(tau**2.0)),2.0)
  def der3(self,A0,Hloc,Hlf,tau,x)
      return -(2.0*(Hloc**2.0)*(gamu**2.0)*((gamu**2.0)*(Hlf**2.0)*(tau**2.0)-1.0))/np.power((1.0+(gamu**2.0)*(Hlf**2.0)*(tau**2.0)),2.0)
    
 
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    Hloc = self.getParameterValue("Hloc");
    Hlf = self.getParameterValue("Hlf");
    tau = self.getParameterValue("tau");
	
    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,Hloc,Hlf,tau,x))
      jacobian.set(i,1, self.der1(A0,Hloc,Hlf,tau,x))
      jacobian.set(i,2, self.der2(A0,Hloc,Hlf,tau,x))
      jacobian.set(i,3, self.der3(A0,Hloc,Hlf,tau,x))


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

  def Hloc(self):
    return self.getParameterValue("Hloc")
  
  def Hlf(self):
    return self.getParameterValue("Hlf")

  def tau(self):
    return self.getParameterValue("tau")


  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)

  def setHloc(self, new_Hloc):
    self.setParameter("Hloc", new_Hloc)

  def setHlf(self, new_Hlf):
    self.setParameter("Hlf", new_Hlf)
    
  def settau(self, new_tau):
    self.setParameter("tau", new_tau)
# Required to have Mantid recognise the new function
FunctionFactory.subscribe(RedfieldCutoff)