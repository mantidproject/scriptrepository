from mantid.api import *
import math
import numpy as np

class SpinGlass(IFunction1D):

  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("Lambda",0.1)
    self.declareParameter("gamma",0.1)
    self.declareParameter("q",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    Lambda = self.getParameterValue("Lambda")
    gamma = self.getParameterValue("gamma")
    q = self.getParameterValue("q")
    ff1 = 4.0*(Lambda**2.0)*(1.0-q)*x/gamma
    ff2 = q*(Lambda**2.0)*(x**2.0)
    return A0*((1.0/3.0)*np.exp(-np.sqrt(ff1)) +(2.0/3.0)*(1.0-(ff2/np.sqrt(ff1+ff2)))*np.exp(-np.sqrt(ff1+ff2))
  
    
  #def functionDerivLocal(self, x, jacobian):
  #  A0 = self.getParameterValue("A0");
  #  Lambda = self.getParameterValue("Lambda");
  #  gamma = self.getParameterValue("gamma");
  #  q = self.getParameterValue("q");
  #  ff1 = 4*(Lambda**2)*(1-q)*x/gamma
  #  ff2 = q*(Lambda**2)*(x**2)
  #  ff3= np.sqrt((1-q)*x*(Lambda**2))
  #  der0 = ((1/3)*np.exp(-np.sqrt(ff1)) +(2/3)*(1-(ff2/np.sqrt(ff1+ff2)))*np.exp(-np.sqrt(ff1+ff2))
  #  der1 = 
  #  der2 = 
  #  der3 = 
    
  #  # X index
  #  i = 0
  #  for x in xvals:
  #    jacobian.set(i,0, der0)
  #    jacobian.set(i,1, der1)
  #    jacobian.set(i,2, der2)
  #    jacobian.set(i,3, der3)
  #    i += 1

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
  def Lambda(self):
    return self.getParameterValue("Lambda")    
  def gamma(self):
    return self.getParameterValue("gamma")
  def q(self):
    return self.getParameterValue("q")
    
  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)
  def setLambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)
  def setgamma(self, new_gamma):
    self.setParameter("gamma", new_gamma)
  def setq(self, new_q):
    self.setParameter("q", new_q)
    
# Required to have Mantid recognise the new function
FunctionFactory.subscribe(SpinGlass)