from mantid.api import *
import math
import numpy as np
from scipy import special as sp

class Bessel(IFunction1D):
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("fi",0.1)
    self.declareParameter("nu",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    fi = self.getParameterValue("fi")
    nu = self.getParameterValue("nu")
    out=A0*(2.0*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)
    return out
  
  def der0(self,A0,fi,nu,x):
      return (2*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)
  def der1(self,A0,fi,nu,x):
      return A0*(np.pi/180)*sp.j0(x)
  def der2(self,A0,fi,nu,x):
      return A0*(2*x*np.pi)*sp.j0(x)
  
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    fi = self.getParameterValue("fi");
    nu = self.getParameterValue("nu");
    #der0 = (2*np.pi*nu*x+(np.pi)/180)*sp.j0(x)
    #der1 = A0*(np.pi/180)*sp.j0(x)
    #der2 = A0*(2*x*np.pi)*sp.j0(x)
    
    # X index
    #i = 0
    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,fi,nu,x))
      jacobian.set(i,1, self.der1(A0,fi,nu,x))
      jacobian.set(i,2, self.der2(A0,fi,nu,x))
      #i += 1

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

  def fi(self):
    return self.getParameterValue("fi")
  
  def nu(self):
    return self.getParameterValue("nu")



  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)

  def setfi(self, new_fi):
    self.setParameter("fi", new_fi)

  def setnu(self, new_nu):
    self.setParameter("nu", new_nu)
# Required to have Mantid recognise the new function
FunctionFactory.subscribe(Bessel)