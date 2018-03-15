from mantid.api import *
import math
import numpy as np
from scipy import special as sp

class InternBessel(IFunction1D):
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("alpha",1.0)
    self.declareParameter("fi",0.1)
    self.declareParameter("nu",0.1)
    self.declareParameter("LamT",0.1)
    self.declareParameter("LamL",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    alpha = self.getParameterValue("alpha")
    fi = self.getParameterValue("fi")
    nu = self.getParameterValue("nu")
    LamT = self.getParameterValue("LamT")
    LamL = self.getParameterValue("LamL")
    return A0*( alpha*(2.0*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)*np.exp(-LamT*x) + (1.0-alpha)*np.exp(-LamL*x) )
	
    
  def der0(self,A0,alpha,fi,nu,LamT,LamL,x):	
      return ( alpha*(2.0*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)*np.exp(-LamT*x) + (1.0-alpha)*np.exp(-LamL*x) )
  def der1(self,A0,alpha,fi,nu,LamT,LamL,x):
      return A0*( (2.0*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)*np.exp(-LamT*x) - np.exp(-LamL*x) )
  def der2(self,A0,alpha,fi,nu,LamT,LamL,x):
      return A0*(np.pi/180.0)*sp.j0(x)*np.exp(-LamT*x)
  def der3(self,A0,alpha,fi,nu,LamT,LamL,x):
      return A0*(2.0*x*np.pi)*sp.j0(x)*np.exp(-LamT*x)
  def der4(self,A0,alpha,fi,nu,LamT,LamL,x):
      return A0*(alpha*(2.0*np.pi*nu*x+(np.pi)/180.0)*sp.j0(x)*(-x)*np.exp(-LamT*x))
  def der5(self,A0,alpha,fi,nu,LamT,LamL,x):
      return A0*(1.0-alpha)*(-x)*np.exp(-LamL*x)
	  
    
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    alpha = self.getParameterValue("alpha");
    fi = self.getParameterValue("fi");
    nu = self.getParameterValue("nu");
    LamT = self.getParameterValue("LamT");
    LamL = self.getParameterValue("LamL");   
    
    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,alpha,fi,nu,LamT,LamL,x))
      jacobian.set(i,1, self.der1(A0,alpha,fi,nu,LamT,LamL,x))
      jacobian.set(i,2, self.der2(A0,alpha,fi,nu,LamT,LamL,x))
      jacobian.set(i,3, self.der3(A0,alpha,fi,nu,LamT,LamL,x))
      jacobian.set(i,4, self.der4(A0,alpha,fi,nu,LamT,LamL,x))
      jacobian.set(i,5, self.der5(A0,alpha,fi,nu,LamT,LamL,x))
      

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
  def alpha(self):
    return self.getParameterValue("alpha")
  def fi(self):
    return self.getParameterValue("fi")
  def nu(self):
    return self.getParameterValue("nu")
  def LamT(self):
    return self.getParameterValue("LamT")
  def LamL(self):
    return self.getParameterValue("LamL")    
    
  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)
  def setalpha(self, new_alpha):
    self.setParameter("alpha",new_alpha)
  def setfi(self, new_fi):
    self.setParameter("fi", new_fi)
  def setnu(self, new_nu):
    self.setParameter("nu", new_nu)
  def setLamT(self, new_LamT):
    self.setParameter("LamT", new_LamT)
  def setLamL(self, new_LamL):
    self.setParameter("LamL", new_LamL)

	# Required to have Mantid recognise the new function
FunctionFactory.subscribe(InternBessel)