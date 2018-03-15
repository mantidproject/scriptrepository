from mantid.api import *
import math
import numpy as np

class CombGaussLorenKT(IFunction1D):
   # 
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("Lambda",0.1)
    self.declareParameter("Sigma",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    Lambda = self.getParameterValue("Lambda")
    Sigma = self.getParameterValue("Sigma")
    out= A0*( (1.0/3.0) + (2.0/3.0) * (1.0-(Sigma**2)*(x**2)-(Lambda*x) ) * np.exp(-0.5*(Sigma**2)*(x**2)-(Lambda*x)) )
    return out
    
  def der0(self,A0,Lambda,Sigma,x):
	  return ( (1.0/3.0) + (2.0/3.0) * (1.0-(Sigma**2)*(x**2)-(Lambda*x) ) * np.exp(-0.5*(Sigma**2)*(x**2)-(Lambda*x)) )
  def der1(self,A0,Lambda,Sigma,x):
      return A0*(2.0/3.0)* ( (Lambda*x + (Sigma**2)*(x**2) - 2.0) * x * np.exp(-0.5*(Sigma**2)*(x**2)-(Lambda*x)) )
  def der2(self,A0,Lambda,Sigma,x):
      return A0*(2.0/3.0)*(Sigma*(x**2)*((Lambda*x)+(Sigma**2)*(x**2) -3.0) * np.exp(-0.5*(Sigma**2)*(x**2)-(Lambda*x)))
    
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    Lambda = self.getParameterValue("Lambda");
    Sigma = self.getParameterValue("Sigma");

    # X index
    #i = 0
    for i, x in enumerate(xvals):
      jacobian.set(i,0,self.der0(A0,Lambda,Sigma,x))
      jacobian.set(i,1,self.der1(A0,Lambda,Sigma,x))
      jacobian.set(i,2,self.der2(A0,Lambda,Sigma,x))
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

  def Lambda(self):
    return self.getParameterValue("Lambda")
    
  def Sigma(self):
    return self.getParameterValue("Sigma")
  
  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)

  def setLambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)
    
  def setSigma(self, new_Sigma):
    self.setParameter("Sigma", new_Sigma)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(CombGaussLorenKT)