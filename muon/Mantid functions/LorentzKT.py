from mantid.api import *
import math
import numpy as np

class StaticLorentzKT(IFunction1D):
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("Lambda",0.1)
       
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    Lambda = self.getParameterValue("Lambda")
    return A0*( (1.0/3.0)+(2.0/3.0)*(1.0-Lambda*x)*np.exp(-Lambda*x) )

  def der0(self,A0,Lambda,x):
      return ( (1.0/3.0)+(2.0/3.0)*(1.0-(Lambda*x))*np.exp(-Lambda*x) )
  def der1(self,A0,Lambda,x):
      return A0*(2.0/3.0)*x*np.exp(-Lambda*x)*((Lambda*x)-2.0)

    
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    Lambda = self.getParameterValue("Lambda");

    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,Lambda,x))
      jacobian.set(i,1, self.der1(A0,Lambda,x))


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
  
  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)
  def setLambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(StaticLorentzKT)