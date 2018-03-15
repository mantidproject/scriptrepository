from mantid.api import *
import math
import numpy as np

class StretchedKT(IFunction1D):
   # 
   
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("A0",1.0)
    self.declareParameter("beta",0.1)
    self.declareParameter("Sigma",0.1)
    
  def function1D(self, x):
    A0 = self.getParameterValue("A0")
    beta = self.getParameterValue("beta")
    Sigma = self.getParameterValue("Sigma")
    return A0*( (1.0/3.0)+(2.0/3.0)*(1.0-( np.power(Sigma*x,beta) )*np.exp(-1.0*np.power(Sigma*x,beta)/beta)) )

	
  def der0(self,A0,beta,Sigma,x):	
    e1 = np.exp( -1.0*( np.power(Sigma*x,beta)/beta ))
    f1 = np.power(Sigma*x, beta) 
    return ( (1.0/3.0)+(2.0/3.0)*(1.0-( np.power(Sigma*x,beta) )*np.exp(-1.0*np.power(Sigma*x,beta)/beta)) )
	
  def der1(self,A0,beta,Sigma,x):	
    e1 = np.exp( -1.0*( np.power(Sigma*x,beta)/beta ))
    f1 = np.power(Sigma*x, beta)	
    return A0*(-2.0/3.0)*e1*f1*( ( (f1/(beta*beta))-( (f1*np.log(Sigma*x))/beta) ) +np.log(Sigma*x) )
	
  def der2(self,A0,beta,Sigma,x):	
    e1 = np.exp( -1.0*( np.power(Sigma*x,beta)/beta ))
    f1 = np.power(Sigma*x, beta)	
    return A0*2.0*e1*f1*(f1-beta)/(3.0*Sigma)
    
	
  def functionDeriv1D(self, xvals, jacobian):
    A0 = self.getParameterValue("A0");
    beta = self.getParameterValue("beta");
    Sigma = self.getParameterValue("Sigma");
 
    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(A0,beta,Sigma,x))
      jacobian.set(i,1, self.der1(A0,beta,Sigma,x))
      jacobian.set(i,2, self.der2(A0,beta,Sigma,x))


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
  def beta(self):
    return self.getParameterValue("beta")
  def Sigma(self):
    return self.getParameterValue("Sigma")
  
  def setA0(self, new_A0):
    self.setParameter("A0",new_A0)
  def setbeta(self, new_beta):
    self.setParameter("beta", new_beta)
  def setSigma(self, new_Sigma):
    self.setParameter("Sigma", new_Sigma)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(StretchedKT)