from mantid.api import *
import math
import numpy as np

class MuH(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("Amp",0.1)
    self.declareParameter("OmegaD",1.0)
    self.declareParameter("Lambda",0.0)
    self.declareParameter("Sigma",0.1)
       
  def function1D(self, x):
    Amp = self.getParameterValue("Amp")
    OmegaD = self.getParameterValue("OmegaD")
    L = self.getParameterValue("Lambda")
    Sigma = self.getParameterValue("Sigma")
    OmegaD = OmegaD*2.0*np.pi
    gau = np.exp(-0.5*Sigma*Sigma*x*x)
    Lor = np.exp(-L*x)
    return Amp*gau*Lor*(1.0+np.cos(OmegaD*x)+2.0*np.cos(0.5*OmegaD*x)+2.0*np.cos(1.5*OmegaD*x))/6.0

  
  def der0(self,Amp,OmegaD,Sigma,Lambda,x):	
      OmegaD = OmegaD*2.0*np.pi
      gau = np.exp(-0.5*Sigma*Sigma*x*x)
      Lor = np.exp(-Lambda*x)
      Posc = (1.0+np.cos(OmegaD*x)+2.0*np.cos(0.5*OmegaD*x)+2.0*np.cos(1.5*OmegaD*x))/6.0
      return gau*Lor*Posc  

  def der1(self,Amp,OmegaD,Sigma,Lambda,x):	
      OmegaD = OmegaD*2.0*np.pi
      gau = np.exp(-0.5*Sigma*Sigma*x*x)
      Lor = np.exp(-Lambda*x)
      Posc = (1.0+np.cos(OmegaD*x)+2.0*np.cos(0.5*OmegaD*x)+2.0*np.cos(1.5*OmegaD*x))/6.0	
      return Amp*gau*Lor*-1.0*(x*np.sin(OmegaD*x)+x*np.sin(0.5*OmegaD*x)+3.0*x*np.sin(1.5*OmegaD*x))/6.0
    
  def der2(self,Amp,OmegaD,Sigma,Lambda,x):	
      OmegaD = OmegaD*2.0*np.pi
      gau = np.exp(-0.5*Sigma*Sigma*x*x)
      Lor = np.exp(-Lambda*x)
      Posc = (1.0+np.cos(OmegaD*x)+2.0*np.cos(0.5*OmegaD*x)+2.0*np.cos(1.5*OmegaD*x))/6.0
      return Amp*(-1.0*Sigma*x*x*gau)*Lor*Posc
	
  def der3(self,Amp,OmegaD,Sigma,Lambda,x):	
      OmegaD = OmegaD*2.0*np.pi
      gau = np.exp(-0.5*Sigma*Sigma*x*x)
      Lor = np.exp(-Lambda*x)
      Posc = (1.0+np.cos(OmegaD*x)+2.0*np.cos(0.5*OmegaD*x)+2.0*np.cos(1.5*OmegaD*x))/6.0		
      return Amp*gau*(-1.0*x*Lor)*Posc	
	
    
  def functionDeriv1D(self, xvals, jacobian):
    Amp = self.getParameterValue("Amp");
    OmegaD = self.getParameterValue("OmegaD");
    Sigma = self.getParameterValue("Sigma")
    Lambda = self.getParameterValue("Lambda") 

    for i, x in enumerate(xvals):
      jacobian.set(i,0, self.der0(Amp,OmegaD,Sigma,Lambda,x))
      jacobian.set(i,1, self.der1(Amp,OmegaD,Sigma,Lambda,x))
      jacobian.set(i,2, self.der2(Amp,OmegaD,Sigma,Lambda,x))
      jacobian.set(i,3, self.der3(Amp,OmegaD,Sigma,Lambda,x))
     

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
  
  def Amp(self):
    return self.getParameterValue("Amp")

  def OmegaD(self):
    return self.getParameterValue("OmegaD")
  
  def Lambda(self):
    return self.getParameterValue("Lambda")
    
  def Sigma(self):
    return self.getParameterValue("Sigma")



  def setAmp(self, new_Amp):
    self.setParameter("Amp",new_Amp)

  def setOmegaD(self, new_OmegaD):
    self.setParameter("OmegaD", new_OmegaD)

  def setlambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)

  def setSigma(self, new_Sigma):
    self.setParameter("Sigma", new_Sigma)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(MuH)