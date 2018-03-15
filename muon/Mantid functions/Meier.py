from mantid.api import *
import math
import numpy as np

class Meier(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("Amp",0.1)
    self.declareParameter("OmegaD",0)
    self.declareParameter("OmegaQ",0)
    self.declareParameter("Spin",0)
    self.declareParameter("Sigma",0)
    self.declareParameter("Lambda",0)
       
  def function1D(self, x):
    Amp = self.getParameterValue("Amp")
    OmegaD = self.getParameterValue("OmegaD")
    OmegaQ = self.getParameterValue("OmegaQ")
    J = self.getParameterValue("Spin")
    L = self.getParameterValue("Lambda")
    sigma = self.getParameterValue("Sigma")
    gau = np.exp(-0.5*sigma*sigma*x*x)
    Lor = np.exp(-L*x)
    
  def W(m):
     return np.sqrt( np.power(2.0*m-1.0,2.0)*np.power(OmegaD+OmegaQ,2.0) + OmegaD*OmegaD*(J*(J+1.0)-m*(m-1.0)) )
     
    
  def lamp(m):
     # it would be none as it does not exist outside of the if
     lampValue = 0
     if m == J+1.0:
         lampValue= (OmegaQ*J*J)-(OmegaD*J)
     else:
         lampValue= 0.5*(OmegaQ*(2.0*m*m-2.0*m+1.0)+OmegaD+W(m))
     return lampValue
    
  def lamm(m):
     lammValue = 0
     if m == -J:
         lammValue= (OmegaQ*J*J)-(OmegaD*J)
     else:
         lammValue= 0.5*(OmegaQ*(2.0*m*m-2.0*m+1.0)+OmegaD-W(m))
     return lammValue
    
  def alpha(m):
    alphaValue = 0
    alphaValue= 0.5*np.arctan( OmegaD*np.sqrt(J*(J+1.0)-m*(m-1.0))/(1.0-2.0*m)*(OmegaD+OmegaQ) )
    return alphaValue
    #you will never get to the following code. The return exits the funnction. Should it be in an if?   
	
    tz=0
    tx=0
    for mm in range(-J+1,J+1):
        tz= tz + np.power(np.cos(2.0*alpha(mm)),2.0) + np.power(np.sin(2.0*alpha(mm)),2.0)*np.cos((lamp(mm)-lamm(mm))*x)
    # I assume this next line is not in the for loop
    Pz=(1.0/(2.0*J+1.0))*(1.0+tz)
    
    for mm in range(-J,J+1):
       a= np.power(np.cos(alpha(mm+1)),2.0) * np.power(np.sin(alpha(mm)),2.0) * np.cos((lamp(mm+1.0)-lamp(mm))*x)
       b= np.power(np.cos(alpha(mm+1)),2.0) * np.power(np.cos(alpha(mm)),2.0) * np.cos((lamp(mm+1.0)-lamm(mm))*x)
       c= np.power(np.sin(alpha(mm+1)),2.0) * np.power(np.sin(alpha(mm)),2.0) * np.cos((lamm(mm+1.0)-lamp(mm))*x)
       d= np.power(np.sin(alpha(mm+1)),2.0) * np.power(np.cos(alpha(mm)),2.0) * np.cos((lamm(mm+1.0)-lamm(mm))*x)
       tx= tx + a + b + c + d 
    Px= (1.0/(2.0*J+1.0))*tx
    
    return Amp*gau*Lor*(1.0/3.0)*(2.0*Px+Pz)
    
    
    
    
  # def functionDerivLocal(self, x, jacobian):
    # Amp = self.getParameterValue("Amp");
    # OmegaD = self.getParameterValue("OmegaD");
    # sigma = self.getParameterValue("Sigma")
    # Lambda = self.getParameterValue("Lambda")
    # gau = np.exp(-0.5*sigma*sigma*x*x)
    # Lor = np.exp(-L*x)
    # Posc = (1+np.cos(Omegad*x)+2*np.cos(0.5*OmegaD.x)+2*np.cos(1.5*OmegaD*x))/6
    # der0 = gau*Lor*Posc    
    # der1 = Amp*gau*Lor*-1*(np.sin(OmegaD*x)+x*np.sin(0.5*OmegaD*x) +3*x*np.sin(1.5*OmegaD*x))/6
    # der2 = Amp* (-sigma*x*x*gau) *Lor*Posc
    # der3 = Amp*gau* (-x*Lor)*Posc
    
    # # X index
    # i = 0
    # for x in xvals:
      # jacobian.set(i,0, der0)
      # jacobian.set(i,1, der1)
      # jacobian.set(i,2, der2)
      # jacobian.set(i,3, der3)
      # i += 1

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
  def OmegaQ(self):
    return self.getParameterValue("OmegaQ")
  def Spin(self):
    return self.getParameterValue("Spin")
  def Lambda(self):
    return self.getParameterValue("Lambda")
  def Sigma(self):
    return self.getParameterValue("Sigma")
    
  def setAmp(self, new_Amp):
    self.setParameter("Amp",new_Amp)
  def setOmegaD(self, new_OmegaD):
    self.setParameter("OmegaD", new_OmegaD)
  def setOmegaQ(self, new_OmegaQ):
    self.setParameter("OmegaQ", new_OmegaQ)    
  def setSpin(self, new_Spin):
    self.setParameter("Spin", new_Spin)
  def setlambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)
  def setSigma(self, new_Sigma):
    self.setParameter("Sigma", new_Sigma)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(Meier)