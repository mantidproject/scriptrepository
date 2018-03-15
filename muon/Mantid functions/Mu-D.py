from mantid.api import *
import math
import numpy as np

class MuD(IFunction1D):
    
  def category(self):
    return "Muon"

  def init(self):
    self.declareParameter("Amp",0.1)
    self.declareParameter("OmegaD",1.0)
    self.declareParameter("OmegaQ",1.0)
    self.declareParameter("beta",1.57)
    self.declareParameter("Lambda",0.2)
    self.declareParameter("Sigma",0.0)
       
  def function1D(self, x):
    Amp = self.getParameterValue("Amp")
    OmegaD = self.getParameterValue("OmegaD")
    OmegaQ = self.getParameterValue("OmegaQ")
    OmegaD = OmegaD*2.0*np.pi
    OmegaQ = OmegaQ*2.0*np.pi
    beta = self.getParameterValue("beta")
    L = self.getParameterValue("Lambda")
    sigma = self.getParameterValue("Sigma")
    gau = np.exp(-0.5*sigma*sigma*x*x)
    Lor = np.exp(-L*x)
    
    k=OmegaQ*(3.0*(np.cos(beta)*np.cos(beta) -1.0))/4.0;
    w1=np.sqrt(6.0*OmegaD*OmegaD+4.0*OmegaD*k+(k*k));
    w2=(k-6.0*OmegaD+w1)/2.0;
    w3=(k-6.0*OmegaD-w1)/2.0;
    A0=5.0*k*k+4.0*OmegaD*OmegaD+20.0*OmegaD*k;
    A1=6.0*OmegaD*OmegaD;
    A2=2.0*(w1*w1-k*w1-2.0*OmegaD*w1);
    A3=2.0*(w1*w1+k*w1+2.0*OmegaD*w1);

    return Amp*gau*Lor*(1.0/(6.0*w1*w1))*(A0+A1*np.cos(w1*x)+A2*np.cos(w2*x)+A3*np.cos(w3*x));
    # NB: jacobian not included!!!
  
  #def functionDerivLocal(self, xvals, jacobian):
  #  Amp = self.getParameterValue("Amp");
  #  OmegaD = self.getParameterValue("OmegaD");
  #  sigma = self.getParameterValue("Sigma")
  #  Lambda = self.getParameterValue("Lambda")
  #  gau = np.exp(-0.5*sigma*sigma*x*x)
  #  Lor = np.exp(-L*x)
  #  Posc = (1+np.cos(Omegad*x)+2*np.cos(0.5*OmegaD.x)+2*np.cos(1.5*OmegaD*x))/6
  #  der0 = gau*Lor*Posc    
  #  der1 = Amp*gau*Lor*-1*(np.sin(OmegaD*x)+x*np.sin(0.5*OmegaD*x) +3*x*np.sin(1.5*OmegaD*x))/6
  #  der2 = Amp* (-sigma*x*x*gau) *Lor*Posc
  #  der3 = Amp*gau* (-x*Lor)*Posc
    
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
  
  def Amp(self):
    return self.getParameterValue("Amp")
  def OmegaD(self):
    return self.getParameterValue("OmegaD")
  def OmegaQ(self):
    return self.getParameterValue("OmegaQ")
  def beta(self):
    return self.getParameterValue("beta")
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
  def setbeta(self, new_beta):
    self.setParameter("beta", new_beta)
  def setlambda(self, new_Lambda):
    self.setParameter("Lambda", new_Lambda)
  def setSigma(self, new_Sigma):
    self.setParameter("Sigma", new_Sigma)

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(MuD)