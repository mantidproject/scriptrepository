import QENSmodels
import matplotlib.pyplot as plt
import numpy as np


#make fake data
dataX = np.linspace(-5,5,200)
Q = np.linspace(0.2,2,10)
chudley_elliot_noisy = QENSmodels.sqwChudleyElliotDiffusion(dataX, Q,scale=1.3, center=0.22, D=0.43,
                              L=0.88)*(1+0.1*np.random.normal(0,1,200)) + 0.01*np.random.normal(0,1,200)
                              
# store in mantid ws
QENS_data = CreateWorkspace(DataX=dataX, DataY=chudley_elliot_noisy, NSpec=10)


#wrap to create mantid fitting function
class sqwChudleyElliotDiffusion(IFunction1D): # or IPeakFunction

   def init(self):
        self.declareParameter("scale", 1.0)
        self.declareParameter("center", 0.0)
        self.declareParameter("D", 0.23)
        self.declareParameter("L", 1.0)
        
        self.declareParameter("Q", 0.1)
        

   def category(self):
       return 'QuasiElastic'
           
   def function1D(self, xvals):
         scale = self.getParameterValue("scale")
         center = self.getParameterValue("center")
         D = self.getParameterValue("D")
         L = self.getParameterValue("L")
         q = self.getParameterValue("Q")
         return QENSmodels.sqwChudleyElliotDiffusion(xvals, q, scale=scale, center=center, D=D, L=L)
# add it to Mantid fitting functions
FunctionFactory.subscribe(sqwChudleyElliotDiffusion)