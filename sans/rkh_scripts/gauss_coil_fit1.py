# RKH's Gaussian coil models Summer 2019
from mantid.api import IFunction1D, FunctionFactory
import numpy as np
import math as math
class monoGaussCoil(IFunction1D):
    """
        Provide a monodisperse Gaussian Coil function for SANS  (- hacked from Guinier & GuinierPorod equation functions, RKH, 23/8/19)

        I(q) = I(0)*2.0*( exp(-X) -1 +X)/(X)^2  where X = q^2.Rg^2 
    """

    def category(self):
        return "SANS"

    def init(self):
        # Active fitting parameters
        self.declareParameter("I0", 50.0, 'I(Q=0)')
        self.declareParameter("Rg", 60.0, 'Radius of gyration')
        self.declareParameter("Background", 0.2, 'Flat background')

    def _monoGaussCoil_core(self, qval):
        """
            Compute the main function for the model
            @param qval: q-value to evaluate at
        """
        Rg = self.getParameterValue('Rg')
        qr = qval*qval*Rg*Rg
        # Taylor expansion from mono_Gauss_Coil in sasview
        # (this is OK in single precision, but there is a better expansion at qr<0.5 for double precision)

        if  qr<0.8:
          C0 = +1.
          C1 = -1./3.
          C2 = +1./12.
          C3 = -1./60.
          C4 = +1./360.
          C5 = -1./2520.
          C6 = +1./20160.
          C7 = -1./181440.
          return ((((((C7*qr + C6)*qr + C5)*qr + C4)*qr + C3)*qr + C2)*qr + C1)*qr + C0

        return 2.0*(math.exp( -qr) -1.0 + qr)/(qr*qr)

    def _first_derivative_rg(self, qval):
        """
            Compute the first derivative dI/d(Rg)
            @param qval: q-value to evaluate at
        """
        Rg = self.getParameterValue('Rg')
        qr = qval*qval*Rg*Rg
        if  qr<0.8:
          C1 = -1./3.
          C2 = +2./12.
          C3 = -3./60.
          C4 = +4./360.
          C5 = -5./2520.
          C6 = +6./20160.
          C7 = -7./181440.
          return 2.0*qval*qval*Rg*((((((C7*qr + C6)*qr + C5)*qr + C4)*qr + C3)*qr +C2)*qr + C1)

        return 4.0*qval*qval*Rg*( (math.exp( -qr) -1.0)*(-2.0/qr -1.0)/(qr*qr) -2.0/qr )


    def function1D(self, xvals):
        """
            Evaluate the model
            @param xvals: numpy array of q-values
        """
        # parameters
        I0 = self.getParameterValue('I0')
        bgd = self.getParameterValue('Background')

        output = np.zeros(len(xvals), dtype=float)
        for i,x in enumerate(xvals):
            output[i] = I0 * self._monoGaussCoil_core(x) + bgd
        return output

    def functionDeriv1D(self, xvals, jacobian):
        """
            Evaluate the first derivatives
            @param xvals: numpy array of q-values
            @param jacobian: Jacobian object
        """
        i = 0
        for x in xvals:
			# assuming these are in order I0, Rg, background (though they come out in a different order in Mantid parameters table!.)
             jacobian.set(i,0, self._monoGaussCoil_core(x))
             jacobian.set(i,1, self._first_derivative_rg(x))
             jacobian.set(i,2, 1.0)
             i += 1
			 
			 
class polyGaussCoil(IFunction1D):
	"""
		Provide a polydisperse Gaussian Coil function for SANS  (python code copied from sasview, RKH 28/8/19)

		I(q) = I(0)*2.0*[ (1+UZ)^(-1/U) +Z-1 ]/[ (1+U)Z^2 ] where Z = q^2.Rg^2/( 1 +2U)  and U = (Mw/Mn) -1 
	"""

	def category(self):
		return "SANS"

	def init(self):
		# Active fitting parameters
		self.declareParameter("I0", 50.0, 'I(Q=0)')
		self.declareParameter("Rg", 60.0, 'Radius of gyration')
		self.declareParameter("Mw_by_Mn", 1.02, 'polydispersity')
		self.declareParameter("Background", 0.2, 'Flat background')

	def _polyGaussCoil_core(self, Rg, polydispersity, q):
		"""
		 Compute the main function for the model
		 """
		
		u = polydispersity - 1.0
		z = q**2 * (Rg**2 / (1.0 + 2.0*u))

		# need to trap the case of the polydispersity being 1 (ie, monodisperse!)
		if polydispersity == 1.0:
			#index = q != 0.
			#result[index] /= z[index]**2
			#result[~index] = 1.0
			if q != 0. : return  2.0 * (expm1(-z) + z)/z**2
			return 1.0
		else:
			#index = z > 1e-4
			#result[index] /= z[index]**2
			#result[~index] = np.polyval(p, z[~index])
			if z > 1e-04 : return  (2.0 * ( np.power(1.0 + u*z, -1.0/u) + z - 1.0) / (1.0 + u) )/z**2
			# Taylor series around z=0 of (2*(1+uz)^(-1/u) + z - 1) / (z^2(u+1))
			p = [
				#(-1 - 20*u - 155*u**2 - 580*u**3 - 1044*u**4 - 720*u**5) / 2520.,
				#(+1 + 14*u + 71*u**2 + 154*u**3 + 120*u**4) / 360.,
				#(-1 - 9*u - 26*u**2 - 24*u**3) / 60.,
				(+1 + 5*u + 6*u**2) / 12.,
				(-1 - 2*u) / 3.,
				(+1),
				]
			return   np.polyval(p, z)


	def _first_derivative_rg(self, qval):
		"""
			Compute the first derivative dI/d(Rg)
			- do this numerically for now, ought to work out the analytic equations ....
		"""
		I0 = self.getParameterValue('I0')
		Rg = self.getParameterValue('Rg')
		polydispersity = self.getParameterValue('Mw_by_Mn')
		Rg_up = min(Rg*1.005,Rg+0.2)
		return I0 * (self._polyGaussCoil_core(Rg_up, polydispersity, qval) - self._polyGaussCoil_core(Rg, polydispersity, qval)  )/(Rg_up - Rg)
		
	def _first_derivative_polydispersity(self, qval):
		"""
			Compute the first derivative dI/d(polydispersity)
			- do this numerically for now, ought to work out the analytic equations ....
		"""
		I0 = self.getParameterValue('I0')
		Rg = self.getParameterValue('Rg')
		polydispersity = self.getParameterValue('Mw_by_Mn')
		polydispersity_up = min(polydispersity*1.005,polydispersity+0.01)
		return I0 * (self._polyGaussCoil_core(Rg, polydispersity_up, qval) - self._polyGaussCoil_core(Rg, polydispersity, qval)  )/(polydispersity_up - polydispersity)

	def function1D(self, xvals):
		"""
			Evaluate the model
			@param xvals: numpy array of q-values
		"""
		# parameters
		I0 = self.getParameterValue('I0')
		Rg = self.getParameterValue('Rg')
		polydispersity = self.getParameterValue('Mw_by_Mn')
		bgd = self.getParameterValue('Background')
		#print I02,Rg2,polydispersity,bgd

		output = np.zeros(len(xvals), dtype=float)
		for i,x in enumerate(xvals):
				output[i] = bgd + I0 * self._polyGaussCoil_core(Rg, polydispersity, x)
		return output

	def functionDeriv1D(self, xvals, jacobian):
		"""
			Evaluate the first derivatives
			@param xvals: numpy array of q-values
			@param jacobian: Jacobian object
		"""
		i = 0
		for x in xvals:
			 jacobian.set(i,0, self._polyGaussCoil_core( self.getParameterValue('Rg'), self.getParameterValue('Mw_by_Mn'), x))
			 jacobian.set(i,1, self._first_derivative_rg(x))
			 jacobian.set(i,2, self._first_derivative_polydispersity(x))
			 jacobian.set(i,3, 1.0)
			 i += 1

# Required to have Mantid recognise the new function
FunctionFactory.subscribe(monoGaussCoil)
FunctionFactory.subscribe(polyGaussCoil)

