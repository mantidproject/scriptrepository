## Module containing lower level functions for Quantum
## Author: James Lord
## Version 1.02, December 2015
import numpy
import math
import time
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
from collections import Counter
import re
import numbers
import bisect
numpy.set_printoptions(linewidth=100)

def normalised(vec):
	# intended for 3 element 1d vectors
	r=numpy.linalg.norm(vec)
	return vec/r

def addZeeman(Ham,spin,B,gamma):
	# Ham= Hamiltonian matrix ( (2I+1)*(2S+1)*(2J+1) , (2I+1)*(2S+1)*(2J+1) ) for spins I,S,J modified in place
	# spin = spin operator matrix for this spin ( (2I+1)*(2S+1)*(2J+1) , (2I+1)*(2S+1)*(2J+1) , 3)
	# gamma = magnetic moment (MHz/T)
	# B (3-vector) not necessarily normalised
	#print "first operand of dot: ",spin
	#print "second operand of dot: ",B
	#print "dot product: ",numpy.tensordot(spin,normalised(B),axes=[[0],[0]])
	# force modification in place!
	if(isinstance(gamma,numbers.Number)):
		Ham[...]=Ham[...]+numpy.tensordot(spin,normalised(B),axes=[[0],[0]])*gamma*math.pi
	else:
		# gamma is a g-value tensor (multiplied by |B|). Calc...?
		eb=numpy.dot(gamma,normalised(B))
		Ham[...]=Ham[...]+numpy.tensordot(spin,eb,axes=[[0],[0]])*math.pi
	#print "modified to ",Ham
#	return Ham

def addZeemanC(HamIn,HamOut,spin,B,gamma):
	# Ham= Hamiltonian matrix ( (2I+1)*(2S+1)*(2J+1) , (2I+1)*(2S+1)*(2J+1) ) for spins I,S,J modified in place
	# spin = spin operator matrix for this spin ( (2I+1)*(2S+1)*(2J+1) , (2I+1)*(2S+1)*(2J+1) , 3)
	# gamma = magnetic moment (MHz/T)
	# B (3-vector) not necessarily normalised
	#print "first operand of dot: ",spin
	#print "second operand of dot: ",B
	#print "dot product: ",numpy.tensordot(spin,normalised(B),axes=[[0],[0]])
	# fill supplied output array
	if(isinstance(gamma,numbers.Number)):
		HamOut[...]=HamIn[...]+numpy.tensordot(spin,normalised(B),axes=[[0],[0]])*gamma*math.pi
	else:
		# gamma is a g-value tensor too. Calc...?
		eb=numpy.dot(gamma,normalised(B))
		HamOut[...]=HamIn[...]+numpy.tensordot(spin,eb,axes=[[0],[0]])*math.pi
	#print "modified to ",Ham
#	return Ham

def createSpinMat(spins):
	# spins = array of 2I+1
	# returns matrix (i,alpha,N,N) where alpha=(x,y,z) and i iterates through spins
	N=numpy.prod(spins)
	Nsp=len(spins)
	if(N>100000):
		raise Exception("Problem almost certainly too big to try")
	if(N>1024):
		print "fairly big problem, prepare for out of memory errors"
	# print ([Nsp,3]+spins+spins)
	spmat=numpy.zeros((Nsp,3,N,N),dtype=complex)
	# fill me (exact copy of Fortran QUANTUM including possible scale errors)
	for thissp in range(Nsp):
		for lowerIz in range(numpy.prod(spins[0:thissp],dtype=numpy.int)):
			for upperIz in range(lowerIz,N,numpy.prod(spins[0:thissp+1])):
				for thisIz in range(spins[thissp]):
					thisIndI=upperIz+thisIz*numpy.prod(spins[0:thissp])
					thisIndJ=upperIz+(thisIz+1)*numpy.prod(spins[0:thissp])
					# element [thissp,alpha,thisIndI,thisIndJ] flips spin from thisIz to thisJz leaving others alone
					spmat[thissp,2,thisIndI,thisIndI]=spins[thissp]-1-2*thisIz # Sz diagonal
					if(thisIz<spins[thissp]-1):
						r=math.sqrt((thisIz+1)*(spins[thissp]-thisIz-1))
						spmat[thissp,0,thisIndI,thisIndJ]=r
						spmat[thissp,0,thisIndJ,thisIndI]=r
						spmat[thissp,1,thisIndI,thisIndJ]=r*1.0j
						spmat[thissp,1,thisIndJ,thisIndI]=-r*1.0j
	return spmat

def AxialHFT(A,D,axis):
	# generate tensor. "axis" need not be normalised
	axisn=normalised(axis)
	rr=numpy.outer(axisn,axisn)
	return numpy.identity(3)*(A-D/2)+rr*(D*3/2)

def addHyperfine(Ham,spin1,spin2,A):
	# Ham=hamiltonian
	# spin1,spin2 (mat) which take part
	# A = 3*3 hyperfine tensor (or scalar for isotropic?) in MHz
	# break down problem?
	A=A*math.pi/2
	spin1x=spin1[0]
	spin1y=spin1[1]
	spin1z=spin1[2]
	spin2x=spin2[0]
	spin2y=spin2[1]
	spin2z=spin2[2]
	hxx=numpy.dot(spin1x,spin2x)
	hyy=numpy.dot(spin1y,spin2y)
	hzz=numpy.dot(spin1z,spin2z)
	if(numpy.isscalar(A)):
		Ham[...]=Ham[...]+(hxx+hyy+hzz)*A
	else:
		hxy=numpy.dot(spin1x,spin2y) # only one of each as operators should commute? spin1 != spin2. No, Sx1Sy2 != Sx2Sy1.
		hxz=numpy.dot(spin1x,spin2z)
		hyx=numpy.dot(spin1y,spin2x)
		hyz=numpy.dot(spin1y,spin2z)
		hzx=numpy.dot(spin1z,spin2x)
		hzy=numpy.dot(spin1z,spin2y)
		Ham[...]=Ham[...]+hxx*A[0,0]+hxy*A[0,1]+hyx*A[1,0]+hxz*A[0,2]+hzx*A[2,0]+hyy*A[1,1]+hyz*A[1,2]+hzy*A[2,1]+hzz*A[2,2]
	# hyper=numpy.einsum('ijk,il,lkm',spin1,A,spin2) # slow?
	# sum over k: get tmpmat(i=3,j=N,l=3,m=N)
	# then sum over i and l
	#tmp=numpy.dot(spin1,spin2)
	#hyper=numpy.tensordot(tmp,A,[[0,2],[0,1]])
	# if A is identity (isotropic) then hyper=A*einsum(ijk,ikm)
	#print "hf="
	#print hyper
	#print "."
	#Ham[...]=Ham[...]+hyper

def addHyperfineC(HamIn,HamOut,spin1,spin2,A):
	# Ham=hamiltonian
	# spin1,spin2 (mat) which take part
	# A = 3*3 hyperfine tensor (or scalar for isotropic?)
	if(numpy.isscalar(A)):
		A=numpy.identity(3)*A
	A=A*math.pi/2
	# hyper=numpy.einsum('ijk,il,lkm',spin1,A,spin2) # slow?
	# sum over k: get tmpmat(i=3,j=N,l=3,m=N)
	# then sum over i and l
	tmp=numpy.dot(spin1,spin2)
	hyper=numpy.tensordot(tmp,A,[[0,2],[0,1]])
	# if A is identity (isotropic) then hyper=A*einsum(ijk,ikm)
	#print "hf="
	#print hyper
	#print "."
	HamOut[...]=HamIn[...]+hyper

def addDipolar(Ham,spin1,spin2,r1,gamma1,r2,gamma2):
	# Ham=hamiltonian
	# spin1 (mat) with location r1 and moment gamma1, likewise spin2
	# units r in A, moments gamma/2pi in MHz/T, time in us
	hmu0=1.054571726E-5 * math.pi**2
	#print "gamma1=",gamma1," and gamma2=",gamma2
	t0=numpy.dot(spin1[0],spin2[0])+numpy.dot(spin1[1],spin2[1])+numpy.dot(spin1[2],spin2[2])
	#print "t0="
	#print t0
	r=r2-r1
	#print "r=",r
	rmag=math.sqrt(numpy.sum(r*r))
	#print "rmag=",rmag
	rn=r/rmag
	#print "rn=",rn
	Dmag=hmu0*gamma1*gamma2/rmag**3 # *3.0?
	#print "Dmag=",Dmag
	t1=rn[0]*spin1[0]+rn[1]*spin1[1]+rn[2]*spin1[2]
	#print "t1="
	#print t1
	t2=rn[0]*spin2[0]+rn[1]*spin2[1]+rn[2]*spin2[2]
	#print "t2="
	#print t2
	t3=numpy.dot(t1,t2)
	#print "t3="
	#print t3
	#print "Ham="
	#print Ham
	Ham[...]=Ham[...]+Dmag*(t0-3*t3)
	#print "now Ham="
	#print Ham
	#print "implemented"

def addQuadrupole(Ham,spin,nutensor):
	# tensor can use same algorithm as HFC.. maybe with scaling or offset?
	addHyperfine(Ham,spin,spin,nutensor)
	
def createInitialDensMat(spin,polvec):
	# fully polarised spin-1/2 in direction polvec, others unpolarised
	# special case "powder" if polvec=None?
	N=spin.shape[1]
	rho=(numpy.tensordot(spin,normalised(polvec),axes=[[0],[0]])+numpy.identity(N,dtype=numpy.complex))/N
	return rho

def createInitialDensMatWithSingletTriplet(spin,polvec,spinE1,spinE2,triplet):
	# fully polarised spin-1/2 in direction polvec
	# spins spinE1 and spinE2 are in singlet or triplet state according to triplet=false,true
	# others unpolarised
	# special case "powder" if polvec=None?
	N=spin.shape[1]
	Edot=numpy.dot(spinE1[0],spinE2[0])+numpy.dot(spinE1[1],spinE2[1])+numpy.dot(spinE1[2],spinE2[2])
	if(triplet):
		Ed2=Edot/3.0+numpy.identity(N)
	else: # singlet
		Ed2=numpy.identity(N)-Edot
	# S1.
	rho=numpy.dot((numpy.tensordot(spin,normalised(polvec),axes=[[0],[0]])+numpy.identity(N,dtype=numpy.complex)),(Ed2/N))
	return rho

def createDetectorOp(spin,polvec):
	# in direction polvec
	# matching special case "powder"?
	rho=numpy.tensordot(spin,normalised(polvec),axes=[[0],[0]])
	return rho

def solveDensityMat(Ham,start,detect,timeend=float("Inf"),tzeros=None):
	# diagonalise Ham
	# start with density matrix "start"
	# returns (omega,cos,sin) each is array of coefficients
	# with optional param timeend, return (omega,cos,sin,rho(timeend)) for use with pulsed mode, and with tzeros too, smooth rho for variable slice length
	eval,evec=numpy.linalg.eigh(Ham)
	evecstar=numpy.conj(evec)
	# transform rho into basis set of evec
	#* rho2=numpy.einsum('ji,jk,kl',evecstar,start,evec)
	rho1a=numpy.dot(start,evec) # start[j,k]*evec[k.l]
	rho2=numpy.tensordot(evecstar,rho1a,(0,0) )
	#print "rho2"
	#print rho2
	# likewise detector operator (CHECK conjg, etc)
	#* detect2=numpy.einsum('ji,jk,kl',evecstar,detect,evec)
	detect1a=numpy.dot(detect,evec,rho1a) # reuse rho1a, it's the right shape
	detect2=numpy.tensordot(evecstar,detect1a,(0,0) )
	#print "detect2"
	#print detect2
	#print numpy.dot(numpy.diag(rho2),numpy.diag(detect2))
	# detectedfreqs Eval(i)-Eval(j) ampl rhoij*Sji
	# detected const dot(diag(rho),diag(s))
	n=Ham.shape[0]
	nc=n*(n-1)/2+1
	omega=numpy.empty(nc,dtype=float)
	ccos=numpy.empty(nc,dtype=float)
	csin=numpy.empty(nc,dtype=float)
	omega[0]=0.0
	ccos[0]=numpy.real(numpy.dot(numpy.diag(rho2),numpy.diag(detect2)))
	csin[0]=0.0
	(ii,jj)=numpy.tril_indices(n,-1)
	omega[1:]=eval[ii]-eval[jj]
	rd1=rho2[ii,jj]*detect2[jj,ii]
	rd2=rho2[jj,ii]*detect2[ii,jj]
	ccos[1:]=numpy.real(rd1)+numpy.real(rd2)
	csin[1:]=numpy.imag(rd1)-numpy.imag(rd2) # check signs!
	#print "Eigenvectors"
	#print evec
	#print "Eigenvalues"
	#print eval
	if(math.isinf(timeend)):
		return (omega,ccos,csin)
	else:
		# calc rho at time=timeend
		# rho2 is in eigenvectors
		# rho2(t) = rho2 * phase factors
		# rho(std basis) by multiplying by Evec^-1
		phases=numpy.exp(1.j*eval*timeend)
		iphases=1./phases
		rho2=rho2*phases*iphases[:,numpy.newaxis]
		if(tzeros is not None):
			freqs=1.j*(eval-eval[:,numpy.newaxis])
			ConvolveTimeResDynListed(freqs,rho2,tzeros) # use Dynamic form to deal with complex coeffs
		#rho2[ii,jj]=rho2[ii,jj]*phases[ii]*iphases[jj]
		#rho2[jj,ii]=rho2[jj,ii]*phases[jj]*iphases[ii]
		rho3a=numpy.dot(evec,rho2)
		rho3=numpy.tensordot(rho3a,evecstar,(1,1))
		return (omega,ccos,csin,rho3)

def evaluateIntoBinsOldLoops(omega,ccos,csin,lam,times):
	# times = array length n+1 of time bin boundaries
	# returns array length n of values, averaged across bins
	# lam = 1/lifetime to weight by, if bins are large (or only one big one for integral counting - times [0,inf] is OK)
	# = integral (cos(w*t))exp(-lam*t)]t1:t2 / integral exp(-lam*t)]t1:t2
	# = (-lam* sin(w*t) exp(-lam*t) + w* cos(w*t) exp*(-lam*t) )/(lam^2+w^2)
	ybins=numpy.zeros(len(times)-1,dtype=numpy.float)
	if(lam>0):
		for i in range(len(times)-1):
			et2=math.exp(-lam*(times[i+1]-times[i]))
			sf=(1-et2)/lam
			for j in range(len(omega)):
				if(omega[j]==0):
					ybins[i]=ybins[i]+ccos[j]
				else:
					c1=math.cos(omega[j]*times[i])
					s1=math.sin(omega[j]*times[i])
					if(et2>0):
						c2=math.cos(omega[j]*times[i+1])*et2
						s2=math.sin(omega[j]*times[i+1])*et2
					else:
						# avoid domain problems with [..,inf] integral counting, would be multiplied by 0 anyway
						c2=0.0
						s2=0.0
					ybins[i]=ybins[i]+(ccos[j]*(omega[j]*(s2-s1 )-lam*(c2-c1 )) +csin[j]*(omega[j]*(-c2+c1 )-lam*(s2-s1 )) )/sf /(lam*lam+omega[j]*omega[j])
	else:
		for i in range(len(times)-1):
			sf=times[i+1]-times[i]
			for j in range(len(omega)):
				if(omega[j]==0):
					ybins[i]=ybins[i]+ccos[j]
				else:
					c1=math.cos(omega[j]*times[i])
					s1=math.sin(omega[j]*times[i])
					c2=math.cos(omega[j]*times[i+1])
					s2=math.sin(omega[j]*times[i+1])
					ybins[i]=ybins[i]+(ccos[j]*(1/omega[j])*(s2-s1 )+csin[j]*(1/omega[j])*(-c2+c1 ))/sf
	return ybins

def evaluateIntoBins(omega,ccos,csin,lam,times):
	# times = array length n+1 of time bin boundaries
	# returns array length n of values, averaged across bins
	# lam = 1/lifetime to weight by, if bins are large (or only one big one for integral counting - times [0,inf] is OK)
	# = integral (cos(w*t))exp(-lam*t)]t1:t2 / integral exp(-lam*t)]t1:t2
	# = (-lam* sin(w*t) exp(-lam*t) + w* cos(w*t) exp*(-lam*t) )/(lam^2+w^2)
	# improved version using time arrays internally (best for time series from fairly simple models)
	ybins=numpy.zeros(len(times)-1,dtype=numpy.float)
	if(lam>0):
#		for i in range(len(times)-1):
		et2=numpy.exp(-lam*(times[1:]-times[:-1]))
		sf=(1-et2)/lam
		for j in range(len(omega)):
			if(omega[j]==0):
				ybins=ybins+ccos[j]
			else:
				c0=numpy.cos(omega[j]*times)
				s0=numpy.sin(omega[j]*times)
				c1=c0[:-1]
				s1=s0[:-1]
				c2=numpy.nan_to_num(c0[1:]*et2)
				s2=numpy.nan_to_num(s0[1:]*et2) # avoid domain problems with [..,inf] integral counting, would be multiplied by 0 anyway
				ybins=ybins+(ccos[j]*(omega[j]*(s2-s1 )-lam*(c2-c1 )) +csin[j]*(omega[j]*(-c2+c1 )-lam*(s2-s1 )) )/sf /(lam*lam+omega[j]*omega[j])
	else:
#		for i in range(len(times)-1):
		sf=times[1:]-times[:-1]
		for j in range(len(omega)):
			if(omega[j]==0):
				ybins=ybins+ccos[j]
			else:
				c0=numpy.cos(omega[j]*times)
				s0=numpy.sin(omega[j]*times)
				c1=c0[:-1]
				s1=s0[:-1]
				c2=c0[1:]
				s2=s0[1:]
				ybins=ybins+(ccos[j]*(1/omega[j])*(s2-s1 )+csin[j]*(1/omega[j])*(-c2+c1 ))/sf
	return ybins

def evaluateIntoBinsNewLoops(omega,ccos,csin,lam,times):
	# times = array length n+1 of time bin boundaries
	# returns array length n of values, averaged across bins
	# lam = 1/lifetime to weight by, if bins are large (or only one big one for integral counting - times [0,inf] is OK)
	# = integral (cos(w*t))exp(-lam*t)]t1:t2 / integral exp(-lam*t)]t1:t2
	# = (-lam* sin(w*t) exp(-lam*t) + w* cos(w*t) exp*(-lam*t) )/(lam^2+w^2)
	# arrays across components: best for integral counting of complex systems
	ybins=numpy.zeros(len(times)-1,dtype=numpy.float)
	if(lam>0):
		for i in range(len(times)-1):
			et2=math.exp(-lam*(times[i+1]-times[i]))
			sf=(1-et2)/lam
			#for j in range(len(omega)):
			#	if(omega[j]==0):
			#		ybins[i]=ybins[i]+ccos[j]
			#	else:
			c1=numpy.cos(omega*times[i])
			s1=numpy.sin(omega*times[i])
			if(et2>0):
				c2=numpy.cos(omega*times[i+1])*et2
				s2=numpy.sin(omega*times[i+1])*et2
			else:
				# avoid domain problems with [..,inf] integral counting, would be multiplied by 0 anyway
				c2=0.0
				s2=0.0
			ybins[i]=ybins[i]+numpy.sum((ccos*(omega*(s2-s1 )-lam*(c2-c1 )) +csin*(omega*(-c2+c1 )-lam*(s2-s1 )) ) /(lam*lam+omega*omega))/sf
	else:
		for i in range(len(times)-1):
			sf=times[i+1]-times[i]
			#for j in range(len(omega)):
			#	if(omega[j]==0):
			#		ybins[i]=ybins[i]+ccos[j]
			#	else:
			ybins[i]=ccos[0] # assume omega[0]=0 i.e. the non oscillating component is first
			c1=numpy.cos(omega[1:]*times[i])
			s1=numpy.sin(omega[1:]*times[i])
			c2=numpy.cos(omega[1:]*times[i+1])
			s2=numpy.sin(omega[1:]*times[i+1])
			ybins[i]=ybins[i]+numpy.sum((ccos[1:]*(s2-s1)+csin[1:]*(c1-c2))/omega[1:])/sf
	return ybins

def ConvolveTimeResolution(omega,ccos,csin,sigma):
	# reduce ccos and csin amplitudes according to Gaussian pulse width / timing resolution FWHM sigma
	# i.e N(t)=exp(-(t/sigma)**2)
	af=numpy.exp(-(omega*sigma)**2/2)
	ccos[:]=ccos[:]*af
	csin[:]=csin[:]*af

def ConvolveTimeResolutionListed(omega,ccos,csin,funcs):
	# funcs is a list/tuple with pairs of values, 1st is key letter and 2nd is "sigma"
	for (code,sigma) in zip(*[iter(funcs)]*2):
		afc=None
		sigma=float(sigma)
		if(sigma==0):
			break
		if(code=="g"):
			# Gaussian
			# f(t)=(1/tau/sqrt2pi) exp(-(t/tau)^2/2) from -inf to +inf
			af=numpy.exp(-(omega*sigma)**2/2)
		elif(code=="l"):
			# Lorentzian
			# f(t)=(1/pi)*tau^2/(t^2+tau^2) from -inf to +inf
			af=numpy.exp(-sigma*numpy.abs(omega))
		elif(code=="p"):
			# parabolic
			# f(t)=(3/4tau)*(1-(t/tau)^2) from -tau to +tau
			so=sigma*omega
			af1=(3/so**3)*(numpy.sin(so)-so*numpy.cos(so))
			af2=1-so**2/10+so**4/280-so**6/15120 # + ... series expansion for omega*sigma ~ 0
			af=numpy.where(numpy.abs(so)>1.E-3,af1,af2)
		elif(code=="c"):
			# half Cosine
			# f(t)=(1/tau) cos(t/tau) from -tau*pi/2 to tau*pi/2
			af=numpy.cos(omega*sigma*math.pi/2)/(1-sigma**2*omega**2)
		elif(code=="u"):
			# Uniform square pulse
			# f(t)=1/tau from -tau/2 to tau/2
			so=omega*sigma/2
			af1=numpy.sin(so)/so
			af2=1-so**2/6+so**4/120 - so**6/5040 # ...
			af=numpy.where(numpy.abs(so)>1.E-3,af1,af2)
		elif(code=="e"):
			# exponential tail
			# f(t)=1/tau *exp(-t/tau) from 0 to +inf
			af=1./(1+sigma**2*omega**2)
			afc=-af*sigma*omega
		if(afc is None):
			ccos[:]=ccos*af
			csin[:]=csin*af
		else:
			tmp=ccos*af+csin*afc
			csin[:]=csin*af-ccos*afc
			ccos[:]=tmp
			
def ConvolveTimeResolutionDynamic(omega,ccossin,sigma):
	# reduce ccos and csin amplitudes according to Gaussian pulse width / timing resolution FWHM sigma
	# for dynamic form, ignore real(omega)=relaxation, use imag(omega)=frequency
	# i.e N(t)=exp(-(t/sigma)**2)
	af=numpy.exp(-(omega.imag*sigma)**2/2)
	ccossin[:]=ccossin[:]*af

def ConvolveTimeResDynListed(omega,ccossin,funcs):
	# funcs is a list/tuple with pairs of values, 1st is key letter and 2nd is "sigma"
	# dynamic, omega complex with Re(omega)<0=relaxation
	# only strictly accurate for t>end-of-pulse so use with care for pulse profiles with tails (Gaussian, exponential).
	# integral of Lorentz(x)*exp(-x/lambda) doesn't even converge for lambda != 0 !
	for (code,sigma) in zip(*[iter(funcs)]*2):
		sigma=float(sigma)
		if(sigma==0):
			break
		if(code=="g"):
			# Gaussian
			# f(t)=(1/tau/sqrt2pi) exp(-(t/tau)^2/2) from -inf to +inf
			af=numpy.exp((omega*sigma)**2/2)
			# warning, diverges for fast relaxing components unless t is large! Need cutoff?
		elif(code=="l"):
			# Lorentzian
			# f(t)=(1/pi)*tau^2/(t^2+tau^2) from -inf to +inf
			#af=numpy.exp(-sigma*numpy.abs(omega))
			if(numpy.any(numpy.real(omega)!=0)):
				raise Exception("Lorentzian does not converge")
			else:
				af=numpy.exp(-sigma*numpy.abs(numpy.imag(omega)))
		elif(code=="p"):
			# parabolic
			# f(t)=(3/4tau)*(1-(t/tau)^2) from -tau to +tau
			so=sigma*omega
			af1=(3/so**3)*(numpy.cosh(so)*so-numpy.sinh(so)) # division 0^3 by 0^3 risk!
			af2=1+so**2/10+so**4/280+so**6/15120 # + ... series expansion for omega*sigma ~ 0
			af=numpy.where(numpy.abs(so)>1.E-3,af1,af2)
		elif(code=="c"):
			# half Cosine
			# f(t)=(1/tau) cos(t/tau) from -tau*pi/2 to tau*pi/2
			af=numpy.cosh(omega*sigma*math.pi/2)/(1+sigma**2*omega**2)
			# = pi/4 when sigma*omega=i
		elif(code=="u"):
			# Uniform square pulse
			# f(t)=1/tau from -tau/2 to tau/2
			so=omega*sigma/2
			af1=numpy.sinh(so)/so
			af2=1+so**2/6+so**4/120+so**6/5040 # ...
			af=numpy.where(numpy.abs(so)>1.E-3,af1,af2)
		elif(code=="e"):
			# exponential tail
			# f(t)=1/tau *exp(-t/tau) from 0 to +inf
			af=1./(1+sigma*omega) # only for slow relaxation omega*Re(sigma)<1
		ccossin[:]=ccossin*af

# generator/iterator for regular orientation choices
def uniformLF(N):
	for i in range(N):
		# field axis, also for beam and detector
		z=(2*i+0.5)/N-1
		phi=i*2*math.pi/math.e
		s=math.sqrt(1-z*z)
		cp=math.cos(phi)
		x=cp*s
		sp=math.sin(phi)
		y=sp*s
		# rf axis rotated faster round field axis (relative to phi2=0 pointing at z axis)
		# "x" axis for rotating vector along tangent:
		z2=s
		x2=-cp*z
		y2=-sp*z
		# "y" axis
		# z3=0
		x3=sp
		y3=-cp
		phi2=i*2*math.pi*math.e
		cp2=math.cos(phi2)
		sp2=math.sin(phi2)
		x4=x2*cp2+x3*sp2
		y4=y2*cp2+y3*sp2
		z4=z2*cp2
		yield ([x,y,z],[x,y,z],[x,y,z],[x4,y4,z4])

def randomLF(N):
	for i in range(N):
		# first axis for B, beam, det
		z=numpy.random.random()*2-1
		phi=numpy.random.random()*2*math.pi
		s=math.sqrt(1-z*z)
		x=math.cos(phi)*s
		y=math.sin(phi)*s
		# second random axis
		z2=numpy.random.random()*2-1
		phi2=numpy.random.random()*2*math.pi
		s2=math.sqrt(1-z2*z2)
		x2=math.cos(phi2)*s2
		y2=math.sin(phi2)*s2
		# third axis perpendicular to 1st for RF
		x3=(y*z2-z*y2)
		y3=(z*x2-x*z2)
		z3=(x*y2-y*x2)
		r=math.sqrt(x3**2+y3**2+z3**2)
		yield ([x,y,z],[x,y,z],[x,y,z],[x3/r,y3/r,z3/r])

def uniformRF(N):
	for i in range(N):
		# field axis, also for beam
		z=(2*i+0.5)/N-1
		phi=i*2*math.pi/math.e
		s=math.sqrt(1-z*z)
		cp=math.cos(phi)
		x=cp*s
		sp=math.sin(phi)
		y=sp*s
		# rf and detector axis rotated faster round field axis (relative to phi2=0 pointing at z axis)
		# "x" axis for rotating vector along tangent:
		z2=s
		x2=-cp*z
		y2=-sp*z
		# "y" axis
		# z3=0
		x3=sp
		y3=-cp
		phi2=i*2*math.pi*math.e
		cp2=math.cos(phi2)
		sp2=math.sin(phi2)
		x4=x2*cp2+x3*sp2
		y4=y2*cp2+y3*sp2
		z4=z2*cp2
		yield ([x,y,z],[x,y,z],[x4,y4,z4],[x4,y4,z4])

def randomRF(N):
	for i in range(N):
		# first axis for B, beam
		z=numpy.random.random()*2-1
		phi=numpy.random.random()*2*math.pi
		s=math.sqrt(1-z*z)
		x=math.cos(phi)*s
		y=math.sin(phi)*s
		# second random axis
		z2=numpy.random.random()*2-1
		phi2=numpy.random.random()*2*math.pi
		s2=math.sqrt(1-z2*z2)
		x2=math.cos(phi2)*s2
		y2=math.sin(phi2)*s2
		# third axis perpendicular to 1st for RF, det
		x3=(y*z2-z*y2)
		y3=(z*x2-x*z2)
		z3=(x*y2-y*x2)
		r=math.sqrt(x3**2+y3**2+z3**2)
		yield ([x,y,z],[x,y,z],[x3/r,y3/r,z3/r],[x3/r,y3/r,z3/r])

def uniformTF(N):
	for i in range(N):
		# field axis as before
		z=(2*i+0.5)/N-1
		phi=i*2*math.pi/math.e
		s=math.sqrt(1-z*z)
		cp=math.cos(phi)
		x=cp*s
		sp=math.sin(phi)
		y=sp*s
		# beam axis rotated faster round field axis (relative to phi2=0 pointing at z axis)
		# "x" axis for rotating vector along tangent:
		z2=s
		x2=-cp*z
		y2=-sp*z
		# "y" axis
		# z3=0
		x3=sp
		y3=-cp
		phi2=i*2*math.pi*math.e
		cp2=math.cos(phi2)
		sp2=math.sin(phi2)
		x4=x2*cp2+x3*sp2
		y4=y2*cp2+y3*sp2
		z4=z2*cp2
		yield ([x,y,z],[x4,y4,z4])

def randomTF(N):
	for i in range(N):
		# first axis for B
		z=numpy.random.random()*2-1
		phi=numpy.random.random()*2*math.pi
		s=math.sqrt(1-z*z)
		x=math.cos(phi)*s
		y=math.sin(phi)*s
		# second random axis
		z2=numpy.random.random()*2-1
		phi2=numpy.random.random()*2*math.pi
		s2=math.sqrt(1-z2*z2)
		x2=math.cos(phi2)*s2
		y2=math.sin(phi2)*s2
		# third axis perpendicular to 1st for Beam
		x3=(y*z2-z*y2)
		y3=(z*x2-x*z2)
		z3=(x*y2-y*x2)
		r=math.sqrt(x3**2+y3**2+z3**2)
		yield ([x,y,z],[x3/r,y3/r,z3/r])

def NullIter(X):
	# for single crystals
	yield X

def ListIter(X):
	# listed orientations such as all inequivalent [111] axes
	for y in X:
		yield y

def GenericOrient(which,N,M0,B0,R0):
	# generic, selectable iterator
	# which= 0:crystal 1:uniformLF 2:randomLF 3:uniformTF 4:randomTF
	# N is number of orientations
	# M0, B0 are suggested axes for field and initial pol
	# returns actual axes for field and initial pol
	if(which==1):
		for (x,y,z) in uniformLF(N):
			yield ([x,y,z],[x,y,z])
	elif(which==2):
		for (x,y,z) in randomLF(N):
			yield ([x,y,z],[x,y,z])
	elif(which==3):
		for (xyz,x2y2z2) in uniformTF(N):
			yield (xyz,x2y2z2)
	elif(which==4):
		for (xyz,x2y2z2) in randomTF(N):
			yield (xyz,x2y2z2)
	else: # (which==0):
		yield M0,B0

def buildDynamicHam(Hams):
	# assemble Hamiltonians (list/tuple = Hams) into Big Ham,for rho, consider taking lower triangle only and unrolling by rows
	# A[0,0] A[0,1] A[0,2]
	# A[1,0] A[1,1] A[1,2]
	# A[2,0] A[2,1] A[2,2]
	# ->
	# [ A[0,0] A[1,0] A[1,1] A[2,0] A[2,1] A[2,2] ]
	# missing above-diagonals: Ham[0,2]=Ham[2,0]* and rho[0,2]=rho[2,0]*
	# need to include in X? How to represent conjugate?
	# simple H=[H00 H01]
	#               [H10 H11]
	# full square: BigHam=i*
	#	ij=	00	01		10		11
	#kl=	00	.	-H01		H10		.
	#	01	-H10	H00-H11	.		H10
	#	10	H01	.		H11-H00	-H01
	#	11	.	H01		-H10		.
	# since H is Hermitian, i*BigHam=Hermitian too. So no need to duplicate?
	# H squashed by triangles on each axis - not so obvious?
	# minimal is:
	#	ij=	00	10		11
	#kl=	00	0
	#	10	H01	H11-H00
	#	11	0	-H10		0
	#
	# (3,3) mat:
	# H=[	H00	H01	H02
	#	H10	H11	H12
	#	H20	H21	H22 ]
	#
	# BigHam=i*
	#	ij=	00	01		02		10		11	12		20		21		22
	#kl=	00	.	-H01		-H02		H10		.	.		H20		.		.
	#	01	-H10	H00-H11	-H12		.		H10	.		.		H20		.
	#	02	-H20	-H21		H00-H22	.		.	H10		.		.		H20
	#	10	H01	.		.		H11-H00	-H01	-H02		H21		.		.
	#	11	.	H01		.		-H10		.	-H12		.		H21		.
	#	12	.	.		H01		-H20		-H21	H11-H22	.		.		H21
	#	20	H02	.		.		H12		.	.		H22-H00	-H01		-H02
	#	21	.	H02		.		.		H12	.		-H10		H22-H11	-H12
	#	22	.	.		H02		.		.	H12		-H20		-H21		.
	#
	# squashed by triangle:
	#	ij=	00	10		11	20		21		22
	#kl=	00	.	H10		,	H20		,		,
	#	10	H01	H11-H00	-H01	H21		.		.
	#	11	.	-H10		.	.		H21		.
	#	20	H02	H12		.	H22-H00	-H01		-H02
	#	21	.	.		H12	-H10		H22-H11	-H12
	#	22	.	.		.	-H20		-H21		.
	#
	nhams=len(Hams)
	nqs=Hams[0].shape[0] # assume square, and all Hams are same size?
	SmallHamSize=nqs*nqs # nqs*(nqs+1)/2
	BigHamSize=nhams*SmallHamSize
	BigHam=numpy.zeros([BigHamSize,BigHamSize],dtype=complex)
	# equivalent Fortran with full array size
	# for i in range(nqs)
	#  for j in range(nqs)
	#   for kl in range(nqs)
	#    BigHam[i,j,kl,j] += Ham[i,kl]
	#    BigHam[i,j,i,kl] += Ham[kl,j]
	h=0
	for HamI in Hams:
		for i in range(nqs):
			for j in range(nqs):
				for kl in range(nqs):
					BigHam[h+i+j*nqs,h+kl+j*nqs] += HamI[i,kl]*1j
					BigHam[h+i+j*nqs,h+i+kl*nqs] -= HamI[kl,j]*1j
		h=h+SmallHamSize
	return (BigHam,SmallHamSize,nhams)
	
def buildDynamicHamTri(Hams):
	# assemble Hamiltonians (list/tuple = Hams) into Big Ham,for rho, consider taking lower triangle only and unrolling by rows
	# A[0,0] A[0,1] A[0,2]
	# A[1,0] A[1,1] A[1,2]
	# A[2,0] A[2,1] A[2,2]
	# ->
	# [ A[0,0] A[1,0] A[1,1] A[2,0] A[2,1] A[2,2] ]
	# missing above-diagonals: Ham[0,2]=Ham[2,0]* and rho[0,2]=rho[2,0]*
	# need to include in X? How to represent conjugate?
	# simple H=[H00 H01]
	#               [H10 H11]
	# full square: BigHam=i*
	#	ij=	00	01		10		11
	#kl=	00	.	-H01		H10		.
	#	01	-H10	H00-H11	.		H10
	#	10	H01	.		H11-H00	-H01
	#	11	.	H01		-H10		.
	# since H is Hermitian, i*BigHam=Hermitian too. So no need to duplicate?
	# H squashed by triangles on each axis - not so obvious?
	# minimal is:
	#	ij=	00	10		11
	#kl=	00	0
	#	10	H01	H11-H00
	#	11	0	-H10		0
	#
	# (3,3) mat:
	# H=[	H00	H01	H02
	#	H10	H11	H12
	#	H20	H21	H22 ]
	#
	# BigHam=i*
	#	ij=	00	01		02		10		11	12		20		21		22
	#kl=	00	.	-H01		-H02		H10		.	.		H20		.		.
	#	01	-H10	H00-H11	-H12		.		H10	.		.		H20		.
	#	02	-H20	-H21		H00-H22	.		.	H10		.		.		H20
	#	10	H01	.		.		H11-H00	-H01	-H02		H21		.		.
	#	11	.	H01		.		-H10		.	-H12		.		H21		.
	#	12	.	.		H01		-H20		-H21	H11-H22	.		.		H21
	#	20	H02	.		.		H12		.	.		H22-H00	-H01		-H02
	#	21	.	H02		.		.		H12	.		-H10		H22-H11	-H12
	#	22	.	.		H02		.		.	H12		-H20		-H21		.
	#
	# squashed by triangle:
	#	ij=	00	10		11	20		21		22
	#kl=	00	.	H10		,	H20		,		,
	#	10	H01	H11-H00	-H01	H21		.		.
	#	11	.	-H10		.	.		H21		.
	#	20	H02	H12		.	H22-H00	-H01		-H02
	#	21	.	.		H12	-H10		H22-H11	-H12
	#	22	.	.		.	-H20		-H21		.
	#
	nhams=len(Hams)
	nqs=Hams[0].shape[0] # assume square, and all Hams are same size?
	SmallHamSize=nqs*(nqs+1)/2 # triangle
	# relate (i,j) of Ham to ij in BigHam:
	# drop all rows/columns
	BigHamSize=nhams*SmallHamSize
	BigHam=numpy.zeros([BigHamSize,BigHamSize],dtype=complex)
	# equivalent Fortran with full array size
	# for i in range(nqs)
	#  for j in range(nqs)
	#   for kl in range(nqs)
	#    BigHam[i,j,kl,j] += Ham[i,kl]
	#    BigHam[i,j,i,kl] += Ham[kl,j]
	h=0
	for HamI in Hams:
		for i in range(nqs):
			for j in range(nqs):
				for kl in range(nqs):
					BigHam[h+i+j*nqs,h+kl+j*nqs] += HamI[i,kl]*1j
					BigHam[h+i+j*nqs,h+i+kl*nqs] -= HamI[kl,j]*1j
		h=h+SmallHamSize
	return (BigHam,SmallHamSize,nhams)

def addSpinRelaxation(ns,BigHam,state,spin,nu,pol,polmag):
	# spin relaxation in specified state, Relaxing to polarisation "pol"
	# spin[] is size (alpha,i,j) representing spin to flip. How to find set of indices i2 which differ from i only by flipping spin s? May have some of spin[alpha,i,i2] non zero..? But operator S only raises or lowers by 1. Should spin flip do the same?
	# pol is 3-vector, or alternatively (2I+1)*(2I+1) matrix for eqm rho for that spin alone (so works for I>1/2 or 0<eqm<1) or 0 for identity rho (unpolarised)
	# scale factor polmag, For 3-vector 1 means fully polarised (I=1/2) or linearly spaced 1 to 0 (I>1/2)
	# state multiplied by SmallHamSize? No.
	# BigHam[I,I]-= lambda (all I)
	# BigHam[I,J] += lambda/2 *beta(i,j) where I=(some others)+i, J=(same others)+j
	# code from CreateSpinMat to scan over all index pairs where only named spin varies?
	nqsq=BigHam.shape[0]/ns
	nq=spin.shape[1] # or numpy.prod(spins)
	if(nq*nq != nqsq):
		raise Exception("Spin matrix and BigHam sizes don't match")
	subset=BigHam[state*nqsq:(state+1)*nqsq,state*nqsq:(state+1)*nqsq].reshape((nq,nq,nq,nq))
	# spin finder
	# scan I until element spin[:,0,I] is non-zero
	# gives modu
	modu=-1
	for I in range(1,nq):
		if(not(numpy.all(spin[:,0,I]==0))):
			modu=I
			break
	if(modu==-1):
		raise Exception("no elements coupling state 0!")
	# scan J in steps of I until element spin[:,J-I,J] is zero
	# gives modv
	modv=-1
	for J in range(2*modu,nq,modu):
		if(numpy.all(spin[:,J-modu,J]==0)):
			modv=J
			break
	if(modv==-1):
		modv=nq # this was the last listed spin
	modw=modv/modu
	if((nq % modv) != 0):
		raise Exception("Spins don't seem to divide Matrix correctly")
	#print "spin ",spin,"stride=",modu,"range=",modv,"states=",modw

	pol2=numpy.array(pol)
	# pol is a numpy array e.g. 2D density matrix
	g=pol2.shape
	if(len(g)==2 and g[0]==modw and g[1]==modw):
		# pol is the density matrix to scale
		beta=pol2*polmag+numpy.identity(modw,dtype=numpy.complex)*(1-polmag)/modw
	elif(len(g)==1 and g[0]==3):
		# pol is an easy axis, find ground state? or linearly scaled populations
		polmag=polmag/math.sqrt(numpy.sum(pol2**2)) # normalise pol, then scale by polmag
		subspin=spin[:,0:modv:modu,0:modv:modu]
		beta=numpy.identity(modw,dtype=numpy.complex)/modw + (subspin[0,:,:]*pol2[0]+subspin[1,:,:]*pol2[1]+subspin[2,:,:]*pol2[2])*polmag/modw
		#print "calc beta=",beta
	else:
		# unpolarised
		beta=numpy.identity(modw,dtype=numpy.complex)/modw

	for I in range(nq):
		for J in range(nq):
			for q in range(modw):
				# I,K coupled by spin op. Need corresponding version of L from J
				K=(I%modu)+q*modu+(I//modv)*modv
				L=(J%modu)+q*modu+(J//modv)*modv
				ttt=nu*beta[(J//modu)%modw,(I//modu)%modw] # transpose or conjg needed?
				subset[I,J,K,L] += ttt
			subset[I,J,I,J] -= nu
	
	
def addConversion(ns,BigHam,state1,state2,nu,rspins,pols):
	# conversion between states preserving coherence except for spins listed, with relaxation to pol listed
	# initially just coherent conversion...
	# BigHam[H1,S1,S2,H2,S1,S2] += ionrate(H1 to H2) (H1 != H2)
	# BigHam[H1,S1,S2,H1,S1,S2] -= sum_H2(ionrate(H1 to H2)) (H1 != H2)
	nqsq=BigHam.shape[0]/ns
	for s1 in range(nqsq):
		BigHam[state2*nqsq+s1,state1*nqsq+s1] += nu # [to, from]?
		BigHam[state1*nqsq+s1,state1*nqsq+s1] -= nu

def createDynamicDensMat(rho,pops):
	# assemble density matrices, weighted
	# normalise pops
	# consider creating triangle
	tpop=numpy.sum(pops)
	npops=numpy.array(pops)/tpop
	# rho is a "typical" density matrix to be duplicated
	ns=npops.shape[0]
	n=rho.shape[0]
	bigrho=numpy.zeros(n*n*ns,dtype=complex)
	for s in range(ns):
		for i in range(n):
			for j in range(n):
				bigrho[s*n*n+i+j*n]=rho[i,j]*npops[s]
	return bigrho

def createDynamicSpinOp(spinop,ns):
	# assemble operators, assuming no state selection
	# consider creating opposite triangle to rho and Ham, and conjugate
	n=spinop.shape[0]
	bigop=numpy.zeros(n*n*ns,dtype=complex)
	for s in range(ns):
		for i in range(n):
			for j in range(n):
				bigop[s*n*n+i+j*n]=spinop[i,j]
	return bigop

def solveDynamicResult(BigHam,BigStart,BigDetect,timeend=float("Inf"),tzeros=None):
	# return tuple (omega,ccos,csin,lamb) -> (Eval=(lamb+j*omega),Eamp=(ccos+j*csin))
	# with timeend, return (Eval,Eamp,BigRho(t=timeend))
	nn=BigHam.shape[0] # not folded, no multiplicity to worry about
	eval,evec=numpy.linalg.eig(BigHam)
	#print "eval"
	#print eval
	#print "evec"
	#print evec
	# rho[i,j](t)=Evecmat[i,j:n]*exp(Eval(n)*t)*InvEvecmat[n:i'',j'']*rho0[i'',j'']
	# S(t) = Trace(rho[i,j](t) * S[i,j]) = sum[i,j](rho(i,j)(t)*S(j,i))
	# but S is Hermitian so S(j,i)=S*(i,j)
	# solve InvEvecmat*rho0 (i.e. convert to basis)
	#print "evec: ",evec.shape
	#print "BigStart: ",BigStart.shape
	Evrho=numpy.linalg.solve(evec,BigStart) # evec*Evrho = BigStart
	# tests>
	#Evrho3=numpy.dot(evec,Evrho) # undo?
	#eee=numpy.amax(numpy.abs(Evrho3-BigStart))
	#if(eee>1.E-6):
	#	print "array re-creation error"
	#	print BigStart
	#	print Evrho3
	#else:
	#	print "probably OK, difference = ",eee
	# <tests
	#print "BigStart in new basis: ",Evrho
	# rho[i,j](t)=Evecmat[i,j:n]*exp(Eval(n)*t)*Evrho[n]
	# S(t)=Trace(rho[i,j](t)*S[i,j]) = Sum(rho[i,j] S[j,i])
	# = Trace(Evecmat[i,j,n]*S*[i,j]) * Eamp(n) * exp(Eval(n)*t)
	Eamp=Evrho* numpy.dot(BigDetect.conjugate(),evec) 
	#Eamp=numpy.zeros(nn,dtype=complex) # complex
	#for ij in range(nn):
	#	for kl in range(nn):
	#		Eamp[ij]+=Evrho[ij]*(BigDetect[kl]).conjugate()*evec[kl,ij]
	if(math.isinf(timeend)):
		return (eval,Eamp)
	else:
		Evrho2=Evrho*numpy.exp(eval*timeend) #phases, relaxation
		if(tzeros is not None):
			ConvolveTimeResDynListed(eval,Evrho2,tzeros) # use Dynamic form to deal with complex coeffs
		Evrho3=numpy.dot(evec,Evrho2) # std basis states..?
		return (eval,Eamp,Evrho3)

def EvaluateDynamicIntoBins(eval,ampl,lamMu,times):
	# fill time bins, lamb=relaxation value lamMu=muon lifetime weighting
	# y(t) = Re(sum(n) ampl(n)*exp(eval(n)*t))
	# bin from t1 to t2:
	# AveAmpl(t1,t2)=sum(n) integral[t1:t2](Re(ampl(n)*exp(eval(n)*t)*exp(-lamMu*t))
	# NEARLY same as EvaluateIntoBins except different lambda for each coefficient = Re(eval(n))-lamMu and omega=Im(eval(n))
	# but don't normalise out all of lam!
	# times = array length n+1 of time bin boundaries
	# returns array length n of values, averaged across bins
	# lam = 1/lifetime to weight by, if bins are large (or only one big one for integral counting - times [0,inf] is OK)
	# = integral (cos(w*t))exp(-lam*t)]t1:t2 / integral exp(-lam*t)]t1:t2
	# = (-lam* sin(w*t) exp(-lam*t) + w* cos(w*t) exp*(-lam*t) )/(lam^2+w^2)
	# improved version using time arrays internally (best for time series from fairly simple models)
	ybins=numpy.zeros(len(times)-1)
#	for i in range(len(times)-1):
	lam=eval-lamMu
	et0=numpy.exp(-lamMu*times)
	sf=-lamMu/(et0[1:]-et0[:-1])
	for j in range(len(lam)):
		cs0=numpy.exp(lam[j]*times)/lam[j]
		cs1=cs0[:-1]
		cs2=numpy.nan_to_num(cs0[1:]) # avoid domain problems with [..,inf] integral counting, would be multiplied by 0 anyway
		ybins=ybins+numpy.real(ampl[j]*(cs2-cs1)*sf)
#		ybins=ybins+numpy.real((ampl[j]*(omega[j]*(s2-s1 )-lam[j]*(c2-c1 )) +csin[j]*(omega[j]*(-c2+c1 )-lam[j]*(s2-s1 )) )/sf /(lam[j]*lam[j]+omega[j]*omega[j]))

	return ybins


def solveRFDensityMat(Ham0,Ham1,Ham1i,omegaRF,start,detect,detecti,RRFharmonic,tsliceIn=None,timeend=float("Inf"),tzeros=None):
	# RF mode (high freq, low B1 limit)
	# create Ham0 and Ham1 (and optionally Ham1i for circular pol)
	# also use Ham1i for RF phase not zero
	# "integrate" up Ham for cycle
	# result is a series of sin/cos terms as for plain solution
	# if rrfharmonic !=0, instantaneous operator is detect*cos(RRFharmonic*omegaRF)+detecti*sin(...)
	# calc time slice for integration:
	# time slice for piecewise approx of RF sufficiently small that harmonics can't cause transitions between any levels
	# and that sampling of detect doesn't cause aliasing
	# therefore need estimate of difference between extreme eigenvalues (roughly 2*magnitude of largest element of H)
	#
	# discrete integral of sin(wt) in N steps * sin(wt) (smooth) has value proportional to N/pi*sin(pi/N) instead of 1
	# enhance Ham1 and Ham1i by factor pi/N/sin(pi/N) to compensate
	# still sometimes leaves a 1/N^4 error (small?) - investigate?
	# nutation?
	#
	tsliceN=7
	# print omegaRF,numpy.amax(numpy.absolute(Ham0)),numpy.amax(numpy.absolute(Ham1))
	Hmax=max(omegaRF,numpy.amax(numpy.absolute(Ham0)),numpy.amax(numpy.absolute(Ham1)))
	if (Ham1i is not None):
		Hmax=max(Hmax,numpy.amax(numpy.absolute(Ham1i)))
	tslices=(tsliceN*Hmax/omegaRF) # detail of sinewave for correct averaging of detector operator
	tslices=int(tsliceIn or tslices) # override on request, force to integer
	boost=math.pi/tslices/math.sin(math.pi/tslices)
	# boost=1.0 # to switch off for testing
	# what to do if omega is small but E is large - tslices gets huge! Reduce tslices and risk harmonics?
	#  was: eslices=int(math.log(1.E9/tslices)/math.log(2.0))+1
	# slice of length 1/omega/tslices has approximately 1/(tsliceN) of phase rotation regardless of B, omega, etc
	# so eslices ought to be constant
	eslices=30 # ~1E9
	eps=1./tslices/2**eslices/omegaRF*2.j*math.pi # time slice for one minimal unit

	#print "tslices=",tslices," eslices=2**",eslices," eps=",eps

	detectacc=numpy.array(detect/2.0)
	#print "da00:",detectacc[0,0]
	Uacc=numpy.eye(Ham0.shape[0],dtype=complex) # complete value here
	Uiacc=numpy.eye(Ham0.shape[0],dtype=complex) # complete value here
	Ubit=numpy.zeros_like(Ham0) # starts as U-1, add one at Tslice stage
	work=numpy.zeros_like(Ham0)
	work2=numpy.zeros_like(Ham0)
	work3=numpy.zeros_like(Ham0)
	for s in range(tslices):
		if(Ham1i is not None):
			Ubit[:,:]=(Ham0+Ham1*math.cos(2.0*math.pi*(s+0.5)/tslices)*boost+Ham1i*math.sin(2.0*math.pi*(s+0.5)/tslices)*boost)*eps
		else:
			Ubit[:,:]=(Ham0+Ham1*math.cos(2.0*math.pi*(s+0.5)/tslices)*boost)*eps
		for ee in range(eslices):
			Ubit[:,:]=Ubit*2+numpy.dot(Ubit,Ubit,out=work)
		# Uacc=Uacc*Ubit = Uacc*(1+U'bit)
		Uacc[:,:]=Uacc+numpy.dot(Uacc,Ubit,out=work)
		Uiacc[:,:]=Uiacc+numpy.dot(numpy.conj(numpy.transpose(Ubit)),Uiacc,out=work)
		#if(not numpy.allclose( numpy.dot(Uacc,Uiacc),numpy.eye(Ham0.shape[0]) )):
		#	print "U inverse failure iter=",s
		# detect(t) = U^-1(t).detect(0).U(t)
		# = (V+1)$ . detect0 . (V+1)
		# = detect + V$.detect + detect.V + V$.detect.V
		# simpler to stick to U...
		# below actually calcs Detect* and assumes U is Hermitian (poor assumption?)
		if(RRFharmonic):
			if(detecti is not None):
				detectacc[:,:]=detectacc+(numpy.dot(numpy.dot(Uiacc,detect,out=work3),Uacc,out=work))*math.cos(2.0*math.pi*(s+1)*RRFharmonic/tslices)+\
								   (numpy.dot(numpy.dot(Uiacc,detecti,out=work3),Uacc,out=work2))*math.sin(2.0*math.pi*(s+1)*RRFharmonic/tslices)
			else:
				detectacc[:,:]=detectacc+(numpy.dot(numpy.dot(Uiacc,detect,out=work3),Uacc,out=work))*math.cos(2.0*math.pi*(s+1)*RRFharmonic/tslices)
		else:
			detectacc[:,:]=detectacc+(numpy.dot(numpy.dot(Uiacc,detect,out=work3),Uacc,out=work))
		#print "da00:",detectacc[0,0]
	detectacc[:,:]=(detectacc-0.5*work) / tslices # from dot(Uacc,detect)
	#print "da00:",detectacc[0,0]
	# then calc exponent
	# exp(t,H) ~ lim(n->inf) (1+(t/n)*H)^n
	# for accuracy define function for (1+H)^2 -> 1+H2, ie. H2=2*H+H*H
	
	# integrate up -> U
	# accumulate new detector op over cycle

	# solving U: ought to have a series of eigenvalues with magnitude 1 and phase E*t_RF
	# no longer needed Uacc=Uacc+numpy.eye(Uacc.shape[0])
	(egvalp,evec)=numpy.linalg.eig(Uacc)
	eval=numpy.angle(egvalp)*omegaRF/2.0/math.pi
	ampls=numpy.absolute(egvalp)
	minamp=numpy.amin(ampls)
	maxamp=numpy.amax(ampls)
	if(minamp<0.999 or maxamp>1.001):
		print "eigenvalues are not normalised? ts=",tslices," e=2**",eslices
		print ampls
	# print "t=",tslices," e=2**",eslices," Egvals=",eval

	# debug tests: eigenvectors should be orthogonal!
	#msg="All good, max="
	#maxone=0
	#for v1 in range(evec.shape[1]-1):
	#	for v2 in range(v1+1,evec.shape[1]):
	#		dot=numpy.dot(evec[:,v1],numpy.conj(evec[:,v2]))
	#		if(abs(dot)>0.01):
	#			print "eigenvecs E1=",eval[v1]," and E2=E1+",eval[v2]-eval[v1]," have dot prod=",dot
	#			msg="remaining ones max="
	#		elif abs(dot)>abs(maxone):
	#			maxone=dot
	#print msg,maxone
	
	# how to "factor out" one of the near degenerate egvals, or to re-normalise the two egvecs?
	# given V1, V2
	# need V2' = (aV1+V2) where V1.V2'=0
	# or does it matter, if sufficient care taken in following code? Provided no linear combinations appear, i.e. evec has an inverse.
	
	#print eval
	#print evec
	#print detect
	#print detectacc

	# evec=numpy.transpose(evec) # in case of row/col mistakes?

	# from here, should be the same as plain SolveDensityMat except detect is delayed by 1/2 a cycle so must adjust phases to compensate
	# halftimes=numpy.exp(1j*math.pi/omegaRF*eval) # relevant phase factors
	# or halftimes = sqrt(egvalp)
	# mhalftimes=numpy.conj(halftimes)
	
	evecstar=numpy.conj(evec)
	evecinv=numpy.linalg.inv(evec)
	# transform rho into basis set of evec
	#* rho2=numpy.einsum('ji,jk,kl',evecstar,start,evec)
	rho1a=numpy.dot(start,evec) # start[j,k]*evec[k.l]
	rho2=numpy.dot(evecinv,rho1a)
	# was rho2=numpy.tensordot(evecstar,rho1a,(0,0) )
	#print "rho2"
	#print rho2
	# likewise detector operator (CHECK conjg, etc)
	#* detect2=numpy.einsum('ji,jk,kl',evecstar,detect,evec)
	detect1a=numpy.dot(detectacc,evec,rho1a) # reuse rho1a, it's the right shape
	detect2=numpy.dot(evecinv,detect1a)
	# was detect2=numpy.tensordot(evecstar,detect1a,(0,0) )
	#print "detect2"
	#print detect2
	#print numpy.dot(numpy.diag(rho2),numpy.diag(detect2))
	# detectedfreqs Eval(i)-Eval(j) ampl rhoij*Sji
	# detected const dot(diag(rho),diag(s))
	n=Ham0.shape[0]
	nc=n*(n-1)/2+1
	omega=numpy.empty(nc,dtype=float)
	ccos=numpy.empty(nc,dtype=float)
	csin=numpy.empty(nc,dtype=float)
	omega[0]=0.0
	ccos[0]=numpy.real(numpy.dot(numpy.diag(rho2),numpy.diag(detect2)))
	csin[0]=0.0
	# note Evec and detect are no longer necessarily Hermitian! Can't short cut with lower triangles any more
	# but OK if all elements used?
	(ii,jj)=numpy.tril_indices(n,-1)
	omega[1:]=numpy.angle(egvalp[ii]/egvalp[jj])*omegaRF/2.0/math.pi # cope with wrap round
	rd1=rho2[ii,jj]*detect2[jj,ii]*numpy.sqrt(egvalp[jj]/egvalp[ii]) # *halftimes[jj]*mhalftimes[ii]
	rd2=rho2[jj,ii]*detect2[ii,jj]*numpy.sqrt(egvalp[ii]/egvalp[jj]) # *halftimes[ii]*mhalftimes[jj]
	ccos[1:]=numpy.real(rd1)+numpy.real(rd2)
	csin[1:]=numpy.imag(rd1)-numpy.imag(rd2) # check signs!
	#print "Eigenvectors"
	#print evec
	#print "Eigenvalues"
	#print eval
	if(math.isinf(timeend)):
		return (omega,ccos,csin)
	else:
		# calc rho at t=timeend
		# exact only if timeend is a whole number of RF cycles (or maybe store fractional Ham from loop stage?)
		# rho2 is in eigenvectors
		# rho2(t) = rho2 * phase factors
		# rho(std basis) by multiplying by Evec^-1
		# print "power for adjustment: ",(timeend*omegaRF/2.0/math.pi)
		egvalp2=egvalp**(timeend*omegaRF/2.0/math.pi)
		rho2=rho2*egvalp2/egvalp2[:,numpy.newaxis]
		if(tzeros is not None):
			# variable time interval (pulse width). Rather approximate here, working with observed freqs <omegaRF/2 only...
			freqs=numpy.log(egvalp/egvalp[:,numpy.newaxis])*omegaRF/2.0/math.pi # i*phase angle over cycle, conv to freq
			ConvolveTimeResDynListed(freqs,rho2,tzeros) # use Dynamic form to deal with complex coeffs
		# care: Evec is not unitary because of RF variations. Need inverse! Or solve eqns directly.
		# already done evecinv=numpy.linalg.inv(evec)
		evecinvstar=numpy.conj(evecinv)
		rho3a=numpy.dot(rho2,evecinv)
		# was rho3=numpy.tensordot(evecinvstar,rho3a,(0,0))
		rho3=numpy.dot(evec,rho3a)
		# print "check1",numpy.dot(evec,evecinv)
		# print "check2",numpy.tensordot(evecstar,evecinvstar,(0,1))
		return (omega,ccos,csin,rho3)

# evaluate sliced
def EvaluateSliced(Hams,start,detect,slices,times,lam=1.0/2.19703,tzeros=None):
	# Hams is array of n Hamiltonians, one per slice
	# slices is array of n-1 end times (absolute, must be increasing). First slice starts at t=0, last slice assumed to go to t=Inf if needed.
	# times is array of m+1 data points wanted (need not coincide with time boundaries in any way! Last can be Inf)
	# return m y-values (averaged within bins of times, weighted by mulife)
	# tzeros = pulse shape rounding factors, apply to 1st slice and to density matrix on 1->2 changeover
	ns=len(Hams)
	yvals=numpy.zeros(len(times)-1,dtype=float)
	weightacc=numpy.zeros(len(times)-1,dtype=float) # debug only
	nextrho=start
	slices2=numpy.array(slices)
	slices2=numpy.insert(slices2,0,0.0)
	slices2=numpy.append(slices2,float("Inf"))
	for i in range(ns):
		# slice with end point given
		# calc start and end points in times which correspond to this slice
		it1=bisect.bisect_right(times,slices2[i])
		it2=bisect.bisect_left(times,slices2[i+1])
		# points from it1 to it2, plus slice ends
		times2=times[it1:it2]
		if(it2<len(times) and slices2[i+1]>times2[-1]+1.E-12):
			times2=numpy.append(times2,slices2[i+1]) # this slice does not go to tmax, so add boundary
			it2=it2+1
		if(it1>0 and slices2[i]<times2[0]-1.E-12):
			times2=numpy.insert(times2,0,slices2[i]) # this slice starts after tmin, add boundary
			it1=it1-1
		weights=numpy.ones(len(times2)-1)
		# broken end points
		# assume that each time bin will be completely filled somehow (may take >2 slices)
		# len=0: pre- (or post-) processing time slice only, no cosines to evaluate
		# len=1, slice times2[0] to times2[1]  within bin times[it1] to times[it1+1]
		# len>1, consider end points only, as above
		# 1st pt of this curve contributes times2[0] to times2[1] out of times[it1-1] to times[it1]
		# weighted, total weight=integral(exp(-t/mulife))_t1^t2 = exp(-t2/mulife)-exp(t1/mulife)
		# sub weights similarly
		if(len(weights)==0):
			# Error? No interval in times() overlaps this slice. Evolve rho anyway, though why bother if now past last time bin?
			pass
		elif(len(weights)==1):
			# one time bin left, may overlap either end?
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
		else:
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
				weights[-1]=(math.exp(-times2[-2]*lam)-math.exp(-times2[-1]*lam))/(math.exp(-times[it2-2]*lam)-math.exp(-times[it2-1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
				weights[-1]=(times2[-2]-times2[-1])/(times[it2-2]-times[it2-1])
				
		te=slices2[i+1]-slices2[i]
		if(math.isinf(te)):
			(omega,ccos,csin) = solveDensityMat(Hams[i],nextrho,detect,timeend=te)
		else:
			if(i==0 and tzeros is not None):
				(omega,ccos,csin,nextrho) = solveDensityMat(Hams[i],nextrho,detect,timeend=te,tzeros=tzeros)
			else:
				(omega,ccos,csin,nextrho) = solveDensityMat(Hams[i],nextrho,detect,timeend=te)
		times3=times2-slices2[i] # relative
		if(len(weights)>0):
			if(i==0 and tzeros is not None):
				ConvolveTimeResolutionListed(omega,ccos,csin,tzeros)
			tmpY=evaluateIntoBins(omega,ccos,csin,lam,times3)
			yvals[it1:it2-1]=yvals[it1:it2-1]+weights*tmpY
			weightacc[it1:it2-1]=weightacc[it1:it2-1]+weights
		#else:
			#print "nothing to do this time"

	return yvals

# evaluate sliced
# solveRFDensityMat(Ham0,Ham1,Ham1i,omegaRF,start,detect,detecti,RRFharmonic,tsliceN=6,timeend=float("Inf"))
def EvaluateSlicedRF(Hams0,Hams1,Hams1i,omegasRF,start,detect,detecti,RRFharmonic,slices,times,lam=1.0/2.19703,tzeros=None):
	# Hams0 is array of n Hamiltonians, one per slice, likewise Hams1 and Hams1i
	# omegasRF is array of freqs, need not be equal. Any zeros mean model that slice in plain mode, Single float value is expanded.
	# RRFharmonic!=0 only makes sense with omegaRF kept constant
	# slice times probably ought to be whole cycles of RF, or alternatively phase of (Ham1+Ham1i) varies to match, to prevent phase jump
	# slices is array of n-1 end times (absolute, must be increasing). First slice starts at t=0, last slice assumed to go to t=Inf if needed.
	# times is array of m+1 data points wanted (need not coincide with time boundaries in any way! Last can be Inf)
	# return m y-values (averaged within bins of times, weighted by mulife)
	ns=len(Hams0)
	if(Hams1i is None):
		Hams1i=[None]*ns
	# omegasRF migt be single number, extend to array
	try:
		if(len(omegasRF) != ns):
			raise Exception("Wrong number of frequencies supplied")
	except:
		omegasRF=[omegasRF]*ns
	yvals=numpy.zeros(len(times)-1,dtype=float)
	weightacc=numpy.zeros(len(times)-1,dtype=float) # debug only
	nextrho=start
	slices2=numpy.array(slices)
	slices2=numpy.insert(slices2,0,0.0)
	slices2=numpy.append(slices2,float("Inf"))
	for i in range(ns):
		# slice with end point given
		# calc start and end points in times which correspond to this slice
		it1=bisect.bisect_right(times,slices2[i])
		it2=bisect.bisect_left(times,slices2[i+1])
		# points from it1 to it2, plus slice ends
		times2=times[it1:it2]
		#print "slicer: slices2=",slices2[i],slices2[i+1]," it1,2=",it1,it2," times2=",times2
		if(it2<len(times) and it2>0 and slices2[i+1]>times[it2-1]+1.E-12):
			times2=numpy.append(times2,slices2[i+1]) # this slice does not go to tmax, so add boundary (if non-negligible overlap)
			it2=it2+1
		if(it1>0 and it1<len(times) and slices2[i]<times[it1]-1.E-12):
			times2=numpy.insert(times2,0,slices2[i]) # this slice starts after tmin, add boundary (if non-negligible overlap)
			it1=it1-1
		if(len(times2)>0):
			weights=numpy.ones(len(times2)-1)
		else:
			weights=[]
		# broken end points
		# assume that each time bin will be completely filled somehow (may take >2 slices)
		# len=0: pre- (or post-) processing time slice only, no cosines to evaluate
		# len=1, slice times2[0] to times2[1]  within bin times[it1] to times[it1+1]
		# len>1, consider end points only, as above
		# 1st pt of this curve contributes times2[0] to times2[1] out of times[it1-1] to times[it1]
		# weighted, total weight=integral(exp(-t/mulife))_t1^t2 = exp(-t2/mulife)-exp(t1/mulife)
		# sub weights similarly
		if(len(weights)==0):
			# Error? No interval in times() overlaps this slice. Evolve rho anyway, though why bother if now past last time bin?
			pass
		elif(len(weights)==1):
			# one time bin left, may overlap either end?
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
		else:
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
				weights[-1]=(math.exp(-times2[-2]*lam)-math.exp(-times2[-1]*lam))/(math.exp(-times[it2-2]*lam)-math.exp(-times[it2-1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
				weights[-1]=(times2[-2]-times2[-1])/(times[it2-2]-times[it2-1])
		#if(numpy.any(numpy.isnan(weights))):
			#print "weights has a NaN value for slice ",i
			#raise Exception("best to crash now")
		#if(numpy.any(numpy.isinf(weights))):
			#print "weights has an infinite value for slice ",i
			#raise Exception("best to crash now")
		#if(numpy.amax(weights)>1.0):
			#print "weights has value(s) exceeding one, i.e.",numpy.amax(weights)
		#if(numpy.amin(weights)<0.0):
			#print "weights has negative value(s), i.e.",numpy.amin(weights)
		te=slices2[i+1]-slices2[i]
		if(omegasRF[i]==0): # plain slice
			if(math.isinf(te)):
				(omega,ccos,csin) = solveDensityMat(Hams0[i],nextrho,detect,timeend=te)
			else:
				if(i==0 and tzeros is not None):
					(omega,ccos,csin,nextrho) = solveDensityMat(Hams0[i],nextrho,detect,timeend=te,tzeros=tzeros)
				else:
					(omega,ccos,csin,nextrho) = solveDensityMat(Hams0[i],nextrho,detect,timeend=te)
		else: # RF slice
			if(math.isinf(te)):
				(omega,ccos,csin) = solveRFDensityMat(Hams0[i],Hams1[i],Hams1i[i],omegasRF[i],nextrho,detect,detecti,RRFharmonic,timeend=te)
			else:
				if(i==0 and tzeros is not None):
					(omega,ccos,csin,nextrho) = solveRFDensityMat(Hams0[i],Hams1[i],Hams1i[i],omegasRF[i],nextrho,detect,detecti,RRFharmonic,timeend=te,tzeros=tzeros)
				else:
					(omega,ccos,csin,nextrho) = solveRFDensityMat(Hams0[i],Hams1[i],Hams1i[i],omegasRF[i],nextrho,detect,detecti,RRFharmonic,timeend=te)
		times3=times2-slices2[i] # relative
		if(len(weights)>0):
			if(i==0 and tzeros is not None):
				ConvolveTimeResolutionListed(omega,ccos,csin,tzeros)
			tmpY=evaluateIntoBins(omega,ccos,csin,lam,times3)
			#if(numpy.any(numpy.isinf(tmpY)) or numpy.any(numpy.isnan(tmpY))):
				#print "tmpY values from ",numpy.amin(tmpY)," to ",numpy.amax(tmpY)
				#print "tmpY has a bad value; times3 was ",times3
				#print "deltas: ",times3[1:]-times3[:-1]
			yvals[it1:it2-1]=yvals[it1:it2-1]+weights*tmpY
			weightacc[it1:it2-1]=weightacc[it1:it2-1]+weights
		#else:
			#print "nothing to do this time"
	#print "final weightacc from ",numpy.amin(weightacc)," to ",numpy.amax(weightacc)
	return yvals

def EvaluateSlicedDynamic(BigHams,BigStart,BigDetect,slices,times,lam=1.0/2.19703,tzeros=None):
	# Hams is array of n Hamiltonians, one per slice
	# slices is array of n-1 end times (absolute, must be increasing). First slice starts at t=0, last slice assumed to go to t=Inf if needed.
	# times is array of m+1 data points wanted (need not coincide with time boundaries in any way! Last can be Inf)
	# return m y-values (averaged within bins of times, weighted by mulife)
	ns=len(BigHams)
	print "Ham length=",ns," and slice times are ",slices
	yvals=numpy.zeros(len(times)-1,dtype=float)
	weightacc=numpy.zeros(len(times)-1,dtype=float) # debug only
	nextrho=BigStart
	slices2=numpy.array(slices)
	slices2=numpy.insert(slices2,0,0.0)
	slices2=numpy.append(slices2,float("Inf"))
	for i in range(ns):
		# slice with end point given
		# calc start and end points in times which correspond to this slice
		it1=bisect.bisect_right(times,slices2[i])
		it2=bisect.bisect_left(times,slices2[i+1])
		# points from it1 to it2, plus slice ends
		times2=times[it1:it2]
		print "slice ",i," times2=",times2
		if(it2<len(times) and (len(times2)==0 or slices2[i+1]>times2[-1]+1.E-12)):
			times2=numpy.append(times2,slices2[i+1]) # this slice does not go to tmax, so add boundary
			it2=it2+1
		if(it1>0 and (len(times2)==0 or slices2[i]<times2[0]-1.E-12)):
			times2=numpy.insert(times2,0,slices2[i]) # this slice starts after tmin, add boundary
			it1=it1-1
		weights=numpy.ones(len(times2)-1)
		# broken end points
		# assume that each time bin will be completely filled somehow (may take >2 slices)
		# len=0: pre- (or post-) processing time slice only, no cosines to evaluate
		# len=1, slice times2[0] to times2[1]  within bin times[it1] to times[it1+1]
		# len>1, consider end points only, as above
		# 1st pt of this curve contributes times2[0] to times2[1] out of times[it1-1] to times[it1]
		# weighted, total weight=integral(exp(-t/mulife))_t1^t2 = exp(-t2/mulife)-exp(t1/mulife)
		# sub weights similarly
		if(len(weights)==0):
			# Error? No interval in times() overlaps this slice. Evolve rho anyway, though why bother if now past last time bin?
			pass
		elif(len(weights)==1):
			# one time bin left, may overlap either end?
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
		else:
			if(lam>0):
				weights[0]=(math.exp(-times2[0]*lam)-math.exp(-times2[1]*lam))/(math.exp(-times[it1]*lam)-math.exp(-times[it1+1]*lam))
				weights[-1]=(math.exp(-times2[-2]*lam)-math.exp(-times2[-1]*lam))/(math.exp(-times[it2-2]*lam)-math.exp(-times[it2-1]*lam))
			else:
				weights[0]=(times2[0]-times2[1])/(times[it1]-times[it1+1])
				weights[-1]=(times2[-2]-times2[-1])/(times[it2-2]-times[it2-1])
				
		te=slices2[i+1]-slices2[i]
		if(math.isinf(te)):
			(omega,ccossin) = solveDynamicResult(BigHams[i],nextrho,BigDetect,timeend=te)
		else:
			if(i==0 and tzeros is not None):
				(omega,ccossin,nextrho) = solveDynamicResult(BigHams[i],nextrho,BigDetect,timeend=te,tzeros=tzeros)
			else:
				(omega,ccossin,nextrho) = solveDynamicResult(BigHams[i],nextrho,BigDetect,timeend=te)
		times3=times2-slices2[i] # relative
		if(len(weights)>0):
			if(i==0 and tzeros is not None):
				ConvolveTimeResDynListed(omega,ccossin,tzeros)
			tmpY=EvaluateDynamicIntoBins(omega,ccossin,lam,times3)
			yvals[it1:it2-1]=yvals[it1:it2-1]+weights*tmpY
			weightacc[it1:it2-1]=weightacc[it1:it2-1]+weights
		#else:
			#print "nothing to do this time"

	return yvals

