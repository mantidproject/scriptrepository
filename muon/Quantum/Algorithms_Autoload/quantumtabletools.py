## Quantum - a program for solving spin evolution of the muon
## Author: James Lord
## Version 1.04, August 2018
# table tools - can be run independently of Mantid except:
#  Point by point fitting
#  Crystal Fields
#  Crystal space group symmetry and coordinates in fractions of unit cell
import numpy
import math

import re
import numbers
#import decimal
#import sys
import traceback
numpy.set_printoptions(linewidth=100)

from quantumtools import *


def ParseStringToDict(st,di0={}):
	# convert String to dict
	# initial part filled dict provided
	# entries separated by semicolon or new line
	# key is type:string
	# value may be tuple, string, etc as necessary
	# "spins=" has comma separated strings
	# "A=" and others have tuple of float
	# "A(1)=A(2)=34" generates two entries for A(1) and A(2)
	# A(*)=44 looks up A() entries and overwrites those which match...No, should expand later.
	# comment lines start with "#"
	
	textlists=("measure","spins","detectspin","initialspin","tzero","recycle","brf","morespin","crystal","goniometer")
	freetext=("fitfunction","loop0par","loop1par","fit0par","fit1par","fit2par","fit3par","fit4par","fit5par","fit6par","fit7par","fit8par","fit9par") # do not parse at all (are on one line)
	di=dict(di0) # copy
	for sl in st.splitlines():
		if(len(sl.strip()) > 0 and sl.strip()[0] != "#"):
			if(sl.split("=")[0] in freetext):
				kvl=sl.split("=",1)
				di[kvl[0].strip()]=kvl[1]
			else:
				for ss in sl.split(";"):
					#print "processing ",ss
					kvl=ss.split("=")
					kvl=map(str.strip,kvl)
					if(len(kvl)<2):
						raise Exception("Bad line in model string")
					val=kvl[-1].split(",")
					kvl0s=kvl[0].split("(")[0]
					if(kvl0s not in textlists):
						#print "would like to float ",val," from ",sl
						val=map(float,val)
					#print "value is ",val
					for key in kvl[0:-1]:
						if(key.find(".*")>=0):
							for ek in di.keys():
								if( re.match(key,ek)):
									di[ek]=val
						else:
							di[key]=val
	return di

def calc_eqm(rates):
	# rates is N*N array (from, to)
	# diagonal elements ignored
	# sum_j(pop(j)*rate(j->i)) = pop(i) * sum_k(rate(i->k))
	# copy from FORTRAN, except replace eqn 0 with the population sum instead of eqn N
	rates2=numpy.array(rates).transpose()
	n=rates2.shape[0]
	for i in range(1,n):
		rates2[i,i]=-numpy.sum(rates2[:i,i])-numpy.sum(rates2[i+1:,i])
	rates2[0,:]=1.0
	vec=numpy.zeros([n],dtype=numpy.float)
	vec[0]=1.0
	#print rates2
	#print vec
	return numpy.linalg.solve(rates2,vec)

def EnumerateSites(p,d,keystr):
	# strip leading ( and trailing ) from keystr
	# split keystr on "," and look for @site
	# return indices of sites to set, and remaining string parts
	# p is number of pulse slices
	# d is number of sites
	kss=((keystr.lstrip("(")).rstrip(")")).split(",") #....
	ds=[]
	ps=[]
	#print kss
	while (len(kss)>0 and len(kss[0])>0 and (kss[0][0]=="@" or kss[0][0]=="|")):
		if(kss[0][0]=="@"):
			ds.append(int(kss[0][1:]))
			kss=kss[1:]
		elif(kss[0][0]=="|"):
			ps.append(int(kss[0][1:]))
			kss=kss[1:]
	#print "ds=",ds,"ps=",ps
	if(d==0):
		d2=1
		if(len(ds)==0):
			ds=(0,)
	else:
		d2=d
		if(len(ds)==0):
			ds=range(d)
	if(len(ps)==0):
		ps=range(p)
	#print "ds=",ds,"ps=",ps
	dps=[(pp*d2+dd) for pp in ps for dd in ds ]
	return dps,kss
	
def ParseAndIterateModel(pars):
	# iterator (for orientations)
	# return (spmat,Ham,rho0,scint) repeatedly until finished
	# first loop to get spins, gamma, etc
	# database, name (Case sensitive!), 2I+1, gamma/2pi (MHz/T)
	# pars with "*" index eg. A(*)=100 are expanded to all possible values of *
	# expand most general (*,*) first then (i,*) (*,i) then do specific ones (i,j)
	#print "pars inside loop is ",pars
	# one 100% isotope with I=0: no need for A. Also default for those with abundance > 90%
	# signs match sign of moment
	PeriodicTable={
	'Mu':(2,135.5342), # Mantid standard value from PhysicalConstants.h
	'e':(2,-28024.95164), # from NIST physical constants tables
	'H':(2,42.5764), # others from CRC handbook pp 9-85 to 9-87
	'D':(3,6.53573),
	'3He':(2,-32.4352),
	'6Li':(3,6.2660),
	'7Li':(4,16.5478),
	'Li':(4,16.5478),
	'Be':(4,-5.986),
	'10B':(7,4.5751),
	'11B':(4,13.6626),
	'13C':(2,10.7081),
	'14N':(3,3.0776),
	'N':(3,3.0776),
	'15N':(2,-4.3172),
	'17O':(6,-5.7741),
	'F':(2,40.0765),
	'21Ne':(4,-3.3630),
	'Na':(4,11.2686),
	'25Mg':(6,-2.6082),
	'Al':(6,11.1028),
	'29Si':(2,-8.4653),
	'P':(2,17.2510),
	'33S':(4,4,3.2716),
	'35Cl':(4,4.1764),
	'37Cl':(4,3.4764),
	'39K':(4,1.9893),
	'K':(4,1.9893),
	'40K':(9,-2.4737),
	'41K':(4,1.0919),
	'43Ca':(8,-2.8688),
	'Sc':(8,10.3588),
	'47Ti':(6,-2.4040),
	'49Ti':(8,-2.4047),
	'50V':(13,4.2504),
	'51V':(8,11.2130),
	'V':(8,11.2130),
	'53Cr':(4,-2.4114),
	'Mn':(6,10.5760),
	'57Fe':(2,1.3815),
	'Co':(8,10.077),
	'61Ni':(4,-3.8113),
	'63Cu':(4,11.2979),
	'65Cu':(4,12.1027),
	'67Zn':(6,2.6693),
	'69Ga':(4,10.2475),
	'71Ga':(4,13.0204),
	'73Ge':(10,-1.4897),
	'As':(4,7.3148),
	'77Se':(2,8.1566),
	'79Br':(4,10.7039),
	'81Br':(4,11.5381),
	'83Kr':(10,-1.6442),
	'85Rb':(6,4.1253),
	'87Rb':(4,13.9807),
	'87Sr':(10,-1.8524),
	'Y':(2,-2.0949),
	'91Zr':(6,-3.9747),
	'Nb':(10,10.4520),
	'95Mo':(6,-2.7874),
	'97Mo':(6,-2.8462),
	'99Ru':(6,-1.9553),
	'101Ru':(6.-2.192),
	'Rh':(2,-1.3476),
	'105Pd':(6,-1.957),
	'107Ag':(2,-1.7330),
	'109Ag':(2,-1.9924),
	'111Cd':(2,-9.0689),
	'113Cd':(2,-9.4868),
	'113In':(10,9.3652),
	'115In':(10,9.3854),
	'In':(10,9.3854),
	'115Sn':(2,-14.0074),
	'117Sn':(2,-15.2606),
	'119Sn':(2,-15.9656),
	'121Sb':(6,10.2549),
	'123Sb':(8,5.5530),
	'123Te':(2,-11.2346),
	'125Te':(2,-13.5451),
	'I':(6,8.5776),
	'129Xe':(2,-11.8601),
	'131Xe':(4,3.5158),
	'Cs':(8,5.6232),
	'135Ba':(4,4.2581),
	'137Ba':(4,4.7633),
	'138La':(11,5.6614),
	'139La':(8,6.0610),
	'La':(8,6.0610),
	'Pr':(6,13.0355),
	'143Nd':(8,-2.319),
	'145Nd':(8,-1.429),
	'147Sm':(8,-1.7747),
	'149Sm':(8,-1.4631),
	'151Eu':(6,10.5854),
	'153Eu':(6,4.6744),
	'155Gd':(4,-1.317),
	'157Gd':(4,-1.727),
	'Tb':(4,10.23),
	'161Dy':(6,-1.4653),
	'163Dy':(6,2.0507),
	'Ho':(8,9.0881),
	'167Er':(8,-1.2281),
	'Tm':(2,-3.531),
	'171Yb':(2,7.5259),
	'173Yb':(6,-2.0730),
	'175Lu':(8,4.8624),
	'Lu':(8,4.8624),
	'176Lu':(15,3.451),
	'177Hf':(8,1.7281),
	'179Hf':(10,-1.0856),
	'180Ta':(19,4.04),
	'181Ta':(8,5.1625),
	'Ta':(8,5.1625),
	'183W':(2,1.7956),
	'185Re':(6,9.717),
	'187Re':(6,9.817),
	'187Os':(2,0.9856),
	'189Os':(4,3.3535),
	'191Ir':(4,0.766),
	'193Ir':(4,0.832),
	'195Pt':(2,9.2920),
	'Au':(4,0.7406),
	'199Hg':(2,7.7121),
	'201Hg':(4,-2.8468),
	'203Tl':(2,24.7310),
	'205Tl':(2,24.9742),
	'207Pb':(2,9.0338),
	'Bi':(10,6.9628),
	'235U':(8,-0.83),
# RE3+ electron shortcuts
	'fCe':(6,6.0/7.0*(-28024.95164)/2.0),
	'fPr':(9,4.0/5.0*(-28024.95164)/2.0),
	'fNd':(10,8.0/11.0*(-28024.95164)/2.0),
	'fPm':(9,3.0/5.0*(-28024.95164)/2.0),
	'fSm':(6,2.0/7.0*(-28024.95164)/2.0),
	# omit Eu (J=0)
	'fGd':(8,2.0*(-28024.95164)/2.0),
	'fTb':(13,3.0/2.0*(-28024.95164)/2.0),
	'fDy':(16,4.0/3.0*(-28024.95164)/2.0),
	'fHo':(17,5.0/4.0*(-28024.95164)/2.0),
	'fEr':(16,6.0/5.0*(-28024.95164)/2.0),
	'fTm':(13,7.0/6.0*(-28024.95164)/2.0),
	'fYb':(8,8.0/7.0*(-28024.95164)/2.0),
	}
	# RE hyperfine constants from Abragam & Bleaney table 5,5 "A_J"
	REHyperfine={
	"Pr":1093.0,
	"143Nd":-220.3,
	"145Nd":-136.9,
	"147Sm":-240.0,
	"149Sm":-194.0,
	"Tb":530.0,
	"161Dy":-109.5,
	"163Dy":152.4,
	"Ho":812.1,
	"167Er":-125.3,
	"Tm":-393.5,
	"171Yb":887.2,
	"173Yb":-243.3,
	}
	spins=[]
	labels={} # name -> index mapper
	slices=1 # until requested

	try:
		slicetimes=map(float,pars["pulsed"])
		slices=len(slicetimes)+1
		#print "sliced mode engaged, time boundaries are ",slicetimes
	except:
		slices=1

	try:
		dynstates=int(pars["dynamic"][0]) # num of states to convert between, 1 for just relaxation
		if(dynstates<1 or dynstates>1000):
			raise Exception("Bad number of sites")
		convtable=numpy.zeros([slices,dynstates,dynstates],dtype=numpy.float)
		poptable=numpy.zeros([dynstates],dtype=numpy.float)
		pops=numpy.zeros([dynstates],dtype=numpy.float)
		dssize=dynstates*slices
	except:
		dynstates=0 # static or RF calc
		dssize=1*slices

	try:
		spinnames=pars["spins"]
	except:
		raise Exception("List of spins must be provided!")
	# spin name forms:
	# H, F, Mu, e : look in table and get I and gamma
	# isotopes 63Cu, 13C, etc: look in table. If >90% abundant then symbol alone accepted as well e.g. N == 14N
	# 1/2, 0.5, 1, 3/2, etc: use I. Expect separate gamma value to follow
	# number identical atoms with suffixes if wanted (H1,H2,H3)
	# refer later with atom name (+suffix if given) or positional param 0..n (required for anonymous I values)
	gammas=[0]*dssize
	for i in range(dssize):
		gammas[i]=[]
	for i in range(len(spinnames)):
		atom = spinnames[i].strip()
		baseatom=atom.rstrip("0123456789")
		try:
			spins.append(PeriodicTable[baseatom][0])
			for d in range(dssize):
				gammas[d].append(PeriodicTable[baseatom][1])
			labels[atom]=i
			labels[str(i)]=i
		except:
			try:
				if(atom[-2:]=="/2"):
					spins.append(int(atom[:-2])+1)
				else:
					spins.append(int(float(atom)*2+1.01))
				for d in range(dssize):
					gammas[d].append(float("NaN"))
				labels[str(i)]=i
			except:
				element=baseatom.lstrip("0123456789")
				ek={x.lstrip("0123456789") for x in PeriodicTable.keys()}
				if(baseatom in ek):
					raise Exception("Atom "+atom+" has several isotopes, choose one")
				elif(element in ek):
					raise Exception("Atom "+atom+" is not a stable isotope with I>0, choose another")
				else:
					raise Exception("Atom "+atom+" is not in the table")

	spmat=createSpinMat(spins)
	if(dynstates>0 or dssize>1):
		Hams=[]
		reltable=numpy.zeros([dssize,len(spins)],dtype=numpy.float)
		reldetails=numpy.zeros([dssize,len(spins),4],dtype=numpy.float)
		for i in range(dssize):
			Hams.append(numpy.zeros_like(spmat[0,0,:,:]) )
	else:
		Ham=numpy.zeros_like(spmat[0,0,:,:])
		Hams=[Ham] # for convenience
	coords=[0]*dssize
	for i in range(dssize):
		coords[i]=[None]*len(labels)
	cfield=[0]*dssize
	for i in range(dssize):
		cfield[i]=[None]*len(labels)
	# second loop to generate H(hyperfine)
	for (key,val) in pars.items():
		#print "2nd loop looking at ",key
		# All have 1st arg @n for specific site or omitted if all or ignored if not dynamic
		if(key[0:2]=="a("):
			spj=key[2:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spj[0])
			for i in ii:
				try:
					(spin1t,spin2t)=sp
				except:
					spin1t=sp[0]
					spin2t="e"
				try:
					spin1=labels[spin1t]
				except:
					raise Exception("Unknown spin "+spin1t)
				try:
					spin2=labels[spin2t]
				except:
					raise Exception("Unknown spin "+spin2t)
				#Aaa=numpy.array(Aval.split(","),dtype=float)
				if(len(val)==1):
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],val[0]) # single A (isotropic)
				elif (len(val)==2):
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],AxialHFT(val[0],val[1],(0.0,0.0,1.0))) # A and D, axis along z by default
				elif (len(val)==3):
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],numpy.array(((val[0]-val[1]/2+val[2]/2,0.0,0.0),(0.0,val[0]-val[1]/2-val[2]/2,0.0),(0.0,0.0,val[0]+val[1]))) ) # A,D,E as in old Quantum, fixed axes for powder use
				elif (len(val)==5):
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],AxialHFT(val[0],val[1],val[2:5])) # A, D, and axis (3)
				elif (len(val)==6):
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],numpy.array(((val[0],val[3],val[4]),(val[3],val[1],val[5]),(val[4],val[5],val[2]))) ) # Tensor for A, excluding duplicate terms. Axx, Ayy, Azz, Axy, Axz, Ayz
				elif(len(val)==9):
					# tensor with principal values and axes: A1,x1,y1,z1,A2,x2,y2,z2,A3. (x3,y3,z3 calculated, x2,y2,z2 adjusted if necessary, all axes will be normalised)
					# put into function for Francis' benefit!
					addHyperfine(Hams[i],spmat[spin1],spmat[spin2],NonAxialHFT(val[0],val[4],val[8],val[1:4],val[5:8]))
					#xyz1=normalised(val[1:4])
					#xyz3=normalised(numpy.cross(val[1:4],val[5:8]))
					#xyz2=normalised(numpy.cross(xyz3,val[1:4]))
					#addHyperfine(Hams[i],spmat[spin1],spmat[spin2],val[0]*numpy.outer(xyz1,xyz1)+val[4]*numpy.outer(xyz2,xyz2)+val[8]*numpy.outer(xyz3,xyz3) )
				elif(len(val)==4):
					# four zeros, try to give "default" e.g. for RE ions
					if(spin1t in REHyperfine and spin2t[0]=="f" and spin1t[-2:]==spin2t[-2:]):
						addHyperfine(Hams[i],spmat[spin1],spmat[spin2],REHyperfine[spin1t])
					else:
						raise KeyError("Don't have a default hyperfine coupling for "+spin1t+" to "+spin2t)
		elif(key[0:2]=="q("):
			spj=key[2:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spj[0])
			for i in ii:
				try:
					spin1=labels[sp[0]]
				except:
					raise Exception("Unknown spin "+spin1t)
			# q(@site,spin)=nuq
				if(len(val)==1):
					addQuadrupole(Hams[i],spmat[spin1],AxialHFT(0.0,val[0],(0.0,0.0,1.0))) # hard coded 001 axis
			# q(@site,spin)=nuq,eta
				elif(len(val)==2):
					addQuadrupole(Hams[i],spmat[spin1],numpy.array([[-val[0]*(1+val[1])/2.0,0.0,0.0],[0.0,-val[0]*(1-val[1])/2.0,0.0],[0.0,0.0,val[0]]]) ) # hard coded 001 and 100 axes
			# q(@site,spin)=nuq,x,y,z axial with arb. axis
				elif(len(val)==4):
					addQuadrupole(Hams[i],spmat[spin1],AxialHFT(0.0,val[0],val[1:4])) # user specified axis, axial
			# q(@site,spin)=6 component tensor hopefully with Tr()=0
				elif(len(val)==6):
					addQuadrupole(Hams[i],spmat[spin1],numpy.array([[val[0],val[3],val[4]],[val[3],val[1],val[5]],[val[4],val[5],val[2]]])) # full tensor
				elif(len(val)==9):
					# tensor with principal values and axes: nuq1,x1,y1,z1,nuq2,x2,y2,z2,nuq3. (x3,y3,z3 calculated, x2,y2,z2 adjusted if necessary, all axes will be normalised)
					xyz1=normalised(val[1:4])
					xyz3=normalised(numpy.cross(val[1:4],val[5:8]))
					xyz2=normalised(numpy.cross(xyz3,val[1:4]))
					addQuadrupole(Hams[i],spmat[spin1],val[0]*numpy.outer(xyz1,xyz1)+val[4]*numpy.outer(xyz2,xyz2)+val[8]*numpy.outer(xyz3,xyz3) )
			
			pass # quadrupole splitting. Specify site?
		elif(key[0:2]=="r("):
			# coords are in crystal lattice units if known
			spint=key[2:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spint[0])
			spin1=labels[sp[0]]
			for i in ii:
				coords[i][spin1]=numpy.dot(val,pars["_unCell"].getBinv())
		elif(key[0:6]=="gamma("):
			spint=key[6:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spint[0])
			for i in ii:
				if(len(val)==1):
					gammas[i][labels[sp[0]]]=float(val[0])
				elif(len(val)==6):
					gammas[i][labels[sp[0]]]=numpy.array(((val[0],val[3],val[4]),(val[3],val[1],val[5]),(val[4],val[5],val[2])))
				elif(len(val)==9):
					xyz1=normalised(val[1:4])
					xyz3=normalised(numpy.cross(val[1:4],val[5:8]))
					xyz2=normalised(numpy.cross(xyz3,val[1:4]))
					gammas[i][labels[sp[0]]]=val[0]*numpy.outer(xyz1,xyz1)+val[4]*numpy.outer(xyz2,xyz2)+val[8]*numpy.outer(xyz3,xyz3)
		elif(key[0:2]=="g("): # as for gamma, but in units of Bohr magneton instead of MHz/T
			muB=13996.24555 # Bohr magneton in MHz/T 
			spint=key[2:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spint[0])
			for i in ii:
				if(len(val)==1):
					gammas[i][labels[sp[0]]]=float(val[0])*muB
				elif(len(val)==6):
					gammas[i][labels[sp[0]]]=numpy.array(((val[0],val[3],val[4]),(val[3],val[1],val[5]),(val[4],val[5],val[2])))*muB
				elif(len(val)==9):
					xyz1=normalised(val[1:4])
					xyz3=normalised(numpy.cross(val[1:4],val[5:8]))
					xyz2=normalised(numpy.cross(xyz3,val[1:4]))
					gammas[i][labels[sp[0]]]=(val[0]*numpy.outer(xyz1,xyz1)+val[4]*numpy.outer(xyz2,xyz2)+val[8]*numpy.outer(xyz3,xyz3))*muB
		elif(key[0]=="b" and key[1] in "246" and key[2] in "0123456" and key[3]=="("): # Crystal field, b<NN>(atom,@site,|slice)=coeff
			cfterm=key[1:3]
			spint=key[4:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spint[0])
			spin1=labels[sp[0]]
			for i in ii:
				if(cfield[i][spin1] is None):
					cfield[i][spin1]={} # dict to store terms
				cfield[i][spin1]["B"+cfterm]=float(val[0])
		elif(key[0:2]=="ib" and key[2] in "246" and key[3] in "123456" and key[4]=="("): # Imaginary part of Crystal field, ib<NN>(atom,@site,|slice)=coeff. Note no ibN0 coefficients.
			cfterm=key[2:4]
			spint=key[5:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,spint[0])
			spin1=labels[sp[0]]
			for i in ii:
				if(cfield[i][spin1] is None):
					cfield[i][spin1]={} # dict to store terms
				cfield[i][spin1]["IB"+cfterm]=float(val[0])
	# now all coords known, process to give dipolar. Sites?
	for d in range(dssize):
		for g in gammas[d]:
			if (numpy.any(numpy.isnan(g))):
				raise Exception("Not all gamma values are known")
	for d in range(dssize):
		for i in range(len(spins)):
			for j in range(i+1,len(spins)):
				if(coords[d][i] is not None and coords[d][j] is not None):
					addDipolar(Hams[d],spmat[i],spmat[j],coords[d][i],gammas[d][i],coords[d][j],gammas[d][j])
	# now all CF terms known, add CF Hamiltonian
	for d in range(dssize):
		for i in range(len(spins)):
			if(cfield[d][i]):
				addSmallHam(Hams[d],spmat[i],calcCFHam(spins[i],cfield[d][i]))
	# now process orientation string
	Bmag=numpy.zeros([slices,3]) # 0.0
	rfomega=numpy.zeros([slices]) # 0.0 # if non-zero, will trigger RF mode
	rfbmag=numpy.zeros([slices])
	rfelip=numpy.zeros([slices])
	rfphase=numpy.zeros([slices])
	rfphaserand=[False]*slices
	rrf=0 # numpy.zeros([slices],dtype=numpy.int)
	iterator=iter(([0.0,0.0,1.0],)) # default zero field, one detector orientation along z
	for (key,val) in pars.items():
		if(key=="lfuniform"):
			iterator=uniformLF(int(val[0]))
		elif(key=="lfrandom"):
			iterator=randomLF(int(val[0]))
		if(key=="rfuniform"):
			iterator=uniformRF(int(val[0]))
		elif(key=="rfrandom"):
			iterator=randomRF(int(val[0]))
		elif(key=="tfuniform"):
			iterator=uniformTF(int(val[0]))
		elif(key=="tfrandom"):
			iterator=randomTF(int(val[0]))
		elif(key=="pquniform"):
			iterator=uniformPQ(int(val[0]))
		elif(key=="pqrandom"):
			iterator=randomPQ(int(val[0]))
		elif(key=="tfphiuniform"):
			iterator=uniformPQ(int(val[0]),val[1])
		elif(key=="tfphirandom"):
			iterator=randomPQ(int(val[0]),val[1])
		elif(key=="lfinh"):
			iterator=inhomogeneousLF(int(val[0]))
		elif(key=="lf"):
			iterator=iter((numpy.array(val),))# 3-vec for orientation
		elif(key=="lfaxes"): # 3-vec, to be permuted except V and -V are equivalent
			a=val[0]
			b=val[1]
			c=val[2]
			axset={} # dict for auto removal of duplicates!
			for permute in ((a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a)):
				for reflectX in (-1,1):
					for reflectY in (-1,1):
						p=(permute[0]*reflectX,permute[1]*reflectY,permute[2])
						if(p[2]<0 or (p[2]==0 and (p[1]<0 or (p[1]==0 and p[0]<0)))):
							p=(-p[0],-p[1],-p[2])
						axset[p]=1
			iterator=iter(axset.keys())
		elif(key=="tf"):
			if(len(val)==6):
				iterator=iter(((val[0:3],val[3:6],val[3:6],val[3:6]),)) # 2 * 3-vec: B, spin0
			elif(len(val)==9):
				iterator=iter(((val[0:3],val[3:6],val[6:9],val[3:6]),)) # 2 * 3-vec: B, spin0
			elif(len(val)==12):
				iterator=iter(((val[0:3],val[3:6],val[6:9],val[9:12]),)) # 2 * 3-vec: B, spin0
			else:
				raise Exception("expected tf B, spin [,detect [,RF]]")
		elif(key=="selectaxis"):
			# use for multiple fits across TF groups, TF & LF, etc (single crystal)
			# args N [,a,b,c [d,e,f]]
			# spin along (a,b,c), default z = 0,0,1
			# field along axis N=(0=x 1=y 2=z)
			# detect along axis N=N+(0=x 3=y 6=z)
			# RF along y (unless d,e,f given)
			# usual cases LF N=8 TF (MuSR) N=6 & Quadrature N=3
			axes1=((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
			if(len(val)>=4):
				spinax1=val[1:4]
			else:
				spinax1=(0.0,0.0,1.0)
			if(len(val)>=7):
				rfax1=val[4:7]
			else:
				rfax1=(0.0,1.0,0.0)
			iterator=iter(((axes1[int(val[0]) % 3],spinax1,axes1[int(val[0]) /3], rfax1),))
		elif(key[0:9]=="bmagGauss"): # must check for this before bmag, which is a substring!
			ii,sp=EnumerateSites(slices,0,key[9:])
			#print "key=",key," to ",ii," and ",sp," into Hams of len ",len(Hams)
			for i in ii:
				if(len(val)==1):
					Bmag[i,:]=[float(val[0])/10000.0,0.0,0.0]
				elif(len(val)==3):
					Bmag[i,:]=map(lambda x: float(x)/10000.0,val)
				else:
					raise Exception("bmagGauss should be a single number or a 3-vector")
		elif(key[0:4]=="bmag"):
			ii,sp=EnumerateSites(slices,0,key[4:])
			#print "key=",key," to ",ii," and ",sp," into Hams of len ",len(Hams)
			for i in ii:
				if(len(val)==1):
					Bmag[i,:]=[float(val[0]),0.0,0.0]
				elif(len(val)==3):
					Bmag[i,:]=map(float,val)
				else:
					raise Exception("bmag should be a single number or a 3-vector")
		# RF mode
		# brf <mag> <freq> [ <phase> [ <ellipticity> [ <rrf> ]]]
		# ellipticity=0: linear, +1: rotating "with" precession about B0, -1: rotating "opposite" to B0 (assuming bmag>0 and gamma>0)
		# requires freq>0
		elif(key[0:3]=="brf"):
			ii,sp=EnumerateSites(slices,0,key[3:])
			#print "key=",key," to ",ii," and ",sp," into Hams of len ",len(Hams)
			for i in ii:
				rfbmag[i]=float(val[0])
				rfomega[i]=2.0*math.pi*float(val[1]) # user specifies freq in MHz
				if(rfomega[i]<=0):
					raise ValueError("Frequency must be positive")
				if(len(val)>2):
					if(val[2]=="r" or val[2]=="R"):
						rfphaserand[i]=True
					else:
						rfphase[i]=float(val[2])*math.pi/180.0 # user specifies degrees
				else:
					rfphase[i]=0.0
				if(len(val)>3):
					rfelip[i]=float(val[3])
				else:
					rfelip[i]=0.0
				if(len(val)>4):
					rrf=int(val[4])
			
		# dynamic params:
		# Relax(@state,spin)=lambda (or Relax(lambda) for all)
		# ..=lambda,pol,x,y,z for relaxation to (partially) polarised state
		# Convert(from,to)=lambda
		# or Convert(from,to)=lambda,relaxthis
		# populations=a,b,c...
		# pop="e" calcs equilibrium,"u"=uniform among those not set
		# default unpopulated state
		elif(key[0:6]=="relax("):
			stsp=key[6:].split(")")
			ii,sp=EnumerateSites(slices,dynstates,stsp[0])
			spin1=labels[sp[0]]
			for i in ii:
				reltable[i,spin1]=float(val[0])
				if(len(val)==5):
					reldetails[i,spin1,:]=val[1:5]
		elif(key[0:8]=="convert("):
			stsp=key[8:].split(")")
			ii,sp=EnumerateSites(slices,0,stsp[0])
			try:
				(cfrom1,cfrom2)=sp[0].split("-")
			except:
				cfrom1=sp[0]
				cfrom2=sp[0]
			try:
				(cto1,cto2)=sp[1].split("-")
			except:
				cto1=sp[1]
				cto2=sp[1]
			for i in ii:
				for cfrom in range(int(cfrom1),int(cfrom2)+1):
					for cto in range(int(cto1),int(cto2)+1):
						if(cfrom!=cto):
							convtable[i,cfrom,cto]=float(val[0])
		elif(key[0:4]=="pop("):
			i=key[4:].split(")")
			i=int(i[0].strip("@"))
			poptable[i]=val[0]
			#if(poptable[i][0]=="e"):
			#	poptable[i]=-1.0
			#if(poptable[i][0]=="u"):
			#	poptable[i]=-2.0

	# process goniometer and crystal last as these iterators must "decorate" the main LF/TF ones
	if("goniometer" in pars):
		val=pars["goniometer"]
		if(len(val)==3):
			gonlist=(("y",val[0]),("z",val[1]),("y",val[2]))
		else:
			gonlist=zip(*[iter(val)]*2) # by twos
		gon=numpy.eye(3)
		for r in gonlist:
			c=math.cos(float(r[1])*math.pi/180.0)
			s=math.sin(float(r[1])*math.pi/180.0)
			if(r[0] in ["x","a"]):
				g2=numpy.array([[1,0,0],[0,c,s],[0,-s,c]])
			if(r[0] in ["y","b"]):
				g2=numpy.array([[c,0,-s],[0,1,0],[s,0,c]])
			if(r[0] in ["z","c"]):
				g2=numpy.array([[c,s,0],[-s,c,0],[0,0,1]])
			if(r[0] in ["x","y","z"]):
				gon=numpy.dot(g2,gon)
			else:
				gon=numpy.dot(gon,g2)
		iterator=GoniometerIterator(iterator,gon)

	if("_crystalEqvs" in pars):
		# "crystalsite" : x y z
		# return transformation matrices for all sites equivalent to representative point xyz
		# no need here to specify unit cell but must know space group
		iterator=CrystalIterator(iterator,pars["_crystalEqvs"])
	
	# debug use, split out averages
	if("selectorient" in pars):
		iterator=SelectOrientIterator(iterator,pars["selectorient"][0])
	
	# now all population data collected, calc eqm pops
	# exactly-specified ones first, if any (scaled down if the sum exceeds 1.0, or up if <1.0 and none unspecified)
	# then equilibrium ones "e" (if any) set from balance within them (scaled) (using 1st of N slices if more than one)
	# all remaining ones "u" equally distributed
	if(dynstates>1):
		ptotal=0.0
		known=0
		for i in range(dynstates):
			if(poptable[i]>=0):
				known=known+1
				pops[i]=poptable[i]
				ptotal=ptotal+poptable[i]
		if(known==dynstates or ptotal>1.0):
			for i in range(dynstates):
				pops[i]=pops[i]/ptotal
			ptotal=1.0
		eqmap=[] # states participating in equilibrium
		sharemap=[] # to sharre remaining
		for i in range(dynstates):
			if(poptable[i]==-1):
				eqmap.append(i)
			if(poptable[i]==-2):
				sharemap.append(i)
		if(len(eqmap)>0 and ptotal<1.0):
			try:
				eqs=calc_eqm(convtable[0,:,:].take(eqmap,axis=0).take(eqmap,axis=1)) # calc_eqm result is normalised
				pops[eqmap]=eqs*((1.0-ptotal))
				ptotal=1.0
			except:
				sharemap.extend[eqmap]
		if(len(sharemap)>0 and ptotal<1.0):
			pops[sharemap]=(1.0-ptotal)/len(sharemap)
			ptotal=1.0
		if(ptotal<1.0):
			raise Exception("Population error")
	else:
		pops=[1.0]

	# pre-check
	if(dynstates>0):
		if(not numpy.all(numpy.isfinite(pops))):
			print "pops=",pops
			raise Exception("Inf or NaN found in population array")
		for i in range(dssize):
			if(not numpy.all(numpy.isfinite(Hams[i]))):
				raise Exception("Inf or NaN found in sub Hamiltonian"+str(i))


	for orient in iterator:
		if(orient[1] is True):
			orient=orient[0]
			inverted=-1.0
		elif(orient[1] is False):
			orient=orient[0]
			inverted=+1.0
		else:
			inverted=+1.0
		if(len(orient)==2): # TF and crystal iterators return 2-tuple of 3-tuples
			fieldax=orient[0]
			beam=orient[1]
			detector=beam
			rfcoil=beam
		elif(len(orient)==4): # general TF and crystal iterators returns 4-tuple of 3-tuples
			fieldax=orient[0]
			beam=orient[1]
			detector=orient[2]
			rfcoil=orient[3]
		elif(len(orient)==3): # LF iterators return one 3-tuple
			fieldax=orient
			beam=orient
			detector=orient
			rfcoil=orient # ideally something perpendicular?
		else:
			raise Exception("Iterator returned "+str(orient))
		rfcoil2=numpy.cross(rfcoil,fieldax) # for rotating RF
		if("initialspin" in pars and pars["initialspin"][0]=="none"):
			rho0spin=None
			rho0 = numpy.identity(spmat.shape[-1])/spmat.shape[-1]
		else:
			try:
				rho0spin=labels[pars["initialspin"][0]]
			except:
				rho0spin=0
			rho0 = createInitialDensMat(spmat[rho0spin],beam)
		try:
			detspin=labels[pars["detectspin"][0]]
		except:
			detspin=0
		if("morespin" in pars):
			morespinlist=pars["morespin"]
			while(len(morespinlist)>=2):
				# morespin <spin>,x,y,z,... or <spin>,<axis> where axis="beam","field","detector","rf"
				# or <spin1>,TS,<spin2>  where TS="T" or "S" (singlet/triplet)
				# add "ground" <e> - populate all states with energy < e, considering all spins except "muon" but with all other interactions included
				# thermal <T> - similar but with thermal populations, temperature T
				if(morespinlist[0]=="thermal"):
					# works from Ham[0], i.e. not taking site variation into account
					# convert from K to Mrad s-1
					if(rho0spin is not None):
						rhomore=BoltzmannDensityMatrix(Hams[0],float(morespinlist[1])*1.38064852E-23/1.0545718E-34/1.E6,spmat[rho0spin]) # exclude the "muon"
					else:
						rhomore=BoltzmannDensityMatrix(Hams[0],float(morespinlist[1])*1.38064852E-23/1.0545718E-34/1.E6,None) # everything in eqm, e.g. for EPR/NMR					
					morespinlist=morespinlist[2:]
				else:
					morespin=labels[morespinlist[0]]
					try:
						moreaxis=numpy.array(map(float,morespinlist[1:4]))*inverted
						morespinlist=morespinlist[4:]
					except:
						if(morespinlist[1]=="beam"):
							moreaxis=beam
						elif(morespinlist[1]=="field"):
							moreaxis=fieldax
						elif(morespinlist[1]=="detector"):
							moreaxis=detector
						elif(morespinlist[1]=="rf"):
							moreaxis=rfcoil
						elif(morespinlist[1]=="t"):
							moreaxis=None
							moretriplet=True
							morespin2=labels[morespinlist[2]]
							morespinlist=morespinlist[1:] # additional trim as well as the following 2
						elif(morespinlist[1]=="s"):
							moreaxis=None
							moretriplet=False
							morespin2=labels[morespinlist[2]]
							morespinlist=morespinlist[1:] # additional trim as well as the following 2
						else:
							raise Exception("Unknown axis for polarising the ",morespinlist[0])
						morespinlist=morespinlist[2:]
					if(moreaxis is not None):
						rhomore=createInitialDensMat(spmat[morespin],moreaxis)
					else:
						Edot=numpy.dot(spmat[morespin,0],spmat[morespin2,0])+numpy.dot(spmat[morespin,1],spmat[morespin2,1])+numpy.dot(spmat[morespin,2],spmat[morespin2,2])
						if(moretriplet):
							rhomore=(Edot/3.0+numpy.identity(spmat.shape[-1]))/spmat.shape[-1]
						else: # singlet
							rhomore=(numpy.identity(spmat.shape[-1])-Edot)/spmat.shape[-1]
				rho0=numpy.dot(rho0,rhomore)*spmat.shape[-1]
		scint = createDetectorOp(spmat[detspin],detector)
		if(numpy.any(rfomega)):
			# generate HamRFs
			HamsRF=[None]*slices
			HamsRFi=[None]*slices
			for k in range(slices):
				HamsRF[k]=numpy.zeros_like(Hams[0])
				if(rfphaserand[k]):
					rfphase[k]=numpy.random.random()*2*math.pi
				if(rfphase[k]!=0 or rfelip[k]!=0):
					HamsRFi[k]=numpy.zeros_like(Hams[0])
				for i in range(len(spins)):
					addZeeman(HamsRF[k],spmat[i],rfcoil,rfbmag[k]*gammas[k][i]*math.cos(rfphase[k]))
					if(rfphase[k]!=0):
						addZeeman(HamsRFi[k],spmat[i],rfcoil,rfbmag[k]*gammas[k][i]*math.sin(rfphase[k]))
					if(rfelip[k]!=0):
						addZeeman(HamsRF[k],spmat[i],rfcoil2,rfbmag[k]*gammas[k][i]*rfelip[k]*math.sin(rfphase[k]))
						addZeeman(HamsRFi[k],spmat[i],rfcoil2,-rfbmag[k]*gammas[k][i]*rfelip[k]*math.cos(rfphase[k]))
		if(rrf!=0 or pars.get("phasequad",(0,))[0]):
			detector2=numpy.cross(detector,fieldax)
			scinti=createDetectorOp(spmat[detspin],detector2)
		else:
			scinti=None
		if((scinti is not None) and not(numpy.any(rfomega))):
			scint=scint+(1.j)*scinti
		if(dynstates==0):
			# process (static)
			Hams1=numpy.array(Hams)
			for k in range(slices):
				for i in range(len(spins)):
					#print "adding Zeeman B=",Bmag," for spin ",i
					addZeeman(Hams1[k],spmat[i],fieldax,gammas[k][i]*Bmag[k][0])
					# off axis components too
					if(Bmag[k][1]!=0):
						addZeeman(Hams1[k],spmat[i],rfcoil,gammas[k][i]*Bmag[k][1])
					if(Bmag[k][2]!=0):
						addZeeman(Hams1[k],spmat[i],rfcoil2,gammas[k][i]*Bmag[k][2])
			if(slices>1):
				if(numpy.any(rfomega)):
					yield (spmat,Hams1,rho0,scint,rfomega,rrf,HamsRF,HamsRFi,scinti,slicetimes)
				else:
					yield (spmat,Hams1,rho0,scint,0.0,0,None,None,None,slicetimes)
			else:
				if(numpy.any(rfomega)):
					yield (spmat,Hams1[0],rho0,scint,rfomega[0],rrf,HamsRF[0],HamsRFi[0],scinti)
				else:
					yield (spmat,Hams1[0],rho0,scint,0.0,0,None,None,None)
		else:
			BigHams=[None]*slices
			BigRho=createDynamicDensMat(rho0,pops)
			if(not numpy.all(numpy.isfinite(BigRho))):
				raise Exception("Inf or NaN found in BigRho")
			BigScint=createDynamicSpinOp(scint,len(pops))
			if(not numpy.all(numpy.isfinite(BigScint))):
				raise Exception("Inf or NaN found in BigScint")
			for k in range(slices):
				Hams1=[0]*dynstates
				for j in range(dynstates):
					jj=j+k*dynstates
					Hams1[j]=numpy.array(Hams[jj],copy=True)
					for i in range(len(spins)):
						addZeeman(Hams1[j],spmat[i],fieldax,gammas[jj][i]*Bmag[k][0])
						if(Bmag[k][1]!=0):
							addZeeman(Hams1[j],spmat[i],rfcoil,gammas[jj][i]*Bmag[k][1])
						if(Bmag[k][2]!=0):
							addZeeman(Hams1[j],spmat[i],rfcoil2,gammas[jj][i]*Bmag[k][2])
					if(not numpy.all(numpy.isfinite(Hams1[j]))):
						raise Exception("Inf or NaN found in Ham"+str(j)+" with Zeeman")
				BigHams[k],SmallHamSize,nhams=buildDynamicHam(Hams1)
				if(not numpy.all(numpy.isfinite(BigHams[k]))):
					raise Exception("Inf or NaN found in BigHam")
				# conversions
				if(dynstates>1):
					for i in range(dynstates):
						for j in range(dynstates):
							if(i!=j and convtable[k,i,j]!=0.0):
								addConversion(len(pops),BigHams[k],i,j,convtable[k,i,j],0.0,0.0) # no relaxing spins or eqm pol for them, yet
				for i in range(dynstates):
					for j in range(len(spins)):
						if(reltable[k*dynstates+i,j] != 0.0):
							if(reldetails[k*dynstates+i,j,0]!=0):
								addSpinRelaxation(len(pops),BigHams[k],i,spmat[j],reltable[k*dynstates+i,j],reldetails[k*dynstates+i,j,1:4],reldetails[k*dynstates+i,j,0]) # with eqm polarisation
							else:
								addSpinRelaxation(len(pops),BigHams[k],i,spmat[j],reltable[k*dynstates+i,j],0.0,0.0) # without eqm polarisation
								
			if(slices>1):
				yield(spmat,BigHams,BigRho,BigScint,[0.0,],0,None,None,None,slicetimes)
			else:
				yield(spmat,BigHams[0],BigRho,BigScint,0.0,0,None,None,None)

# debug helper non-iterator. Sticks at last value if loop too long
def SelectOrientIterator(iterator,which):
	for i in range(int(which)+1):
		try:
			this=next(iterator)
		except:
			pass
	yield this

# preprocessor options
def processor_TimeSpectra(pars,ybins,ebins,dest,asym):
	# copy spectrum. "dest" is an integer
	ybins[dest,:]=asym[:]
	ebins[dest,:]=0.0

def processor_TimeSpectraPQ(pars,ybins,ebins,dest,asym):
	# copy spectra. "dest" is an integer
	ybins[dest,:]=asym[:].real
	ebins[dest,:]=0.0
	ybins[dest+1,:]=asym[:].imag
	ebins[dest+1,:]=0.0

def processor_IntegralAsym(pars,ybins,ebins,dest,asym):
	# just copy integral asymmetry. "dest" must be a tuple
	# print "dest is ",dest,"in array shape ",ybins.shape
	ybins[dest]=asym[0]
	ebins[dest]=0.0

def processor_FittedCurve(pars,ybins,ebins,dest,asym):
	# fit
	# needs Mantid, imports specific functions
	from mantid.simpleapi import Fit,DeleteWorkspace
	# copy pars into different spectra
	fitfn=pars["fitfunction"]
	# temporary workspace (reused per point) for fitting
	tmpws = pars["tmpws"]
	tmpws.dataY(0)[:]=asym
	fitstuff=Fit(Function=fitfn,InputWorkspace=tmpws,CreateOutput=True)
	# store in bin. "dest" is an integer,
	#print "fit pars ",fitstuff[3].column("Value")," to fit in ",(ybins[:,dest]).shape
	ybins[:,dest]=fitstuff[3].column("Value")
	ebins[:,dest]=fitstuff[3].column("Error")
	if("recycle" in pars):
		pff=pars["fitfunction"]
		# recycle named fit pars as initial value next time
		for (i,parname) in enumerate(fitstuff[3].column("Name")):
			if (parname in pars["recycle"]):
				matc=re.match("f([0-9]+)\.(.+)",parname)
				if(matc):
					baseparname=matc.group(2)
					funcindex=int(matc.group(1))
				else:
					baseparname=parname
					funcindex=0
				finder="(.*?name=){"+str(funcindex+1)+"}.*?"+baseparname+"=([0-9eE.+-]+)([;,]|$)"
				matc2=re.match(finder,pff)
				if(matc2):
					pff=pff[:matc2.start(2)]+str(fitstuff[3].column("Value")[i])+pff[matc2.end(2):]
				else:
					print "warning, couldn't find anywhere to recycle ",parname," into ",pff," using finder=",finder
		nextpars={"fitfunction":pff}
	else:
		nextpars=None
	#print fitfn," & ",pff
	DeleteWorkspace(fitstuff[2])
	DeleteWorkspace(fitstuff[3])
	DeleteWorkspace(fitstuff[4])
	return nextpars
	
def tidyup_FittedCurve(pars,x0,x1):
	from mantid.simpleapi import DeleteWorkspace
	try:
		DeleteWorkspace(pars["tmpws"])
	except:
		print "workspace '",pars["tmpws"],"' won't die!!!"

def processor_moments(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):
	# first or second moment, or amplitude, within freq range
	# ignore phases and relaxation
	f1s=pars["minfreq"]
	if(len(f1s)==2):
		f1=(float(f1s[0])+float(pars["bmag"][0])*float(f1s[1])) * 2.0*math.pi
	else:
		f1=float(f1s[0]) * 2.0*math.pi
	f2s=pars["maxfreq"]
	if(len(f2s)==2):
		f2=(float(f2s[0])+float(pars["bmag"][0])*float(f2s[1])) * 2.0*math.pi
	else:
		f2=float(f2s[0]) * 2.0*math.pi
	mt=pars["measure"][0]
	sum0=0.0
	sum1=0.0
	sum2=0.0
	#sum00=0
	if(len(bigcsin)==0):
		# results from Dynamic calc or PQ mode
		for i in range(len(bigomega)):
			f=(bigomega[i]).imag
			if(f>f1 and f<f2):
				a=abs(bigccos[i])
				#sum00=sum00+1
				sum0=sum0+a
				sum1=sum1+a*f
				sum2=sum2+a*f*f
	else:
		# results from plain or RF
		for i in range(len(bigomega)):
			f=abs(bigomega[i])
			if(f>f1 and f<f2):
				a=math.sqrt(bigccos[i]**2+bigcsin[i]**2)
				#if(a>1.0):
				#	print "coeff ",i," at f=",f," has ampl=",a
				#sum00=sum00+1
				sum0=sum0+a
				sum1=sum1+a*f
				sum2=sum2+a*f*f
	#print "moment sum complete, pts=",len(bigomega)," used=",sum00," numave=",numave," sum012=",sum0,sum1,sum2
	if(mt=="m0"):
		ybins[dest]=sum0/numave
		ebins[dest]=0.0
	elif(mt=="m1"):
		if(sum0>0):
			ybins[dest]=sum1/sum0/2/math.pi
			ebins[dest]=0.0
		else:
			ybins[dest]=0.0
			ebins[dest]=1.0
	elif(mt=="m2"):
		try:
			ybins[dest]=math.sqrt(sum2*sum0-sum1*sum1)/sum0/2/math.pi
			ebins[dest]=0.0
		except MathError:
			ybins[dest]=0.0
			ebins[dest]=1.0

def processor_freqspec(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):
	# frequency spectrum (binned) from coefficients.
	# ignore phases and relaxation for now!
	# freq limits as for moments. Shift versus B if requested, if so the X axis for display is not shifted (0 means Larmor freq)
	# respect sign of omega in PQ mode
	f1s=pars["minfreq"]
	if(len(f1s)==2):
		f1=(float(f1s[0])+float(pars["bmag"][0])*float(f1s[1])) * 2.0*math.pi
	else:
		f1=float(f1s[0]) * 2.0*math.pi
	f2s=pars["maxfreq"]
	if(len(f2s)==2):
		f2=(float(f2s[0])+float(pars["bmag"][0])*float(f2s[1])) * 2.0*math.pi
	else:
		f2=float(f2s[0]) * 2.0*math.pi
	nfbins=int(pars["ntbins"][0])
	yacc=numpy.zeros(nfbins)
	fs=nfbins/(f2-f1)
	if(len(bigcsin)==0):
		# results from Dynamic calc, or plain calc in PQ mode
		# note Dynamic single det mode gives paired equal and opposite freqs
		for i in range(len(bigomega)):
			f=(bigomega[i]).imag
			fi=int(math.floor((f-f1)*fs))
			if(fi>=0 and fi<nfbins):
				yacc[fi]=yacc[fi]+abs(bigccos[i])
	else:
		# results from plain or RF, non-PQ
		for i in range(len(bigomega)):
			f=abs(bigomega[i])
			fi=int(math.floor((f-f1)*fs))
			if(f>f1 and f<f2):
				yacc[fi]=yacc[fi]+math.sqrt(bigccos[i]**2+bigcsin[i]**2)	
	ybins[dest,:]=yacc/numave
	ebins[dest,:]=0.0

def processor_breitrabi(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):
	# freqs generated by code
	#(ii,jj)=numpy.tril_indices(n,-1)
	#omega[1:]=eval[ii]-eval[jj]
	# last n-1 elements are differences between 0..n-2 and n-1
	# uncertainty in absolute energy: but all interactions leave average e=0
	ns=len(pars["axis1"])
	lbo=len(bigomega)
	tmpybins=numpy.array(bigomega[lbo-ns:])
	tmpybins[0]=0.0
	offset=numpy.average(tmpybins)
	ybins[:,dest]=(offset-numpy.sort(tmpybins))/(numave*2.0*math.pi) # should never actually use averaging with B-R diag!
	ebins[:,dest]=0.0

def PreParseLoop(pars,hadaxis0,prog=None):
	# identify loops and measure lengths
	# keywords: loop0par = param,
	# loop0range=start,end,nsteps
	# rewrites param=xi into pars
	# multiple scans in parallel
	# use loop0par=param1;param2;param3 and loop0range=start1,end1,start2,end2,Nsteps
	# first loop is used for axis coords
	# Log scale if Npoints is negative. Can scan completely-negative region
	# loop1 likewise
	# only 1 loop max if Fitting or if Time Spectra
	# no loops at all is valid - single point / single bin / single spectrum WS may result
	# specials for angle, or A or D with anisotropic HFC, or 1 coord of position? put param[i], will find line and insert
	# what if both refer to same line eg. scanning X and Y components of r?
	# normally line need not exist before.
	# fills pars["axis0"]
	axctr=-1
	progLen=1
	for loopctr in range(2):
		if("loop"+str(loopctr)+"par" in pars and "loop"+str(loopctr)+"range" in pars):
			axctr=axctr+1
			while((loopctr>0 or not(hadaxis0)) and len(pars["axis"+str(axctr)])>1):
				axctr=axctr+1
				if(axctr>=2):
					raise Exception("Too many nested loops")
			lr=pars["loop"+str(loopctr)+"range"]
			lpn=pars["loop"+str(loopctr)+"par"].split(";")
			lpn=map(str.strip,lpn)
			nv=len(lpn)
			if(len(lr) != 2*nv+1):
				raise Exception("mismatch between loop parameters and ranges")
			for i in range(nv):
				pars["_loopvar"+str(axctr+i*100)]=lpn[i]
				if("axis"+str(axctr+i*100) in pars and hadaxis0):
					pass # in use as fit function, values already set
				elif(lr[-1]<0):
					pars["axis"+str(axctr+i*100)] = numpy.logspace(0.0, math.log10(lr[2*i+1]/lr[2*i]),-int(lr[-1]),endpoint=True,base=10.0)*lr[2*i]
				else:
					pars["axis"+str(axctr+i*100)] = numpy.linspace(lr[2*i],lr[2*i+1],int(lr[-1]))
			pars["axis"+str(axctr)+"name"]=(lpn[0],) # pars["loop"+str(loopctr)+"par"]
			pars["_axis"+str(axctr)+"ID"]=(axctr,) # was loopctr. Now redundant param?
			progLen=progLen*abs(int(lr[-1]))
	if(prog is not None):
		prog.resetNumSteps(progLen,0.0,1.0)

def ParseMeasureType(pars,prog=None):
	# extract loop params and measurements
	# measure can return N values (time spectrum, different detectors, etc)
	# time dependence or multi detectors bumps loop indices to next dimension
	# Fitted data option: create tmp WS for time spectrum, do dummy fit on it to get parameter names
	# do crystal structure identification here too (needs Mantid)
	method=0
	processor=None
	tidyup=None

	if("crystal" in pars):
		from mantid.geometry import SpaceGroupFactory,UnitCell
		vars=pars["crystal"]
		# Mantid tables might have difficulty with spaces. Allow compressed form (but match on full string first in case of duplicates)
		sglist=SpaceGroupFactory.getAllSpaceGroupSymbols()
		sgname=vars[0].strip()
		if sgname in sglist:
			spGrp=SpaceGroupFactory.createSpaceGroup(vars[0])
		else:
			sgname=sgname.translate(None," ")
			sgtrans={s.translate(None," "):s for s in sglist}
			if sgname in sgtrans:
				spGrp=SpaceGroupFactory.createSpaceGroup(sgtrans[sgname])
			else:
				raise ValueError("Space group "+vars[0]+" unknown")
		unCell=UnitCell(*map(float,vars[1:]))
		if not (spGrp.isAllowedUnitCell(unCell)):
			raise ValueError("Unit Cell and space group are incompatible")
	else: # 1:1 scale, no symmetry. Dummy functions independent of Mantid
		#spGrp=SpaceGroupFactory.createSpaceGroup("P 1")
		spGrp=None # not used anyway
		#unCell=UnitCell(1,1,1)
		class DummyUnitCell:
			def getBinv(self):
				return numpy.identity(3)
		unCell=DummyUnitCell()
	pars["_spGrp"]=spGrp
	pars["_unCell"]=unCell
	if(("crystal" in pars) and ("crystalsite" in pars)):
		# getCrystalEquivalents is in quantumtools (but doesn't import Mantid geometry stuff directly)
		pars["_crystalEqvs"]=getCrystalEquivalents(spGrp, unCell, pars["crystalsite"],"morespin" not in pars) # morespin might have fixed axis polarisation of e; add other commands here
		# for loop length adjustment (debug use):
		if("loop0par" in pars and "selectorient" in pars["loop0par"]):
			pars["loop0range"]=(0,len(pars["_crystalEqvs"])-1,len(pars["_crystalEqvs"]))
		if("loop1par" in pars and "selectorient" in pars["loop1par"]):
			pars["loop1range"]=(0,len(pars["_crystalEqvs"])-1,len(pars["_crystalEqvs"]))

	hadaxis0 = ("axis0" in pars)
	if(not(hadaxis0)):
		pars["axis0"]=[0.0]
	if("axis0extra" not in pars):
		pars["axis0extra"]=(0,) # point data unless stated elsewhere
	pars["axis0name"]=("-",)
	pars["_axis0ID"]=(-1,) # not linked to a loop
	pars["axis1"]=(0.0,)
	pars["axis1name"]=("-",) # not yet used
	pars["_axis1ID"]=(-1,)
	# measure types
	if "starttime" not in pars:
		pars["starttime"]=(0.0,)
	mtype=pars["measure"][0]
	if(mtype=="integral"):
		# expect par mulife=2.197, starttime=0, endtime=Inf
		if "endtime" not in pars:
			pars["endtime"]=(float("Inf"),)
		if "ntbins" not in pars:
			pars["ntbins"]=(1,)
		method=2 # accumulate average P(t)
		processor=processor_IntegralAsym
	elif (mtype=="timespectra"):
		# expect pars starttime=0, endtime=20,ntbins=1000
		if "endtime" not in pars:
			pars["endtime"]=(20.0,)
		if "ntbins" not in pars:
			pars["ntbins"]=(1000,)
		pars["axis0name"]=("Time","microsecond") # ("Time","\xB5s")
		method=2 # accumulate average P(t)
		processor=processor_TimeSpectra
	elif (mtype=="phasequad"):
		# expect pars starttime=0, endtime=20,ntbins=1000
		if "endtime" not in pars:
			pars["endtime"]=(20.0,)
		if "ntbins" not in pars:
			pars["ntbins"]=(1000,)
		pars["axis0name"]=("Time","microsecond") # ("Time","\xB5s")
		pars["axis1name"]=("Detectors",)
		pars["axis1"]=["0","90"]
		method=2 # accumulate average P(t)
		pars["phasequad"]=(1,)
		processor=processor_TimeSpectraPQ
	elif (mtype=="fit"):
		from mantid.api import WorkspaceFactory
		from mantid.simpleapi import Fit
		# expect par starttime=0, endtime=20,ntbins=1000
		if "endtime" not in pars:
			pars["endtime"]=(20.0,)
		if "ntbins" not in pars:
			pars["ntbins"]=(1000,)
		# expect FitFn=... as for Mantid
		# create temporary workspace
		pars["axis1name"]=("Parameters",)
		tmpws=WorkspaceFactory.create("Workspace2D",NVectors=1,XLength=int(pars["ntbins"][0]+1),YLength=int(pars["ntbins"][0]))
		pars["tmpws"]=tmpws
		fitfn=pars["fitfunction"]
		fitstuff=Fit(Function=fitfn,InputWorkspace=tmpws,CreateOutput=True)
		pars["axis1"]=fitstuff[3].column("Name")
		method=2 # accumulate average P(t)
		processor=processor_FittedCurve
		tidyup=tidyup_FittedCurve
	elif(mtype=="breitrabi"):
		# one spectrum per eigenstate
		# haven't yet parsed enough of model to know how many states! Do it...
		if "endtime" not in pars:
			pars["endtime"]=(float("Inf"),)
		if "ntbins" not in pars:
			pars["ntbins"]=(1,)
		method=1
		gen=ParseAndIterateModel(pars)
		(spmat,ham,rho0,scint,rffrq,rrf,j1,j2,j3)=gen.next()
		gen.close()
		# (spmat,Hams1[0],rho0,scint,0.0,0,None,None,None)
		nlevels=ham.shape[0]
		pars["axis1name"]=("Levels",)
		pars["axis1"]=numpy.arange(nlevels*1.0)
		processor=processor_breitrabi		
	elif(mtype=="m0" or mtype=="m1" or mtype=="m2"):
		# amplitude, first or second moment of line (within window)
		# minfreq and maxfreq can be 1 number (fixed) or 2 (A+B*bmag), in MHz and MHz/T
		if("minfreq" not in pars):
			pars["minfreq"]=(0.0,) # small but exclude zero
		if("maxfreq" not in pars):
			pars["maxfreq"]=(float("Inf"),) # anything
		if "endtime" not in pars:
			pars["endtime"]=(20.0,) # anything
		if "ntbins" not in pars:
			pars["ntbins"]=(1,)
		method=1
		processor=processor_moments
	elif(mtype=="freqspec"):
		# frequencies, binned
		# minfreq, maxfreq assumed set
		if "ntbins" not in pars:
			pars["ntbins"]=(1000,)
		if "endtime" not in pars:
			pars["endtime"]=(20.0,) # anything
		method=1
		processor=processor_freqspec
		pars["axis0"]=numpy.linspace(pars["minfreq"][0],pars["maxfreq"][0],pars["ntbins"][0]+1,endpoint=True) # these will be bin boundaries, maybe ref to Larmor freq
		pars["axis0name"]=("Frequency","MHz")
		pars["axis0extra"]=(1,) # binned data
	if((mtype=="timespectra" or mtype=="phasequad") and hadaxis0):
		timebins=pars["axis0"] # passed in from fit function
	elif(pars["ntbins"][0]==1):
		timebins=numpy.array([pars["starttime"][0],pars["endtime"][0]],dtype=numpy.float) # exact ends
	else:
		timebins=numpy.linspace(pars["starttime"][0],pars["endtime"][0],pars["ntbins"][0]+1,endpoint=True) # these will be bin boundaries
	if((mtype=="timespectra" or mtype=="phasequad")and not(hadaxis0)):
		pars["axis0"]=timebins
		pars["axis0extra"]=(1,) # X has one more point than Y and E; WS to be created later.
	elif (mtype=="fit"):
		tmpws.dataX(0)[:]=timebins
		# if("mulife" in pars):
		eevents=float(pars.get("fitstats",(10.0,))[0])*1.E6
		#tmpws.dataE(0)[:]=numpy.exp((timebins[:-1]+timebins[1:])/4.0/float(pars["mulife"]))/1000.0 # arbitrary scale factor!
		if("mulife" in pars): # exponential decay weighting, including events outside window
			tmpws.dataE(0)[:]=4.0/numpy.sqrt(eevents*(numpy.exp(-timebins[:-1]/float(pars["mulife"][0]))-numpy.exp(-timebins[1:]/float(pars["mulife"][0]))))
		else: # uniformly distributed events over whole range
			tmpws.dataE(0)[:]=4.0/numpy.sqrt(eevents*(timebins[1:]-timebins[:-1])/(timebins[-1]-timebins[0]))
		if (numpy.any(numpy.isnan(tmpws.dataE(0)))):
			print "E has at least one NaN"
			print tmpws.dataE(0)[0],"...",tmpws.dataE(0)[-1]
		if (numpy.any(numpy.isinf(tmpws.dataE(0)))):
			print "E has at least one Inf"
			print tmpws.dataE(0)[0],"...",tmpws.dataE(0)[-1]
			
	PreParseLoop(pars,hadaxis0,prog) # get axis lengths
	
	ybins=numpy.zeros([len(pars["axis1"]),len(pars["axis0"])-pars["axis0extra"][0]])
	ebins=numpy.zeros([len(pars["axis1"]),len(pars["axis0"])-pars["axis0extra"][0]])
	x0axis=pars["axis0"]
	x1axis=pars["axis1"]

	if("selectorient" in pars["axis1name"]):
		x1axis=[m[2] for m in pars["_crystalEqvs"]]
		print "added names to axis1"

	return (ybins,ebins,timebins,x0axis,x1axis,method,processor,tidyup)
		
def ParseAndIterateLoop(pars):
	# return addition to pars which sets loop position (as dict)
	# may have 0, 1 or 2 nested loops
	# also coords in "ybins" for storage of result (tuple for 2 loops)
	# no loops: emits dest=0 and dict={}
	if("_loopvar0" in pars):
		N0=len(pars["axis0"])
	else:
		N0=1
	if("_loopvar1" in pars):
		N1=len(pars["axis1"])
	else:
		N1=1
	#print "starting loop with N0=",N0,"and N1=",N1
	#print "pars=",pars
	for L0 in range(N0):
		for L1 in range(N1):
			loopy=[L0,L1]
			dest=(L1,L0) # reversed, outer loop over spectra and inner loop over bins
			d={}
			for axctr in range(2):
				axl=pars["_axis"+str(axctr)+"ID"][0]
				if(axl>=0):
					#print "found axis: axctr=",axctr,"axl=",axl
					i=0
					while("_loopvar"+str(axctr+i*100) in pars):
						# found a loop variable
						#print "found a loop variable i=",i
						x=pars["axis"+str(axctr+i*100)][loopy[axl]]
						lpnp=pars["_loopvar"+str(axctr+i*100)].split("[")
						if(len(lpnp)==1):
							# simple par=value, always overwrite
							#print "setting simple par ",lpnp[0],"=",x
							d[lpnp[0]]=(x,)
						else:
							# par[..]=value, need to loop up previous one and adjust
							lpi=lpnp[1].strip("[]")
							if lpnp[0] in d:
								d2=d[lpnp[0]][:]
							elif lpnp[0] in pars:
								d2=pars[lpnp[0]][:]
							else:
								raise Exception("Scanning one component of "+lpnp[0]+" but rest is undefined")
							lpi2=lpi.split(",")
							if(len(lpi2)==1):
								d2[int(lpi2[0])]=x
							elif(len(lpi2)==2): # 2 of 3 components allowing scanning of angle. Use for A axes, LF axis, etc
								d2[int(lpi2[0])]=math.cos(x*math.pi/180.0)
								d2[int(lpi2[1])]=math.sin(x*math.pi/180.0)
							else:
								raise Exception("Loop varying more than 2 components")
							#print "setting complex par ",lpnp[0],"=",d2
							d[lpnp[0]]=d2						
						i=i+1
				else:
					try:
						if(len(pars["axis"+str(axctr)])>1):
							dest=L0+L1 # the axis was not a loop, but does have non-trivial size, so overall index is not a tuple. Rethink if 3D data...
					except:
						pass
			yield(d,dest)

def RunModelledSystem(pars0,prog=None):
	# loop through Lstring
	# model from Mstring, oriented by Ostring (sum these)
	#  method 1: accumulate all the coefficients and evaluate at end (e.g. second moment of line), also for non-averageable such as dominant freq
	#  method 2: evaluate asym and then average those (most cases)
	# process result by Astring
	# return results array(s): N*sweep var[s] (1D), average[s] (m*ND), errors similarly

	(ybins,ebins,timebins,x0axis,x1axis,method,processor,tidyup)=ParseMeasureType(pars0,prog)
	# ybins is simulated final result(s) - numpy 2D array or tuple of them, zeroed and ready to fill
	# timebins for evaluation, if method 2
	# x1axis for Mantid workspace "X" (numpy 1D array)
	# x2axis is numeric value for "spectrum" axis of Mantid workspace ("None" means keep the detector mapping in a blank workspace?) (String axis for parameter names from fit)
	# call processor(pars,ybins,freqs,coeffs,method=1) or processor(pars,ybins,asym,method=2) at end
	#print "pars0=",pars0
	#print "ybins has shape",ybins.shape
	#print "ebins has shape",ebins.shape

	Dynamic = ("dynamic" in pars0)
	PQ=pars0.get("phasequad",(0,))[0]
	if("mulife" in pars0):
		mulife=float(pars0["mulife"][0])
	else:
		mulife=2.19703
	for (loopvar,dest) in ParseAndIterateLoop(pars0):
		#print "dest for this iteration will be ",dest
		#print "pars to set for this iteration: ",loopvar
		#print "pars0 insode loop: ",pars0
		pars=pars0.copy()
		for (k,v) in loopvar.items():
			pars[k]=v
		if(method==1):
			bigomega=None # numpy.array([],dtype=numpy.float)
			bigccos=None # numpy.array([],dtype=numpy.float)
			bigcsin=None # numpy.array([],dtype=numpy.float)
			numave=0
		else:
			yvals=None # don't yet know size
		for foo in ParseAndIterateModel(pars):
			# single time slice is (spmat,Ham,rho0,scint,rfomega,rrf,HamRF,HamRFi,scinti)
			# multi slices (spmat,Ham[],rho0,scint,rfomega[],rrf,HamRF[],HamRFi[],scinti,slicetimes)
			if(len(foo)==9):
				(spmat,Ham,rho0,scint,rfomega,rrf,HamRF,HamRFi,scinti)=foo
				slices=1
				# solve
				if(Dynamic):
					omega,ccos=solveDynamicResult(Ham,rho0,scint) # (BigHam,BigRho,BigScint)
					csin=[] # unused as ccos is complex
				else:
					if(rfomega==0):
						if(PQ):
							omega,ccos,csin=solveDensityMatComplex(Ham,rho0,scint)
						else:
							omega,ccos,csin=solveDensityMat(Ham,rho0,scint)
					else:
						omega,ccos,csin=solveRFDensityMat(Ham,HamRF,HamRFi,rfomega,rho0,scint,scinti,rrf,tsliceIn=pars.get("rfslices",(None,))[0])
						
				#omega,ccossin=solveDynamicResult(BigHam,BigRho,BigScint)
				if(method==1):
					# accumulate coeffs
					if(bigomega is None):
						bigomega=numpy.array(omega,copy=True)
						bigccos=numpy.array(ccos,copy=True)
						bigcsin=numpy.array(csin,copy=True)
					else:
						bigomega=numpy.concatenate((bigomega,omega))
						bigccos=numpy.concatenate((bigccos,ccos))
						if(len(csin)>0):
							bigcsin=numpy.concatenate((bigcsin,csin))
					numave=numave+1
				else:
					# evaluate asym
					if(Dynamic):
						#print "y=["
						#for (om,cc) in zip(omega,ccos):
						#	print "(",cc,")*exp((",om,")*t)"
						#print "]"
						if("tzero" in pars):
							ConvolveTimeResDynListed(omega,ccos,pars["tzero"])
						if(PQ):
							yvals1=EvaluateDynamicIntoBinsComplex(omega,ccos,1.0/mulife,timebins) # lam=1.0/2.197
						else:
							yvals1=EvaluateDynamicIntoBins(omega,ccos,1.0/mulife,timebins) # lam=1.0/2.197
					else:
						if("tzero" in pars):
							ConvolveTimeResolutionListed(omega,ccos,csin,pars["tzero"])
						if(PQ):
							yvals1=EvaluateDynamicIntoBinsComplex(omega,ccos,1.0/mulife,timebins) # lam=1.0/2.197
						else:
							yvals1=evaluateIntoBins(omega,ccos,csin,1.0/mulife,timebins) # lam=1.0/2.197
					if (yvals is None):
						yvals=yvals1
						numave=1
					else:
						yvals=yvals+yvals1
						numave += 1
			else:
				(spmat,Hams,rho0,scint,rfomegas,rrf,HamsRF,HamsRFi,scinti,slicetimes)=foo
				slices=len(slicetimes)+1
				# sliced mode ONLY with method=2 (processing the time spectrum)
				if("tzero" in pars):
					tzeros=pars["tzero"]
				else:
					tzeros=None
				if(Dynamic):
					# print "evaluating dynamic, Hams=",Hams," and slices=",slicetimes
					yvals1=EvaluateSlicedDynamic(Hams,rho0,scint,slicetimes,timebins,1.0/mulife,tzeros)
					# pass # implement EvaluateSlicedDynamic
				elif(numpy.any(rfomegas)):
					yvals1=EvaluateSlicedRF(Hams,HamsRF,HamsRFi,rfomegas,rho0,scint,scinti,rrf,slicetimes,timebins,1.0/mulife,tzeros)
				else:
					yvals1=EvaluateSliced(Hams,rho0,scint,slicetimes,timebins,1.0/mulife,tzeros)
				if (yvals is None):
					yvals=yvals1
					numave=1
				else:
					yvals=yvals+yvals1
					numave += 1
				
		if(method==1):
			recycle=processor(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave)
		else:
			yvals=yvals/numave
			recycle=processor(pars,ybins,ebins,dest,yvals)
		if(recycle is not None):
			for (k,v) in recycle.items():
				pars0[k]=v
		if (prog is not None):
			prog.report()
	if(tidyup):
		tidyup(pars0,x0axis,x1axis)
	
	if(len(x1axis)>1):
		return ([x0axis,x1axis],ybins, ebins)
	else:
		return ([x0axis],ybins, ebins)

def GetUserAxisName(axname):
	# given Axis Name from pars(), format it for user use
	# return 2 parts, variable and units
	# "" and "unit" means use a Mantid standard one, see print UnitFactoryImpl.Instance().getKeys()
	if(axname is None):
		return ("-","")
	if(len(axname)==2):
		return axname
	# is axname a loop variable?
	matches=(
		(r"^bmag","Magnetic Field","T"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)$","\\1 hyperfine constant","MHz"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[0\]$","\\1 isotropic hyperfine constant","MHz"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[1\]$","\\1 dipolar hyperfine constant","MHz"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+),(\w+)\)$","J coupling between \\1 and \\2","MHz"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+),(\w+)\)\[0\]$","Isotropic coupling between \\1 and \\2","MHz"),
		(r"^a\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+),(\w+)\)\[1\]$","Dipolar coupling between \\1 and \\2","MHz"),
		(r"^q\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)$","\\1 quadrupole splitting","MHz"),
		(r"^q\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[0\]$","\\1 axial quadrupole","MHz"),
		(r"^q\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[1\]$","\\1 quadrupole eta",""),
		(r"^r\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[0\]$","\\1 x-position","\xC5"),
		(r"^r\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[1\]$","\\1 y-position","\xC5"),
		(r"^r\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)\[2\]$","\\1 z-position","\xC5"),
		(r"^relax\((?:(?:\@[0-9]+,)|(?:\|[0-9]+,))*(\w+)\)$","Relaxation of \\1","\xB5s-1"),
		(r"^convert\(([0-9]+),([0-9]+)\)$","Conversion rate from state \\1 to \\2","\xB5s-1"),
		(r"^pulsed\[([0-9]+)\]$","Pulse time","\xB5s"),
		(r"^lf\[[0-9]+,[0-9]+\]$","Field angle","degrees"),
		(r"^tf\[[0-9]+,[0-9]+\]$","Field angle","degrees"),
		(r"^brf(?:\(.+\))*\[0\]$","RF magnetic field","T"),
		(r"^brf(?:\(.+\))*\[1\]$","RF frequency","MHz"),
		(r"^brf(?:\(.+\))*\[2\]$","RF phase","degrees"),
	)
	for m in matches:
		ma=re.match(m[0],axname[0])
		if(ma):
			uservar=ma.expand(m[1])
			userunit=ma.expand(m[2])
			return (uservar,userunit)
	return (axname[0],"") # not found, use raw value

# tester
#examples=("bmag","a(Mu)","r(H)[1]","brf(|1)[0]","convert(2,3)","q(@1,Cu)","notknown","lf[1,2]")
#for e in examples:
#	v,u=GetUserAxisName([e])
#	print e," -> ", v,"/",u

def GetUserYAxisName(axname):
	#label the Y (counts) axis if known
	if(axname=="timespectra"):
		return "Polarisation"
	elif(axname=="integral"):
		return "Integral Polarisation"
	elif(axname=="fit"):
		return "Fit parameters"
	elif(axname=="breitrabi"):
		return "Energy (MHz)"
	else:
		return axname
