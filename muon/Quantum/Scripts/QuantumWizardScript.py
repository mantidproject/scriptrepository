## QuantumWizardScript - a user friendly interface to Quantum, a program for solving spin evolution of the muon
## Author: James Lord
## Version 1.02, December 2015
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
from collections import OrderedDict

# GUI Quantum table generator
# stage 0 (own dialog?): enter name for table and for results workspace
# stage 1, select spins. Dialog with element names
# stage 2, select whether to use Dynamic mode (and if so how many states), RF mode, and/or Pulsed mode (if so how many slices)
# stage 2a: select pulse division times
# stage 3, ask what parameters vary with site or with pulse slice. Ask if HFCs cover all possible combinations or just to electron (or none, if no electron)
# stage 4, input coordinates
# stage 5, input hyperfine constants
# stage 6, input quadrupole splittings
# stage 7, if dynamic: input spin relaxation rates
# stage 8, if dynamic: input conversion rates
# stage 9: dynamic, enter initial populations
# stage 10, if RF, input freq, B1, etc
# stage 11: input averaging options, LF/TF
# stage 12: input magnetic field
# stage 13: select "measure" type and params
# stage 14: select what to vary in loops. (Restricted list?) and ranges.
# stage 15: generate table. Add spare rows.
# stage 16: Run it!
# stage 17: Plot it. If 1D, plot line spectrum. If 2D with few spectra, plot them. If 2D with many, do Colour Plot.
# Stage 1 dialog "algorithm"

PeriodicTable={
		'Mu':(2,135.5),
		'e':(2,-28025.0),
		'H':(2,42.5764),
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
		'49Ti':(8,2.4047),
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
		'109Ag':(2,1.9924),
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
		'235U':(8,-0.83) }
		
wiz1=QuantumWizard1Dialog()
s0=wiz1.getProperty("Spin0").value
s1=wiz1.getProperty("Spin1").value
s2=wiz1.getProperty("Spin2").value
s3=wiz1.getProperty("Spin3").value
s4=wiz1.getProperty("Spin4").value
s5=wiz1.getProperty("Spin5").value
spins=[s0,s1,s2,s3,s4,s5]
try:
	ns=spins.index('')
except:
	ns=len(spins)
if(spins.count('') + ns != len(spins)):
	raise Exception("Please don't leave gaps, I'm confused!")
spins=spins[:ns]
quadmap=[PeriodicTable[x][0]>2 for x in spins]
# unique names for multiple identical nuclei
for i in range(ns-1):
	if spins[i] in spins[i+1:]:
		wa=spins[i]
		k=1
		for j in range(len(spins)):
			if(spins[j]==wa):
				spins[j]=spins[j]+str(k)
				k=k+1
if(len(spins)==0):
	raise Exception("You need at least one spin!")
hypopt=wiz1.getProperty('HyperfineOptions').value

pars=OrderedDict({"spins":','.join(spins)})

wiz2=QuantumWizard2Dialog()
do_RF=False
do_dynamic=False
ct=wiz2.getProperty('Calculation type').value
sfc=wiz2.getProperty('States_for_conversion').value
if(ct=='Relaxation and Conversion'):
	pars['dynamic']=str(sfc)
	do_dynamic=True
if(ct=='RF resonance'):
	do_RF=True

is_anisotropic=False

if(do_dynamic and sfc>1):
	opts={}
	if(not(max(quadmap))):
		opts['Vary_quadrupole_between_sites']=False
	if(hypopt=='None (nuclear)'):
		opts['Vary_hyperfine_between_sites']=False
	wiz3=QuantumWizard3Dialog(**opts)
	differ_r=wiz3.getProperty('Vary_position_between_sites').value
	differ_a=wiz3.getProperty('Vary_hyperfine_between_sites').value
	differ_q=wiz3.getProperty('Vary_quadrupole_between_sites').value
else:
	differ_r=False
	differ_a=False
	differ_q=False
	
if(differ_r):
	for s in range(sfc):
		for i in range(len(spins)):
			wiz4=QuantumWizard4Dialog(Message="Enter coordinates for spin "+str(i)+"("+spins[i]+") in state "+str(s))
			used=wiz4.getParameters('Use_dipolar').value
			if(used):
				pars["r(@"+str(s)+","+spins[i]+")"]=str(wiz4.getProperty('x').value)+","+str(wiz4.getProperty('y').value)+","+str(wiz4.getProperty('z').value)
				is_anisotropic=True
else:
	for i in range(len(spins)):
		wiz4=QuantumWizard4Dialog(Message="Enter coordinates for spin "+str(i)+"("+spins[i]+")")
		used=wiz4.getProperty('Use_dipolar').value
		if(used):
			pars["r("+spins[i]+")"]=str(wiz4.getProperty('x').value)+","+str(wiz4.getProperty('y').value)+","+str(wiz4.getProperty('z').value)
			is_anisotropic=True

if(differ_a):
	for s in range(sfc):
		if(hypopt=='To One Electron Only' and ('e' in spins)):
			for i in range(len(spins)):
				e=spins.index('e')
				if(i!=e):
					wiz5=QuantumWizard5Dialog(Message="Enter hyperfine constant for spin "+str(i)+"("+spins[i]+") in state "+str(s))
					a=wiz5.getProperty("A").value
					d=wiz5.getProperty("D").value
					ax=wiz5.getProperty("AxisX").value
					ay=wiz5.getProperty("AxisY").value
					az=wiz5.getProperty("AxisZ").value
					if(d!=0):
						pars["a(@"+str(s)+","+spins[i]+")"]=",".join(map(str,[a,d,ax,ay,az]))
						is_anisotropic=True
					else:
						pars["a(@"+str(s)+","+spins[i]+")"]=str(a)					
		elif(hypopt=='All Pairs'):
			for i in range(len(spins)-1):
				for j in range(i+1,len(spins)):
					wiz5=QuantumWizard5Dialog(Message="Enter hyperfine constant between spins "+str(i)+"("+spins[i]+") and "+str(j)+"("+spins[j]+") in state "+str(s))
					a=wiz5.getProperty("A").value
					d=wiz5.getProperty("D").value
					ax=wiz5.getProperty("AxisX").value
					ay=wiz5.getProperty("AxisY").value
					az=wiz5.getProperty("AxisZ").value
					if(d!=0):
						pars["a(@"+str(s)+","+spins[i]+","+spins[j]+")"]=",".join(map(str,[a,d,ax,ay,az]))
						is_anisotropic=True
					else:
						pars["a(@"+str(s)+","+spins[i]+","+spins[j]+")"]=str(a)					
else:
	if(hypopt=='To One Electron Only' and ('e' in spins)):
		e=spins.index('e')
		for i in range(len(spins)):
			if(i!=e):
				wiz5=QuantumWizard5Dialog(Message="Enter hyperfine constant for spin "+str(i)+"("+spins[i]+")")
				a=wiz5.getProperty("A").value
				d=wiz5.getProperty("D").value
				ax=wiz5.getProperty("AxisX").value
				ay=wiz5.getProperty("AxisY").value
				az=wiz5.getProperty("AxisZ").value
				if(d!=0):
					pars["a("+spins[i]+")"]=",".join(map(str,[a,d,ax,ay,az]))
					is_anisotropic=True
				else:
					pars["a("+spins[i]+")"]=str(a)
	elif(hypopt=='All Pairs'):
		for i in range(len(spins)-1):
			for j in range(i+1,len(spins)):
				wiz5=QuantumWizard5Dialog(Message="Enter hyperfine constant between spins "+str(i)+"("+spins[i]+") and "+str(j)+"("+spins[j]+")")
				a=wiz5.getProperty("A").value
				d=wiz5.getProperty("D").value
				ax=wiz5.getProperty("AxisX").value
				ay=wiz5.getProperty("AxisY").value
				az=wiz5.getProperty("AxisZ").value
				if(d!=0):
					pars["a("+spins[i]+","+spins[j]+")"]=",".join(map(str,[a,d,ax,ay,az]))
					is_anisotropic=True
				else:
					pars["a("+spins[i]+","+spins[j]+")"]=str(a)
if(differ_q):
	for s in range(sfc):
		for i in range(len(spins)):
			if(quadmap[i]):
				wiz6=QuantumWizard6Dialog(Message="Enter quadrupole splitting for spin "+str(i)+"("+spins[i]+") in state "+str(s))
				if(wiz6.getProperty('nuQ').value!=0):
					pars["q(@"+str(s)+","+spins[i]+")"]=str(wiz6.getProperty('nuQ').value)+","+str(wiz6.getProperty('AxisX').value)+","+str(wiz6.getProperty('AxisY').value)+","+str(wiz6.getProperty('AxisZ').value)
					is_anisotropic=True
else:
	for i in range(len(spins)):
		if(quadmap[i]):
			wiz6=QuantumWizard6Dialog(Message="Enter quadrupole splitting for spin "+str(i)+"("+spins[i]+")")
			if(wiz6.getProperty('nuQ').value!=0):
				pars["q("+spins[i]+")"]=str(wiz6.getProperty('nuQ').value)+","+str(wiz6.getProperty('AxisX').value)+","+str(wiz6.getProperty('AxisY').value)+","+str(wiz6.getProperty('AxisZ').value)
				is_anisotropic=True
			
if(do_dynamic):
	opts={}
	for i in range(len(spins),6):
		opts["RelaxSpin"+str(i)]=0.0
	for s in range(sfc):
		wiz7=QuantumWizard7Dialog(Message="Enter spin relaxation rates in state "+str(s),**opts)
		for i in range(len(spins)):
			rat=wiz7.getProperty("RelaxSpin"+str(i)).value
			if(rat!=0):
				pars["relax(@"+str(s)+","+spins[i]+")"]=str(rat)

if(do_dynamic and sfc>1):
	if(sfc==2):
		wiz8=QuantumWizard8Dialog()
	else:
		wiz8=QuantumWizard8Dialog(ConversionOption='AllEqual',From_0_to_1=0.0,From_1_to_0=0.0)
	cvopt=wiz8.getProperty('ConversionOption').value
	if(cvopt=='AllEqual'):
		rat=wiz8.getProperty('BetweenAllStates').value
		pars["convert(0-"+str(sfc-1)+",0-"+str(sfc-1)+")"]=str(rat)
		pops=[1.0/sfc]*sfc
	elif(cvopt=='TwoLevels'):
		rat01=wiz8.getProperty('From_0_to_1').value
		rat10=wiz8.getProperty('From_1_to_0').value
		if(rat01==0 and rat10==0):
			pops=[0.5,0.5] # pointless, but...
		elif(rat01==0 and rat10!=0): # one way must start in non-equilibrium level
			pops=[0.0,1.0]
			pars["convert(1,0)"]=str(rat10)
		elif(rat01!=0 and rat10==0):
			pars["convert(0,1)"]=str(rat01)
			pops=[1.0,0.0]
		else:
			pops=[rat10/(rat01+rat10),rat01/(rat01+rat10)] # equilibrium.
			pars["convert(0,1)"]=str(rat01)
			pars["convert(1,0)"]=str(rat10)
	opts={}
	enab=[]
	for i in range(len(pops)):
		opts['PopulationState'+str(i)]=str(pops[i])
		enab.append('PopulationState'+str(i))
	opts['Enable']=",".join(enab)
	for i in range(len(pops),6):
		opts['PopulationState'+str(i)]=0.0
	wiz9=QuantumWizard9Dialog(**opts)
	for i in range(len(spins)):
		pars["pop("+str(i)+")"]=str(wiz9.getProperty('PopulationState'+str(i)).value)
		
if(do_RF):
	wiz10=QuantumWizard10Dialog()
	pars["brf"]=str(wiz10.getProperty('RF_Field_strength_G').value / 10000.0)+","+str(wiz10.getProperty('RF_Frequency_MHz').value)+","+str(wiz10.getProperty('RF_Phase_Deg').value)

if(is_anisotropic):
	opts={"NumberOfAverages":100,"Enable":"NumberOfAverages"}
else:
	opts={"NumberOfAverages":1,"Enable":"NumberOfAverages"}	
wiz11=QuantumWizard11Dialog(**opts)
atype=wiz11.getProperty('Averaging').value
if(atype=='Fixed axis LF'):
	pars["lf"]=str(wiz11.getProperty('AxisX').value)+","+str(wiz11.getProperty('AxisY').value)+","+str(wiz11.getProperty('AxisZ').value)
else:
	pars[atype.lower()]=str(wiz11.getProperty('NumberOfAverages').value)

wiz12=QuantumWizard12Dialog()
pars["bmag"]=str(wiz12.getProperty('MagneticFieldTesla').value)
		
wiz13=QuantumWizard13Dialog()
meas=wiz13.getProperty('MeasurementType').value
if(meas=='Time Spectra'):
	allowedLoops=1
	pars["measure"]="timespectra"
elif(meas=="Integral Asymmetry"):
	allowedLoops=2
	pars["measure"]="integral"

for i in range(allowedLoops):
	defa="Magnetic Field"
	if(allowedLoops==1):
		mess="Enter the loop variable (optional, for Spectra axis)"
	elif(i==0):
		mess="Enter the loop variable for the X (bins) axis"
	else:
		mess="Enter the optional loop variable for the Y (spectra) axis"
		defa="<none>"
	wiz14=QuantumWizard14Dialog(Message=mess,LoopVariable=defa,Enable="LoopVariable")
	lv=wiz14.getProperty("LoopVariable").value
	#print "["+lv+"]"
	if(lv != "<none>"):
		if(lv=="Magnetic Field"):
			pars["loop"+str(i)+"par"]="bmag"
		elif(lv=='Relaxation rate of Electron'):
			pars["loop"+str(i)+"par"]="relax(e)"
		elif(lv=='Conversion rates'):
			keys=list(pars.keys())
			for j in range(len(keys)):
				if(keys[j][0:7]=="convert"):
					pars["loop"+str(i)+"par"]=keys[j]
		elif(lv=='RF Frequency'):
			pars["loop"+str(i)+"par"]="brf[1]"
		elif(lv=='Muon isotropic hyperfine constant'):
			pars["loop"+str(i)+"par"]="a(Mu)[0]"
		start=wiz14.getProperty('StartValue').value
		endv=wiz14.getProperty('EndValue').value
		npts=wiz14.getProperty('NumberOfPoints').value
		logscale=wiz14.getProperty('LogScale').value
		if(logscale):
			npts=-npts
		pars["loop"+str(i)+"range"]=str(start)+","+str(endv)+","+str(npts)
		
# ready to generate table!
#cetw=self.createChildAlgorithm("CreateEmptyTableWorkspace")
#cetw.execute()
#Tab3=cetw.getProperty("OutputWorkspace").value
Table=CreateEmptyTableWorkspace()
Table.addColumn("str","Code")
Table.addColumn("str","Value")

for (cod,valu) in list(pars.items()):
	Table.addRow([cod,valu])

for i in range(20):
	Table.addRow(["",""])

# and GO!!!!

QuantumTableDrivenSimulation(ModelTable="Table",Results="Results")

# plot it
nh=mtd["Results"].getNumberHistograms()
if(nh ==1):
	plotSpectrum("Results",0)
elif(nh <= 5):
	plotSpectrum("Results",list(range(nh)))
else:
	plot2D("results")
	