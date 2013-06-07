#8922 DEPOLARISED

pl=6.1  #enter cell pressure length

#Choose which runs to look at

#X
#AFPDOWN=8912
#AFPUP=8913

#Y
#AFPUP=8914
#AFPDOWN=8915

#Z
#AFPUP=8910
#AFPDOWN=8911

##########################################################

nl=0.0733*pl

##Load Depolarised Run in units lamda
Load(Filename=r'\\britannic\3He\LET Data\0213\LET00008922.raw',OutputWorkspace='DepolarisedRun',SpectrumList='40966')
NormaliseByCurrent(InputWorkspace='DepolarisedRun',OutputWorkspace='DepolarisedRun')
ConvertUnits(InputWorkspace='DepolarisedRun',OutputWorkspace='DepolarisedRun',Target='Wavelength')

##Load UP Polarised Run
Load(Filename=r'\\britannic\3He\LET Data\0213\LET0000'+str(AFPUP)+'.raw',OutputWorkspace='polup',SpectrumList='40966')
NormaliseByCurrent(InputWorkspace='polup',OutputWorkspace='polup')
ConvertUnits(InputWorkspace='polup',OutputWorkspace='polup_lam',Target='Wavelength')

##Pn=sqrt(1-(Tn0/Tn)^2) for up data
Divide(LHSWorkspace='polup_lam',RHSWorkspace='polup_lam',OutputWorkspace='one')
Divide(LHSWorkspace='DepolarisedRun',RHSWorkspace='polup_lam',OutputWorkspace='Tratio')
Power(InputWorkspace='Tratio',OutputWorkspace='Tratio2',Exponent='2')
Minus(LHSWorkspace='one',RHSWorkspace='Tratio2',OutputWorkspace='Pn UP')
Power(InputWorkspace='Pn UP',OutputWorkspace='Pn UP',Exponent='0.5')

##Load DOWN Polarised Run
Load(Filename=r'\\britannic\3He\LET Data\0213\LET0000'+str(AFPDOWN)+'.raw',OutputWorkspace='poldown',SpectrumList='40966')
NormaliseByCurrent(InputWorkspace='poldown',OutputWorkspace='poldown')
ConvertUnits(InputWorkspace='poldown',OutputWorkspace='poldown_lam',Target='Wavelength')

##Pn=sqrt(1-(Tn0/Tn)^2) for down data
Divide(LHSWorkspace='poldown_lam',RHSWorkspace='poldown_lam',OutputWorkspace='one')
Divide(LHSWorkspace='DepolarisedRun',RHSWorkspace='poldown_lam',OutputWorkspace='Tratio')
Power(InputWorkspace='Tratio',OutputWorkspace='Tratio^2',Exponent='2')
Minus(LHSWorkspace='one',RHSWorkspace='Tratio2',OutputWorkspace='Pn DOWN')
Power(InputWorkspace='Pn DOWN',OutputWorkspace='Pn DOWN',Exponent='0.5')

##Find 3He Polarisation
Fit(Function='name=UserFunction,Formula=tanh('+str(nl)+'*x*PHe),PHe=0.5,constraints=(1.00>=PHe>=0.00)',InputWorkspace='Pn UP'  ,Output='Pn UP',StartX='1.5',EndX='12.5')
Fit(Function='name=UserFunction,Formula=tanh('+str(nl)+'*x*PHe),PHe=0.5,constraints=(1.00>=PHe>=0.00)',InputWorkspace='Pn DOWN'  ,Output='Pn DOWN',StartX='1.5',EndX='12.5')
