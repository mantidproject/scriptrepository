from numpy import *
from pylab import *

##Optimum Time calculator
pl=6.1
wavelength=8  #time optimisation wavelength
lmda_min=0       #Range of interest to minimis error
lmda_max=15
Final_Polarisation=0.8
start_time=1  #in terms of spin up
run_length=0.25  #hours
Tup=10
run=8922
#################################################################################################
#CLOSING GRAPHS OUT OF SCRIPT WILL CAUSE A CRASH! DO NOT!
#################################################################################################
##Relable params
nl=0.0733*float(pl)
lmda=float(wavelength)
Pol=float(Final_Polarisation)
a=float(start_time)
r=float(run_length)
##close graphs from previous runs
gui_cmd(close)
gui_cmd(close)
##create time error graph
figure(0)

title(r'Error Dependence on Time')

xlabel('Time / Hours')

ylabel('Relative Error')
ax=gca()
ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
##plot DeltaS0
t=arange(0,r,0.001)
E=abs( (cosh(nl*lmda*Pol*(1-exp(-a/Tup)))+cosh(nl*lmda*Pol*(1-exp(-(a+r)/Tup)))+4*cosh(nl*lmda*Pol*(1-exp(-(2*a+r)/(2*Tup)))))/( 6*cosh(nl*lmda*Pol*(1-exp(-(a+t)/10))) ) -1 )
gui_cmd(plot,t,E)
##find optimum time(point of intersect)
Y=-a-Tup*log(1-(1/(nl*lmda*Pol))*arccosh((cosh(nl*lmda*Pol*(1-exp(-a/Tup)))+cosh(nl*lmda*Pol*(1-exp(-(a+r)/Tup)))+4*cosh(nl*lmda*Pol*(1-exp(-(2*a+r)/(2*Tup)))) )/6 ))

print 'Optimum time for this wavelength is '+str(Y)+' hours'

##create wavelength error graph
x=arange(0,15.05,0.05)
y=abs( (   cosh(nl*x*Pol*(1-exp(-a/Tup) ) )+cosh(nl*x*Pol*(1-exp(-(a+r)/Tup) ) )+4*cosh(nl*x*Pol*(1-exp(-(2*a+r)/(2*Tup)) ) )   ) / (   6*cosh(nl*x*Pol*(1-exp(-(a+Y)/10) ) ) ) -1  )

figure(1)
gui_cmd(plot,x,y,'b-')

title(r'Error Dependence on Wavelength')

xlabel('Wavelength/Angstrom')

ylabel('Relative Error')
ax=gca()
ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y')

##write wavelength error to a data file
Error=open(r'\\Britannic\3he\LET Data\\Error','w')
Error.write('')
Error=open(r'\\Britannic\3he\LET Data\\Error','a')
for i in range(len(x)):
	Error.write(str(x[i])+'	' +str(y[i]) +'\n')
Error.close()
##load 'Error' as workspace and multiply by the monitor to obtain 'AbsError'
LoadAscii(Filename=r'//Britannic/3He/LET Data/Error',OutputWorkspace='Error',Unit='Wavelength')
Load(Filename=r'\\isis\inst$\cycle_12_5\NDXLET\LET0000'+str(run)+'.raw',OutputWorkspace='Monitor',SpectrumList='40966')
ConvertUnits(InputWorkspace='Monitor',OutputWorkspace='Monitor',Target='Wavelength')
Rebin(InputWorkspace='Monitor',OutputWorkspace='Monitor',Params='0,0.05,15',PreserveEvents='0')
Rebin(InputWorkspace='Error',OutputWorkspace='Error',Params='0,0.05,15',PreserveEvents='0')
ConvertToHistogram(InputWorkspace='Error',OutputWorkspace='Error')
ConvertToHistogram(InputWorkspace='Monitor',OutputWorkspace='Monitor')
Multiply(LHSWorkspace='Error',RHSWorkspace='Monitor',OutputWorkspace='AbsError',AllowDifferentNumberSpectra='1')
Integration(InputWorkspace='AbsError',OutputWorkspace='Integral of AbsError',RangeLower=str(lmda_min),RangeUpper=str(lmda_max))

##Open graphs
gui_cmd(show)


