from mantid.simpleapi import *
from numpy import *
try:
	from scipy import signal
except:
	print 'Scipy module not in scope'
	
#uses the new api all input workpsaces are expected to have |Q| along X and Energy Transfer along Y
def yscale(wkspin,mass):
	#a routine that converts mari data in energy transfer to Y in momentum transfer.
	
	dat=mtd[wkspin]
	
x=dat.extractX()
	
y=dat.extractY()
Err=dat.extractE()
	
	
ei=data.efixed
M=4.0
thimax=2.44
q2=(2.0*ei-2.0*np.cos(thimax)*((ei*ei)**0.5))/2.072


#get q grid the energy grid is the x data histograms from the workspace
axis=dat.getAxis(1)
qq=axis.extractValues()
q=(qq[1:len(qq)]+qq[0:len(qq)-1])/2


#for each value of Q renormalise the energy grid to Y =(M/Q)(w-wr) :wr=Q^2/2M

for i in range(len(x)):
	wr= (q[i]**2)/2*M
	x[i]= (M/q[i]) * (x[i]-wr)


e1=2.072.*q1.^2./m;
energyGrid=(energyGrid-e2)./((2.072.*q2)./(m/2))



	
	
	
	
axis=dat.getAxis(0)
qq=axis.extractValues()
	q=(qq[1:len(qq)]+qq[0:len(qq)-1])/2
	
	qgrid=ones_like(y)
	bosegrid=ones_like(y)
	for i in range(0,shape(y)[0]):
		
		qgrid[i,:]=q
	
	Q2=qgrid**2
	DW=(Q2*dbwfac)
	factor2=exp((-2*DW));
	
	
	y=y/Q2
	Err=Err/Q2
	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	
	bose=1-exp(-energy*11.604/(T))
	
	for i in range(0,shape(y)[1]):
		
		y[:,i]=y[:,i]*energy
		Err[:,i]=Err[:,i]*energy
		bosegrid[:,i]=bose
	
	y=y*(bosegrid/factor2)	
	Err=Err*(bosegrid/factor2)	
	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='pdos',WorkSpaceTitle='Density of States')
	return wkspOut