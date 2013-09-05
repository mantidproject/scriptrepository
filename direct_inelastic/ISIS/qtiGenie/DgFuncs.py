from mantid.simpleapi import *
from numpy import *
try:
	from scipy import signal
except:
	print 'Scipy module not in scope'
	
#uses the new api all input workpsaces are expected to have |Q| along X and Energy Transfer along Y
def gofe(wkspin,T,dbwfac):
	dat=mtd[wkspin]
	print 'correcting workspace' , dat.name(),' to be the pdos with T=',T,'K and |U^2| as ',dbwfac 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
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

def byE(wkspin):
	dat=mtd[wkspin]
	print 'Flattening workspace' , dat.name() 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	qq=axis.extractValues()
	q=(qq[1:len(qq)]+qq[0:len(qq)-1])/2

	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2

	for i in range(0,shape(y)[1]):
		
		y[:,i]=y[:,i]*energy
		Err[:,i]=Err[:,i]*energy

	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='byE',WorkSpaceTitle='Flattened')
	return wkspOut
	
def gaussSmooth(wkspin,n):
	ReplaceSpecialValues(InputWorkspace=wkspin,OutputWorkspace=wkspin,NaNValue='0',InfinityValue='0')
	dat=mtd[wkspin]
	print 'Smoothing workspace' , dat.name(), 'by Gaussian convolution ',n,' times'
	
	x=dat.extractX()	
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	qq=axis.extractValues()
	q=(qq[1:len(qq)]+qq[0:len(qq)-1])/2

	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	
	y=blur_image(y,n)
	try:	
		wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='byE',WorkSpaceTitle='smoothed '+str(n)+'times')
	except:

		wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en),VerticalAxisUnit='DeltaE',VerticalAxisValues=en,UnitX='|Q|',YUnitLabel='byE',WorkSpaceTitle='smoothed '+str(n)+'times')

	return wkspOut
	
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='same')
    return(improc)

def bose(wkspin,T):
	dat=mtd[wkspin]
	print 'correcting workspace ' , dat.name(),' by the Bose Factor with T=',T,'K' 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	
	bosegrid=ones_like(y)
	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	#bose=1-exp(-energy*11.604/(T))
	bose=exp(-energy*11.604/(T))
	
	for i in range(0,shape(y)[1]):
		
		bosegrid[:,i]=bose
	y=y*bosegrid
	Err=Err*bosegrid
	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='T corrected',WorkSpaceTitle='Bose factor corrected')
	return wkspOut
	
def DetailBalance(wkspin,T):
	dat=mtd[wkspin]
	print 'correcting workspace for DB' , dat.name(),' with T=',T,'K' 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	
	bosegrid=ones_like(y)
	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	bose=exp(-energy*11.604/(T))
	
	for i in range(0,shape(y)[1]):
		
		bosegrid[:,i]=bose
	y=y*bosegrid
	Err=Err*bosegrid
	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='T corrected',WorkSpaceTitle='Bose factor corrected')
	return wkspOut

def perCm(wkspin):
	dat=mtd[wkspin]
	print 'converting energy units to per cm' 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	
	
	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	energy=energy*8.0654
	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='cm^-1')
	return wkspOut
def meV(wkspin):
	dat=mtd[wkspin]
	print 'converting energy units to meV' 
	x=dat.extractX()
		
	y=dat.extractY()
	Err=dat.extractE()
	axis=dat.getAxis(0)
	
	
	#deal with the energy axis
	axis=dat.getAxis(1)
	en=axis.extractValues()
	energy=(en[1:len(en)]+en[0:len(en)-1])/2
	energy=energy*0.124
	wkspOut=CreateWorkspace(x,y,Err,Nspec=len(en)-1,VerticalAxisUnit='DeltaE',VerticalAxisValues=energy,UnitX='|Q|',YUnitLabel='meV')
	return wkspOut
#ww=gofe('w2',1,1)
