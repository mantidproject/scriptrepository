import offspec_red as nr
reload(nr)
from math import *
import numpy as np
import os,sys,re, time
from mantid.simpleapi import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace 
import matplotlib.pyplot as plt
try:
  from mantidplot import *
except ImportError:
  pass
  

def getLog(w,log_name):
        # Get handle to the workspace
        try:
          h=mtd[w]
        except:
          print "Can't get Workspace handle"
        #
        # Get access to SampleDetails
	s=h.getRun().getLogData(log_name).value
	
	return s


def writemap_csv(wksp,times,fname):
	dir=os.path.dirname(fname+'Qscale.csv')
	try:
	  os.stat(dir)
	except:
	  os.mkdir(dir)
	f=open(fname+'Qscale.csv','w')
	w1=mtd[wksp]
	xarray=w1.readX(0)
	npts=len(xarray)
	nhist=w1.getNumberHistograms()
	x1=np.zeros(npts-1)
	for i in range(npts-1):
		x1[i]=(xarray[i]+xarray[i+1])/2.0
	s=""
	for i in range(npts-2):
		s+="%g,"  %  (x1[i])
	s+="%g\n" % (x1[npts-2])
	f.write(s)
	f.close()
	f=open(fname+'timeScale.csv','w')
	s=""
	for i in range(len(times)-1):
		s+="%g,"  %  (times[i])
	s+="%g\n" % (times[len(times)-1])
	f.write(s)
	f.close()
	f=open(fname+'ZData.csv','w')
	s=""
	for i in range(nhist):
		yarray=w1.readY(i)
		s=""
		for j in range(npts-2):
			s+="%g,"  %  (yarray[j])
		s+="%g\n" % (yarray[npts-2])
		f.write(s)
	f.close()
	f=open(fname+'Errors.csv','w')
	for i in range(nhist):
		earray=w1.readE(i)
		s=""
		for j in range(npts-2):
			s+="%g,"  %  (earray[j])
		s+="%g\n" % (earray[npts-2])
		f.write(s)
	f.close()

def writeXYE(wksp,fname):
	a1=Mantid.getMatrixWorkspace(wksp)
	x1=a1.readX(0)
	X1=n.zeros((len(x1)-1))
	for i in range(0,len(x1)-1):
		X1[i]=(x1[i]+x1[i+1])/2.0
	y1=a1.readY(0)
	e1=a1.readE(0)
	f=open(fname,'w')
	for i in range(len(X1)):
		s=""
		s+="%g " % X1[i]
		s+="%g " % y1[i]
		s+="%g\n" % e1[i]  
		f.write(s)
	f.close()

def loadlatest(currentrun=None):
    path = "Z:/"
    runfiles = os.listdir(path)
    targetname = "OFFSPEC000"+currentrun
    endings = []
    for file in runfiles:
        filelist = file.split('.n',1)
        if targetname in filelist[0]:
            try:
                endings.append(filelist[1])
            except: pass   
    sortedendings = sorted(endings)
    print targetname+'.n'+sortedendings[-1]
    return targetname+'.n'+sortedendings[-1]

def loaddata(rnum, path = '//ndloffspec1/L/RawData/',loadcrpt=0):
    try:
        Load(str(rnum), OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
    except:
        #try:
            if loadcrpt == 0:
              updatefile=loadlatest(str(rnum))   
              Load(Filename='z:/'+updatefile, OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
            else:
              print 'trying to load crpt snapshot'
              Load('z:/snapshot_crpt.nxs', OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
        #except: 
        #    raise Exception('Could not find data')
    return str(rnum)        

def doDSSCorrections(wksp,angle=1.2,nper=0,Nqx=200,Nqz=200):
    # get the lambda and theta arrays from the original data
    thetaf=[]
    if nper == 0:
        suffix1='detnorm'
        suffix2='qxqz'
        suffix3='qxlam'
    else:
        suffix1='detnorm_'+str(nper)
        suffix2='qxqz_'+str(nper)
        suffix3='qxlam_'+str(nper)
    a1=mtd[wksp+suffix1]
    nth=a1.getNumberHistograms()
    ntc=len(a1.dataY(0))
    thetaf=a1.getAxis(1).extractValues()
    thetaf=thetaf*pi/180.0
    lambda0=a1.getAxis(0).extractValues()
    lambda1=[]
    for i in range(len(lambda0)-1):
        lambda1.append(0.5*(lambda0[i]+lambda0[i+1]))
    dthf=float(nth-1)/(thetaf[-1]-thetaf[0])
    dlam=float(ntc-1)/(lambda1[-1]-lambda1[0])
    # get the qx and qz arrays from the data we just created
    a2=mtd[wksp+suffix2]
    lmin=lambda0[0]
    lmax=lambda0[-1]
    lamstep=(lmax-lmin)/(Nqz-1)
    lam2=[]
    for i in range(Nqz):
        lam2.append(lmin+i*lamstep)
    qz=a2.getAxis(1).extractValues()
    qx=a2.getAxis(0).extractValues()
    cthetai=cos(angle*pi/180.0)
    sthetai=sin(angle*pi/180.0)
    thetai=angle*pi/180.0
    thetaf0=thetaf[0]
    lambda0=lambda1[0]
    for i in range(Nqz):
        qzi=qz[i]
        #qzi=lam2[i]
        pi2qz=2.0*pi/qzi
        for j in range(Nqx):
            qxj=qx[j]
            ang=(qzi*cthetai-qxj*sthetai)/sqrt(qzi*qzi+qxj*qxj)
            #ang=cthetai-(qxj*qzi/(2.0*pi))
            ang=min(1.0,ang)
            ang=asin(ang)
            #ang=acos(ang)
            if qxj==0.0:
                ang=thetai
            else:
                ang=pi-ang-atan(qzi/qxj)
                if ang > pi:
                    ang=ang-pi
            lam=pi2qz*(sthetai+sin(ang))
            #lam=qzi
            xind=(ang-thetaf0)*dthf
            yind=(lam-lambda0)*dlam
            indy=int(yind)
            indx=int(xind)
            if indy >=  0 and indy <= ntc-2 and indx >= 0 and indx <= nth-2:
                dyind=yind-float(indy)
                dxind=xind-float(indx)
                ofsp00=a1.dataY(indx)[indy]
                ofsp01=a1.dataY(indx)[indy+1]
                ofsp10=a1.dataY(indx+1)[indy]
                ofsp11=a1.dataY(indx+1)[indy+1]
                offsp1=(1.0-dxind)*ofsp00+dxind*ofsp10
                offsp2=(1.0-dxind)*ofsp01+dxind*ofsp11
                a2.dataY(i)[j]=(1.0-dyind)*offsp1+dyind*offsp2
                ofsp00=a1.dataE(indx)[indy]
                ofsp00=ofsp00*ofsp00
                ofsp01=a1.dataE(indx)[indy+1]
                ofsp01=ofsp01*ofsp01
                ofsp10=a1.dataE(indx+1)[indy]
                ofsp10=ofsp10*ofsp10
                ofsp11=a1.dataE(indx+1)[indy+1]
                ofsp11=ofsp11*ofsp11
                offsp1=((1.0-dxind)*(1.0-dxind))*ofsp00+(dxind*dxind)*ofsp10
                offsp2=((1.0-dxind)*(1.0-dxind))*ofsp01+(dxind*dxind)*ofsp11
                a2.dataE(i)[j]=sqrt(abs((1.0-dyind)*(1.0-dyind)*offsp1+dyind*dyind*offsp2))
            else:
                a2.dataY(i)[j]=0.0
                a2.dataE(i)[j]=0.0
    w1=mtd[wksp+suffix2]*1.0
    w2=mtd[wksp+suffix2]
    a2=w1.getAxis(1)

def DSqxqz(run,wksp,angle=1.2,qxqzlimits='-5e-4,5e-4,0.02,0.1',binning=["1.5","0.02","14.0","2"],Nqx=200,Nqz=200,withpol=False):
    halftheta = angle/2.0
    if withpol=='pnr':
        nr.nrPNRFn(run,wksp,str(halftheta),'none',nr.current_detector.specular,nr.current_detector.specular,nr.current_detector.specular,binning,"",PNRwithPA=False,pnums=['2','1'],doCorrs=False,dofloodnorm=True)
        periods = ("_1", "_2")
    elif withpol=='pa':
        nr.nrPNRFn(run,wksp,str(halftheta),'none',nr.current_detector.specular,nr.current_detector.specular,nr.current_detector.specular,binning,floodfile="",PNRwithPA=True,pnums =["1","2","3","4"],doCorrs=False)
        periods = ("_1", "_2", "_3", "_4")
    else:
        nr.nrNRFn(run,wksp,str(halftheta),'none',nr.current_detector.specular,nr.current_detector.specular,snr.current_detector.specular,binning,"",dofloodnorm=True)
        periods = ("",)
    # Delete the norm and RvQ workspaces as they have the wrong angle
    DeleteWorkspace(wksp+'RvQ')
    DeleteWorkspace(wksp+'norm')
    #CloneWorkspace(wksp+"detnorm", OutputWorkspace = wksp+"det")
    ConvertSpectrumAxis(InputWorkspace=wksp+"detnorm", OutputWorkspace=wksp+"detnorm", Target='SignedTheta')
    wksplist = []
    for period in periods:
        ConvertToReflectometryQ(InputWorkspace=wksp+'detnorm'+period, OverrideIncidentTheta=True, IncidentTheta=angle, Extents=qxqzlimits, OutputAsMDWorkspace=False, OutputWorkspace=wksp+"qxqz"+period, NumberBinsQx=Nqx, NumberBinsQz=Nqz, OutputVertexes="somevertices")
        wksplist.append(wksp+'qxqz'+period)
        
    if withpol == 'pnr' or withpol == 'pa':
        doDSSCorrections(wksp,angle,1,Nqx=Nqx,Nqz=Nqz)
        doDSSCorrections(wksp,angle,2,Nqx=Nqx,Nqz=Nqz)
        try:
            doDSSCorrections(wksp,angle,3,Nqx=Nqx,Nqz=Nqz)
            doDSSCorrections(wksp,angle,4,Nqx=Nqx,Nqz=Nqz)
        except: pass
        GroupWorkspaces(wksplist, Outputworkspace = wksp+"qxqz")
    else:
        doDSSCorrections(wksp,angle,0,Nqx=Nqx,Nqz=Nqz)
        
      
def QxQzcuts(name, qzmin=None, qzmax=None,plot=False):
    if qzmin:
        outputname = name+'_cut_'+str(qzmin)+'-'+str(qzmax)
    else:
        outputname = name+'_cut_all'
    Transpose(InputWorkspace=name, OutputWorkspace=outputname)
    Integration(InputWorkspace=outputname, RangeLower=qzmin, RangeUpper=qzmax, OutputWorkspace=outputname)
    Transpose(InputWorkspace=outputname, OutputWorkspace=outputname)
   
    if plot:
        plot(outputname,0, tool='plot_spectrum', error_bars=True)
        yscale('log')


def offspecPlot(wksp, (xmin, xmax), (ymin,ymax), (zmin,zmax),logscale='z'):
    
    p = plot2D(wksp)
    l=p.activeLayer()
    l.setScale(0, ymin, ymax)
    l.setScale(1, zmin, zmax)
    l.setScale(2, xmin, xmax)
    
    if 'z' in logscale:
        l.setAxisScale(1, zmin, zmax, Layer.Log10)
    elif 'x' in logscale:
        l.setAxisScale(2, xmin, xmax, Layer.Log10)
    elif 'y' in logscale:   
        l.setAxisScale(0, ymin, ymax, Layer.Log10)
              
