from __future__ import print_function
import offspec_offset2 as nr
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
reload(nr)
nr.current_detector = nr.old_detector

from math import *
import numpy as np
import os,sys,re, time

def getLog(w,log_name):
        # Get handle to the workspace
        try:
          h=mtd[w]
        except:
          print("Can't get Workspace handle")
        #
        # Get access to SampleDetails
	s=h.getSampleDetails().getLogData(log_name).value
	
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
    print(targetname+'.n'+sortedendings[-1])
    return targetname+'.n'+sortedendings[-1]

def loaddata(rnum, path = 'L://RawData/cycle_15_3/',loadcrpt=0):
    try:
        Load(Filename=path+'OFFSPEC000'+str(rnum)+'.nxs', OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoaderVersion=1, LoadMonitors=True)
    except:
        try:
            if loadcrpt == 0:
              updatefile=loadlatest(str(rnum))   
              Load(Filename='z:/'+updatefile, OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoaderVersion=1, LoadMonitors=True)
            else:
              print('trying to load crpt snapshot')
              Load(Filename='z:/snapshot_crpt.nxs', OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoaderVersion=1, LoadMonitors=True)
        except: 
            raise Exception('Could not find data')
    return str(rnum)        

def timeslice(rnum,btime,etime,output,loadcrpt=0):
    loaddata(rnum,loadcrpt=loadcrpt)
    try:
        FilterByTime(InputWorkspace=str(rnum), OutputWorkspace=str(rnum)+'_slice', StartTime=btime,StopTime=etime)
    except:
        raise Exception('Error in slicing')
    Rebin(str(rnum)+'_slice','5.0,20.0,100000.0',PreserveEvents=False,OutputWorkspace=str(rnum)+'_slicereb')
    a1=mtd[str(rnum)]
    gr=a1.getRun()
    tamps=gr.getProtonCharge()
    print('tamps=',str(tamps))
    a2=mtd[str(rnum)+'_slice']
    ua=a2.getRun().getProtonCharge()
    print('ua=',str(ua))
    monnorm=mtd[str(rnum)+'_monitors']*ua/tamps
    Rebin(monnorm,'5.0,20.0,100000.0',OutputWorkspace=str(rnum)+'monreb')
    ConjoinWorkspaces(str(rnum)+'monreb',str(rnum)+'_slicereb',CheckOverlapping=False)
    RenameWorkspace(str(rnum)+'monreb',OutputWorkspace=output+'_'+str(btime)+'-'+str(etime))
    DeleteWorkspace(str(rnum))
    DeleteWorkspace(str(rnum)+'_monitors')
    DeleteWorkspace(str(rnum)+'_slice')
    return output+'_'+str(btime)+'-'+str(etime)


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

def DSqxqz(run1,wksp,angle=1.2,qxqzlimits='-5e-4,5e-4,0.02,0.1',binning1=["1.5","0.02","14.0","2"],Nqx=200,Nqz=200,withpol=1):
    halftheta = angle/2.0
    #binning1=["1.0","0.05","14.0","2"]
    if withpol==1:
        nr.nrPNRFn(run1,wksp,str(halftheta),'none',"114","112","116",binning1,"",'0',['2','1'],'0',dofloodnorm=2)
    else:
        nr.nrNRFn(run1,wksp,str(halftheta),'none',"114","112","116",binning1,"",dofloodnorm=2)
    # Delete the norm and RvQ workspaces as they have the wrong angle
    DeleteWorkspace(wksp+'RvQ')
    DeleteWorkspace(wksp+'norm')
    ConvertSpectrumAxis(InputWorkspace=wksp+"detnorm", OutputWorkspace=wksp+"detnorm", Target='SignedTheta')
    ConvertToReflectometryQ(InputWorkspace=wksp+'detnorm', OverrideIncidentTheta=True, IncidentTheta=angle, Extents=qxqzlimits, OutputAsMDWorkspace=False, OutputWorkspace=wksp+"qxqz", NumberBinsQx=Nqx, NumberBinsQz=Nqz)
    if withpol == 1:
        doDSSCorrections(wksp,angle,1,Nqx=Nqx,Nqz=Nqz)
        doDSSCorrections(wksp,angle,2,Nqx=Nqx,Nqz=Nqz)
    else:
        doDSSCorrections(wksp,angle,0,Nqx=Nqx,Nqz=Nqz)

def _offspecslice_simple(rnum,btime,etime,qmin,qmax,output, binning,theta=0.7, DB="LDDB05k",spec=114,loadcrpt=0):
    wksp =timeslice(rnum,btime,etime,output,loadcrpt=loadcrpt)
    nr.nrNRFn("",wksp,str(theta),DB,spec,"105","122",binning,"",usewkspname=1)
    DeleteWorkspace(wksp)
    DeleteWorkspace(wksp+'norm')

   
def offspecslice2(rnum,qmin,qmax,output,start = 0, tslice=None,nslices = None,sarray=[], theta=0.7, binning=["1.5","0.02","14.0","2"],spec=114,loadcrpt=0):
    slicearray = sarray[:]
    slicenames=[] #this will be a list of workspace names for all the slices created
    datatimes = [] # this will contain the meantime for each dataset
    if tslice or nslices: # if tslice or nslices exist they will take precedence over slicearray
        testws = loaddata(rnum,loadcrpt=loadcrpt)
        runtotaltime = getLog(testws, 'duration')
        print("Total runtime in seconds: " + str(runtotaltime))
        DeleteWorkspace(testws)
        if nslices: 
            tslice = ceil((runtotaltime - start)/(nslices))
        slicearray.append(start)
        while slicearray[-1] < runtotaltime:
            slicearray.append(slicearray[-1]+tslice)
        slicearray[-1] = runtotaltime # lastentry is some random number > than total runduration, set equal to runduration, this means the last slice has a different length to the others
        print("Time boundaries:\n") 
        print(slicearray)
        print("Start making slices:\n")
    for idx in range(len(slicearray)):
        try:
            start = slicearray[idx]; end = slicearray[idx+1]
            datatimes.append(0.5*(start+end)) # calculate the time for this dataset for saving later
            print("\nCreated slice "+str(datatimes[-1]))
            _offspecslice_simple(rnum, start, end, qmin, qmax, output, binning=binning, theta=theta,spec=spec,loadcrpt=loadcrpt)
        except:
            print(datatimes)
            break
       

def offspecPlot(wksp, xrange, yrange, zrange,logscale='z'):
    
    (xmin, xmax) = xrange
    (ymin,ymax) = yrange
    (zmin,zmax) = zrange
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

def offspecQplot(rnum,qmin,qmax,output, nslices=None,sarray = [], angle=0.7,Nqx=50,Nqz=50,qxqzlimits='-2e-4,2e-4,0.01,0.05',  zmin=1e-4, zmax=0.01,qzmin=None,qzmax=None,binning=["1.5","0.02","14.0","2"],spec=114,loadcrpt=0):
    limitlist = qxqzlimits.split(',')
    xmin=float(limitlist[0]); xmax=float(limitlist[1])
    ymin=float(limitlist[2]);ymax=float(limitlist[3])
   
    offspecslice2(rnum,qmin,qmax,'wrong',nslices = nslices,sarray=sarray, theta=angle/2.0,spec=spec,loadcrpt=loadcrpt)

    names=mtd.getObjectNames()
    for name in names:
        m = re.search('^wrong',name)
        if m:
            n = re.search('^wrong{1}(.*)detnorm{1}$', name)
            if n:
                print(name)
                newname = re.sub('wrong',output, name)
                newname = re.sub('detnorm','',newname)
                print("newname: "+newname)
                ConvertSpectrumAxis(InputWorkspace=name, OutputWorkspace=name, Target='SignedTheta')
                try:
                    ConvertToReflectometryQ(InputWorkspace=name, OverrideIncidentTheta=True, IncidentTheta=angle, Extents=qxqzlimits, OutputAsMDWorkspace=False, OutputWorkspace=newname+"qxqz", NumberBinsQx=Nqx, NumberBinsQz=Nqz)
                except:
                    ConvertToReflectometryQ(InputWorkspace=name, OverrideIncidentTheta=True, IncidentTheta=angle, Extents=qxqzlimits, OutputAsMDWorkspace=False, OutputWorkspace=newname+"qxqz", NumberBinsQx=Nqx, NumberBinsQz=Nqz,OutputVertexes='somevertexes')
                        
                CloneWorkspace(name,OutputWorkspace=newname+'detnorm')
                doDSSCorrections(newname,angle,0,Nqx=Nqx,Nqz=Nqz)
                DeleteWorkspace(newname+'detnorm')
                offspecPlot(newname+'qxqz', (xmin, xmax), (ymin,ymax), (zmin,zmax),logscale='z')
                QxQzcuts(newname+'qxqz',qzmin,qzmax)
            DeleteWorkspace(mtd[name])
    offspecslice2(rnum,qmin,qmax,output,nslices = nslices,sarray=sarray, theta=angle,spec=spec,loadcrpt=loadcrpt)
   

    
def chopit(rnum,btime,etime,tslice,output,slicearray=None,usearray=0,sf=1.0, save = True,binning=["1.5","0.02","14.0","2"],loadcrpt=0):
    
    nslice=int((etime*1.0-btime)/(tslice*1.0))
    slicenames=[]
    print('nslice=',str(nslice))
    if usearray==0:
        slicearray=[]
        slicearray.append(btime)
        for i in range(1,nslice+1):
            slicearray.append(btime+(i*tslice))
        if slicearray[-1] < etime:
            slicearray.append(etime)
            nslice=nslice+1
    for i in range(nslice-1):
        btime2=slicearray[i]
        etime2=slicearray[i+1]
        try:
            wksp=timeslice(rnum,btime2,etime2,output,loadcrpt=loadcrpt)
            slicenames.append(wksp)
            print(slicenames)
        except:
            print('time slicing failed')
            break
        nr.nrNRFn("",wksp,"0.700","LDDB05k","114","110","120",binning,"",usewkspname=1,sf=sf)
        #Rebin(wksp+"RvQ","0.011,-0.025,0.09",OutputWorkspace=wksp+"RvQ")
        Rebin(wksp+"RvQ","0.011,-0.02,0.09",OutputWorkspace=wksp+"RvQ")
        DeleteWorkspace(wksp)
        DeleteWorkspace(wksp+'detnorm')
        DeleteWorkspace(wksp+'norm')
    CloneWorkspace(slicenames[0]+'RvQ',OutputWorkspace=output+'_allslices'    )
    for i in range(1,len(slicenames)):
        ConjoinWorkspaces(output+'_allslices',slicenames[i]+'RvQ',CheckOverlapping=0)
    DeleteWorkspace(slicenames[0]+'RvQ')
    datatimes=[]
    for i in range(nslice-1):
        datatimes.append(0.5*(slicearray[i]+slicearray[i+1]))
    writemap_csv(output+'_allslices',datatimes,'C:/everything/userthings/'+output+'/'+output)
    if save:
        print("\n Trying to save the following slices: \n")
        saveslices(output+'_allslices','C:/everything/userthings/'+output+'/')

def saveslices(inputwksp, dir = None):
    if dir:
        userdirectory = dir
    else:
        userdirectory = "C:/everything/userthings/"
    spectrum = 0
    print(spectrum)
    while True:
        try:
            filename = userdirectory + inputwksp + "_" + str(spectrum) + ".dat"
            SaveAscii(inputwksp, filename, SpectrumList = [spectrum], WriteSpectrumID = False, CommentIndicator = "#", Separator = "Tab", ColumnHeader = False)
            spectrum += 1
            print(spectrum)
        except: 
            print("End of slices reached, this one does not exist: " + str(spectrum))
            break


def chopit2(rnum,output,start = 0, tslice=None,nslices = None ,sarray=[],usearray=0,sf=1.0,userdirectory = 'U:/vanWell/April_2015/',binning=["1.5","0.02","14.0","2"],loadcrpt=0):
    
    slicearray = sarray[:]
    slicenames=[] #this will be a list of workspace names for all the slices created
    datatimes = [] # this will contain the meantime for each dataset
    if tslice or nslices: # if tslice or nslices exist they will take precedence over slicearray
        testws = loaddata(rnum,loadcrpt=loadcrpt)
        runtotaltime = getLog(testws, 'duration')
        print("Total runtime in seconds: " + str(runtotaltime))
        DeleteWorkspace(testws)
        if nslices: 
            tslice = ceil((runtotaltime - start)/(nslices))
        slicearray.append(start)
        while slicearray[-1] < runtotaltime:
            slicearray.append(slicearray[-1]+tslice)
        slicearray[-1] = runtotaltime # lastentry is some random number > than total runduration, set equal to runduration, this means the last slice has a different length to the others
        print("Time boundaries:\n") 
        print(slicearray)
        print("Start making slices:\n")
    for idx in range(len(slicearray)):
        try:
            start = slicearray[idx]; end = slicearray[idx+1]
            datatimes.append(0.5*(start+end)) # calculate the time for this dataset for saving later
            wksp = timeslice(rnum,start,end,output,loadcrpt=loadcrpt) 
        except:
            break
        slicenames.append(wksp)    
        nr.nrNRFn("",wksp,"0.7","LDDB05k","114","110","118",binning,"",usewkspname=1,sf=sf)
        Rebin(wksp+"RvQ","0.011,-0.01,0.09",OutputWorkspace=wksp+"RvQ")
        DeleteWorkspace(wksp)
        DeleteWorkspace(wksp+'detnorm')
        DeleteWorkspace(wksp+'norm')
        if idx == 0:
            CloneWorkspace(slicenames[0]+'RvQ',OutputWorkspace=output+'_allslices')
        else:
            ConjoinWorkspaces(output+'_allslices',slicenames[-1]+'RvQ',CheckOverlapping=0)
    DeleteWorkspace(slicenames[0]+'RvQ')  
    writemap_csv(output+'_allslices',datatimes,userdirectory + output + '/'+output)
    print("\n Trying to save the following slices: \n")
    saveslices(output+'_allslices', userdirectory+ output + '/')


    
def offspecslice(rnum,btime,etime,qmin,qmax,output,spec=114,loadcrpt=0):
    wksp=timeslice(rnum,btime,etime,output,loadcrpt=loadcrpt)
    nr.nrNRFn("",wksp,"0.7","LDDB05k",spec,"110","118",binning,"",usewkspname=1)
    ConvertUnits(wksp+'detnorm','MomentumTransfer',OutputWorkspace=wksp+'detnormQ')
    Rebin(wksp+'detnormQ','0.011,-0.01,0.09',OutputWorkspace=wksp+'detnormQ')
    Integration(wksp+'detnormQ',qmin,qmax,OutputWorkspace=wksp+'detnormQ_Int')
    Transpose(wksp+'detnormQ_Int',OutputWorkspace=wksp+'detnormQ_'+str(qmin)+'_'+str(qmax))
