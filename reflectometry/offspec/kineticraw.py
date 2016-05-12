import offspec_offset2 as nr
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
          print "Can't get Workspace handle"
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
    print targetname+'.n'+sortedendings[-1]
    return targetname+'.n'+sortedendings[-1]

def loaddata(rnum, path = 'u://',loadcrpt=0):
    try:
        Load(Filename=path+'OFFSPEC000'+str(rnum)+'.nxs', OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoaderVersion=1, LoadMonitors=True)
    except:
        try:
            if loadcrpt == 0:
              updatefile=loadlatest(str(rnum))   
              Load(Filename='z:/'+updatefile, OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoaderVersion=1, LoadMonitors=True)
            else:
              print 'trying to load crpt snapshot'
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
    print 'tamps=',str(tamps)
    a2=mtd[str(rnum)+'_slice']
    ua=a2.getRun().getProtonCharge()
    print 'ua=',str(ua)
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
        print "Total runtime in seconds: " + str(runtotaltime)
        DeleteWorkspace(testws)
        if nslices: 
            tslice = ceil((runtotaltime - start)/(nslices))
        slicearray.append(start)
        while slicearray[-1] < runtotaltime:
            slicearray.append(slicearray[-1]+tslice)
        slicearray[-1] = runtotaltime # lastentry is some random number > than total runduration, set equal to runduration, this means the last slice has a different length to the others
        print "Time boundaries:\n" 
        print slicearray
        print "Start making slices:\n"
    for idx in range(len(slicearray)):
        try:
            start = slicearray[idx]; end = slicearray[idx+1]
            datatimes.append(0.5*(start+end)) # calculate the time for this dataset for saving later
            print "\nCreated slice "+str(datatimes[-1])
            _offspecslice_simple(rnum, start, end, qmin, qmax, output, binning=binning, theta=theta,spec=spec,loadcrpt=loadcrpt)
        except:
            print datatimes
            break
       

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
                print name
                newname = re.sub('wrong',output, name)
                newname = re.sub('detnorm','',newname)
                print "newname: "+newname
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
    print 'nslice=',str(nslice)
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
            print slicenames
        except:
            print 'time slicing failed'
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
        print "\n Trying to save the following slices: \n"
        saveslices(output+'_allslices','C:/everything/userthings/'+output+'/')

def saveslices(inputwksp, dir = None):
    if dir:
        userdirectory = dir
    else:
        userdirectory = "C:/everything/userthings/"
    spectrum = 0
    print spectrum
    while True:
        try:
            filename = userdirectory + inputwksp + "_" + str(spectrum) + ".dat"
            SaveAscii(inputwksp, filename, SpectrumList = [spectrum], WriteSpectrumID = False, CommentIndicator = "#", Separator = "Tab", ColumnHeader = False)
            spectrum += 1
            print spectrum
        except: 
            print "End of slices reached, this one does not exist: " + str(spectrum)
            break
     

def slice_the_data(rnum,output,start = 0, tslice=None,nslices = None ,sarray=[],usearray=0,sf=1.0,userdirectory = 'U://',binning=["1.5","0.02","14.0","2"],loadcrpt=0):
    
    slicearray = sarray[:]
    slicenames=[] #this will be a list of workspace names for all the slices created
    datatimes = [] # this will contain the meantime for each dataset
    if tslice or nslices: # if tslice or nslices exist they will take precedence over slicearray
        testws = loaddata(rnum,loadcrpt=loadcrpt)
        runtotaltime = getLog(testws, 'duration')
        print "Total runtime in seconds: " + str(runtotaltime)
        DeleteWorkspace(testws)
        if nslices: 
            tslice = ceil((runtotaltime - start)/(nslices))
        slicearray.append(start)
        while slicearray[-1] < runtotaltime:
            slicearray.append(slicearray[-1]+tslice)
        slicearray[-1] = runtotaltime # lastentry is some random number > than total runduration, set equal to runduration, this means the last slice has a different length to the others
        print "Time boundaries:\n" 
        print slicearray
        print "Start making slices:\n"
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
    print "\n Trying to save the following slices: \n"
    saveslices(output+'_allslices', userdirectory+ output + '/')


    
def offspecslice(rnum,btime,etime,qmin,qmax,output,spec=114,loadcrpt=0):
    wksp=timeslice(rnum,btime,etime,output,loadcrpt=loadcrpt)
    nr.nrNRFn("",wksp,"0.7","LDDB05k",spec,"110","118",binning,"",usewkspname=1)
    ConvertUnits(wksp+'detnorm','MomentumTransfer',OutputWorkspace=wksp+'detnormQ')
    Rebin(wksp+'detnormQ','0.011,-0.01,0.09',OutputWorkspace=wksp+'detnormQ')
    Integration(wksp+'detnormQ',qmin,qmax,OutputWorkspace=wksp+'detnormQ_Int')
    Transpose(wksp+'detnormQ_Int',OutputWorkspace=wksp+'detnormQ_'+str(qmin)+'_'+str(qmax))
    
binning=["1.5","0.02","14.0","2"]
combine_binning=["0.0085","-0.015","0.3"]

times = [i for i in range(0,19800,1800)]
times2= [i for i in range(25200,36000,7200)]
timearray = times+times2+[36644]
#[0, 1800, , 900, 1200, 1500,1800,2100,2400,2700,3000,3300,3600,7200,10800,14400,18000, 36644.0]
#offspecQplot('35723',0.01,0.06,'test2',nslices=20,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

QxQzcuts('test2_3665.0-7330.0qxqz', qzmin=0.03, qzmax=0.033)

nr.nrDBFn("34193+34195","w93","34194","w94","LDDB05","108","120","6.0",binning,"",fitspline=10,diagnostics="0")

nr.nrDBFn("34220+34222","w20","34221","w21","LDDB05piezo","108","120","8.0",binning,"",fitspline=10,diagnostics="0")

saveslices("Mg1_loading1_1160mbar_chopit_60sec_allslices")


#Old Direct Beam
#nr.nrDBFn("35737","w37","35738","w38","LDDB05k","108","120","4.5",binning,"",fitspline=10,diagnostics="0")
#nr.nrDBFn("35739","w39","35740","w40","LDDB05s","108","120","5.3",binning,"",fitspline=10,diagnostics="0")

#New Direct Beam
nr.nrDBFn("35816+35818+35820+35822+35824+35826+35828+35830","w37","35817+35819+35821+35823+35825+35827+35829","w38","LDDB05k","108","120","10.0",binning,"",fitspline=10,diagnostics="0")
nr.nrDBFn("35785+35787+35789+35791+35793+35795+35797+35799+35801+35803+35805+35807+35809+35811+35813+35815","w39","35786+35788+35790+35792+35794+35796+35798+35800+35802+35804+35806+35808+35810+35812+35814","w40","LDDB05s","108","120","8.0",binning,"",fitspline=10,diagnostics="0")
#nr.nrDBFn("35816","w37","35817","w38","LDDB05k","108","120","4.5",binning,"",fitspline=10,diagnostics="0")



#Copy and paste the list called 'j' here from the output bar below,
j=[0, 3601.0, 7202.0, 10803.0, 14404.0, 18005.0, 21606.0, 25207.0, 28808.0, 32409.0, 36010.0, 39611.0, 43212.0, 46813.0, 50414.0, 54015.0, 57616.0, 61217.0, 64818.0, 68419.0, 72019.0]
qzmin = 0.022
qzmax=0.035
for i in range(len(j)-1):
    
    QxQzcuts('test2_'+str(j[i])+'-'+str(j[i+1])+'qxqz', qzmin=qzmin, qzmax=qzmax,plot=False)
    SaveAscii('test2_'+str(j[i])+'-'+str(j[i+1])+'qxqz'+'_cut_'+str(qzmin)+'-'+str(qzmax), 'U:/VanWell/July_2015/Cuts/test2_'+str(j[i])+'-'+str(j[i+1])+'qxqz'+'_cut_'+str(qzmin)+'-'+str(qzmax)+'.dat', WriteSpectrumID = False, CommentIndicator = "#", Separator = "Tab", ColumnHeader = False)



print os.getcwd()
#
#=====================================================================================================================================================
#

###########################
#Mg-1
###########################

#Virgin Sample in Air 
#Start Time: 14/07/2015 14:53
#Folder Name Pictures: ISIS/July2015/ NA

nr.nrNRFn("35715","Mg1_VirginAir_1000mbar_030C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35716","Mg1_VirginAir_1000mbar_030C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_VirginAir_1000mbar_030C_th=0.5RvQ,Mg1_VirginAir_1000mbar_030C_th=2.0RvQ","Mg1_VirginAir_1000mbar_030C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#Run 35717 Not usefull.
#Increased temperature to T=80 C in about 5 min before start of the run.
#    <T<81.6    <T_heater<91 C
#Some dirt on the O-ring of the top  part of the caused a failed attempt to vacuum pump the sample. 

#Virgin Sample in Vacuum 
#Start Time: 14/07/2015 15:22
#Folder Name Pictures: ISIS/July2015/ NA

nr.nrNRFn("35718","Mg1_Virgin_0000mbar_080C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35719","Mg1_Virgin_0000mbar_080C_th=1.7","1.7","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_Virgin_0000mbar_080C_th=0.5RvQ,Mg1_Virgin_0000mbar_080C_th=1.7RvQ","Mg1_Virgin_0000mbar_080C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#Virgin Sample in Vacuum: Kinetic run
#Start Time: 14/07/2015 18:03
#Folder Name Pictures: ISIS/July2015/ NA
chopit2(35720,'Mg1_Virgin_0000mbar_080C_th0.7_600', nslices=1) 
#chopit(35720,0,601,60,'Mg1_Virgin_0000mbar_080C_th0.7_60') 
#noticed from plots that graphs are different from static measurements

loaddata(35720)


#Scanned realized theta, th=0.707
#Start Time: 14/07/2015 18:50
#Folder Name Pictures: ISIS/July2015/ NA
chopit(35721,0,601,600,'Mg1_Virgin_0000mbar_080C_th0.7_600') 

#offset in time considering the pictures: 2minutes earlier than Time at Offspec
#Start Time: 14/07/2015 19:39:50
#Stop Time: 15/07/2015 05:50:50
#Total length 10:10:14
#Folder Name Pictures: ISIS/July2015/Mg1/loading1_1200mbar_080C
chopit(35722,120,1321,60,'Mg1_loading1_1200mbar_080C_60s') 
chopit(35722,120,3721,300,'Mg1_loading1_1200mbar_080C_300s') 
chopit(35722,120,7321,600,'Mg1_loading1_1200mbar_080C_600s') 
chopit(35722,120,36600,1800,'Mg1_loading1_1200mbar_080C_1800s') 
chopit(35722,120,36600,3600,'Mg1_loading1_1200mbar_080C_3600s') 


#2:00 start to increase the pressure (97230 at Project_X pressure software) Flow=100 sscm, Vout=0V
#2:30 P=150 mbar
#3:00 P=250
#4:00 P=400
#4:45 P=500
#5:10 P=600
#6:49 P=800
#7:50 P=900
#9:06 P=1000
#11:21 P=1100 flow to 10 sccm
#14:30 P=1185
#14:55 P=1200
#16:18 P=1220
#19:46 P=1200
#24:55 P=1150 mbar


#Change of Temperature to 70 degrees.


#Start Time: 15/07/2015 05:50:53
#Stop Time: 16/07/2015 01:51:12
#Total length:20:00:19
#Folder Name Pictures: ISIS/July2015/Mg1/loading1_1200mbar_080C (Unchanged as compared with T=80 C)
chopit(35723,0,1201,60,'Mg1_loading1_1200mbar_070C_60s') 
chopit(35723,0,3601,600,'Mg1_loading1_1200mbar_070C_600s') 
chopit(35723,0,72001,1800,'Mg1_loading1_1200mbar_070C_1800s') 
chopit(35723,0,72001,3600,'Mg1_loading1_1200mbar_070C_3600s') 

offspecslice(35723,120,3701,0.015,0.04,"Mg1_loading1_1200mbar_070C_th0.7_1800_offspec",spec=116.5)

#Static Measurement Mg1_loaded1_1200mbar_070C
#Start Time: 16/07/2015 02:00:00
#Stop Time: 16/07/2015 ????
#Longer Since beam was down for about 16 minutes.
#Folder Name Pictures: NA
nr.nrNRFn("35724","Mg1_loaded1_1200mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35725","Mg1_loaded1_1200mbar_070C_th=1.7","1.7","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_loaded1_1200mbar_070C_th=0.5RvQ,Mg1_loaded1_1200mbar_070C_th=1.7RvQ","Mg1_loaded1_1200mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#Unloading1
#Start Time: 16/07/2015 03:40:54
#Stop Time: 16/07/2015 09:54:51
#Total length: 06:13:57
#Folder Name Pictures: ISIS/July2015/Mg1/unloading1_0000mbar_070C
chopit(35726,0,1201,60,'Mg1_unloading1_0000mbar_070C_60s') 
chopit(35726,0,3601,300,'Mg1_unloading1_0000mbar_070C_300s') 
chopit(35726,0,18001,600,'Mg1_unloading1_0000mbar_070C_600s') 
chopit(35726,0,21601,1800,'Mg1_unloading1_0000mbar_070C_1800s') 
#01:00: Started to decrease pressure from 1180mbar to 0mbar.
#01:30: Final pressure of 0000mbar reached.


reload(nr)
nr.current_detector = nr.old_detector
offspecQplot('35726',0.01,0.06,'Mg1_unloading1_0000mbar_070C_OFFSPEC',spec="117",nslices=2,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)



#Static Measurement Mg1_unloaded1_0000mbar_070C
#Start Time: 16/07/2015 09:59:05
#Stop Time: 16/07/2015 11:00:42
#Folder Name Pictures: NA
nr.nrNRFn("35727","Mg1_unloaded1_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35728","Mg1_unloaded1_0000mbar_070C_th=1.7","1.7","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_unloaded1_0000mbar_070C_th=0.5RvQ,Mg1_unloaded1_0000mbar_070C_th=1.7RvQ","Mg1_unloaded1_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")


#loading2
#Start Time: 16/07/2015 11:19:36
#Stop Time: 16/07/2015 15:20:01
#Total length: 4:00:25
#Folder Name Pictures: ISIS/July2015/Mg1/loading2_1200mbar_070C
chopit(35729,0,601,60,'Mg1_loading2_1200mbar_070C_60s')
chopit(35729,0,3601,300,'Mg1_loading2_1200mbar_070C_300s') 
chopit(35729,0,3601,600,'Mg1_loading2_1200mbar_070C_600s') 
chopit(35729,0,14401,1800,'Mg1_loading2_1200mbar_070C_1800s') 
#240056 at Project_X pressure software corresponds with 00:00
#00:30: Started to increase pressure from 0mbar to 1200mbar. Flow=10 sscm, Vout=0V
#01:10: 100 mbar
#01:42: 200 mbar
#02:58: 400 mbar
#03:30  500 mbar
#04:10 600 mbar
#05:20 800 mbar
#06:05 900 mbar
#07:00 1000 mbar
#08:00 1100 mbar
#10:00 1200 mbar

offspecQplot('35729',0.01,0.06,'Mg1_loading2_1200mbar_070C_OFFSPEC',nslices=4,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)


#Static Measurement Mg1_loaded2_1200mbar_070C
#Start Time: 16/07/2015 09:59:05
#Stop Time: 16/07/2015 10:59:34
#Folder Name Pictures: NA
nr.nrNRFn("35730","Mg1_loaded2_1200mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35731","Mg1_loaded2_1200mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_loaded2_1200mbar_070C_th=0.5RvQ,Mg1_loaded2_1200mbar_070C_th=2.0RvQ","Mg1_loaded2_1200mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#kinetic measurement Mg1_unloading2_0000mbar_070C

#Start Time: 16/07/2015 16:51:14
#Stop Time: 16/07/2015 21:57:36
#Total length: 5:06:22
#Folder Name Pictures: ISIS/July2015/Mg1/unloading2_0000mbar_070C
chopit(35732,0,1201,120,'Mg1_unloading2_0000mbar_070C_th0.7_120sec') 
chopit(35732,0,3601,600,'Mg1_unloading2_0000mbar_070C_th0.7_600sec') 
chopit(35732,0,18001,1800,'Mg1_unloading2_0000mbar_070C_1800s') 
#01:00: Started to decrease pressure from 1180mbar to 0mbar.
#02:00: Final pressure of 0000mbar reached.
#03:00: Pictures started.

offspecQplot('35732',0.01,0.06,'Mg1_unloading2_0000mbar_070C_OFFSPEC',nslices=10,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Static Measurement Mg1_unloaded2_0000mbar_070C
#Start Time: 16/07/2015 22:03:39
#Stop Time: 16/07/2015 23:04:05
#Folder Name Pictures: NA
nr.nrNRFn("35733","Mg1_unloaded2_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35734","Mg1_unloaded2_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_unloaded2_0000mbar_070C_th=0.5RvQ,Mg1_unloaded2_0000mbar_070C_th=2.0RvQ","Mg1_unloaded2_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#kinetic measurement Mg1_unloaded2_000mbar_cooling070to030C

#in run title: Mg1_unloading_0000mbar_070C !!!!!!!!!!!!!!
#no piezoslit installed !!!!!!! This means that effective 'sample slit' = 70*sim(0.7) = 0.85 mm

#Start Time: 16/07/2015 23:18:54
#Stop Time: 16/07/2015 00:15
#Total length: 57:10
#Folder Name Pictures: NA
chopit(35735,0,1201,120,'Mg1_unloading2_0000mbar_cooling070to030C_th0.7_120sec') 
chopit(35735,0,3001,300,'Mg1_unloading2_0000mbar_cooling070to030C_th0.7_300sec') 
#23:15: Started to decrease temperature: set value -> 30C, sample T starts with 80C !
#23:21 T_sample=75C
#23:28 65C
#23:33 60C
#23:48 50C
#00:00 45C
#00:15 40C
#during cooling no changes visible

#Install piezo slit
#Vent sample cell at 00:33

#kinetic measurement Mg1_unloaded2_air_030C

#Start Time: 17/07/2015 00:36:27
#Stop Time: 17/07/2015 01:36
#Total length: 1:00
#Folder Name Pictures: NA
chopit(35736,0,1201,120,'Mg1_unloading2_air_030C_th0.7_120sec') 
chopit(35736,0,3001,300,'Mg1_unloading2_air_030C_th0.7_300sec') 
#sample changes as a result of the air! after 50 min no changes visible

offspecQplot('35736',0.01,0.06,'Mg1_unloading2_air_070C_OFFSPEC',nslices=2,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)


#######################
#direct beam measurements
#######################
#for kinetic measurements:
#piezo slit should be in the beam. This slit was removed and then the sample taken out. 
#Piezo slit replaced and aligned without sample in the straight beam (theta=0)
# intensity was too high. width 1st slit changed from 30 to 3 mm, then coutrate 0.4 kHz 
# runs 35737 and 3538: Start Tine: 17/07/2015 02:39:01; Stop time: 17/07/2015 3:42:27 

#for static measurements:
#no piezo slit
# runs 35739 and 35740: Start Tine: 17/07/2015 3:46:03; Stop time: 17/07/2015 4:46


nr.nrDBFn("35737","w37","35738","w38","LDDB05k","108","120","4.5",binning,"",fitspline=10,diagnostics="0")

nr.nrDBFn("35739","w39","35740","w40","LDDB05s","108","120","5.3",binning,"",fitspline=10,diagnostics="0")


##################
#Mg-1
#################
#Static Measurement Mg1_unloaded2Air_1000mbar_030C
#Start Time: 17/07/2015 05;30:49
#Stop Time: 17/07/2015  06:33:55
#Total length: 1:03:06
#Folder Name Pictures: NA
nr.nrNRFn("35741","Mg1_unloaded2_air_030C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35742","Mg1_unloaded2_air_030C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_unloaded2_air_030C_th=0.5RvQ,Mg1_unloaded2_air_030C_th=2.0RvQ","Mg1_unloaded2_air_030C_anglesCombined","0","","","0",combine_binning,1.0,"2")


#Mg1_unloaded2_Air_21072015_13:00_30C
nr.nrNRFn("35779","Mg1_unloaded2_Air_21072015_13:00_30C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35780","Mg1_unloaded2_Air_21072015_13:00_30C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg1_unloaded2_Air_21072015_13:00_30C_th=0.5RvQ,Mg1_unloaded2_Air_21072015_13:00_30C_th=2.0RvQ","Mg1_unloaded2_Air_21072015_13:00_30C_anglesCombined","0","","","0",combine_binning,1.0,"2")


######################################
#Hf-1
######################################
#Static Measurement Hf1_Virgin_0000mbar_120C
#Start Time: 17/07/2015 07;02:39
#Stop Time: 17/07/2015  08:04;46
#Total length: 1:03:06
#Folder Name Pictures: NA
nr.nrNRFn("35743","Hf1_Virgin_0000mbar_120C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35744","Hf1_Virgin_0000mbar_120C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Hf1_Virgin_0000mbar_120C_th=0.5RvQ,Hf1_Virgin_0000mbar_120C_th=2.0RvQ","Hf1_Virgin_0000mbar_120C_anglesCombined","0","","","0",combine_binning,1.0,"2")

offspecQplot('35744',0.01,0.06,'Hf1_Virgin_0000mbar_120C_OFFSPEC',nslices=1,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)


#Kinetic Measurement Hf1_loading_0010mbar_120C
#Start Time: 17/07/2015 08:22:54
#Stop Time: 17/07/2015 16:12:47
#Total length: 07:49;53
#Folder Name Pictures: NA
chopit(35745,480,1681,300,'Hf1_loading1_0010mbar_120C_th0.7_300sec') 
chopit(35745,480,6481,600,'Hf1_loading1_0010mbar_120C_th0.7_600sec') 
chopit(35745,480,27481,1800,'Hf1_loading1_0010mbar_120C_th0.7_1800sec') 
chopit(35745,480,27481,3600,'Hf1_loading1_0010mbar_120C_th0.7_3600sec') 

#5520 at Project_X pressure software corresponds with 00:00
#08:00 Started to increase pressure from 0mbar to 10mbar. Flow=10 sscm, Vout=0V. Ppump=2.94
#08:15 10 mbar
#2:20:00 increased Vout to 9.9V (auto) Ppum 3.69

offspecQplot('35745',0.01,0.06,'Hf1_loading1_0010mbar_120C_OFFSPEC',nslices=7,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Kinetic Measurement Hf1_loading1_1000mbar_120C
#Start Time: 17/07/2015 16:16:59
#Stop Time: 17/07/2015 19:17:00
#Total length: 03:14:01
#Folder Name Pictures: NA
chopit(35746,0,1201,60,'Hf1_loading1_1000mbar_120C_th0.7_60sec') 
chopit(35746,0,10801,600,'Hf1_loading1_1000mbar_120C_th0.7_600sec') 
#01:00 150 mbar
#02:00 300 mbar
#03;30 600 mbar
#05:00 800 mbar
#06:00 900 mbar
#07:30 1000 mbar

offspecQplot('35746',0.01,0.06,'Hf1_loading1_1000mbar_120C_OFFSPEC',nslices=6,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Kinetic Measurement Hf1_unloading1_0000mbar_120C
#Start Time: 17/07/2015 19:34:23
#Stop Time: 17/07/2015 20:35:00
#Total length: 1:00:37
#Folder Name Pictures: NA
chopit(35747,0,601,60,'Hf1_unloading1_0000mbar_120C_60sec') 
chopit(35747,0,3601,300,'Hf1_unloading1_0000mbar_120C_300sec') 
#chopit(35747,0,3601,600,'Hf1_unloading1_0000mbar_120C_th0.7_600sec') 
#01:00 0 mbar

offspecQplot('35747',0.01,0.06,'Hf1_unloading1_0000mbar_120C_OFFSPEC',nslices=4,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Kinetic Measurement Hf1_unloading1_Air_120C
#Start Time: 17/07/2015 20:37:34
#Stop Time: 17/07/2015  21:08:14
#Total length: 30:40
#Folder Name Pictures: NA
chopit(35748,0,1801,300,'Hf1_unloading1_Air_120C_300sec') 
#Delay between allowing air to enter the cell and the start of the measurement of about 2 min.
#Relatively Large initial effect seen. (almost) nothing afterwards)

offspecQplot('35748',0.01,0.06,'Hf1_unloading1_Air_120C_OFFSPEC',nslices=2,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

###########################################################################################
#Mg-2
###########################################################################################
#Virgin Sample in Vacuum 
#Start Time: 17/07/2015 22:12:25
#Stop Time: 17/07/2015  23:12:54
#Total length: 1:00:29
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35749","Mg2_Virgin_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35750","Mg2_Virgin_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_Virgin_0000mbar_070C_th=0.5RvQ,Mg2_Virgin_0000mbar_070C_th=2.0RvQ","Mg2_Virgin_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#loading1
#Start Time: 17/07/2015 23:29:13
#Stop Time: 18/07/2015 22:56:30
#Total length: 23:33:17
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_loading1_1200mbar_070C
chopit(35751,120,721,60,'Mg2_loading1_1200mbar_th0.7_60s') 
chopit(35751,120,3721,300,'Mg2_loading1_1200mbar_070C_300s') 
chopit(35751,120,84721,3600,'Mg2_loading1_1200mbar_070C_3600s') 
#stopped since DAQ did not work properly.

offspecQplot('35751',0.01,0.06,'Mg2_loading1_1200mbar_070C_OFFSPEC',nslices=10,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Start Time: 18/07/2015 23:04:11
#Stop Time: 19/07/2015 08:04:29
#Total length: 9:00;18
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_loading1_1200mbar_070C
chopit(35752,0,32401,3600,'Mg2_loading1_1200mbar_3600sec-2') 

offspecQplot('35752',0.01,0.06,'Mg2_loading1_1200mbar_070C_OFFSPEC-2',nslices=5,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#Loaded1 Sample
#Start Time: 19/07/2015 08:07:55
#Stop Time: 19/07/2015  09:38:25
#Total length: 1;00:30
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35753","Mg2_loaded1_1200mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35754","Mg2_loaded1_1200mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_loaded1_1200mbar_070C_th=0.5RvQ,Mg2_loaded1_1200mbar_070C_th=2.0RvQ","Mg2_loaded1_1200mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#unloading1 Sample
#Start Time: 19/07/2015 09:27:23
#Stop Time: 19/07/2015 17:12:39
#Total length: 07:45:16
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_unloading1_0000mbar_070C
chopit(35755,0,601,60,'Mg2_unloading1_0000mbar_070C_th0.7_60sec') 
chopit(35755,0,3601,300,'Mg2_unloading1_0000mbar_070C_th0.7_300sec') 
chopit(35755,0,27001,1800,'Mg2_unloading1_0000mbar_070C_th0.7_1800sec') 
#1:00 Final pressure of 0 mbar reached.

offspecQplot('35755',0.01,0.06,'Mg2_unloading1_0000mbar_070C_OFFSPEC',nslices=7,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#unloaded1 Sample
#Start Time: 19/07/2015 17:15:45
#Stop Time: 19/07/2015  18:19:41
#Total length: 1:03:56
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35756","Mg2_unloaded1_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35757","Mg2_unloaded1_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_unloaded1_0000mbar_070C_th=0.5RvQ,Mg2_unloaded1_0000mbar_070C_th=2.0RvQ","Mg2_unloaded1_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")


#Start Time: 19/07/2015 18:32:46
#Stop Time: 19/07/2015 23:10:01
#Total length: 4:37:15
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_loading2_1200mbar_070C
chopit(35758,0,601,60,'Mg2_loading2_1200mbar_070C_th0.7_60sec') 
chopit(35758,0,3601,300,'Mg2_loading2_1200mbar_070C_th0.7_300sec') 
chopit(35758,0,16201,1800,'Mg2_loading2_1200mbar_070C_th0.7_1800sec') 
#214550 at Project_X pressure software corresponds with 00:00
#00:39 Started to increase pressure from 0mbar to 1200mbar. Flow=10 sscm, Vout=0V. Ppump=2.94
#08:30 Final pressure reached

offspecQplot('35758',0.01,0.06,'Mg2_loading2_0000mbar_070C_OFFSPEC',nslices=4,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#loaded2 Sample
#Start Time: 19/07/2015 23:13:53
#Stop Time: 20/07/2015  00:14:22
#Total length: 1:00:29
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35759","Mg2_loaded2_1200mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35760","Mg2_loaded2_1200mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_loaded2_1200mbar_070C_th=0.5RvQ,Mg2_loaded2_1200mbar_070C_th=2.0RvQ","Mg2_loaded2_1200mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#Start Time: 20/07/2015 00:28:42 
#Stop Time: 20/07/2015 06:45:27
#Total length: 6:16:45
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_unloading2_0100mbar_070C
chopit(35761,0,601,60,'Mg2_unloading2_0100mbar_070C_th0.7_60sec') 
chopit(35761,0,3601,300,'Mg2_unloading2_0100mbar_070C_th0.7_300sec') 
chopit(35761,0,21601,1800,'Mg2_unloading2_0100mbar_070C_th0.7_1800sec') 
#236000 at Project_X pressure software corresponds with 00:00
#00:25 Started to decrease pressure from 1200mbar to 100mbar. Flow=10 sscm, Vout=5V. Ppump=
#03:00 Camara Switched on.
#03:00 200 mbar reached

offspecQplot('35761',0.01,0.06,'Mg2_unloading2_0100mbar_070C_OFFSPEC',nslices=6,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#unloaded2 Sample @100 mbar
#Start Time: 20/07/2015 
#Stop Time: 20/07/2015 
#Total length: 1:00:29
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35762","Mg2_unloaded2_0100mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35763","Mg2_unloaded2_0100mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_unloaded2_0100mbar_070C_th=0.5RvQ,Mg2_unloaded2_0100mbar_070C_th=2.0RvQ","Mg2_unloaded2_0100mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#Start Time: 20/07/2015 08:02:28
#Stop Time: 20/07/2015 12:09:55
#Total length: 04:07:28
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_unloading2_0040mbar_070C
chopit(35764,0,301,60,'Mg2_unloading2_0040mbar_070C_th0.7_60sec') 
chopit(35764,0,7201,600,'Mg2_unloading2_0040mbar_070C_th0.7_600sec') 
chopit(35764,0,14401,1800,'Mg2_unloading2_0040mbar_070C_th0.7_1800sec') 
#263400 at Project_X pressure software corresponds with 00:00
#00:15 Started to decrease pressure from 1200mbar to 100mbar. Flow=10 sscm, Vout=5V. Ppump=
#01:30 70 mbar
#03:00 50 mbar
#05:30 40 mbar

offspecQplot('35764',0.01,0.06,'Mg2_unloading2_0040mbar_070C_OFFSPEC',nslices=6,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#unloaded2 Sample @ 40 mbar
#Start Time:20/07/2015 12:12:45
#Stop Time: 20/07/2015 13:13:14
#Total length: 1:00:29
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35765","Mg2_unloaded2_0040mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35766","Mg2_unloaded2_0040mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_unloaded2_0040mbar_070C_th=0.5RvQ,Mg2_unloaded2_0040mbar_070C_th=2.0RvQ","Mg2_unloaded2_0040mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#unloading2 Sample @ 0 mbar
#Start Time: 20/07/2015 14:04:22
#Stop Time: 21/07/2015 01:09:09
#Total length: 11:04:47
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_unloading2_0000mbar_070C
chopit(35767,0,601,60,'Mg2_unloading2_0000mbar_070C_th0.7_60sec') 
chopit(35767,0,2401,300,'Mg2_unloading2_0000mbar_070C_th0.7_300sec') 
chopit(35767,0,14401,600,'Mg2_unloading2_0000mbar_070C_th0.7_600sec') 
chopit(35767,0,39601,1800,'Mg2_unloading2_0000mbar_070C_th0.7_1800sec') 
#instantaniously set to vacuum 

offspecQplot('35767',0.01,0.06,'Mg2_unloading2_0000mbar_070C_OFFSPEC',nslices=11,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#unloaded2 Sample @ 0 mbar
#Start Time:21/07/2015 01:12:23
#Stop Time: 21/07/2015 02:12:53 
#Total length: 01:00:30
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35768","Mg2_unloaded2_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35769","Mg2_unloaded2_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_unloaded2_0000mbar_070C_th=0.5RvQ,Mg2_unloaded2_0000mbar_070C_th=2.0RvQ","Mg2_unloaded2_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#unloading2 Sample @Air
#Start Time: 21/07/2015 02:23:21
#Stop Time: 21/07/2015 03:23:45
#Total length: 1:00:24
#Folder Name Pictures: ISIS/July2015/Mg2/Mg2_unloading2_0000mbar_070C
chopit(35770,0,601,60,'Mg2_unloading2_air_030C_th0.7_60sec') 
chopit(35770,0,3601,300,'Mg2_unloading2_air_030C_th0.7_300sec') 
#First minute of unloading not captured since valve in blockhouse had to be opened.
#Gradual decrease of temperature during the run.

offspecQplot('35770',0.01,0.06,'Mg2_unloading2_air_030C_OFFSPEC',nslices=2,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)


##############
#Mg2_unloaded2_Air_21072015_13:00_30C
nr.nrNRFn("35781","Mg2_unloaded2_Air_21072015_13:00_30C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35782","Mg2_unloaded2_Air_21072015_13:00_30C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg2_unloaded2_Air_21072015_13:00_30C_th=0.5RvQ,Mg2_unloaded2_Air_21072015_13:00_30C_th=2.0RvQ","Mg2_unloaded2_Air_21072015_13:00_30C_anglesCombined","0","","","0",combine_binning,1.0,"2")
##############

############################################################################################################
#Mg-3
############################################################################################################
#General Remark: Sample looks extremely dirty!
#Virgin state looks reasonably similar to Mg-1 and Mg-2

#Virgin Sample @ 0 mbar
#Start Time:21/07/2015  03:43:40
#Stop Time: 21/07/2015 
#Total length: 
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35771","Mg3_Virgin_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35772","Mg3_Virgin_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg3_Virgin_0000mbar_070C_th=0.5RvQ,Mg3_Virgin_0000mbar_070C_th=2.0RvQ","Mg3_Virgin_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")
#Gradual increase of Temperature from 30 to 70C During the run. (reasonably fast)

#loading1 @ 300 mbar
#Start Time:21/07/2015 04:56:38  
#Stop Time: 21/07/2015 06:47:22
#Total length: 1:50:44 
#Folder Name Pictures: ISIS/July2015/Mg3/Mg3_loading1_0300mbar_070C
chopit(35773,0,1801,60,'Mg3_loading1_0300mbar_070C_60s') 
chopit(35773,0,6001,300,'Mg3_loading1_0300mbar_070C_300s')
chopit(35773,0,6001,600,'Mg3_loading1_0300mbar_070C_600s') 
# 338627 t Project_X pressure software corresponds with 00:00
# 01:00 150 mbar
# 02:00

offspecQplot('35773',0.01,0.06,'Mg3_loading1_0300mbar_070C_OFFSPEC',nslices=5,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#loaded1 @ 300 mbar
#Start Time:21/07/2015  06:52:25
#Stop Time: 21/07/2015 
#Total length: 
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35774","Mg3_loaded1_0300mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35775","Mg3_loaded1_0300mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg3_loaded1_0300mbar_070C_th=0.5RvQ,Mg3_loaded1_0300mbar_070C_th=2.0RvQ","Mg3_loaded1_0300mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

#unloading1 @ 000 mbar
#Start Time:21/07/2015  8:07:22
#Stop Time: 21/07/2015  09:20:54
#Total length:  1:13:32
#Folder Name Pictures: ISIS/July2015/Mg3/Mg3_unloading1_0000mbar_070C
chopit(35776,0,601,60,'Mg3_unloading1_0000mbar_070C_60s') 
chopit(35776,0,4201,300,'Mg3_unloading1_0000mbar_070C_300s')
chopit(35776,0,4201,600,'Mg3_unloading1_0000mbar_070C_600s') 
# instantaneous vacuum at start of run

offspecQplot('35776',0.01,0.06,'Mg3_unloading1_0000mbar_070C_OFFSPEC',nslices=5,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)

#unloaded1 @ 000 mbar 
#Start Time:21/07/2015  9:24:54
#Stop Time: 21/07/2015 
#Total length: 
#Folder Name Pictures: ISIS/July2015/ NA
nr.nrNRFn("35777","Mg3_unloaded1_0000mbar_070C_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35778","Mg3_unloaded1_0000mbar_070C_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("Mg3_unloaded1_0000mbar_070C_th=0.5RvQ,Mg3_unloaded1_0000mbar_070C_th=2.0RvQ","Mg3_unloaded1_0000mbar_070C_anglesCombined","0","","","0",combine_binning,1.0,"2")

###########################################################################################################
#cleaned substrate

nr.nrNRFn("35783","cleaned_substrate_th=0.5","0.5","LDDB05s","114","110","118",binning,"")
nr.nrNRFn("35784","cleaned_substrate_th=2.0","2.0","LDDB05s","114","110","118",binning,"") 
nr.NRCombineDatafn("cleaned_substrate_th=0.5RvQ,cleaned_substrate_th=2.0RvQ","cleaned_substrate_anglesCombined","0","","","0",combine_binning,1.0,"2")

#################################################################################3
#looking at off-spec intensities

offspecQplot('35726',0.01,0.06,'Mg1_unloading1_0000mbar_070C_OFFSPEC',nslices=6,sarray = [] , Nqx=150, Nqz=150, zmin=5e-7, zmax=0.01)
#j=[0, 3601.0, 7202.0, 10803.0, 14404.0, 18005.0, 21606.0, 25207.0, 28808.0, 32409.0, 36010.0, 39611.0, 43212.0, 46813.0, 50414.0, 54015.0, 57616.0, 61217.0, 64818.0, 68419.0, 72019.0]
j=[0, 3740.0, 7480.0, 11220.0, 14960.0, 18700.0, 22437.0]
qzmin = 0.025
qzmax=0.030
for i in range(len(j)-1):
    
    QxQzcuts('Mg1_unloading1_0000mbar_070C_OFFSPEC_'+str(j[i])+'-'+str(j[i+1])+'qxqz', qzmin=qzmin, qzmax=qzmax,plot=False)
    SaveAscii('Mg1_unloading1_0000mbar_070C_OFFSPEC_'+str(j[i])+'-'+str(j[i+1])+'qxqz'+'_cut_'+str(qzmin)+'-'+str(qzmax), 'U:/VanWell/July_2015/Cuts/Mg1_unloading1_0000mbar_070C_OFFSPEC_'+str(j[i])+'-'+str(j[i+1])+'qxqz'+'_cut_'+str(qzmin)+'-'+str(qzmax)+'.dat', WriteSpectrumID = False, CommentIndicator = "#", Separator = "Tab", ColumnHeader = False)















New function: can be used like this:
chopit2(34253, 'test', tslice = 30000, userdirectory = "U:/vanWell/April_2015/savetest/" )    
or like this:
chopit2(34253, 'test', nslices = 5, userdirectory = "U:/vanWell/April_2015/savetest/" )    
or like this:
chopit2(35755, 'test', start = 10000, nslices = 4, userdirectory = "U:/vanWell/April_2015/savetest/",loadcrpt=1 )    
saves the individual data slices as dat files  as well in the folders already created.

help(chopit2)




############################################




