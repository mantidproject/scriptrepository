from __future__ import print_function
import offspec_red as nr
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
reload(nr)
#last modified: 12/05/2016
#by: njs
from math import *
import numpy as np
import os,sys,re, time
from mantid.simpleapi import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace 
try:
  from mantidplot import *
except ImportError:
  pass

def getLog(w,log_name):
        # Get handle to the workspace
        try:
          h=mtd[w]
        except:
          print("Can't get Workspace handle")
        #
        # Get access to SampleDetails
	s=h.getRun().getLogData(log_name).value
	
	return s

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
    
def loaddata(rnum, path = '//ndloffspec1/L/RawData/',loadcrpt=0):
    try:
        Load(str(rnum), OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
    except:
        #try:
            if loadcrpt == 0:
              updatefile=loadlatest(str(rnum))   
              Load(Filename='z:/'+updatefile, OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
            else:
              print('trying to load crpt snapshot')
              Load('z:/snapshot_crpt.nxs', OutputWorkspace=str(rnum), LoaderName='LoadEventNexus', LoadMonitors=True)
        #except: 
        #    raise Exception('Could not find data')
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
     
#k.slice_the_data(rnum=run,output=name,angle = angle, start = start, tslice = tslice, nslices = nslices ,sarray=sarray,usearray=usearray, DB = kindb, sf=kinetic_scalefactor, qlimits = qlimits, dlimits = dlimits, userdirectory = savepath,binning=binningpars,loadcrpt=0, diagnostics = False)
def slice_the_data(rnum,output,angle = None,start = 0, tslice=None,nslices = None ,sarray=[],usearray=0, DB = None, sf=1.0,qlimits = None, dlimits = None,userdirectory = None,binning=["1.5","0.02","14.0","2"],loadcrpt=0, diagnostics = False):
    
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
        nr.nrNRFn("",wksp,str(angle),DB,dlimits['specular'],dlimits['low'],dlimits['high'],binning,"",usewkspname=1,sf=sf)
        Rebin(wksp+"RvQ",qlimits[0]+","+qlimits[1]+","+qlimits[2],OutputWorkspace=wksp+"RvQ")
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

def _offspecslice_simple(rnum,btime,etime,qmin,qmax,output, binning,theta, DB,loadcrpt=0):
    wksp =timeslice(rnum,btime,etime,output,loadcrpt=loadcrpt)
    nr.nrNRFn("",wksp,str(theta),DB,nr.current_detector.specular,nr.current_detector.minspec,nr.current_detector.maxspec,binning,usewkspname=True)
    DeleteWorkspace(wksp)
    DeleteWorkspace(wksp+'norm')

def offspecslice(rnum,qmin,qmax,output, theta,binning, DB,start = 0, tslice=None,nslices = None,sarray=[] ,loadcrpt=0):
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
            _offspecslice_simple(rnum, start, end, qmin, qmax, output, binning=binning, theta=theta,DB=DB,loadcrpt=loadcrpt)
        except:
            print(datatimes)
            break

