''' This module should be imported into a user reduction script, for instance: import offspec_pnruser as p
The following variables should/can then be defined:

p.dbnr = 'none' # direct beam workspace for pnr
savepath = "u:/" # The directory in which to save the reduced data

p.nrangles = (0.6, 1.6, 3.0) # The angles at which you measure!
p.nrqlimits = ["0.008","-0.02","0.15"] # The q ranges and binning resolution
p.nrscalefactor =1.0                   # If you want to normalise to 1.0

p.nrqlimits = ["0.008","-0.02","0.15"] # Same for PA
p.nrscalefactor =1.0
p.nrangles = (0.6, 1.6, 3.0)

p.detectorlimits = {'low': 390, 'specular': 404, 'high': 410} # Low is lower limit for region of interest on the
                detector, high upper limit for ROI and specular the position of the reflection

p.binningpars = ["2.0","0.02","12.0","2"]  # wavelength binning, for historical reasons this includes the choice of monitor
                                           # default is usually ok

lmin, lmax = 2.2, 12.0 #Determine at what wavelength to crop the workspace, should be an ok default

The functions provided in here are:

NRreduce(runs, name, angles=pnrangles, qlimits = pnrqlimits, scalefactor=pnrscalefactor, diagnostics = False)
OFSreduce(runs, name, withpol = False, angle = None, qlimits = None, Nqx = None, Nqz = None)

where runs is a list of strings containing the runnumbers for each angle, for example:
runs = ["36422","36423+36424","36425-36430"]
name is the name of your resulting workspace/file

'''
import offspec_red as nr
reload(nr)
#last modified: 12/05/2016
#by: njs
import numpy as np
from math import *
from mantid.simpleapi import *
import itertools, re
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace 
import matplotlib.pyplot as plt
try:
  from mantidplot import *
except ImportError:
  pass
import offspec_offspec as of
reload(of)
import nrkinetic as k
reload(k)

detectorlimits={'low': 395, 'specular': 404, 'high': 414}
binningpars = ["1.0","0.02","14.0","2"]
dbnr = 'none'
savepath = "u:/randomfiles/"

nrqlimits = ["0.008","-0.02","0.2"]
nrscalefactor =1.0
nrangles = [0.5, 1.0, 2.0]
nrqgroup=False
nrbackground=False

ofangle = 0.7
ofqlimits = '-5e-4,5e-4,0.02,0.1'
ofNqx = 500
ofNqz = 500

kinetic_scalefactor = nrscalefactor
kinetic_angle = 0.7
kinetic_qlimits = nrqlimits
dbkin = dbnr

plotzmin, plotzmax = 1e-6, 0.01

lmin, lmax = 1.5 , 12.0
fourpi = 4*np.pi
deg2rad = np.pi/180

def NRreduce(runs, name, angles=None, qlimits = None, scalefactor=None, bckg= None, specular = None, diagnostics = False):
    if not angles:
        angles=nrangles
    if not qlimits:
        qlimits = nrqlimits
    if not scalefactor:
        scalefactor = nrscalefactor
    if not bckg:
        bckg = nrbackground
    if not specular:
        specular = itertools.cycle([detectorlimits['specular']])
    wkspnames = ('low','mid', 'high', 'veryhigh', 'extremelyhigh', 'toomany', 'waytoomany', '8', '9', '10')
    tocombine = []
    for run, angle, wname, thisspecular in zip(runs,angles,wkspnames,specular):
        try:
            print "Runs: "+run+", angle: "+ str(angle)+", specular: "+str(thisspecular)+", name:"+name+wname
            nr.nrNRFn(run,name+wname,str(angle),dbnr,thisspecular,detectorlimits['low'],detectorlimits['high'],binningpars,qgroup=nrqgroup,floodfile="",subbgd=nrbackground, diagnostics=diagnostics)
            qmin = fourpi*np.sin(angle*deg2rad)/lmax
            qmax = fourpi*np.sin(angle*deg2rad)/lmin
            CropWorkspace(name+wname+'RvQ',qmin,qmax,OutputWorkspace=name+wname+"normcorrRvQ")
            tocombine.append(name+wname+"normcorrRvQ")
        except:
            continue
        if not diagnostics:
            DeleteWorkspace(name+wname+'detnorm')
            DeleteWorkspace(name+wname+'norm')
            DeleteWorkspace(name+wname+'RvQ')
                
    tocombine =  ','.join(tocombine)
    print tocombine
    try:
        nr.NRCombineDatafn(tocombine,name,"0","","","0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        plotSpectrum([name],0,1,2)

    except:
        Rebin(tocombine,qlimits[0]+","+qlimits[1]+","+qlimits[2],OutputWorkspace=name)
        
    SaveAscii(name, savepath+name+".dat", CustomSeparator="    ", ColumnHeader = False, WriteSpectrumID=False)
    

    
def OFSreduce(runs, name, withpol = False, angle = None, qlimits = None, Nqx = None, Nqz = None):
    if not angle:    angle = ofangle
    if not qlimits: qlimits = ofqlimits
    if not Nqx: Nqx = ofNqx
    if not Nqz: Nqz = ofNqz
    print angle
    #of.DSqxqz(runs,name,angle=angle,qxqzlimits=qlimits,binning=binningpars,Nqx=Nqx,Nqz=Nqz,withpol=withpol, dspec=detectorlimits['specular'], dmin=detectorlimits['low'], dmax=detectorlimits['high'])    
    of.DSqxqz(runs,name,angle=angle,qxqzlimits=qlimits,binning=binningpars,Nqx=Nqx,Nqz=Nqz,withpol=withpol) 
 
def time_resolve(run, name, start = 0, angle=None, tslice=None,nslices = None ,sarray=[], dlimits = None, diagnostics = False):  
    if not angle:
        angle=kinetic_angle
    if not dlimits:
        dlimits = detectorlimits
    if len(sarray): usearray = True
    else: usearray = False
    k.slice_the_data(rnum=run,output=name,angle = angle, start = start, tslice = tslice, nslices = nslices ,sarray=sarray,usearray=usearray, DB = dbkin, sf=kinetic_scalefactor, qlimits = kinetic_qlimits, dlimits = dlimits, userdirectory = savepath,binning=binningpars,loadcrpt=True, diagnostics = False)
 
def tr_stitch(wslist, name, qlimits = None, scalefactor=None, diagnostics=False):
    if not qlimits:
        qlimits = nrqlimits
    if not scalefactor:
        scalefactor = nrscalefactor
    tocombine = ','.join(wslist)
    nr.NRCombineDatafn(tocombine,name,"0","","","0",qlimits,str(scalefactor),"2",diagnostics=True)
    plotSpectrum([name],0,1,2)
    if not diagnostics:
        [DeleteWorkspace(ws+"reb") for ws in wslist]
    

def tr_offspec(rnum,output, nslices=None, tslice= None, sarray = [], angle=None, qlimits=None, zmin=plotzmin, zmax=plotzmax,qzmin=None,qzmax=None,loadcrpt=True):
    if not angle:
        angle = kinetic_angle
    if not qlimits:
        qlimits = ofqlimits
    limitlist = qlimits.split(',')
    xmin=float(limitlist[0]); xmax=float(limitlist[1])
    ymin=float(limitlist[2]);ymax=float(limitlist[3])
    qmin = fourpi*np.sin(angle*deg2rad)/lmax
    qmax = fourpi*np.sin(angle*deg2rad)/lmin
    
    k.offspecslice(rnum,qmin,qmax,'wrong',nslices = nslices,tslice=tslice,sarray=sarray,DB=dbkin, theta=angle/2.0,binning = binningpars,loadcrpt=loadcrpt)

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
                ConvertToReflectometryQ(InputWorkspace=name, OverrideIncidentTheta=True, IncidentTheta=angle, Extents=qlimits, OutputAsMDWorkspace=False, OutputWorkspace=newname+"qxqz", NumberBinsQx=ofNqx, NumberBinsQz=ofNqz,OutputVertexes="somevertices")
                CloneWorkspace(name,OutputWorkspace=newname+'detnorm')
                of.doDSSCorrections(newname,angle,0,Nqx=ofNqx,Nqz=ofNqz)
                DeleteWorkspace(newname+'detnorm')
                of.offspecPlot(newname+'qxqz', (xmin, xmax), (ymin,ymax), (zmin,zmax),logscale='z')
                of.QxQzcuts(newname+'qxqz',qzmin,qzmax)
            DeleteWorkspace(mtd[name])
    k.offspecslice(rnum,qmin,qmax,output,nslices = nslices,tslice=tslice,sarray=sarray,DB=dbkin, binning=binningpars,theta=angle,loadcrpt=loadcrpt)
   
 
 
 
 
 
 
 
 
 
 
 