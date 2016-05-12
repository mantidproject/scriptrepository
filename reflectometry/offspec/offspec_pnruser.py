''' This module should be imported into a user reduction script, for instance: import offspec_pnruser as p
The following variables should/can then be defined:

p.dbpnr = 'none' # direct beam workspace for pnr
p.dbpa = 'none'  # directbeamworkspace for pa once those have been measured
savepath = "u:/" # The directory in which to save the reduced data

p.pnrangles = (0.5, 1.0, 2.0) # The angles at which you measure!
p.pnrqlimits = ["0.008","-0.02","0.15"] # The q ranges and binning resolution
p.pnrscalefactor =1.0                   # If you want to apply a scalefactor

paangles = [0.5, 1.0, 2.0]  # Same for PA
p.paqlimits = ["0.008","-0.02","0.15"]
p.pascalefactor =1.0


Only adjust the following if you know what they mean:

p.detectorlimits = {'low': 390, 'specular': 404, 'high': 410} # Low is lower limit for region of interest on the
                detector, high upper limit for ROI and specular the position of the reflection

p.binningpars = ["2.0","0.02","12.0","2"]  # wavelength binning, for historical reasons this includes the choice of monitor
                                           # default is usually ok

lmin, lmax = 2.2, 12.0 #Determine at what wavelength to crop the workspace, should be an ok default

pnrqgroup = True, pnrbckg = True # and their equivalents in PA. Decide if you sum in Q or lambda, 
                                 # subtract background or not.

calibration: default is set to current polarisation corrections
                                 
-------------------------------------------------                                 
For offspecular data:            
ofangle = 0.7
ofqlimits = '-5e-4,5e-4,0.02,0.1'
ofsNqx = 200
ofsNqz = 200
                      
The functions provided in this module are:

PNRreduce(runs, name, angles=pnrangles, qlimits = pnrqlimits, scalefactor=pnrscalefactor, diagnostics = False)
PAreduce(runs, name, angles=paangles, qlimits = paqlimits, scalefactor=pascalefactor, diagnostics = False)

where runs is a list of strings containing the runnumbers for each angle, for example:
runs = ["36422","36423+36424","36425-36430"]
name is the name of your resulting workspace/file

print_calibrations() 
prints a list of all the available polarisation corrections

OFSreduce(runs, name, withpol = 'pnr')
Reduces offspcular data

'''
import offspec_red as nr
reload(nr)
#last modified: 12/05/2016
#by: njs
import numpy as np
import itertools
#from math import *
from mantid.simpleapi import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace 
import matplotlib.pyplot as plt
try:
  from mantidplot import *
except ImportError:
  pass
import offspec_offspec as of
reload(of)

detectorlimits={'low': 390, 'specular': 404, 'high': 410}
binningpars = ["2.0","0.02","13.0","2"]
dbpnr = 'none'
dbpa = 'none'
savepath = "u:/"

pnrqlimits = ["0.008","-0.02","0.15"]
pnrscalefactor =1.0
pnrangles = [0.5, 1.0, 2.0]
pnrqgroup = False
pnrbckg = False


paqlimits = ["0.008","-0.02","0.15"]
pascalefactor =1.0
paangles = [0.5, 1.0, 2.0]
paqgroup = False
pabckg = False

calibration = nr.PolarisationCalibration.current

ofangle = 0.7
ofqlimits = '-5e-4,5e-4,0.02,0.1'
ofsNqx = 200
ofsNqz = 200

lmin, lmax = 2.2, 12.0
fourpi = 4*np.pi
deg2rad = np.pi/180

def print_calibrations():
    calibrations = nr.PolarisationCalibration.get_calibrations()
    print "========================================\n"
    for each in calibrations.keys():
        print each + ":  "+calibrations[each].comment+"\n"
    print "========================================\n"
    print "Current calibration: "+nr.PolarisationCalibration.current.name
    print "========================================\n"
    

def PNRreduce(runs, name, angles=None, qlimits = None, scalefactor=None, bckg = None, specular = None, diagnostics = False):
    if not angles:
        angles=pnrangles
    if not qlimits:
        qlimits = pnrqlimits
    if not scalefactor:
        scalefactor = pnrscalefactor
    if not bckg:
        bckg = pnrbckg
    if not specular:
        specular = itertools.cycle([detectorlimits['specular']])
    wkspnames = ('low','mid', 'high', 'veryhigh', 'extremelyhigh', 'toomany', 'waytoomany', '8', '9', '10')
    tocombineup = ""
    tocombinedown = ""
    for run, angle, wname, thisspecular in zip(runs,angles,wkspnames,specular):
        print "Runs: "+run+", angle: "+ str(angle)+", specular:"+str(thisspecular)+", name:"+name+str(wname)
        nr.nrPNRFn(run,name+wname,str(angle),dbpnr,thisspecular,detectorlimits['low'],detectorlimits['high'],binningpars,floodfile="",qgroup=pnrqgroup,PNRwithPA=False,pnums=["1","2"],doCorrs=True, calibration = calibration, subbgd=bckg, diagnostics=diagnostics)
        qmin = fourpi*np.sin(angle*deg2rad)/lmax
        qmax = fourpi*np.sin(angle*deg2rad)/lmin
        CropWorkspace(name+wname+'RvQ',qmin,qmax,OutputWorkspace=name+wname+"normcorrRvQ")
        tocombineup = tocombineup+name+wname+"normcorrRvQ_1,"
        tocombinedown = tocombinedown+name+wname+"normcorrRvQ_2,"
    tocombineup = tocombineup[:-1]
    tocombinedown = tocombinedown[:-1]
    print tocombineup, tocombinedown
    try:
        a1=nr.NRCombineDatafn(tocombineup,name+"U","0","","","0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        nr.NRCombineDatafn(tocombinedown,name+"D","2",a1[0],a1[1],"0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        plotSpectrum([name+"U",name+"D"],0,1)
        GroupWorkspaces(name+"U,"+name+"D", Outputworkspace = name)
    except:
       Rebin(tocombineup,pnrqlimits[0]+","+pnrqlimits[1]+","+pnrqlimits[2],OutputWorkspace=name+"U")
       Rebin(tocombinedown,pnrqlimits[0]+","+pnrqlimits[1]+","+pnrqlimits[2],OutputWorkspace=name+"D")
       GroupWorkspaces(name+"U,"+name+"D", Outputworkspace = name)
    SaveAscii(name+"U", savepath+name+".u", CustomSeparator="    ", ColumnHeader = False,WriteSpectrumID=False)
    SaveAscii(name+"D", savepath+name+".d",CustomSeparator="    ",ColumnHeader = False, WriteSpectrumID=False)
    

def PAreduce(runs, name, angles=None, qlimits = None, scalefactor=None, bckg = None, specular = None, diagnostics = False):
    if not angles:
        angles=paangles
    if not qlimits:
        qlimits = paqlimits
    if not scalefactor:
        scalefactor = pascalefactor
    if not bckg:
        bckg = pabckg
    if not specular:
        specular = itertools.cycle([detectorlimits['specular']])

    wkspnames = ('low','mid', 'high', 'veryhigh', 'extremelyhigh', 'toomany', 'waytoomany', '8', '9', '10')
    tocombineuu, tocombineud, tocombinedu, tocombinedd  = "", "", "", ""
    print "Angles used: "+str(angles)
    for run, angle, wname, thisspecular in zip(runs,angles,wkspnames, specular):
        print "Runs: "+run+", angle: "+ str(angle)+", specular:"+str(thisspecular)+", name:"+name+str(wname)
        nr.nrPNRFn(run,name+wname,str(angle),dbpa,thisspecular,str(detectorlimits['low']),str(detectorlimits['high']),binningpars,floodfile="",qgroup=paqgroup,PNRwithPA=True,pnums=["1","2","3","4"],doCorrs=True, calibration = calibration, subbgd=bckg, diagnostics=diagnostics)
        qmin = fourpi*np.sin(angle*deg2rad)/lmax
        qmax = fourpi*np.sin(angle*deg2rad)/lmin
        extension = "normcorrRvQ"
        try:
            CropWorkspace(name+wname+extension,qmin,qmax,OutputWorkspace=name+wname+extension)
        except:
            extension = "RvQ"
            CropWorkspace(name+wname+extension,qmin,qmax,OutputWorkspace=name+wname+extension)
        tocombineuu = tocombineuu +name+wname+extension+"_1,"
        tocombineud = tocombineud +name+wname+extension+"_2,"
        tocombinedu = tocombinedu +name+wname+extension+"_3,"
        tocombinedd = tocombinedd +name+wname+extension+"_4,"
        if not diagnostics:
            DeleteWorkspace(name+wname+"RvQ")
    tocombineuu = tocombineuu[:-1]
    tocombineud = tocombineud[:-1]
    tocombinedu = tocombinedu[:-1]
    tocombinedd = tocombinedd[:-1]
    print tocombineuu+tocombineud+tocombinedu+tocombinedd
    try:
        a1=nr.NRCombineDatafn(tocombineuu,name+"UU","0","","","0",qlimits,str(scalefactor),"2", diagnostics=diagnostics)
        nr.NRCombineDatafn(tocombineud,name+"UD","2",a1[0],a1[1],"0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        nr.NRCombineDatafn(tocombinedu,name+"DU","2",a1[0],a1[1],"0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        nr.NRCombineDatafn(tocombinedd,name+"DD","2",a1[0],a1[1],"0",qlimits,str(scalefactor),"2",diagnostics=diagnostics)
        Plus(name+"UD", name+"DU", OutputWorkspace = name+"FL")
        GroupWorkspaces(name+"UU,"+name+"UD,"+name+"DU,"+name+"DD,"+name+"FL", OutputWorkspace = name)
        plotSpectrum([name+"UU",name+"DD"],0,1)
        for dataset, ending in zip((name+'UU', name+'UD', name+'DU', name+'DD', name+'FL'), ('.uu', '.ud', '.du', '.dd', '.flip')):
            SaveAscii(dataset, savepath+name+ending, CustomSeparator="    ", ColumnHeader = False)
       
    except:
        Rebin(tocombineuu,paqlimits[0]+","+paqlimits[1]+","+paqlimits[2],OutputWorkspace=name+"UU")
        Rebin(tocombineud,paqlimits[0]+","+paqlimits[1]+","+paqlimits[2],OutputWorkspace=name+"UD")
        Rebin(tocombinedu,paqlimits[0]+","+paqlimits[1]+","+paqlimits[2],OutputWorkspace=name+"DU")
        Rebin(tocombinedd,paqlimits[0]+","+paqlimits[1]+","+paqlimits[2],OutputWorkspace=name+"DD")
        Plus(name+"UD", name+"DU", OutputWorkspace = name+"FL")
        mtd[name+'FL'] = mtd[name+'FL']/2.0
        GroupWorkspaces(name+"UU,"+name+"UD,"+name+"DU,"+name+"DD,"+name+"FL", OutputWorkspace = name)
        plotSpectrum([name+"UU",name+"DD"],0,1)
    for dataset, ending in zip((name+'UU', name+'UD', name+'DU', name+'DD', name+'FL'), ('.uu', '.ud', '.du', '.dd', '.flip')):
        SaveAscii(dataset, savepath+name+ending, CustomSeparator="    ", ColumnHeader = False)
        
        


    
def OFSreduce(runs, name, withpol = 'pnr', angle = None, qlimits = None, Nqx = None, Nqz = None):
    if not angle:    angle = ofangle
    if not qlimits: qlimits = ofqlimits
    if not Nqx: Nqx = ofsNqx
    if not Nqz: Nqz = ofsNqz
    print angle
    of.DSqxqz(runs,name,angle=angle,qxqzlimits=qlimits,binning=binningpars,Nqx=Nqx,Nqz=Nqz,withpol=withpol)    
    
 
