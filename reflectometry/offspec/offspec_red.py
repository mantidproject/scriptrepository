from math import *
#27/06 njs: add check on latest run files and change data loader frfom Load to LoadISISNexus
from mantid.simpleapi import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace 
import numpy as np
import matplotlib.pyplot as plt
try:
  from mantidplot import *
except ImportError:
  pass

import numpy as n
import os, socket, re
from weakref import WeakValueDictionary

class LinearDetector():
    def __init__(self,specular, minspec, maxspec, floodfile, detectorposition, name = 'lineardetector', nickname = 'LD', btm_background = False, top_background = False, pixelsize=0, idf = None):
        self.minspec =  minspec
        self.maxspec = maxspec
        self.nspec = maxspec - minspec+1
        self.floodfile = os.path.join(os.path.dirname(__file__), floodfile) #default flat in same directory
        self.detectorposition = detectorposition #sample to detector position in metres
        self.name = name
        self.nickname = nickname
        self.btm_background = btm_background
        self.top_background = top_background
        try:
            self.nbackgroundspectra = np.absolute(self.btm_background[1]-self.btm_background[0]) + np.absolute(self.top_background[1]-self.top_background[0])
        except:
            self.nbackgroundspectra = 0
        self.pixelsize = pixelsize
        self.specular = specular
        self.detectortype = "Linear"
    def croptodetector(self, inputworkspace, suffix = 'det'):
        CropWorkspace(InputWorkspace=inputworkspace,OutputWorkspace=inputworkspace+suffix,StartWorkspaceIndex=self.minspec-1,EndWorkspaceIndex=self.maxspec-1)
        return inputworkspace+suffix #Workspace name of cropped workspace


class CurrentRun():
    def __init__(self, localdata):
        self.localdata = localdata
        self.currentcycle = self.find_current_cycle()
        self.lastRun = None
    
    def find_current_cycle(self):
        cycles = os.listdir(self.localdata)
        cycles.sort()
        cycles.reverse()
        reCycle = re.compile('cycle', re.IGNORECASE)
        for element in cycles:
            m = reCycle.search(element)
            if m: 
                return element
    
    def _extractRunFromName(self,Filename):
        reNumber = re.compile('[0-9]+')
        reMatch = reNumber.search(Filename)
        runnumber = int(reMatch.group())
        return runnumber    
        
    def _check_runnumber(self,runnumber):   
        if self.lastRun >= runnumber: # we know about those already
            return "old"
        else: 
            if runnumber == (self.lastRun+1): # The current runnumber is +1 the last saved one
                return "current"
            else:
                return "DoesNotExistYet" # User has given a runnumber that hasn't been measured yet
    
    def _update_lastrun_if_required(self,runnumber):
        if self.lastRun < runnumber:
            self.set_lastRun()
            
    def get_lastRun(self):
        currentCycleData = os.listdir(self.localdata+self.currentcycle)
        while currentCycleData:
            lastRunFilename = currentCycleData.pop()
            if lastRunFilename.endswith('.nxs'):
                lastRun = self._extractRunFromName(lastRunFilename)
                return lastRun

    
    def set_lastRun(self):
        self.lastRun = self.get_lastRun()
        return self.lastRun        
        
    def get_currentRun(self, runnumber):
        self._update_lastrun_if_required(runnumber)
        runnumberIs = self._check_runnumber(runnumber)
        return runnumberIs, runnumber
                    

                    
class NotAtHome():
	def get_currentRun(self,run):
		return ("old",run)

                    
def check_is_host_at_isis():
    hostname = socket.getfqdn()
    reobject = re.compile('isis', re.IGNORECASE)
    matchobject = reobject.search(hostname.lower()) 
    if matchobject: 
        print "This is "+thisinstrument+" at ISIS"
        return True
    else: 
        print "You are on a local machine"
        return False
    
old_detector = LinearDetector(specular = 114, minspec = 5, maxspec = 244, floodfile = "LD240flood_premarch2012.nxs", detectorposition = 3.63, name = 'Old_LD', nickname = 'oldLD', btm_background = (10, 50), top_background = (160, 230), pixelsize = 1.2e-3)        
wsf_detector = LinearDetector(specular = 404, minspec = 6, maxspec = 772, floodfile = "WSF_Flood.nxs", detectorposition = 3.46, name = 'WSD_LD', nickname = 'wsf', btm_background = (440,500), top_background = (440, 500), pixelsize = 0.5e-3)
# position= 3.46 new best guess after Tommy's LC experiment December 2015
wsf_detector.idf = "C:/mantidinstall/instrument/OFFSPEC_Definition.xml"
old_detector.idf = "C:/MantidInstall/Instrument/Offspec_Definition.xml"
babylarmor = LinearDetector(specular = 120, minspec = 4, maxspec = 125, floodfile = '', name = 'Baby_Larmor', detectorposition = 3.53, nickname = 'bLar', pixelsize = 8e-3)       

current_detector = wsf_detector
thishost = "ndxoffspec"
thisinstrument = "OFFSPEC"
print "Current default detector is: "+current_detector.name
print "Change by setting the current_detector in this module (ask your local contact)."

rad2deg=180.0/pi
pulsestart = 0.5 # be generous here, this is just a rough chop to avoid 0 errors in the unit conversion
pulseend = 20.0


if check_is_host_at_isis():
    try:
        currentRun = CurrentRun("//ndloffspec1/L/RawData/")
    except:
        currentRun = CurrentRun("//isis/inst$/NDXOFFSPEC/Instrument/data/")
        print "ndloffspec is missing, using the archive."
else:
	currentRun = NotAtHome()

 
               
class PolarisationCalibration():
    calibrations = WeakValueDictionary()
    
    def __init__(self, name, rhocoefficients, alphacoefficients, Apcoefficients, Ppcoefficients, runs = None, comment = None):
        self.rhocoefficients = rhocoefficients
        self.alphacoefficients = alphacoefficients
        self.Apcoefficients = Apcoefficients
        self.Ppcoefficients = Ppcoefficients
        self.name = name
        self.comment = comment
        PolarisationCalibration.calibrations[name] = self
        self.current = False
        
    @classmethod
    def get_calibrations(cls):
        return cls.calibrations
    @classmethod
    def set_current_calibration(cls, instancename):
        thecalibration = cls.calibrations[instancename]
        thecalibration.current = True
        cls.current = thecalibration
        return thecalibration
        
    def correctPNR(self, UpWksp, DownWksp):  
        print "Applying polarisation corrections "+str(self.name)
        Ip = mtd[UpWksp]
        Ia = mtd[DownWksp]
        CloneWorkspace(Ip,OutputWorkspace="PCalpha")
        CropWorkspace(InputWorkspace="PCalpha",OutputWorkspace="PCalpha",StartWorkspaceIndex="0",EndWorkspaceIndex="0")
        PCalpha=(mtd['PCalpha']*0.0)+1.0
        alpha=mtd['PCalpha']
        CloneWorkspace("PCalpha",OutputWorkspace="PCrho")
        CloneWorkspace("PCalpha",OutputWorkspace="PCAp")
        CloneWorkspace("PCalpha",OutputWorkspace="PCPp")
        rho=mtd['PCrho']
        Ap=mtd['PCAp']
        Pp=mtd['PCPp']
        PolynomialCorrection(InputWorkspace="PCalpha",OutputWorkspace="PCalpha",Coefficients=self.alphacoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCrho",OutputWorkspace="PCrho",Coefficients=self.rhocoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCAp",OutputWorkspace="PCAp",Coefficients=self.Apcoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCPp",OutputWorkspace="PCPp",Coefficients=self.Ppcoefficients,Operation="Multiply")
        D=Pp*(1.0+rho)
        nIp=(Ip*(rho*Pp+1.0)+Ia*(Pp-1.0))/D
        nIa=(Ip*(rho*Pp-1.0)+Ia*(Pp+1.0))/D
        RenameWorkspace(nIp,OutputWorkspace=str(Ip)+"corr")
        RenameWorkspace(nIa,OutputWorkspace=str(Ia)+"corr")
        iwksp=mtd.getObjectNames()
        list_n=[str(Ip),str(Ia),"PCalpha","PCrho","PCAp","PCPp","1_p"]
        for i in range(len(iwksp)):
            for j in list_n:
                lname=len(j)
                if iwksp[i] [0:lname+1] == j+"_":
                    DeleteWorkspace(iwksp[i])
        DeleteWorkspace("PCalpha")
        DeleteWorkspace("PCrho")
        DeleteWorkspace("PCAp")
        DeleteWorkspace("PCPp")
        DeleteWorkspace("D")
    
    def correctPA(self, UpUpWksp, UpDownWksp, DownUpWksp, DownDownWksp):
        print "Applying polarisation corrections "+str(self.name)
        Ipp = mtd[UpUpWksp]
        Ipa = mtd[UpDownWksp]
        Iap = mtd[DownUpWksp]
        Iaa = mtd[DownDownWksp]
        CloneWorkspace(Ipp,OutputWorkspace="PCalpha")
        CropWorkspace(InputWorkspace="PCalpha",OutputWorkspace="PCalpha",StartWorkspaceIndex="0",EndWorkspaceIndex="0")
        PCalpha=(mtd['PCalpha']*0.0)+1.0
        alpha=mtd['PCalpha']
        
        CloneWorkspace("PCalpha",OutputWorkspace="PCrho")
        CloneWorkspace("PCalpha",OutputWorkspace="PCAp")
        CloneWorkspace("PCalpha",OutputWorkspace="PCPp")
        rho=mtd['PCrho']
        Ap=mtd['PCAp']
        Pp=mtd['PCPp']
    
        # Use the polynomial corretion fn instead
        PolynomialCorrection(InputWorkspace="PCalpha",OutputWorkspace="PCalpha",Coefficients=self.alphacoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCrho",OutputWorkspace="PCrho",Coefficients=self.rhocoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCAp",OutputWorkspace="PCAp",Coefficients=self.Apcoefficients,Operation="Multiply")
        PolynomialCorrection(InputWorkspace="PCPp",OutputWorkspace="PCPp",Coefficients=self.Ppcoefficients,Operation="Multiply")
    
        A0 = (Iaa * Pp * Ap) + (Ap * Ipa * rho * Pp) + (Ap * Iap * Pp * alpha) + (Ipp * Ap * alpha * rho * Pp)
        A1 = Pp * Iaa
        A2 = Pp * Iap
        A3 = Ap * Iaa
        A4 = Ap * Ipa
        A5 = Ap * alpha * Ipp
        A6 = Ap * alpha * Iap
        A7 = Pp * rho  * Ipp
        A8 = Pp * rho  * Ipa
        D = Pp * Ap *( 1.0 + rho + alpha + (rho * alpha) )  
        nIpp = (A0 - A1 + A2 - A3 + A4 + A5 - A6 + A7 - A8 + Ipp + Iaa - Ipa - Iap) / D
        nIaa = (A0 + A1 - A2 + A3 - A4 - A5 + A6 - A7 + A8 + Ipp + Iaa - Ipa - Iap) / D
        nIpa = (A0 - A1 + A2 + A3 - A4 - A5 + A6 + A7 - A8 - Ipp - Iaa + Ipa + Iap) / D
        nIap = (A0 + A1 - A2 - A3 + A4 + A5 - A6 - A7 + A8 - Ipp - Iaa + Ipa + Iap) / D
        RenameWorkspace(nIpp,OutputWorkspace=str(Ipp)+"corr")
        RenameWorkspace(nIpa,OutputWorkspace=str(Ipa)+"corr")
        RenameWorkspace(nIap,OutputWorkspace=str(Iap)+"corr")
        RenameWorkspace(nIaa,OutputWorkspace=str(Iaa)+"corr")
        ReplaceSpecialValues(str(Ipp)+"corr",OutputWorkspace=str(Ipp)+"corr",NaNValue="0.0",NaNError="0.0",InfinityValue="0.0",InfinityError="0.0")
        ReplaceSpecialValues(str(Ipp)+"corr",OutputWorkspace=str(Ipp)+"corr",NaNValue="0.0",NaNError="0.0",InfinityValue="0.0",InfinityError="0.0")
        ReplaceSpecialValues(str(Ipp)+"corr",OutputWorkspace=str(Ipp)+"corr",NaNValue="0.0",NaNError="0.0",InfinityValue="0.0",InfinityError="0.0")
        ReplaceSpecialValues(str(Ipp)+"corr",OutputWorkspace=str(Ipp)+"corr",NaNValue="0.0",NaNError="0.0",InfinityValue="0.0",InfinityError="0.0")
        iwksp=mtd.getObjectNames()
        list_n=[str(Ipp),str(Ipa),str(Iap),str(Iaa),"PCalpha","PCrho","PCAp","PCPp","1_p"]
        for i in range(len(iwksp)):
            for j in list_n:
                lname=len(j)
                if iwksp[i] [0:lname+1] == j+"_":
                    DeleteWorkspace(iwksp[i])
        DeleteWorkspace("PCalpha")
        DeleteWorkspace("PCrho")
        DeleteWorkspace("PCAp")
        DeleteWorkspace("PCPp")
        DeleteWorkspace('A0')
        DeleteWorkspace('A1')
        DeleteWorkspace('A2')
        DeleteWorkspace('A3')
        DeleteWorkspace('A4')
        DeleteWorkspace('A5')
        DeleteWorkspace('A6')
        DeleteWorkspace('A7')
        DeleteWorkspace('A8')
        DeleteWorkspace('D')
    

   
nocorr = PolarisationCalibration('nocorr',  rhocoefficients=[1.0],   alphacoefficients=[1.0], Apcoefficients=[1.0], Ppcoefficients=[1.0], runs = [00000], comment = "all are = 1.0")

    
oldest = PolarisationCalibration('oldest',rhocoefficients=[1.006831,-0.011467,0.002244,-0.000095], alphacoefficients=[1.017526,-0.017183,0.003136,-0.000140], Apcoefficients=[0.917940,0.038265,-0.006645,0.000282], Ppcoefficients=[0.972762,0.001828,-0.000261,0.0], runs = [1850,1855,1851,1856], comment = "Constants Based on Runs 18350+18355 and 18351+18356 analyser theta at -0.1deg 2 RF Flippers as the polarising system")   
   
december2012 = PolarisationCalibration('december2012',  rhocoefficients=[0.970257,0.016127,-0.002318,0.000090],   alphacoefficients=[0.975722,0.012464,-0.002408,0.000105], Apcoefficients=[1.030894,-0.040847,0.006069,-0.000247], Ppcoefficients=[0.961900,-0.003722,0.001094,-0.000057], runs = [i for i in range(19438,19460)], comment = "Drabkin on incident side RF Flipper on analyser side Dec 2012")

march2016 = PolarisationCalibration('march2016', rhocoefficients=[1.00002594,-0.004658,-0.000336,0.000023], alphacoefficients=[1.01626,-0.006871,0.000509,-0.00002], Apcoefficients=[0.971155,-0.00469,0.000601,-0.000071], Ppcoefficients=[1.133720,-0.152533,0.050698,-0.008468,0.000699,-0.000022], runs=[i for i in range(37660,37709)], comment= "RF system removed, Drabkin incident between MagA and B and double Vs between C and D, analyser at nominal -0.1")

february2013_2 = PolarisationCalibration('february2013_2', rhocoefficients=[0.945927,0.025421,-0.003647,0.000156], alphacoefficients=[0.940769,0.027250,-0.003848,0.000164],  Apcoefficients=[0.974374,-0.005334,0.001313,-0.000115], Ppcoefficients=[1.023141,-0.024548,0.003398,-0.000134], runs =  [i for i in range(19628,19671)], comment = "analyser -0.2deg RF Flippers on polariser and analyser side Feb 2013")

february2013_1 = PolarisationCalibration('february2013_1', rhocoefficients=[0.955384,0.021501,-0.002962,0.000112], alphacoefficients=[0.944022,0.025639,-0.0035,0.000145],  Apcoefficients=[0.980605,-0.002201,-0.00097,0.000014], Ppcoefficients=[0.992212,-0.009161,0.000601,0.000025], runs =  [i for i in range(19628,19671)], comment = "analyser -0.1deg RF Flippers on polariser and analyser side Feb 2013")

november2015 = PolarisationCalibration('november2015', rhocoefficients=[0.944514,0.026548,-0.003837,0.000168], alphacoefficients=[0.944022,0.025639,-0.0035,0.000145], Apcoefficients=[0.980605,-0.002201,-0.00097,0.000014], Ppcoefficients=[0.992212,-0.009161,0.000601,0.000025], runs =  [i for i in range(36680,36695)], comment = "analyser -0.1deg RF Flippers on polariser and analyser side 0.5MHz November 2015")

november2015_2 = PolarisationCalibration('november2015_2', rhocoefficients=[0.702221,0.229436,-0.067308,0.009353,-0.000619,0.000016], alphacoefficients=[0.697757,0.230859,-0.067212,0.009276,-0.000610,0.000015], Apcoefficients=[1.065455,-0.067522,0.017810,-0.002430,0.000147,-3.298886e-6], Ppcoefficients=[1.030972,-0.046076,0.014934,-0.002600,0.000226,-7.240769e-06], runs =  [i for i in range(36680,36695)], comment = "analyser -0.1deg RF Flippers on polariser and analyser side 0.5MHz November 2015")


february2016 = PolarisationCalibration('february2016', rhocoefficients=[0.998515,-0.001878,-0.000902,0.000057], alphacoefficients=[1.020642,-0.008133,0.000511,-0.000008346674], Apcoefficients=[0.959164,0.000919,-0.000169,-0.00004], Ppcoefficients=[-0.470963,2.055378,-1.193303,0.363635,-0.0063189,0.006286,-0.000333], runs=[i for i in range(37660,37709)], comment= "RF system removed, Drabkin incident between MagA and B and double Vs between C and D, analyser at nominal -0.1")


alongbeamspecial2016 = PolarisationCalibration('alongbeamspecial2016', rhocoefficients=[1.00002594,-0.004658,-0.000336,0.000023], alphacoefficients=[0.999856,0.001013,-0.000645,0.000017], Apcoefficients=[-0.276464,1.16719,-0.440508,0.084601,-0.008725,0.000458153,-0.000009597066], Ppcoefficients=[1.133720,-0.152533,0.050698,-0.008468,0.000699,-0.000022], runs=[i for i in range(37660,37709)], comment= "RF system removed, Drabkin incident between MagA and B and double Vs between C and D, analyser at nominal -0.1. Guide field along beam at sample position. This measurement is done at max field ~110mT. Kept incident coefficients from march2016 (better stats) and adjusted scattered beam side.")


# 7th order fit [1.510208,0.653148,-0.374399,0.110249,-0.018274,0.001717,-0.000085,0.000001703707



PolarisationCalibration.set_current_calibration('march2016')
    
    
    
def load_eventmode(run, outputworkspace="wtemp"):
    LoadISISNexus(str(run),OutputWorkspace=outputworkspace,LoadMonitors=1)
    Rebin(outputworkspace,'5.0,20.0,100000.0',PreserveEvents=False,OutputWorkspace=outputworkspace+'reb')
    Rebin(outputworkspace+'_monitors','5.0,20.0,100000.0',OutputWorkspace=outputworkspace+'monreb')
    ConjoinWorkspaces(outputworkspace+'monreb',outputworkspace+'reb',CheckOverlapping=True)
    RenameWorkspace(outputworkspace+'monreb',OutputWorkspace=outputworkspace)
    DeleteWorkspace(outputworkspace+'_monitors')

def load_livedata_once(run, outputworkspace):
    StartLiveData(thisinstrument,AccumulationMethod="Replace",UpdateEvery=0.0,OutputWorkspace=outputworkspace)
            
    
def addRuns(runlist,wname):
  output=str(wname)
  for run in runlist:
    run = int(run)
    runis, run = currentRun.get_currentRun(run)
    if runis == "old":
        try:
            LoadISISNexus(str(run),OutputWorkspace="wtemp",LoadMonitors="Include")
        except:
            load_eventmode(run, outputworkspace="wtemp")
    elif runis == "current":
        # Live data doesn't return an event mode data set.. only histograms
        load_livedata_once(run, outputworkspace="wtemp")
    elif runis == "DoesNotExistYet":
        print "This is not a valid Run. I will ignore it."
        return "Run Does Not exist"
    # why are we doing this? need to find out    
    if isinstance(mtd["wtemp"], WorkspaceGroup):
        for k in mtd["wtemp"].getNames():
            mtd[k].setYUnit('Counts')
    else:
        mtd["wtemp"].setYUnit('Counts')
    try:
        Plus(output,"wtemp",OutputWorkspace=output)
    except:
        CloneWorkspace("wtemp", OutputWorkspace=output)
    DeleteWorkspace("wtemp")

def parseRunList(istring):
    if len(istring) >0: 
        s1=istring.split(',')
        rlist1=[]
        for i in range(len(s1)):
            tstr=s1[i].strip()
            if len(tstr) > 0:
                rlist1.append(tstr)
        rlist=[]
        for i in range(len(rlist1)):
            rlist2=[]
            if rlist1[i].find('+') >= 0:
                tstr=rlist1[i].split('+')
                for j in range(len(tstr)):
                    if tstr[j].find(':') >=0 and tstr[j].find('-') >=0:
                        tstr[j].strip()
                        tstr2=tstr[j].split('-')
                        tstr3=tstr2[1].split(':')
                        r1=range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1]))
                        for k in r1:
                            rlist2.append(str(k))
                    elif tstr[j].find('-') >=0:
                        tstr[j].strip()
                        tstr2=tstr[j].split('-')
                        r1=range(int(tstr2[0]),int(tstr2[1])+1)
                        for k in r1:
                            rlist2.append(str(k))
                    else:
                        rlist2.append(tstr[j])
            else:
                if rlist1[i].find(':') >=0 and rlist1[i].find('-')>=0:
                    rlist1[i].strip()
                    tstr2=rlist1[i].split('-')
                    tstr3=tstr2[1].split(':')
                    r1=range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1]))
                    for k in r1:
                        rlist2.append(str(k))
                elif rlist1[i].find('-')>=0:
                    rlist1[i].strip()
                    tstr2=rlist1[i].split('-')
                    r1=range(int(tstr2[0]),int(tstr2[1])+1)
                    for k in r1:
                        rlist2.append(str(k))
                else:
                    rlist2.append(rlist1[i])
            rlist.append(rlist2)
    return rlist


def parseNameList(istring):
    s1=istring.split(',')
    namelist=[]
    for i in range(len(s1)):
        tstr=s1[i].strip()
        namelist.append(tstr)   
    return namelist

def floodnorm(wkspName,floodfile='',floodopt='default'):
   #
   # pixel by pixel efficiency correction for the linear detector
   #

   if len(floodfile)==0:
        if not floodopt:
            return
        else:
            flood_file = current_detector.floodfile
   else:
      flood_file = floodfile
   
   #print "using flood normalisation file "+flood_file
    
   flood_wksp = current_detector.nickname + "_flood"
   if  flood_wksp not in mtd:
       LoadNexusProcessed(Filename=flood_file,OutputWorkspace=flood_wksp)
       if current_detector.idf:
        print "change idf"
        LoadInstrument(flood_wksp, InstrumentName=thisinstrument, RewriteSpectraMap=False)
       ConvertUnits(flood_wksp, Target='Wavelength', OutputWorkspace =  flood_wksp)

   Divide(LHSWorkspace=wkspName, RHSWorkspace=flood_wksp, OutputWorkspace=wkspName)

def nrNRFn(runList,nameList,incidentAngles,DBList,specChan,minSpec,maxSpec,gparams,floodfile='',subbgd=False,btmbgd=current_detector.btm_background,topbgd=current_detector.top_background,qgroup=False,dofloodnorm=True,usewkspname=False,sf=1.0,diagnostics=False):
    nlist=parseNameList(nameList)
    #logger.notice("This is the sample nameslist:"+str(nlist))
    if not usewkspname:
        rlist=parseRunList(runList)
    else:
        rlist=parseNameList(runList)
        #print rlist
    logger.notice("This is the sample runlist:"+str(rlist))
    dblist=parseNameList(DBList)
    #logger.notice("This is the Direct Beam nameslist:"+str(dblist))
    incAngles=parseNameList(incidentAngles)
    #logger.notice("This incident Angles are:"+str(incAngles))

    if not usewkspname:
        for i in range(len(rlist)):
            try:
                addRuns(rlist[i],nlist[i])
            except:
                break
    mon_spec=int(gparams[3])-1
    reb=gparams[0]+","+gparams[1]+","+gparams[2]
    
    k=0
    for i in nlist:
      try:  
        if isinstance(mtd[i], WorkspaceGroup):
            snames=mtd[i].getNames()
            Plus(LHSWorkspace=i+"_1",RHSWorkspace=i+"_2",OutputWorkspace="wtemp")
            if len(snames) > 2:
                for j in range(2,len(snames)-1):
                    Plus(LHSWorkspace="wtemp",RHSWorkspace=snames[j],OutputWorkspace="wtemp")
            for j in snames:
                if not diagnostics:
                   DeleteWorkspace(j)
            RenameWorkspace(InputWorkspace="wtemp",OutputWorkspace=i)
        ConvertUnits(InputWorkspace=i,OutputWorkspace=i,Target="Wavelength",AlignBins=True)
        #print i
        a1=mtd[i]
        nspec=a1.getNumberHistograms()
        if current_detector.__class__.__name__ == "PointDetector":
            Rebin(InputWorkspace=i,OutputWorkspace=i,Params=reb)
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"mon",StartWorkspaceIndex=mon_spec,EndWorkspaceIndex=mon_spec)
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"det",StartWorkspaceIndex=3,EndWorkspaceIndex=3)
            RotateInstrumentComponent(i+"det","DetectorBench",X="-1.0",Angle=str(2.0*float(incAngles[k])))
            Divide(LHSWorkspace=i+"det",RHSWorkspace=i+"mon",OutputWorkspace=i+"norm")
            if dblist[k] != "none":
                Divide(LHSWorkspace=i+"norm",RHSWorkspace=dblist[k],OutputWorkspace=i+"norm")
                ReplaceSpecialValues(InputWorkspace=i+"norm",OutputWorkspace=i+"norm",NanValue=0.0,NaNError=0.0,InfinityValue=0.0,InfinityError=0.0)
            ConvertUnits(InputWorkspace=i+"norm",OutputWorkspace=i+"RvQ",Target="MomentumTransfer")
        else:
            minSp=int(minSpec)
            maxSp=int(maxSpec)
            current_detector.croptodetector(i)
            if dofloodnorm:
                floodnorm(i+"det",floodfile,floodopt=dofloodnorm)
            twotheta=2.0*float(incAngles[k])+atan((-float(specChan)+current_detector.specular)*current_detector.pixelsize/current_detector.detectorposition)*rad2deg 
            RotateInstrumentComponent(i+"det","DetectorBench",X="-1.0",Angle=str(twotheta))
            #print str(2.0*float(incAngles[k])+atan((-float(specChan)+current_detector.specular)*current_detector.pixelsize/current_detector.detectorposition)*rad2deg)
            if qgroup:
                print "Summing detectors in Q."
                CropWorkspace(InputWorkspace=i+"det",OutputWorkspace=i+"detQ",XMin = pulsestart, XMax= pulseend, StartWorkspaceIndex=int(minSpec)-current_detector.minspec,EndWorkspaceIndex=int(maxSpec)-current_detector.minspec)
                ConvertUnits(InputWorkspace=i+"detQ",OutputWorkspace=i+"detQ",Target="MomentumTransfer",AlignBins=True)
                GroupDetectors(i+"detQ",OutputWorkspace=i+"sum",WorkspaceIndexList=range(int(maxSpec)-int(minSpec)+1),KeepUngroupedSpectra=False)
                ConvertUnits(InputWorkspace=i+"sum",OutputWorkspace=i+"sum",Target="Wavelength",AlignBins=True)

                if not diagnostics:
                    DeleteWorkspace(i+"detQ")
            else:
                print "Summing detectors in lambda."
                GroupDetectors(i+"det",OutputWorkspace=i+"sum",WorkspaceIndexList=range(int(minSpec)-current_detector.minspec,int(maxSpec)-current_detector.minspec+1),KeepUngroupedSpectra=False)
                
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"mon",StartWorkspaceIndex=mon_spec,EndWorkspaceIndex=mon_spec)
            Rebin(InputWorkspace=i+"mon",OutputWorkspace=i+"mon",Params=reb)
            RebinToWorkspace(WorkspaceToRebin=i+"sum",WorkspaceToMatch=i+"mon",OutputWorkspace=i+"sum")
            RebinToWorkspace(WorkspaceToRebin=i+"det",WorkspaceToMatch=i+"mon",OutputWorkspace=i+"det")
            if subbgd:
                print "Applying very light background correction. This is not necessarily all(or what) you want to do."
                # Calculate a background correction
                GroupDetectors(i+"det",OutputWorkspace=i+"bgd2",WorkspaceIndexList=range(*btmbgd),KeepUngroupedSpectra="0")
                GroupDetectors(i+"det",OutputWorkspace=i+"bgd1",WorkspaceIndexList=range(*topbgd),KeepUngroupedSpectra="0")
                Plus(i+"bgd1",i+"bgd2",OutputWorkspace=i+"bgd")
                wbgdtemp=mtd[i+"bgd"]/((btmbgd[1]-btmbgd[0])+(topbgd[1]-topbgd[0]))
                DeleteWorkspace(i+"bgd1")
                DeleteWorkspace(i+"bgd2")
                DeleteWorkspace(i+"bgd")
                # Subract a per spectrum background
                Minus(i+"det",wbgdtemp,OutputWorkspace=i+"det")
                if not diagnostics:
                    DeleteWorkspace("wbgdtemp")
                
            Divide(LHSWorkspace=i+"sum",RHSWorkspace=i+"mon",OutputWorkspace=i+"norm")
            Divide(LHSWorkspace=i+"det",RHSWorkspace=i+"mon",OutputWorkspace=i+"detnorm")
            if dblist[k]  == "none":
                a1=0
            elif dblist[k] == "function":
                # polynomial + power law corrections based on Run numbers 8291 and 8292
                Divide(LHSWorkspace=i+'norm',RHSWorkspace=i+'norm',OutputWorkspace=i+'normt1')
                PolynomialCorrection(InputWorkspace=i+'normt1',OutputWorkspace=i+'normPC',Coefficients='-0.0177398,0.00101695,0.0',Operation='Multiply')
                PowerLawCorrection(InputWorkspace=i+'normt1',OutputWorkspace=i+'normPLC',C0='2.01332',C1='-1.8188')
                Plus(LHSWorkspace=i+'normPC',RHSWorkspace=i+'normPLC',OutputWorkspace=i+'normt1')
                Divide(LHSWorkspace=i+'norm',RHSWorkspace=i+'normt1',OutputWorkspace=i+'norm')
                ReplaceSpecialValues(InputWorkspace=i+'norm',OutputWorkspace=i+'norm',NanValue=0.0,NaNError=0.0,InfinityValue=0.0,InfinityError=0.0)
                DeleteWorkspace(i+'normPC')
                DeleteWorkspace(i+'normPLC')
                DeleteWorkspace(i+'normt1')
            else:
                Divide(LHSWorkspace=i+"norm",RHSWorkspace=dblist[k],OutputWorkspace=i+"norm")
                ReplaceSpecialValues(InputWorkspace=i+"norm",OutputWorkspace=i+"norm",NanValue=0.0,NaNError=0.0,InfinityValue=0.0,InfinityError=0.0)
                Divide(LHSWorkspace=i+"detnorm",RHSWorkspace=dblist[k],OutputWorkspace=i+"detnorm")
                ReplaceSpecialValues(InputWorkspace=i+"detnorm",OutputWorkspace=i+"detnorm",NanValue=0.0,NaNError=0.0,InfinityValue=0.0,InfinityError=0.0)
            ConvertUnits(InputWorkspace=i+"norm",OutputWorkspace=i+"RvQ",Target="MomentumTransfer")
            DeleteWorkspace(i+"sum")

            
        CreateSingleValuedWorkspace(sf,0.0,OutputWorkspace='sf')
        Multiply(i+"RvQ",'sf',OutputWorkspace=i+'RvQ')
        
        k=k+1
        if not usewkspname:
            DeleteWorkspace(i)
        if not diagnostics:
            DeleteWorkspace(i+"mon")
            DeleteWorkspace(i+"det")
      except KeyError:
        print str(i)+"Does Not exist"  

def numberofbins(wksp):
    a1=mtd[wksp]
    y1=a1.readY(0)
    return len(y1)-1
#
#===========================================================
#
def maskbin(wksp,val):
    a1=mtd[wksp]
    x1=a1.readX(0)
    for i in range(len(x1)-1):
        if x1[i] > val:
            break
    a1.dataY(0)[i-1]=0.0
    a1.dataE(0)[i-1]=0.0
#
#===========================================================
#
def arr2list(iarray):
    # convert array of strings to a single string with commas
    res=""
    for i in range(len(iarray)-1):
        res=res+iarray[i]+","
    res=res+iarray[len(iarray)-1]
    return res
#
#===========================================================
#
def findbin(wksp,val):
    a1=mtd[wksp]
    x1=a1.readX(0)
    bnum=-1
    for i in range(len(x1)-1):
        if x1[i] > val:
            break
    return i-1
#
def nrDBFn(runListShort,nameListShort,runListLong,nameListLong,nameListComb,minSpec,maxSpec,minWavelength,gparams,floodfile="",dofloodnorm=True,fitspline=0,diagnostics=False):
    nlistS=parseNameList(nameListShort)
    rlistS=parseRunList(runListShort)
    nlistL=parseNameList(nameListLong)
    rlistL=parseRunList(runListLong)
    nlistComb=parseNameList(nameListComb)

    for i in range(len(rlistS)):
        addRuns(rlistS[i],nlistS[i])
    for i in range(len(rlistL)):
        addRuns(rlistL[i],nlistL[i])
    
    mon_spec=int(gparams[3])-1
    minSp=int(minSpec)-1
    maxSp=int(maxSpec)-1
    reb=gparams[0]+","+gparams[1]+","+gparams[2]
    
    for i in nlistS:
        ConvertUnits(InputWorkspace=i,OutputWorkspace=i,Target="Wavelength",AlignBins="1")
        Rebin(InputWorkspace=i,OutputWorkspace=i,Params=reb)
        CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"mon",StartWorkspaceIndex=mon_spec,EndWorkspaceIndex=mon_spec)
        if isinstance(mtd[i], WorkspaceGroup):
            snames=mtd[i].getNames()
            a1=mtd[snames[0]]
        else:
            a1=mtd[i]
        
        nspec=a1.getNumberHistograms()
        
        if nspec == 4:
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"det",StartWorkspaceIndex=3,EndWorkspaceIndex=3)
            Divide(i+"det",i+"mon",OutputWorkspace=i+"norm")
            ReplaceSpecialValues(i+"norm","0.0","0.0","0.0","0.0",OutputWorkspace=i+"norm")
        else:
            #CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"det",StartWorkspaceIndex=current_detector.minspec-1,EndWorkspaceIndex=current_detector.maxspec-1) #old ld start=4 end=243, WorkSpaceIndex = Spectrum-1 
            current_detector.croptodetector(i)
            print "short wavelengths: "+i+"det"
            floodnorm(i+"det",floodfile,floodopt=dofloodnorm)
            GroupDetectors(i+"det",OutputWorkspace=i+"sum",WorkspaceIndexList=range(int(minSpec)-current_detector.minspec,int(maxSpec)-current_detector.minspec+1),KeepUngroupedSpectra="0")
            Divide(i+"sum",i+"mon",OutputWorkspace=i+"norm")
            ReplaceSpecialValues(i+"norm","0.0","0.0","0.0","0.0", OutputWorkspace=i+"norm")
            
    for i in nlistL:
        ConvertUnits(InputWorkspace=i,OutputWorkspace=i,Target="Wavelength",AlignBins="1")
        Rebin(InputWorkspace=i,OutputWorkspace=i,Params=reb)
        CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"mon",StartWorkspaceIndex=mon_spec,EndWorkspaceIndex=mon_spec)
        if isinstance(mtd[i], WorkspaceGroup):
            lnames=mtd[i].getNames()
            a1=mtd[lnames[0]]
        else:
            a1=mtd[i]

        nspec=a1.getNumberHistograms()
        # duplicates above, needs work
        if nspec == 4:
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"det",StartWorkspaceIndex=3,EndWorkspaceIndex=3)
            Divide(i+"det",i+"mon",OutputWorkspace=i+"norm")
            ReplaceSpecialValues(i+"norm","0.0","0.0","0.0","0.0",OutputWorkspace=i+"norm")
        else: 
            #CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"det",StartWorkspaceIndex=current_detector.minspec-1,EndWorkspaceIndex=current_detector.maxspec-1) 
            current_detector.croptodetector(i)
            print "long wavelengths: "+i+"det"
            floodnorm(i+"det",floodfile,floodopt=dofloodnorm)
            GroupDetectors(i+"det",OutputWorkspace=i+"sum",WorkspaceIndexList=range(int(minSpec)-current_detector.minspec,int(maxSpec)-current_detector.minspec+1),KeepUngroupedSpectra="0") #+1 for range function
            Divide(i+"sum",i+"mon",OutputWorkspace=i+"norm")
            ReplaceSpecialValues(i+"norm","0.0","0.0","0.0","0.0",OutputWorkspace=i+"norm")
        
    for i in range(len(nlistS)):
        if isinstance(mtd[nlistS[i]+"norm"], WorkspaceGroup):
            snames=mtd[nlistS[i]+"norm"].getNames()
            lnames=mtd[nlistL[i]+"norm"].getNames()
            for k in range(len(snames)):
                Integration(snames[k],minWavelength,gparams[2], OutputWorkspace=snames[k]+"int")
                Integration(lnames[k],minWavelength,gparams[2], OutputWorkspace=lnames[k]+"int")
                Multiply(snames[k],lnames[k]+"int",OutputWorkspace=snames[k])
                Divide(snames[k],snames[k]+"int",OutputWorkspace=snames[k])
                a1=findbin(lnames[k],float(minWavelength))
                MultiplyRange(lnames[k],"0",str(a1),"0.0",OutputWorkspace=lnames[k])
                WeightedMean(snames[k],lnames[k],OutputWorkspace=nlistComb[i]+'_'+str(k+1))
                if(fitspline > 0):
                    SmoothData(nlistComb[i]+'_'+str(k+1),5,OutputWorkspace=nlistComb[i]+'_'+str(k+1))
                    SplineBackground(InputWorkspace=nlistComb[i]+'_'+str(k+1),NCoeff=fitspline,OutputWorkspace=nlistComb[i]+'_spline_'+str(k+1))
                if not diagnostics:
                    DeleteWorkspace(snames[k]+"int")
                    DeleteWorkspace(lnames[k]+"int")
        else:
            Integration(nlistS[i]+"norm",minWavelength,gparams[2],OutputWorkspace=nlistS[i]+"int")
            Integration(nlistL[i]+"norm",minWavelength,gparams[2],OutputWorkspace=nlistL[i]+"int")
            Multiply(nlistS[i]+"norm",nlistL[i]+"int",OutputWorkspace=nlistS[i]+"norm")
            Divide(nlistS[i]+"norm",nlistS[i]+"int",OutputWorkspace=nlistS[i]+"norm")
            a1=findbin(nlistL[i]+"norm",float(minWavelength))
            MultiplyRange(nlistL[i]+"norm","0",str(a1),"0.0", OutputWorkspace=nlistL[i]+"norm")
            WeightedMean(nlistS[i]+"norm",nlistL[i]+"norm", OutputWorkspace=nlistComb[i])
            if(fitspline > 0):
                SmoothData(nlistComb[i],5,OutputWorkspace=nlistComb[i])
                SplineBackground(InputWorkspace=nlistComb[i],NCoeff=fitspline,OutputWorkspace=nlistComb[i]+'_spline')
            if not diagnostics:
                DeleteWorkspace(nlistS[i]+"int")
                DeleteWorkspace(nlistL[i]+"int")

            if not diagnostics:
                DeleteWorkspace(nlistS[i]+"mon")
                DeleteWorkspace(nlistS[i]+"det")
                if nspec != 4:
                    DeleteWorkspace(nlistS[i]+"sum")
                DeleteWorkspace(nlistS[i]+"norm")
                DeleteWorkspace(nlistS[i])
                DeleteWorkspace(nlistL[i]+"mon")
                DeleteWorkspace(nlistL[i]+"det")
                if nspec != 4:
                    DeleteWorkspace(nlistL[i]+"sum")
                DeleteWorkspace(nlistL[i]+"norm")
                DeleteWorkspace(nlistL[i])




def NRCombineDatafn(RunsNameList,CombNameList,applySFs,SFList,SFError,scaleOption,bparams,globalSF,applyGlobalSF,diagnostics=False):
    qmin=bparams[0]
    bin=bparams[1]
    qmax=bparams[2]
    rlist=parseNameList(RunsNameList)
    listsfs=parseNameList(SFList)
    listsfserr=parseNameList(SFError)
    sfs=[]
    sferrs=[]
    for i in rlist:
        Rebin(i,qmin+","+bin+","+qmax,OutputWorkspace=i+"reb")
    # find the overlap ranges
    bol=[] #beginning of overlaps
    eol=[] #end of overlaps
    for i in range(len(rlist)-1):
        a1=mtd[rlist[i+1]]
        x=a1.readX(0)
        bol.append(x[0])
        a1=mtd[rlist[i]]
        x=a1.readX(0)
        eol.append(x[len(x)-1])
    # set the edges of the rebinned data to 0.0 to avoid partial bin problems
    maskbin(rlist[0]+"reb",eol[0])
    if len(rlist) > 2:
        for i in range(1,len(rlist)-1):
            maskbin(rlist[i]+"reb",bol[i-1])
            maskbin(rlist[i]+"reb",eol[i])
    maskbin(rlist[len(rlist)-1]+"reb",bol[len(rlist)-2])
    # Now find the various scale factors and store in temp workspaces
    for i in range(len(rlist)-1):
        Integration(rlist[i]+"reb",str(bol[i]),str(eol[i]), OutputWorkspace="i"+str(i)+"1temp")
        Integration(rlist[i+1]+"reb",str(bol[i]),str(eol[i]), OutputWorkspace="i"+str(i)+"2temp")
        if scaleOption != "2":
            Divide("i"+str(i)+"1temp","i"+str(i)+"2temp",OutputWorkspace="sf"+str(i))
            a1=mtd["sf"+str(i)]
            print "sf"+str(i)+"="+str(a1.readY(0))+" +/- "+str(a1.readE(0))
            sfs.append(str(a1.readY(0)[0]))
            sferrs.append(str(a1.readE(0)[0]))
        else:
            Divide("i"+str(i)+"2temp","i"+str(i)+"1temp",OutputWorkspace="sf"+str(i))
            a1=mtd["sf"+str(i)]
            print "sf"+str(i)+"="+str(a1.readY(0))+" +/- "+str(a1.readE(0))
            sfs.append(str(a1.readY(0)[0]))
            sferrs.append(str(a1.readE(0)[0]))
        DeleteWorkspace("i"+str(i)+"1temp")
        DeleteWorkspace("i"+str(i)+"2temp")
    # if applying pre-defined scale factors substitute the given values now 
    # Note the errors are now set to 0
    if applySFs == "2":
        for i in range(len(rlist)-1):
            a1=mtd["sf"+str(i)]
            a1.dataY(0)[0]=float(listsfs[i])
            a1.dataE(0)[0]=float(listsfserr[i])
            #a1.setY(0, float(listsfs[i]))
            #a1.setE(0, float(listsfserr[i]))
    # Now scale the various data sets in the correct order
    if scaleOption != "2":
        for i in range(len(rlist)-1):
            for j in range(i+1,len(rlist)):
                Multiply(rlist[j]+"reb","sf"+str(i),OutputWorkspace=rlist[j]+"reb")
    else:
        for i in range(len(rlist)-1,0,-1):
            for j in range(i,0,-1):
               
                Multiply(rlist[j-1]+"reb","sf"+str(i-1),OutputWorkspace=rlist[j-1]+"reb")

    WeightedMean(rlist[0]+"reb",rlist[1]+"reb",OutputWorkspace="currentSum")
    if len(rlist) > 2:
        for i in range(2,len(rlist)):
            WeightedMean("currentSum",rlist[i]+"reb",OutputWorkspace="currentSum")
    
    # if applying a global scale factor do it here
    if applyGlobalSF == "2":
        scaledData=mtd['currentSum']/float(globalSF)
        RenameWorkspace('scaledData',OutputWorkspace=CombNameList)
        DeleteWorkspace('currentSum')
    else:
        RenameWorkspace('currentSum',OutputWorkspace=CombNameList)
    for i in range(len(rlist)-1):
        DeleteWorkspace("sf"+str(i))
    if not diagnostics:
        for i in range(len(rlist)):
            DeleteWorkspace(rlist[i]+"reb")
            DeleteWorkspace(rlist[i])
    return [arr2list(sfs),arr2list(sferrs)] 
    
def nrPNRFn(runList,nameList,incidentAngles,DBList,specChan,minSpec,maxSpec,gparams,floodfile,PNRwithPA,pnums,doCorrs=True,doLDCorrs="0",subbgd=False,qgroup=False,calibration=PolarisationCalibration.current,usewkspname=False,dofloodnorm=True,diagnostics=False):
    nlist=parseNameList(nameList)
    logger.notice("This is the sample nameslist:"+str(nlist))
    if not usewkspname:
        rlist=parseRunList(runList)
    else:
        rlist=parseNameList(runList)
    logger.notice("This is the sample runlist:"+str(rlist))
    dblist=parseNameList(DBList)
    incAngles=parseNameList(incidentAngles)
    
    if PNRwithPA:
        nper=4
        logger.notice("PNRwithPA = "+str(PNRwithPA))
    else:
        nper=2
    
    if not usewkspname:
        for i in range(len(rlist)):
            addRuns(rlist[i],nlist[i])
    
    mon_spec=int(gparams[3])-1
    minSp=int(minSpec)
    maxSp=int(maxSpec)
    reb=gparams[0]+","+gparams[1]+","+gparams[2]
    
    k=0
    for i in nlist:
        a1=mtd[i+"_1"]
        nspec=a1.getNumberHistograms()
        
        wksp=i
        ConvertUnits(InputWorkspace=wksp,OutputWorkspace=wksp,Target="Wavelength",AlignBins=True)
        
        if current_detector.__class__.__name__ == "PointDetector":
            Rebin(InputWorkspace=wksp,OutputWorkspace=wksp,Params=reb)
            Rebin(InputWorkspace=wksp+'mon',OutputWorkspace=wksp+'mon',Params=reb)
            CropWorkspace(InputWorkspace=wksp,OutputWorkspace=wksp+"det",StartWorkspaceIndex=3,EndWorkspaceIndex=3)
            RotateInstrumentComponent(wksp+"det","DetectorBench",X="-1.0",Angle=str(2.0*float(incAngles[k])))
            Divide(LHSWorkspace=wksp+"det",RHSWorkspace=wksp+"mon",OutputWorkspace=wksp+"norm")
            if dblist[k] != "none":
                Divide(LHSWorkspace=wksp+"norm",RHSWorkspace=dblist[k],OutputWorkspace=wksp+"norm")
                ReplaceSpecialValues(wksp+"norm","0.0","0.0","0.0","0.0",OutputWorkspace=wksp+"norm")
                ConvertUnits(wksp+"norm",Target="MomentumTransfer",OutputWorkspace=wksp+"RvQ")
            if not diagnostics:
                DeleteWorkspace(wksp+"mon")
                DeleteWorkspace(wksp+"det")
        else:
            current_detector.croptodetector(wksp)
            if dofloodnorm:
                floodnorm(wksp+"det",floodfile,floodopt=dofloodnorm)
            a1=2.0*(float(incAngles[k]))+(atan((float(specChan)-current_detector.specular)*current_detector.pixelsize/current_detector.detectorposition)*rad2deg)
            RotateInstrumentComponent(wksp+"det","DetectorBench",X="-1.0",Angle=str(a1))
            
            # Experimental convert to Q before summing 
            if qgroup:
                print "Summing detectors in Q."
                CropWorkspace(InputWorkspace=i+"det",OutputWorkspace=i+"detQ",XMin=pulsestart, XMax=pulseend, StartWorkspaceIndex=int(minSpec)-current_detector.minspec,EndWorkspaceIndex=int(maxSpec)-current_detector.minspec)
                ConvertUnits(InputWorkspace=i+"detQ",OutputWorkspace=i+"detQ",Target="MomentumTransfer", AlignBins = True)
                GroupDetectors(i+"detQ",OutputWorkspace=i+"sum",WorkspaceIndexList=range(int(maxSpec)-int(minSpec)+1),KeepUngroupedSpectra=False)
                ConvertUnits(InputWorkspace=i+"sum",OutputWorkspace=i+"sum",Target="Wavelength",AlignBins=True)
                if not diagnostics:
                    DeleteWorkspace(i+"detQ")
            
            else:
                print "Summing detectors in lambda"
                GroupDetectors(wksp+"det",OutputWorkspace=wksp+"sum",WorkspaceIndexList=range(int(minSpec)-current_detector.minspec,int(maxSpec)-current_detector.minspec+1),KeepUngroupedSpectra=False)
            
            CropWorkspace(InputWorkspace=i,OutputWorkspace=i+"mon",StartWorkspaceIndex=mon_spec,EndWorkspaceIndex=mon_spec)
            Rebin(InputWorkspace=i+"mon",OutputWorkspace=i+"mon",Params=reb)
            RebinToWorkspace(WorkspaceToRebin=wksp+"sum",WorkspaceToMatch=wksp+"mon",OutputWorkspace=wksp+"sum")
            RebinToWorkspace(WorkspaceToRebin=wksp+"det",WorkspaceToMatch=wksp+"mon",OutputWorkspace=wksp+"det")
             
            if subbgd:
                print "Applying very light background correction. This is not necessarily all (or what) you want to do."
                CloneWorkspace(i, OutputWorkspace="bgdtemp")
                ConvertUnits(InputWorkspace="bgdtemp",OutputWorkspace="bgdtemp",Target="Wavelength",AlignBins=True)
                Rebin(InputWorkspace="bgdtemp",OutputWorkspace="bgdtemp",Params=reb)
                current_detector.croptodetector("bgdtemp")
                Plus("bgdtemp"+"_"+pnums[0],"bgdtemp"+"_"+pnums[1],OutputWorkspace="wbgdsum")
                for j in range(2,nper):
                    Plus("wbgdsum","bgdtemp"+"_"+pnums[j],OutputWorkspace="wbgdsum")
                GroupDetectors("wbgdsum",OutputWorkspace="bgd2",WorkspaceIndexList=range(*current_detector.btm_background),KeepUngroupedSpectra=False)
                GroupDetectors("wbgdsum",OutputWorkspace="bgd1",WorkspaceIndexList=range(*current_detector.top_background),KeepUngroupedSpectra=False)
                Plus("bgd1","bgd2",OutputWorkspace="bgd")
                wbgdtemp=mtd["bgd"]/(current_detector.nbackgroundspectra*nper)
                DeleteWorkspace("bgdtemp")
                DeleteWorkspace("wbgdsum")
                DeleteWorkspace("bgd1")
                DeleteWorkspace("bgd2")
                DeleteWorkspace("bgd")
               
                # Subract a per spectrum background
                Minus(wksp+"det",wbgdtemp,OutputWorkspace=wksp+"det")
                ResetNegatives(InputWorkspace=wksp+"det",OutputWorkspace=wksp+"det",AddMinimum='0',ResetValue="0.0")
                numberofpixelstosum = int(maxSpec)-int(minSpec)+1
                Minus(wksp+"sum",wbgdtemp*numberofpixelstosum,OutputWorkspace=wksp+"sum")
                ResetNegatives(InputWorkspace=wksp+"sum",OutputWorkspace=wksp+"sum",AddMinimum='0',ResetValue="0.0")
                     
            Divide(LHSWorkspace=wksp+"sum",RHSWorkspace=wksp+"mon",OutputWorkspace=wksp+"norm")
            Divide(LHSWorkspace=wksp+"det",RHSWorkspace=wksp+"mon",OutputWorkspace=wksp+"detnorm")
            if dblist[k]  != "none":
                Divide(LHSWorkspace=wksp+"norm",RHSWorkspace=dblist[k],OutputWorkspace=wksp+"norm")
                ReplaceSpecialValues(wksp+"norm","0.0","0.0","0.0","0.0",OutputWorkspace=wksp+"norm")
                Divide(LHSWorkspace=wksp+"detnorm",RHSWorkspace=dblist[k],OutputWorkspace=wksp+"detnorm")
                ReplaceSpecialValues(wksp+"detnorm","0.0","0.0","0.0","0.0",OutputWorkspace=wksp+"detnorm")
            ConvertUnits(wksp+"norm",OutputWorkspace=wksp+"RvQ",Target="MomentumTransfer")
            if not diagnostics:
                DeleteWorkspace(wksp+"sum")
                DeleteWorkspace(wksp+"mon")
                DeleteWorkspace(wksp+"det")
        if doCorrs:
            if nper == 2:
                calibration.correctPNR(i+"norm_"+pnums[0],i+"norm_"+pnums[1])
                for j in range(2):
                    RenameWorkspace(InputWorkspace=i+"norm"+"_"+pnums[j]+"corr",OutputWorkspace=i+"normcorr"+"_"+pnums[j])
                GroupWorkspaces(InputWorkspaces=i+"normcorr_"+pnums[0]+","+i+"normcorr_"+pnums[1],OutputWorkspace=i+"normcorr")
                ConvertUnits(InputWorkspace=i+"normcorr",OutputWorkspace=i+"normcorrRvQ",Target="MomentumTransfer")
                if (nspec > 4 and doLDCorrs != "0"):
                    calibration.correctPNR(i+"detnorm_"+pnums[0],i+"detnorm_"+pnums[1],calibration=calibration)
                    for j in range(4):
                        RenameWorkspace(InputWorkspace=i+"detnorm"+"_"+pnums[j]+"corr",OutputWorkspace=i+"detnormcorr"+"_"+pnums[j])
                    GroupWorkspaces(InputWorkspaces=i+"detnormcorr_"+pnums[0]+","+i+"detnormcorr_"+pnums[1],OutputWorkspace=i+"detnormcorr")
            else:
                calibration.correctPA(i+"norm"+"_"+pnums[0],i+"norm"+"_"+pnums[1],i+"norm"+"_"+pnums[2],i+"norm"+"_"+pnums[3])
                for j in range(4):
                    RenameWorkspace(InputWorkspace=i+"norm"+"_"+pnums[j]+"corr",OutputWorkspace=i+"normcorr"+"_"+pnums[j])
                GroupWorkspaces(InputWorkspaces=i+"normcorr_"+pnums[0]+","+i+"normcorr_"+pnums[1]+","+i+"normcorr_"+pnums[2]+","+i+"normcorr_"+pnums[3]+"",OutputWorkspace=i+"normcorr")
                ConvertUnits(InputWorkspace=i+"normcorr",OutputWorkspace=i+"normcorrRvQ",Target="MomentumTransfer")
                if (nspec > 4 and doLDCorrs != "0"):
                    calibration.correctPA(i+"detnorm_"+pnums[0],i+"detnorm_"+pnums[1],i+"detnorm_"+pnums[2],i+"detnorm_"+pnums[3])
                    for j in range(4):
                        RenameWorkspace(InputWorkspace=i+"detnorm"+"_"+pnums[j]+"corr",OutputWorkspace=i+"detnormcorr"+"_"+pnums[j])
                    GroupWorkspaces(InputWorkspaces=i+"detnormcorr_"+pnums[0]+","+i+"detnormcorr_"+pnums[1]+","+i+"detnormcorr_"+pnums[2]+","+i+"detnormcorr_"+pnums[3]+"",OutputWorkspace=i+"detnormcorr")
            if (diagnostics == 0 and doCorrs != "0"):
                DeleteWorkspace(i+"norm")
                DeleteWorkspace(i+'normcorr')
            #if (diagnostics == 0 and doLDCorrs != "0"):
            #    DeleteWorkspace(i+"detnorm")
        k=k+1
        if not usewkspname:
            DeleteWorkspace(i)
        if subbgd:
            DeleteWorkspace("wbgdtemp")
#
#===========================================================
#
    
        
if __name__ == '__main__':

    print socket.getfqdn()
    print check_is_host_at_isis()
