"""
Reduction including generation of an experiment log 
"""
#!/usr/bin/env python
from __future__ import print_function

import sys,os
sys.path.append("/opt/Mantid/bin")
from mantid.simpleapi import *
from numpy import *
from string import *
from MaskBTP import *

# Logs at: /var/log/SNS_applications/autoreduce.log

class ExperimentLog(object):
    def __init__(self):
        self.log_list=['vChTrans','Speed1','Phase1','Speed2','Phase2','Speed3','Phase3','EnergyRequest','s1t','s1r','s1l','s1b','vAttenuator2','vAttenuator1','svpressure','dvpressure']
        self.cols=[1,4,4,4,4,4,4,1,1,1,1,1,1,1,4,4,4]
        self.SERotOptions=['CCR13VRot','SEOCRot']
        self.SETempOptions=['SampleTemp']

    def log_line_gen(self,IWSName):
        """
        IWSName is a string of the workspace name
        """
        h=mtd[IWSName]
        outstr=str(h.getRunNumber())+','+str(h.getTitle()).replace(' ','_').replace(',','_')+','+ h.getRun()['start_time'].value+','+ h.getRun()['end_time'].value+','+str(h.getRun()['duration'].value)+','+str(h.getRun().getProtonCharge())+','
        for cn,idx in zip(self.cols,self.log_list):
            try:
                jk=h.getRun()[idx].getStatistics()
                valave=jk.mean                
                valmax=jk.maximum
                valmin=jk.minimum
                valstdev=jk.standard_deviation
                if cn==1:
                    outstr=outstr+'%s,'%(valave)
                else:
                    outstr=outstr+'%s,%s,%s,%s,'%(valave,valmax,valmin,valstdev)
            except:
                if cn==1:
                    outstr=outstr+'N/A,'
                else:
                    outstr=outstr+'N/A, N/A, N/A, N/A,'
  
        #check for sample environment temperature reading
        for SET in self.SETempOptions:
            if self.log_list.count(SET)==0:
                if h.getRun().hasProperty(SET):
                    SEVals=h.getRun().getProperty(SET).getStatistics()
                    outstr=outstr+"%s,%s,%s,%s,"%(SEVals.mean,SEVals.maximum,SEVals.minimum,SEVals.standard_deviation)
                    self.log_list.append(SET)
                    self.cols.append(4)
        #check sample environment rotation stage
        angle=0.
        for SE in self.SERotOptions:
            if self.log_list.count(SE)==1:
                jk=h.getRun()[SE].getStatistics()
                angle=jk.mean 
            else:
                if h.getRun().hasProperty(SE):
                    SEVals=h.getRun().getProperty(SE).getStatistics()
                    outstr=outstr+"%s,%s,%s,%s,"%(SEVals.mean,SEVals.maximum,SEVals.minimum,SEVals.standard_deviation)
                    self.log_list.append(SE)
                    self.cols.append(4)
                    angle=SEVals.mean
        return [outstr,angle]

    def save_line(self,fname,IWSName,Ei,T0):
        """
        fname is a file name
        IWSName is a string of the workspace name
        Ei,T0 are outputs from GetEiT0
        """
        # if the file does not exist, set a flga to add a header line
        create_h_flag=False
        try:
            h_ft=open(fname,'r')
        except IOError:
            create_h_flag=True
            pass
        else:   
            h_ft.close()   
        
        [str_data,angle]=self.log_line_gen(IWSName)
        str_data=str_data+'%s, %s,\n'%(Ei,T0)  

        h_f=open(fname,'a')  # append to the file unless it does not exist 
        if create_h_flag:
            strtmp='run_number,title,start_time,end_time,duration,proton_charge,'
            for idx,nc in zip(self.log_list,self.cols):
                if nc==1:
	                strtmp=strtmp+idx+','
                else:
                    strtmp=strtmp+'%s mean,%s maximum,%s minimum,%s stddev,'%(idx,idx,idx,idx)
            strtmp=strtmp+'Ei, T0,\n'
            h_f.write(strtmp)	   
           
        h_f.write(str_data)
        h_f.close()
        return angle



def GetEiT0(ws_name,EiGuess):
    try:
        alg=GetEi(InputWorkspace=ws_name,Monitor1Spec="1",Monitor2Spec="2",EnergyEstimate=float(EiGuess))				#Run GetEi algorithm
        Ei=alg[0]
        Tzero=-alg[3]					#Extract incident energy and T0
    except:
        Ei='N/A'
        Tzero='N/A'
    return [Ei,Tzero]


def autoloader(filename,outdir):
    """
    function for autoloader
    """

    __MonWS=LoadNexusMonitors(Filename=filename)
    Eguess=array(__MonWS.getRun()['EnergyRequest'].value).mean()
    [Efixed,T0]=GetEiT0("__MonWS",Eguess)
    elog=ExperimentLog()
    angle=elog.save_line(outdir+'experiment_summary_rered.csv','__MonWS',Efixed,T0)

    runnum=str(__MonWS.getRun()['run_number'].value)     
    Estep=Eguess*0.005
    Elow=-Eguess*0.5
    Ehigh=Eguess*0.95
    Erange='%g,%g,%g'%(Elow,Estep,Ehigh)   
    DeleteWorkspace('__MonWS')
 
    if Efixed!='N/A':
        LoadEventNexus(Filename=filename,OutputWorkspace="__IWS") #Load an event Nexus file
        #Fix that all time series log values start at the same time as the proton_charge
        r=mtd['__IWS'].getRun()
        for x in list(r.keys()):
        	if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
        		try:
        			ShiftTime('__IWS',x)
        		except:
    			    pass

	#Filter chopper 3 bad events
	valC3=r['Phase3'].getStatistics().median
	FilterByLogValue(InputWorkspace='__IWS',OutputWorkspace='__IWS',LogName='Phase3',MinimumValue=valC3-0.15,MaximumValue=valC3+0.15)
        #FilterBadPulses(InputWorkspace="__IWS",OutputWorkspace = "__IWS",LowerCutoff = 50)
        return [runnum,Efixed,T0,Erange,angle]
    else:
        #do not load data if we cannot process it
        return None 

def quick_process(IWS,Erange,Efixed,T0):
    """
    normalize by current, correct for Tube efficiency, compensate for ki/kf, divide by bin width
    """
    ChangeBinOffset(InputWorkspace=IWS,OutputWorkspace="__OWS",Offset=500,IndexMin=54272,IndexMax=55295) # adjust time for pack C17 wired backward
    ChangeBinOffset(InputWorkspace="__OWS",OutputWorkspace="__OWS",Offset=T0)
    ConvertUnits(InputWorkspace="__OWS",OutputWorkspace="__OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)	#The algorithm for He3 tube efficiency requires wavelength units
    He3TubeEfficiency(InputWorkspace="__OWS",OutputWorkspace="__OWS")												#Apply correction due to absorption in He3
    ConvertUnits(InputWorkspace="__OWS",OutputWorkspace="__OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)		#Switch  to energy transfer
    CorrectKiKf(InputWorkspace="__OWS",OutputWorkspace="__OWS")													    #Apply k_i/k_f factor
    Rebin(InputWorkspace="__OWS",OutputWorkspace="__OWS",Params=Erange,PreserveEvents=True)						#Make sure the bins are correct
#    ConvertToDistribution(Workspace="__OWS") 		                                                                #Divide by bin width
    
def WS_clean():
    DeleteWorkspace('__IWS')
    DeleteWorkspace('__OWS')
    DeleteWorkspace('__VAN')



class V_norm_obj(object):
    def __init__(self,vanfile,wlstr,outdir,maskfile='',ld_saved_fl=True):
        self.outdir=outdir
        self.vanfile=vanfile
        self.wlstr=wlstr
        self.ld_saved_fl=ld_saved_fl
        self.maskpars=[]
        self.maskfile=maskfile

    def MaskBTP(self,**kwargs):
        self.maskpars.append(kwargs)

    def CreateMasksAndVanadiumNormalization(self):
        if (os.path.isfile(self.outdir+"van.nx5") and (self.ld_saved_fl)):
            LoadNexus(Filename=outdir+"van.nx5",OutputWorkspace="__VAN")
        else:   
            LoadEventNexus(Filename=self.vanfile,OutputWorkspace="__VAN")
            ChangeBinOffset(InputWorkspace="__VAN",OutputWorkspace="__VAN",Offset=500,IndexMin=54272,IndexMax=55295) # adjust time for pack C17 wired backward
	    ConvertUnits(InputWorkspace="__VAN",OutputWorkspace="__VAN",Target="Wavelength",EMode="Elastic")
            Rebin(InputWorkspace="__VAN",OutputWorkspace="__VAN",Params=self.wlstr,PreserveEvents=False)			#integrate all events in the range given by wlstr
            ConvertToDistribution("__VAN")
            NormaliseByCurrent(InputWorkspace="__VAN",OutputWorkspace="__VAN")									#normalize by proton charge
            for d in self.maskpars:
                MaskBTP(Workspace="__VAN",**d)
            MedianDetectorTest(InputWorkspace="__VAN",OutputWorkspace="__MASK")			#determine which detectors to mask, and store them in the "MASK" workspace
	  
            if len(self.maskfile)>0:
                LoadNexus(Filename=self.maskfile,OutputWorkspace="__temp_mask")
                MaskDetectors(Workspace="__MASK",MaskedWorkspace="__temp_mask")		    #add detectors masked in "temp_mask" to "MASK"
                DeleteWorkspace(Workspace="__temp_mask")
            MaskDetectors(Workspace="__VAN",MaskedWorkspace="__MASK")												#Mask "VAN". This prevents dividing by 0		
            DeleteWorkspace(Workspace="__MASK")																	#Mask is carried by VAN workspace
            SaveNexus(InputWorkspace="__VAN",Filename=self.outdir+"van.nx5")


def ShiftTime(WName,lg_name):
	"""
	shift the time in a given log to match the time in the proton charge log"
	"""
	H_IN = mtd[WName]
	PC =  H_IN.getRun()['proton_charge'].firstTime()
	#print "P="+str(PC)+"\n"
	P =  H_IN.getRun()[lg_name].firstTime()
	#print "P="+str(P)+"\n"
	Tdiff = PC-P
	Tdiff_num = Tdiff.total_milliseconds()*1E-3
	#print "Tdiff="+str(Tdiff_num)+"\n"
	ChangeLogTime(InputWorkspace=WName, OutputWorkspace = WName, LogName = lg_name, TimeOffset = Tdiff_num)




if __name__ == "__main__":
    #check number of arguments
    if (len(sys.argv) != 3): 
        print("autoreduction code requires a filename and an output directory")
        sys.exit()
    if not(os.path.isfile(sys.argv[1])):
        print("data file ", sys.argv[1], " not found")
        sys.exit()
    else:
        filename = sys.argv[1]
        outdir = sys.argv[2]
    #----------------------------------------------------------------------------------
    # changes


    clean=True
    NXSPE_flag=True
    outpre="reduced"
    #Vanadium and masking    
    Vanadium="/SNS/SEQ/shared/2013_A/V_files/SEQ_31279_event.nxs"
    maskfile=''
    Norm=V_norm_obj(Vanadium,"0.3,0.9,1.2",outdir,maskfile=maskfile,ld_saved_fl=True)
    Norm.MaskBTP(Bank="70,99,100,101,102,110")
    Norm.MaskBTP(Pixel="1,2,3,4,5,6,7,8,121,122,123,124,125,126,127,128")
    Norm.MaskBTP(Bank="96",Tube="8")
    Norm.CreateMasksAndVanadiumNormalization()

    #end changes
    #-----------------------------------------------------------------------------------

    
    al=autoloader(filename,outdir) #this line also records the information in outdir+'experiment_summary.csv'

    if al!=None:
        [runnum,Efixed,T0,Erange,angle]=al   
        quick_process('__IWS',Erange,Efixed,T0)
	outfile=outpre+'_'+runnum
	# save nexus file for combining data later before V normalization
	AddSampleLog(Workspace="__OWS",LogName="psi",LogText=str(angle),LogType="Number")
	SaveNexus(InputWorkspace="__OWS", Filename= outdir+outfile+".nxs")
    #FilterBadPulses(InputWorkspace="__OWS",OutputWorkspace = "__OWS",LowerCutoff = 50)
    NormaliseByCurrent(InputWorkspace="__OWS",OutputWorkspace="__OWS")
    ConvertToPointData(InputWorkspace="__OWS",OutputWorkspace="__OWS") 
    ConvertToHistogram(InputWorkspace="__OWS",OutputWorkspace="__OWS") 
    ConvertToDistribution(Workspace="__OWS") 		                                                                #Divide by bin width
    Divide(LHSWorkspace="__OWS", RHSWorkspace="__VAN",OutputWorkspace="__OWS")
    if NXSPE_flag:            
        SaveNXSPE(InputWorkspace="__OWS", Filename= outdir+outfile+".nxspe",Efixed=Efixed,Psi=angle,KiOverKfScaling=True) 
    if clean:
        WS_clean()

    

