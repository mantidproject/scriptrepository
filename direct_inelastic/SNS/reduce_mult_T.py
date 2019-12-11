"""
Reduce temperature scan. Going as far in events as possible, split, then extract useful information
"""
import numpy as np

def GetEiT0(ws_name,EiGuess):
    try:
        alg=GetEi(InputWorkspace=ws_name,Monitor1Spec="1",Monitor2Spec="2",EnergyEstimate=float(EiGuess))				#Run GetEi algorithm
        Ei=alg[0]
        Tzero=-alg[3]					#Extract incident energy and T0
    except:
        Ei='N/A'
        Tzero='N/A'
    return [Ei,Tzero]

def quick_process(IWS,OWS,Erange,Efixed,T0):
    """
    normalize by current, correct for Tube efficiency, compensate for ki/kf, divide by bin width
    """
    ChangeBinOffset(InputWorkspace=IWS,OutputWorkspace="__OWS",Offset=T0)
    ConvertUnits(InputWorkspace="__OWS",OutputWorkspace="__OWS",Target="Wavelength",EMode="Direct",EFixed=Efixed)	#The algorithm for He3 tube efficiency requires wavelength units
    He3TubeEfficiency(InputWorkspace="__OWS",OutputWorkspace="__OWS")												#Apply correction due to absorption in He3
    ConvertUnits(InputWorkspace="__OWS",OutputWorkspace="__OWS",Target="DeltaE",EMode="Direct",EFixed=Efixed)		#Switch  to energy transfer
    CorrectKiKf(InputWorkspace="__OWS",OutputWorkspace="__OWS")
    GroupDetectors(InputWorkspace="__OWS",OutputWorkspace=OWS,MapFile=r'/SNS/users/19g/granroth/MCCL/p5deg_group.xml',Behaviour='Average')

    #Apply k_i/k_f factor
    #Rebin(InputWorkspace="__OWS",OutputWorkspace="__OWS",Params=Erange,PreserveEvents=True)



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


run=31276
dir='/SNS/SEQ/IPTS-8807/data/'
fname='SEQ_'+str(run)+'_event.nxs'
ws_name='SEQ_'+str(run)+'_event'
mon_name=ws_name+'_monitors'
Erange1='-200.0,5,490.0'
Erange2='-20.0,0.5,32.0'
outdir='/SNS/SEQ/IPTS-8807/shared/'
Vanadium="/SNS/SEQ/shared/2012_B/V_files/SEQ_30675_event.nxs"
Norm=V_norm_obj(Vanadium,"0.3,0.9,1.2",outdir,ld_saved_fl=True)
Norm.MaskBTP(Bank="38,75,76,99,100,101,102,114,115,120")
Norm.MaskBTP(Pixel="1,2,3,4,5,6,7,8,121,122,123,124,125,126,127,128")
Norm.MaskBTP(Bank="74",Tube="8")
Norm.CreateMasksAndVanadiumNormalization()
Load(Filename=fname,OutputWorkspace=ws_name,LoadMonitors='1')
ChangeBinOffset(InputWorkspace=ws_name,OutputWorkspace=ws_name,Offset=500,IndexMin=54272,IndexMax=55295) # adjust time for pack C17 wired backward
[Ei1,Tzero1]=GetEiT0(mon_name,500.0)
[Ei2,Tzero2]=GetEiT0(mon_name,36.0)
GenerateEventsFilter(InputWorkspace=ws_name,OutputWorkspace='split_ws',SplittersInformationWorkspace='Sp_Info',LogName='SampleTemp',MinimumLogValue='276',MaximumLogValue='288',LogValueInterval='2')
quick_process(ws_name,'E_500',Erange1,Ei1,Tzero1)
quick_process(ws_name,'E_37',Erange2,Ei2,Tzero2)
OWSBN='E_500_ev_'
FilterEvents(InputWorkspace='E_500',OutputWorkspaceBaseName=OWSBN,SplittersInformationWorkspace='Sp_Info',InputSplittersWorkspace='split_ws')
h_Sp_Info=mtd['Sp_Info']
Sp_index=np.array(h_Sp_Info.column(0))-1
ky_lst=list(mantid.keys())
for idx in Sp_index[1:]:
   WS_name=OWSBN+'_'+str(idx)
   skip_flag=False
   try:
	   ky_lst.index(WS_name)	   
   except:
     skip_flag=True
     pass    
   if  not skip_flag:
	   NormaliseByCurrent(WS_name,OutputWorkspace=WS_name)
	   Rebin(InputWorkspace=WS_name,OutputWorkspace=WS_name+'_h',Params=Erange1,PreserveEvents=True)
	   SofQW3(InputWorkspace=WS_name+'_h',OutputWorkspace=WS_name+'_QW',QAxisBinning="0.1,0.1,15.0",Emode='Direct',Efixed=Ei1)
	   SumSpectra(InputWorkspace=WS_name+'_QW',OutputWorkspace=WS_name+'_sum')


