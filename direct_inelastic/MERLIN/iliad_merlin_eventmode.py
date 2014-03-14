def find_binning_range(energy,ebin,ls,lm2,mult,dt_DAE):
    """ function finds the binning range used in multirep mode 
        for merlin ls=11.8,lm2=10. mult=2.8868 dt_DAE=1;
        for LET    ls=25,lm2=23.5 mult=4.1     dt_DAE=1.6;
        all these values have to be already present in IDF and should be taken from there

        # THIS FUNCTION SHOULD BE MADE GENERIG AND MOVED OUT OF HERE
    """
    energy=float(energy)

    emin=(1.0-ebin[2])*energy   #minimum energy is with 80% energy loss
    lam=(81.81/energy)**0.5
    lam_max=(81.81/emin)**0.5
    tsam=252.82*lam*ls   #time at sample
    tmon2=252.82*lam*lm2 #time to monitor 6 on LET
    tmax=tsam+(252.82*lam_max*mult) #maximum time to measure inelastic signal to
    t_elastic=tsam+(252.82*lam*mult)   #maximum time of elastic signal
    tbin=[int(tmon2),1.6,int(tmax)]				
    energybin=[float("{0: 6.4f}".format(elem*energy)) for elem in ebin]

    return (energybin,tbin,t_elastic);

#--------------------------------------------------------------------------------------------------------
def find_background(ws_name,bg_range,dt_DAE):
    """ Function to find background from multirep event workspace
    dt_DAE = 1 for MERLIN and 1.6 for LET
     should be precalculated or taken from IDF

        # THIS FUNCTION SHOULD BE MADE GENERIC AND MOVED OUT OF HERE
    """
    delta=bg_range[1]-bg_range[0]
    Rebin(InputWorkspace='w1',OutputWorkspace='bg',Params=[bg_range[0],delta,bg_range[1]],PreserveEvents=False)	
    v=(delta)/dt_DAE
    CreateSingleValuedWorkspace(OutputWorkspace='d',DataValue=v)
    Divide(LHSWorkspace='bg',RHSWorkspace='d',OutputWorkspace='bg')

"""      
MERLIN TRANSIENT REDUCTION SCRIPT MARCH 2014;
"""
from qtiGenie import *
#from PySlice2 import *

# Calculates 

inst='MER'
iliad_setup(inst)
##############################################################################################################
# mulitple PC convenience section. parameters from here can be usually set up from GUI
##############################################################################################################
config['defaultsave.directory']=r"d:\Data\Mantid_Testing\14_03_14" 
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    save_dir = config.getString('defaultsave.directory')
    
#print "Data will be saved into: ",save_dir
# map mask and cal file, again the values from Mantid, data search directories can be modified here
#config.appendDataSearchDir('/home/merlin/mprogs/InstrumentFiles/merlin') 
config.appendDataSearchDir(r'c:\Users\wkc26243\Documents\work\Libisis\InstrumentFiles\merlin');
# data (raw or nxs) run files -- values from data search directories can be modified here
#config.appendDataSearchDir('/isisdatar55/NDXMERLIN/Instrument/data/cycle_13_3') 
#config.appendDataSearchDir('/archive/NDXMERLIN/Instrument/data/cycle_13_3') 
config.appendDataSearchDir(r'd:\Data\Mantid_Testing\14_03_14') 

##############################################################################################################
# USERS SECTION -- particular run parameters
##############################################################################################################
white_run = 17186
run_no=[17335] #[17314]           #event mode run numbers here or use next line for a continous sequence of runs i.e range(first run, last run +1)
ei=[10,17,35,103]        # incident energies you want analysed
ebin=[-0.3,0.005,0.950]  #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
#   Other positional  parameters
# key-coded parameters
params={}
params['norm_method']='current'



# should go out of here -- there is internal bkg reduction procedure 
remove_background = False  #if true then will subtract a flat background in time from the time range given below otherwise put False
bg_range=[92000,98000] # range of times to take background in
params['bkgd_range']=[bg_range[0],bg_range[1]] # these two are currently not used as we remove background externaly but can be used 
params['check_background']=False;
##############################################################################################################
# ABSOLUTE UNITS RELATED
##############################################################################################################
MonoVanRun=17335  #  vanadium run in the same configuration as your sample used for absolute notmalization (None for no absolute nomalization)
params['sample_mass']=7.85     
params['sample_rmm']=50.9415 # 50.9415 The rmm of Vanadium. Put correct value for sample here
params['monovan_mapfile']='rings_125.map'  
##############################################################################################################
# ILLIAD MULTIREP RUN: Instrument scientist specified parameters
##############################################################################################################
# Instrument scientist specified parameters
# map file to combime instrument spectra
mapping ='rings_125.map'  # ring map file is used for powder.  if absend idf file value is used instead
params['det_cal_file']='det_corr_125.dat'  #det_cal_file must be specified if the reduction sends out put to a workpsace
#params['det_cal_file']='det_corrected7.nxs' # ASCII correction file provides different results on different OS for LET. Nexus solves this proble,
params['hardmaskPlus']='hard_mask_sp.msk'
#params['hardmaskOnly']=mask_file   # diag does not work well on LET. At present only use a hard mask RIB has created
# New vanadium mass. This should move to IDF
params['vanadium-mass']=7.85
# does not work in event mode TODO: investigate
params['diag_bleed_test']=False;

# White beam
if 'wb_wksp' in mtd:
        wb_wksp=mtd['wb_wksp']
else:  #only load whitebeam if not already there
    wb_wksp=LoadRaw(Filename='MER00'+str(white_run),OutputWorkspace="wb_wksp",LoadMonitors='Exclude')
    

for sample_run in run_no:
      
    fname='MER00'+str(sample_run)+'.nxs'
    print ' processing file ', fname
    #w1 = dgreduce.getReducer().load_data(run,'w1')
    Load(Filename=fname,OutputWorkspace='w1',LoadMonitors='1');

    
    if remove_background:
        find_background('w1',bg_range,1);

     #############################################################################################
     # this section finds all the transmitted incident energies if you have not provided them
     #if len(ei) == 0:  -- not tested here -- should be unit test for that. 
           #ei = find_chopper_peaks('w1_monitors');       
    print 'Energies transmitted are:'
    print (ei)

    RenameWorkspace(InputWorkspace = 'w1',OutputWorkspace='w1_storage');
    RenameWorkspace(InputWorkspace = 'w1_monitors',OutputWorkspace='w1_mon_storage');
                    
    #now loop around all energies for the run
    for ind,energy in enumerate(ei):
        print "Reducting around energy: {0}".format(float(energy))
        (energybin,tbin,t_elastic) = find_binning_range(energy,ebin,11.8,10,2.8868,1);
        print " Rebinning will be performed in the range: ",energybin
        # if we calculate more then one energy, initial workspace will be used more then once
        if ind <len(ei)-1:
              CloneWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1')
              CloneWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors')
        else:
              RenameWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1');
              RenameWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors');

        if remove_background:
              w1=Rebin(InputWorkspace='w1',OutputWorkspace='w1',Params=tbin,PreserveEvents=False)            
              Minus(LHSWorkspace='w1',RHSWorkspace='bg',OutputWorkspace='w1')
               

       ######################################################################
       # ensure correct round-off procedure
        params['monovan_integr_range']=[round(ebin[0]*energy,4),round(ebin[2]*energy,4)]; # integration range of the vanadium 
        # range to check detector's efficiency using white beam vanadium . If wb is not in multirep mode, here should be optimal WB integration range, 
        params['wb_integr_range']   = [-0.3*45,1.2*45.];

       # absolute unit reduction -- if you provided MonoVan run or relative units if monoVan is not present
        out=iliad(wb_wksp,"w1",energy,energybin,mapping,MonoVanRun,**params)

        ws_name = 'MER00{0}_absEi{1:4.1f}'.format(sample_run,energy);
        RenameWorkspace(InputWorkspace=out,OutputWorkspace=ws_name);
        SaveNXSPE(InputWorkspace=ws_name,Filename=ws_name+'.nxspe');
