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
#config.appendDataSearchDir(r'd:\Data\Mantid_Testing\14_03_14') 

##############################################################################################################
# USERS SECTION -- particular run parameters
##############################################################################################################
white_run = 17186
run_no=[17335] #[17314]           #event mode run numbers here or use next line for a continuous sequence of runs i.e range(first run, last run +1)
ei=[10,17,35,103]        # incident energies you want analyzed
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
MonoVanRun=17335  #  vanadium run in the same configuration as your sample used for absolute normalization (None for no absolute normalization)
params['sample_mass']=7.85     
params['sample_rmm']=50.9415 # 50.9415 The rmm of Vanadium. Put correct value for sample here
# New vanadium mass. This should move to IDF
params['vanadium-mass']=7.85
params['monovan_mapfile']='rings_125.map'  
##############################################################################################################
# ILLIAD MULTIREP RUN: Instrument scientist specified parameters
##############################################################################################################
# Instrument scientist specified parameters
# map file to combine instrument spectra
mapping ='rings_125.map'  # ring map file is used for powder.  if absent idf file value is used instead
params['det_cal_file']='det_corr_125.dat'  #det_cal_file must be specified if the reduction sends out put to a workpsace
#params['det_cal_file']='det_corrected7.nxs' # ASCII correction file provides different results on different OS for LET. Nexus solves this proble,
params['hardmaskPlus']='hard_mask_sp.msk'
#params['hardmaskOnly']=mask_file   # diag does not work well on LET. At present only use a hard mask RIB has created
# does not work in event mode TODO: investigate
params['diag_bleed_test']=False;
# this parameter may need checking for an instrument with guides
params['detector_van_range']=[0.5,200]
params['norm_method']='current'



#white beam
loadFreshWB=True;
if 'wb_wksp' in mtd:
    wb_wksp=mtd['wb_wksp']
    if wb_wksp.getRunNumber()==wb:
        loadFreshWB = False;
        
#only load whitebeam if not already there
if loadFreshWB:  
    wb_wksp=LoadRaw(Filename=wbFile,OutputWorkspace="wb_wksp",LoadMonitors='Exclude')
    

for sample_run in run_no:
     
    fname='MER00'+str(sample_run)+'.nxs'
    print ' processing file ', fname
    #w1 = dgreduce.getReducer().load_data(run,'w1')
    Load(Filename=fname,OutputWorkspace='w1',LoadMonitors='1');
    

    
    if remove_background:
        bg_ws_name=find_background('w1',bg_range);

     #############################################################################################
     # this section finds all the transmitted incident energies if you have not provided them
    if len(ei) == 0:
        ei = find_chopper_peaks('w1_monitors');
        print ei
    print 'Energies transmitted are:'
    print (ei)

    RenameWorkspace(InputWorkspace = 'w1',OutputWorkspace='w1_storage');
    RenameWorkspace(InputWorkspace = 'w1_monitors',OutputWorkspace='w1_mon_storage');
                    
    #now loop around all energies for the run
    for ind,energy in enumerate(ei):
        print "Reducing around energy: {0}".format(float(energy))
        (energybin,tbin,t_elastic) = find_binning_range(energy,ebin);
        print " Rebinning will be performed in the range: ",energybin
        
        # if we calculate more then one energy, initial workspace will be used more then once. Is this enough for that?    
        #Rebin(InputWorkspace='w1',OutputWorkspace='w1reb',Params=tbin,PreserveEvents='1')                        
        # safe option would be:
        if ind <len(ei)-1:
              CloneWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1')
              CloneWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors')
        else:
              RenameWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1');
              RenameWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors');

        if remove_background:
              w1=Rebin(InputWorkspace='w1',OutputWorkspace='w1',Params=tbin,PreserveEvents=False)            
              Minus(LHSWorkspace='w1',RHSWorkspace=bg_ws_name,OutputWorkspace='w1')
               

       ######################################################################
        params['monovan_integr_range']=[energybin[0],energybin[2]]; # integration range of the absolute units vanadium 

       # absolute unit reduction -- if you provided MonoVan run or relative units if monoVan is not present
        out=iliad(wb_wksp,"w1",energy,energybin,mapping,MonoVanRun,**params)

        ws_name = 'MER00{0}_absEi{1:4.1f}'.format(sample_run,energy);
        RenameWorkspace(InputWorkspace=out,OutputWorkspace=ws_name);
        SaveNXSPE(InputWorkspace=ws_name,Filename=ws_name+'.nxspe');


