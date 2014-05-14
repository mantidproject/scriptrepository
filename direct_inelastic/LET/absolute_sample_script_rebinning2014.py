"""      
LET TRANSIENT REDUCTION SCRIPT MARCH 2014;
"""

# program to crunch down event mode from LET and produce output SPE files. 
# Program can automatically find all incident energies in rep rate mode and write out spe files in the # form of LET'run no: +ei'.spe
from qtiGenie import *
iliad_setup('let')
##############################################################################################################
# mulitple PC convenience section. parameters from here can be usually set up from GUI
##############################################################################################################
#config['defaultsave.directory'] = 'd:/Data/isis/Let/2013_12'   
# map mask and cal file, again the values from Mantid, data search directories can be modified here
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    
# map mask and cal file, again the values from Mantid, data search directories can be modified here
config.appendDataSearchDir('/home/let/Desktop/LET_maps') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/isisdatar55/ndxlet/Instrument/data/cycle_14_1') 
config.appendDataSearchDir('/isisdatar55/ndxlet/Instrument/data/cycle_13_5') 
config.appendDataSearchDir('/isisdatar55/ndxlet/Instrument/data/cycle_13_4') 

##############################################################################################################
# USERS SECTION -- particular run parameters
##############################################################################################################
white_run = 15961   # enter whitebeam run number here
run_no=[15987] #[17314]  #event mode run numbers here or use next line for a continuous sequence of runs i.e range(first run, last run +1)
ei=[1.6,2.8,5.9,19.9]#[8]          # incident energies you want analyzed
ebin=[-4,0.002,0.8]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
#   Other positional  parameters
# key-coded parameters

params={}
params['norm_method']='current'
# External background reduction should go out of here -- there is internal bkg reduction procedure 
remove_background = False  #if true then will subtract a flat background in time from the time range given below otherwise put False
bg_range=[92000,98000] # range of times to take background in
params['bkgd_range']=[bg_range[0],bg_range[1]] # these two are currently not used as we remove background externaly but can be used 
# we use external background removal. One needs to compare it with internal and move to the internal one. Meanwhile, disable internal
params['check_background']=False;
##############################################################################################################
# ABSOLUTE UNITS RELATED
##############################################################################################################
MonoVanRun=None  #  vanadium run in the same configuration as your sample used for absolute normalization (None for no absolute normalization)
params['sample_mass']=7.85     
params['sample_rmm']=50.9415 # 50.9415 The rmm of Vanadium. Put correct value for sample here
params['monovan_mapfile']='LET_rings_141.map'  
##############################################################################################################
# ILLIAD MULTIREP RUN: Instrument scientist specified parameters
##############################################################################################################
# Instrument scientist specified parameters
# map file to combine instrument spectra
mapping ='LET_rings_141'  # ring map file is used for powder.  if absent idf file value is used instead
params['det_cal_file']='det_LET_cycle141.dat'  #det_cal_file must be specified if the reduction sends out put to a workpsace
#params['det_cal_file']='det_corrected7.nxs' # ASCII correction file provides different results on different OS for LET. Nexus solves this proble,
params['hardmaskOnly']='hard_2014_1.msk'   # diag does not work well on LET. At present only use a hard mask RIB has created
# does not work in event mode TODO: investigate
params['diag_bleed_test']=False;
# this parameter  need carefull checking and fine tunning for an instrument with guides
params['detector_van_range']=[2,7]
params['norm_method']='current'

# loads the whitebeam (or rather the long monovan ). Does it as a raw file to save time as the event mode is very large
#white beam
loadFreshWB=True;
if 'wb_wksp' in mtd:
    wb_wksp=mtd['wb_wksp']
    if wb_wksp.getRunNumber()==white_run:
        loadFreshWB = False;        
#only load whitebeam if not already there
if loadFreshWB:  
    wb_wksp=LoadRaw(Filename=str(white_run),OutputWorkspace="wb_wksp",LoadMonitors='Exclude')    
    #iliad_reducer().det_cal_file = 'det_LET_cycle133.dat'
    #wb_wksp = iliad_reducer().load_data('LET000'+str(wb)+'.raw','wb_wksp')
    #iliad_reducer().det_cal_file = wb_wksp;
    #LoadRaw(Filename=str(wb),OutputWorkspace="wb_wksp") # load whitebeam
                        
for sample_run in run_no:     
    fname='LET000'+str(sample_run)+'.nxs'
    print ' processing file ', fname
    #print ' Det Cal file:  ', iliad_reducer().det_cal_file;
    #w1 = iliad_reducer().load_data(run,'w1')
    w1 =Load(Filename=fname,OutputWorkspace='w1',LoadMonitors='1');

    start = time.clock()   

    
    if remove_background:
        bg_ws_name=find_background('w1',bg_range);

     #############################################################################################
     # this section finds all the transmitted incident energies if you have not provided them
    if len(ei) == 0: 
        # does not currently work properly due to generig loader problem in Mabtid. Needs changes in qtiGenie and Mantid fixes
        ei = find_chopper_peaks('w1_monitors');
        print 'Found energies: ',ei
    print 'Energies reduced are:'
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

        ws_name = 'LET{0}_absEi{1:5.2f}'.format(sample_run,energy);
        RenameWorkspace(InputWorkspace=out,OutputWorkspace=ws_name);
        SaveNXSPE(InputWorkspace=ws_name,Filename=ws_name+'.nxspe');

elapsed = (time.clock() - start)	
print "Time to complete reduction : {0: 5.2f} min ".format(elapsed/60)

