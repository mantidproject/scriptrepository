"""
The sample script for MERLIN absolute units reduction. Copied from LET system tests. 
"""
from __future__ import print_function
from qtiGenie import *
iliad_setup('MER')


# program to crunch down event mode from LET and produce output SPE files. Program can automatically find all incident energies in rep rate mode and write out spe files in the # form of LET'run no: +ei'.spe


#############################################
# where to save results (usually specified in Mantid, data search directories but can be redefined here (comment two following rows)
save_dir = '/home/merlin/users/Adroja/NewScriptFromAlex'
config['defaultsave.directory']=save_dir

# this is the user input section
wb=23012   # enter whitebeam run number here 
#run_no=[8570,8581] # event mode run numbers here or use next line for a continuous sequence of runs i.e range(first run, last run +1)
#run_no=range(23408,23409)
run_no=[23408, 23409] #
# how many files to sum together. 
num_files2sum=2; #1 -- nothing is summed. If >1, the number of runs above should be divisible by this number without an reminder
ei=[75.4, 20,9.08]
#ei=[191,80.1, 43.8]  #80meV 550Hz
#ei = [102,38.7,20.2, 12.3]  #, 100meV 350Hz
#ei=[5.8,15]           # incident energies you want analysed
ebin=[-0.8,0.002,0.95]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
mapping='rings_125.map'  # rings mapping file for powders
file = 'Bjorn_mask.msk'    # standard hard mask file  for MDRLIN
##############################################################
# Background
remove_background=False;
bg_range=[12000,19000] # range of times to take background in
############################################
MonoVanRun=None #23206 # None vanadium run in the same configuration as your sample
sampleMass=17.1;
sampleRMM=612.4109;
vanadium_mass=7.85;  #added by dta
monovan_mapfile='rings_125.map'
 
 
##############################################
# map mask and cal file, again the values from Mantid, data search directories can be modified here
config.appendDataSearchDir('/home/merlin/mprogs/InstrumentFiles/merlin') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/archive/ndxmerlin/Instrument/data/cycle_14_2') 

##########################################
present_run=-1
if 'wb_wksp' in mtd:
	m=mtd['wb_wksp']
	present_run= m.getRunNumber()
if present_run !=wb:  #only load whitebeam if not already there
	LoadRaw(Filename=str(wb),OutputWorkspace="wb_wksp") # load whitebeam

 
 
######################################################################
present_run=-1
for irun in range(0,len(run_no),num_files2sum):    #loop around runs
       run = run_no[irun];
       fname='MER0000'+str(run)+'.nxs'
       print(' processing file ', fname)
       #if 'w1' in mtd:  # Uncomment this to save time and not to load existing workspace twice
       #     m=mtd['w1']
       #     present_run= m.getRunNumber()
       #print run, present_run
       # Load single run or sum runs 
       if  present_run != run:
            w1, w1_monitors=Load(Filename=fname,LoadMonitors='1',MonitorsAsEvents='0')  
            w1_mon=RenameWorkspace(w1_monitors)
            for ns in range(num_files2sum-1):
                print(" Adding run N {0} to run N {1}".format(run_no[irun+ns+1],run_no[irun]), end=' ')
                run = run_no[irun+ns+1]
                ws,ws_monitors= Load(Filename='MER0000'+str(run)+'.nxs',SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='0')
                w1 += ws
                w1_mon += ws_monitors
                print("got {0} events".format(w1.getNumberEvents()))
                AddSampleLog(Workspace='w1',LogName='run_number_'+str(run),LogText=str(run),LogType='Number') #only have to do this as run_number is passed to workspace for event files       
              
       if remove_background:
                find_background('w1',bg_range);

        #############################################################################################
        # this section finds all the transmitted incident energies if you have not provided them
       if len(ei) == 0:  #-- not tested here -- should be unit test for that. 
            ei = find_chopper_peaks('w1_monitors');       
       print('Energies transmitted are:', end=' ')
       print (ei)

       RenameWorkspace(InputWorkspace = 'w1',OutputWorkspace='w1_storage');
       RenameWorkspace(InputWorkspace = 'w1_mon',OutputWorkspace='w1_mon_storage');
                    
       #now loop around all energies for the run
       result =[];
       for ind,energy in enumerate(ei):
                print(float(energy))
                (energybin,tbin,t_elastic) = find_binning_range(energy,ebin);
                print(" Rebinning will be performed in the range: ",energybin)
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
                argi={};
                
                argi['bleed']=False
                argi['norm_method']='current'
                argi['det_cal_file']='det_corr_125.dat'
                argi['detector_van_range']=[40,55]
                argi['bkgd_range']=[int(t_elastic),int(tbin[2])]
                argi['hardmaskOnly']=file		
                argi['check_background']=False;

                # abs units
                argi['sample_mass']=sampleMass;
                argi['sample_rmm'] =sampleRMM;
                argi['monovan_mapfile']=monovan_mapfile;
                argi['vanadium-mass']=vanadium_mass
                argi['monovan_integr_range']=[round(ebin[0]*energy,4),round(ebin[2]*energy,4)]; # integration range of the vanadium 
                #MonoVanWSName = None;

                # absolute unit reduction -- if you provided MonoVan run or relative units if monoVan is not present
                out=iliad("wb_wksp","w1",energy,energybin,mapping,MonoVanRun,**argi)

                ws_name = 'MER{0}_ei{1:2.1f}mev'.format(run,energy);
                #RenameWorkspace(InputWorkspace=out,OutputWorkspace=ws_name);
                fileName=ws_name+'_Rings_125abs.nxspe';
                print('Saving results into file {0}'.format(fileName))
                SaveNXSPE(InputWorkspace=out,Filename=fileName)		
