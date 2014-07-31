"""
The sample script for MERLIN absolute units reduction. 
Should be used for non_abs units too by setting MonoVanRun number to None
File copied from current (31/07.2014) Mantid System tests for LET and modified with MERLIN settings
"""
from qtiGenie import *
import datetime
iliad_setup('merlin')

# this is the user input section
#############################################
# where to save results (overwrites specified in Mantid, data search directories, if you want to use those one, comment two lines below here)
#save_dir = '/home/xfz39506/RB1320468/'
#config['defaultsave.directory']=save_dir

wb=23013   # enter whitebeam run number here 
#run_no=[8570,8581] # event mode run numbers here or use next line for a continuous sequence of runs i.e range(first run, last run +1)
#run_no=range(11376,11379)
run_no=[23047,23048] #
# this string allows to sum all or part of the files from above
num_files2sum=1; #1 -- nothing is summed the number of runs above should be divisible by this number without an reminder
#
ei = [40,22,14] #,1.46]# # incident energies you want analysed, 
#ei=[5.8,15]           
ebin=[-0.2,0.002,0.8]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
mapping='one2one_125.map'  # or provide rings mapping file for powders (see the map used for WB grouping)
file = 'Bjorn_mask.msk'    # standard hard mask file  for LET
##############################################################
# Background  -- should be used for powders only
remove_background=False;
bg_range=[12000,19000] # range of times to take background in
###################################################
# Absolute Units Related
###################################################
argi={}
MonoVanRun=None #10433 # vanadium run in the same configuration as your sample. If none, relative units will be used

argi['sample_mass']=1  #mass of your sample in g
argi['sample_rmm'] =1 #Enter the correct relative molecular mass of your sample
#argi['vanadium-mass']=7.85
argi['monovan_mapfile']='rings_125.map'

##############################################
# map mask and cal file, again the values from Mantid, data search directories can be modified here or in Mantid DataSearch directories
config.appendDataSearchDir('/home/merlin/mprogs/InstrumentFiles/merlin') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/archive/ndxmerlin/Instrument/data/cycle_14_2') 

config.appendDataSearchDir('/home/merlin/merlin_share')


# loads the whitebeam (or rather the long monovan ). Does it as a raw file to save time as the e/Bjorn_maskvent mode is very large
present_run=-1
if 'wb_wksp' in mtd:
    m=mtd['wb_wksp']
    present_run= m.getRunNumber()
if present_run !=wb:  #only load whitebeam if not already there
    LoadRaw(Filename=str(wb),OutputWorkspace="wb_wksp") # load whitebeam

 
 
######################################################################
start = time.clock()
present_run=-1
for irun in xrange(0,len(run_no),num_files2sum):    #loop around runs
       run = run_no[irun];
       fname='MER0000'+str(run)+'.nxs'
       print ' processing file ', fname
       if 'w1' in mtd:
            m=mtd['w1']
            present_run= m.getRunNumber()
       print run, present_run
       # Load single run or sum runs 
       if  present_run != run:
            w1, w1_monitors=Load(Filename=fname,LoadMonitors='1',MonitorsAsEvents='0')  
            w1_mon=RenameWorkspace(w1_monitors)
            for ns in xrange(num_files2sum-1):
                print " Adding run N {0} to run N {1}".format(runno[i+ns+1],runno[i]),
                run = runno[irun+ns+1]
                ws,ws_monitors= Load(Filename='MER0000'+str(run)+'.nxs',SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='0')		
                w1 += ws
                w1_mon += ws_monitors
                print "got {0} events".format(w1.getNumberEvents())
                AddSampleLog(Workspace='w1',LogName='run_number_'+str(run),LogText=str(run),LogType='Number') #only have to do this as run_number is passed to workspace for event files       
       
       
       if remove_background:
                find_background('w1',bg_range);

        #############################################################################################
        # this section finds all the transmitted incident energies if you have not provided them
       if len(ei) == 0:  #-- not tested here -- should be unit test for that. 
            ei = find_chopper_peaks('w1_monitors');       
       print 'Energies transmitted are:',
       print (ei)

       RenameWorkspace(InputWorkspace = 'w1',OutputWorkspace='w1_storage');
       RenameWorkspace(InputWorkspace = 'w1_mon',OutputWorkspace='w1_mon_storage');
                    
       #now loop around all energies for the run
       result =[];
       for ind,energy in enumerate(ei):
                print float(energy)
                (energybin,tbin,t_elastic) = find_binning_range(energy,ebin);
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
                
                argi['bleed']=False
                argi['norm_method']='current'
                argi['det_cal_file']='det_corr_125.dat'
                argi['detector_van_range']=[40,55]
                argi['bkgd_range']=[int(t_elastic),int(tbin[2])]
                argi['hardmaskOnly']=file		
                argi['check_background']=False;

                # abs units
                argi['monovan_integr_range']=[round(ebin[0]*energy,4),round(ebin[2]*energy,4)]; # integration range of the vanadium 

                #MonoVanWSName = None;

                # absolute unit reduction -- if you provided MonoVan run or relative units if monoVan is not present
                out=iliad("wb_wksp","w1",energy,energybin,mapping,MonoVanRun,**argi)

                ws_name = 'MER{0}_ei{1:3.1f}mev'.format(run,energy);
                #RenameWorkspace(InputWorkspace=out,OutputWorkspace=ws_name);
                fileName=ws_name+'.nxspe';
                print 'Saving results into file {0}'.format(fileName)
                SaveNXSPE(InputWorkspace=out,Filename=fileName)
##############################
elapsed = (time.clock() - start)
print "Time to run the whole reduction",(elapsed)
