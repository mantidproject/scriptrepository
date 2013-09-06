"""
The sample script for LET absolute units reduction
"""
from qtiGenie import *
iliad_setup('let')


# program to crunch down event mode from LET and produce output SPE files. Program can automatically find all incident energies in rep rate mode and write out spe files in the # form of LET'run no: +ei'.spe


#############################################
# where to save resutls (usually specified in Mantid, data search directories)
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    save_dir = config.getString('defaultsave.directory')
    
print "Data will be saved into: ",save_dir
# map mask and cal file, again the values from Mantid, data search directories can be modified here
config.appendDataSearchDir('/home/let/Desktop/LAT_maps') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/isisdatar55/NDXLET/Instrument/data/cycle_12_3') 


# this is the user input section
wb=10431   # enter whitebeam run number here (cycle 2013/1)
#run_no=[8570,8581] # event mode run numbers here or use next line for a continous sequence of runs i.e range(first run, last run +1)
#run_no=range(11376,11379)
run_no=[11398]#range(11440,11516)
ei = [2.3]#[0.943] 
#ei=[5.8,15]           # incident energies you want analysed, or leave as ei=[]  if you want all incident energies analysed
ebin=[-0.2,0.002,0.8]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
mapping='LET_one2one_123'  # rings mapping file for powders
file = 'hard_2013_2.msk'    # standard hard mask file  for LET
#file = 'magnet_9T_pm45_hard.msk'  #mask for pm 45 9T magnet orientation
#file = '9Tmagnet_0to90_hard.msk'  #mask for pm 45 9T magnet orientation
#file = 'hard_14Tmagnet.msk'  #14T mask
#file = 'magnet7T_hard.msk'    #7T mask
############################################
 
 
##########################
 
LoadRaw(Filename=str(wb),OutputWorkspace="wb_wksp") # load whitebeam
 
######################################################################
 
 
for run in run_no:     #loop around runs
        fname='LET0000'+str(run)+'.nxs'
        print ' processing file ', fname
        LoadEventNexus(Filename=fname,OutputWorkspace='w1',SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='1')
        Rebin(InputWorkspace='w1_monitors',OutputWorkspace='mon',Params=[1000,100,90000])
        ExtractSingleSpectrum(InputWorkspace='mon',OutputWorkspace='mon',WorkspaceIndex='5')  #extract monitor 6
        ConvertToMatrixWorkspace(InputWorkspace='mon',OutputWorkspace='mon')
        ConvertUnits(InputWorkspace='mon',OutputWorkspace='mon',Target='Energy')
        NormaliseByCurrent(InputWorkspace='mon',OutputWorkspace='mon')     #monitor 6 converted to energy and normalised
        ConjoinWorkspaces(InputWorkspace1='w1',InputWorkspace2='w1_monitors')
        ##################################300
        # this section finds all the transmitted incident energies
        if len(ei) == 0:
           for x in range(0,15):
               Max(InputWorkspace='mon',OutputWorkspace='maxval')
               mv=mtd['maxval']
               if mv.dataY(0)[0] >= 250:
                     min=mv.dataX(0)[0] -0.02
                     max=mv.dataX(0)[1] +0.02
                     RemoveBins(InputWorkspace='mon',OutputWorkspace='mon',XMin=min,XMax=max)
                     ei.append(mv.dataX(0)[0])
        ei.sort()     #sorts energies into order
        print ei
        if run == run_no[0]:
                        ei = [ '%.2f' % elem for elem in ei ]
        print 'energies transmitted are:'
        print (ei)
 
        for energy in ei:
            energy=float(energy)
            print (energy)
            emin=0.2*energy   #minimum energy is with 80% energy loss
            lam=(81.81/energy)**0.5
            lam_max=(81.81/emin)**0.5
            tsam=252.82*lam*25   #time at sample
            tmon2=252.82*lam*23.5 #time to monitor 6 on LET
            tmax=tsam+(252.82*lam_max*4.1) #maximum time to measure inelastic signal to
            t_elastic=tsam+(252.82*lam*4.1)   #maximum time of elastic signal
            tbin=[int(tmon2),1.6,int(tmax)]
            Rebin(InputWorkspace='w1',OutputWorkspace='w1reb',Params=tbin,PreserveEvents='1')        
 
                    
            energybin=[ebin[0]*energy,ebin[1]*energy,ebin[2]*energy]
            energybin = [ '%.4f' % elem for elem in energybin ] 
            ebinstring=str(energybin[0])+','+str(energybin[1])+','+str(energybin[2])
            print ebinstring
            argi={};
            argi['save_format']=''
            argi['fixei']=False
            argi['bleed']=False
            argi['norm_method']='current'
            argi['det_cal_file']='det_LET_cycle12-3.dat'
            argi['detector_van_range']=[0.5,200]
            argi['bkgd_range']=[int(t_elastic),int(tmax)]
            argi['hardmaskOnly']=file

            for kk,val in argi.iteritems():
                print "arguments :" ,kk,val
            out=iliad("wb_wksp","w1reb",energy,ebinstring,mapping,**argi)
            SaveNXSPE(out,'LET'+str(run)+'_'+str(energy)+'mask_mev.nxspe') #
