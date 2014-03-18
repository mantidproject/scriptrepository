import os
import sys
import re
sys.path.append("/opt/Mantid/bin")
from mantid.simpleapi import *
import mantid

# Crab information from autoreduction call
#############################################
# first argument is full path of the data file 
dataFile = sys.argv[1]
# second argument is the output directory
outputDir = sys.argv[2]

# extract information from dataFile path and filename
dataFileName = os.path.split(dataFile)[-1]
dataFilePath = dataFile.replace(dataFileName, '')
dataFileNameMinuxExt = dataFileName.split('.')[0]
runNumber = int(re.findall('\d+', dataFileNameMinuxExt)[0])
instrument = dataFilePath.lower().split('/ndx')[1].split('/')[0]
#############################################

print 'data file path: ' + dataFilePath
sys.path.append(dataFilePath)

from qtiGenie import *
iliad_setup(instrument)

# program to crunch down event mode from LET and produce output SPE files. Program can automatically find all incident energies in rep rate mode and write out spe files in the # form of LET'run no: +ei'.spe

#############################################
# this is the user input section
wb=14392   # enter whitebeam run number here (cycle 2013/3)
#ei = [10.2,4.0,2.13,1.31]  #ei=[5.8,15]           # incident energies you want analysed, or leave as ei=[]  if you want all incident energies analysed
ei = [4.0]  #ei=[5.8,15]           # incident energies you want analysed, or leave as ei=[]  if you want all incident energies analysed
#ei = [30.6,10.2,5.01,2.97]  #ei=[5.8,15]           # incident energies you want analysed, or leave as ei=[]  if you want all incident energies analysed
ebin=[-0.8,0.002,0.9]    #binning of the energy for the spe file. The numbers are as a fraction of ei [from ,step, to ]
mapping='LET_rings_133'  # rings mapping file for powders
file = 'hard_133.msk'    # standard hard mask file  for LET
#file = '/home/let/Desktop/LET_maps/hard_133.msk'    # standard hard mask file  for LET
#file = '/home/let/Desktop/LET_maps/magnet_9T_pm45_hard.msk'  #mask for pm 45 9T magnet orientation
#file = '/home/let/Desktop/LET_maps/9Tmagnet_0to90_hard.msk'  #mask for pm 45 9T magnet orientation
#file = '/home/let/Desktop/LET_maps/hard_14Tmagnet.msk'  #14T mask
#file = '/home/let/Desktop/LET_maps/magnet7T_hard.msk'    #7T mask
############################################
# currently done here on-
remove_background = False  #if true then will subtract a flat background in time from the time range given below otherwise put False
bg_range=[92000,98000] # range of times to take background in

# Absolute units reduction MonoVanRun=None disables it. Sample and vanadium mass parameters have to be right otherwise
# Vanadium labelled Dec 2011 - flat plate of dimensions: 40.5x41x2.0# volume = 3404.025 mm**3 mass= 20.79
sampleMass=20.79 # 17.25  # mass of your sample (PrAl3)
sampleRMM= 50.9415 # 221.854  # molecular weight of your sample
MonoVanRun=None # vanadium run in the same configuration as your sample
monovan_mapfile='rings_103.map'  #

##########################

wbFile= dataFilePath + 'LET000'+str(wb)+'.raw'
# White beam
if 'wb_wksp' in mtd and wb.getRunNumber()==wb:
    wb_wksp=mtd['wb_wksp']
else:  #only load whitebeam if not already there
    wb_wksp=LoadRaw(Filename=wbFile,OutputWorkspace="wb_wksp",LoadMonitors='Exclude')

######################################################################

print ' processing file ', dataFile
#w1 = dgreduce.getReducer().load_data(run,'w1')
Load(Filename=dataFile,OutputWorkspace='w1',LoadMonitors='1');

    
if remove_background:
    bg_ws_name=find_background('w1',bg_range,1.6);

	
##################################300
# this section finds all the transmitted incident energies
if len(ei) == 0:
   ei = find_chopper_peaks('w1_monitors');
   print ei
    
print 'energies to process are:'
print (ei)

RenameWorkspace(InputWorkspace = 'w1',OutputWorkspace='w1_storage');
RenameWorkspace(InputWorkspace = 'w1_monitors',OutputWorkspace='w1_mon_storage');

# generic iliad parameters:
argi={};
argi['bleed']=False
argi['norm_method']='current'
argi['det_cal_file']='det_LET_cycle133.dat'
argi['detector_van_range']=[0.5,200]
argi['hardmaskOnly']=file
# abs units
argi['sample_mass']=sampleMass;
argi['sample_rmm'] =sampleRMM;
argi['monovan_mapfile']=monovan_mapfile;


for ind,energy in enumerate(ei):
    print "Reducing around energy: {0}".format(float(energy))
    (energybin,tbin,t_elastic) = find_binning_range(energy,ebin,25.,23.5,4.1,1.6);
                       
    energybin=[ebin[0]*energy,ebin[1]*energy,ebin[2]*energy]
    energybin = [ '%.4f' % elem for elem in energybin ] 
    ebinstring=str(energybin[0])+','+str(energybin[1])+','+str(energybin[2])
    print ebinstring

    # background removal range  -- used for standard background removal
    argi['bkgd_range']=[int(t_elastic),int(tmax)]
    # abs units
    argi['monovan_integr_range']=[round(ebin[0]*energy,4),round(ebin[2]*energy,4)]; # integration range of the vanadium 

    # if we calculate more then one energy, initial workspace will be used more then once. Is this enough for that?    
    #Rebin(InputWorkspace='w1',OutputWorkspace='w1reb',Params=tbin,PreserveEvents='1')                
    # safe option
    if ind <len(ei)-1:
        CloneWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1')
        CloneWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors')
    else:
        RenameWorkspace(InputWorkspace = 'w1_storage',OutputWorkspace='w1');
        RenameWorkspace(InputWorkspace = 'w1_mon_storage',OutputWorkspace='w1_monitors');
    
    if remove_background:
        Minus(LHSWorkspace='w1',RHSWorkspace=bg_ws_name,OutputWorkspace='w1')
        
		
    # absolute unit reduction -- if you provided MonoVan run or relative units if monoVan is not present
    out=iliad("wb_wksp","w1reb",energy,energybin,mapping,MonoVanRun,**argi);
    file_name = 'LET{0}_reducedEi{1:4.2f}_mev.nxspe'.format(runNumber,energy);
    SaveNXSPE(InputWorkspace=out,Filename=file_name);

