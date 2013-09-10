"""
Script to perform absolute units data reduction for MAPS
"""
from qtiGenie import *
from mantid.simpleapi import *
from mantid import config

import time

#instrument name:
inst='map'
iliad_setup(inst)
ext='.raw'

# where to save resutls (usually specified in Mantid, data search directories)
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    save_dir = config.getString('defaultsave.directory')
    
print "Data will be saved into: ",save_dir
# map mask and cal file, again the values from Mantid, data search directories can be modified here
config.appendDataSearchDir('/home/maps/mprogs/InstrumentFiles/maps') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/isisdatar55/NDXMAPS/Instrument/data/cycle_12_3') 

maskfile='4to1_022.msk' #'testMask2.msk'#hard mask out the edges of detectors, which tend to be noisy

#map file
mapfile='4to1' #single crystal mapping file
#mapfile='/opt/Mantid/instrument/mapfiles/maps/parker_rings' #powder mapping file
mv_mapfile='4to1_mid_lowang'

#If run number is 00000 (from updatestore) delete existing workspace so that new raw data file is loaded
try:
    map00000=CloneWorkspace('MAP00000')
    DeleteWorkspace('MAP00000')
    DeleteWorkspace('map00000')
except:
    print('Workspace zero did not exist anyway')


try:
    map00000=CloneWorkspace('MAP00000.raw')
    DeleteWorkspace('MAP00000.raw')
    DeleteWorkspace('map00000')
except:
    print('Workspace zero did not exist anyway')
    
argi = {};
argi['bkgd_range'] = [13000,19000]
argi['hardmaskPlus']=maskfile 
#argi['hardmaskOnly']=maskfile 
argi['diag_remove_zero']=False
argi['abs_units_van_range']=[-40,40]   
argi['wb_integr_range'] = [20,100] 
argi['save_format']   = 'none'
## dgREDUCE old
#argi['diag_van_median_rate_limit_hi'] = 100
#argi['diag_van_median_rate_limit_lo'] = 0.01 
#argi['diag_van_median_sigma_lo']=0.1
#argi['diag_van_median_sigma_hi']=1.5
#argi['diag_samp_median_sigma_lo']=0.0
#argi['diag_samp_median_sigma_hi']=2.0
#argi['diag_samp_median_sigma']=3.0
#argi['diag_variation']=1.1

#argi['save_and_reuse_masks']=False


####################################
### DO NOT EDIT ABOVE THIS LINE #############
####################################

#########################################
### BELOW THIS LINE ARE EDITABLE PARAMETERS ##########
#########################################


# latest white beam vanadium file for bad detector diagnosis
wbvan=19327

#Run numbers can be specified as a list:
#runno=[17422,17423, etc]
runno=[19403] #[19399] #,[00004] 19402]

#Incident energy list e.g. ei=[20,30,40]
ei=[60,60]

#Sample info
#sam_rmm=350.653
sam_rmm = 50.9415
#sam_mass= 2.465
sam_mass = 30.1

if ( sam_rmm==0 or sam_mass==0 ) :
	abs_units=1
	argi['sample_mass'] = sam_mass   
	argi['sample_rmm']   =sam_rmm
else:
	abs_units=0
	



for i in range(len(runno)):
    if ei[i]==60:
        if abs_units==1:
		rebin_pars=[-6,0.6,54]
		monovan=19403
		w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],rebin_pars,mapfile,mv_mapfile,**argi)
	else:
		rebin_pars=[-6,0.6,54]
		w1=iliad_abs(wbvan,runno[i],ei[i],rebin_pars,mapfile,**argi)
		
	#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],rebin_pars,mapfile,mv_mapfile,**argi)
        # this does absolute units normalization as far as monovan is not None. Uses default map file, provided in MAPS_Parameters.xml file. Any changes from defaults should be provided here or above as parameters
        
	#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,**argi)
        #w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,bkgd_range=[13000,19000],\
        #                     hardmaskPlus=maskfile,diag_remove_zero=False,save_format='none')

        
    #Alternative (abs units):
    #w1=iliad_abs(wbvan,runno[i],monovan[i],wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars[i]).strip('[]'),mapfile,mapfile,bkgd_range=[14000,19000],hardmaskPlus=maskfile,diag_remove_zero=False)
    save_file=inst+str(runno[i])+'_ei'+str(ei[i])  
    SaveNXSPE(w1,save_file+'Abs_DgrdOld.nxspe')
    #SaveNexus(w1,save_file+'newDgrd_NewQTG_NewDirectConv.nxs')	
    DeleteWorkspace(w1)
    if runno[i]==0:
        DeleteWorkspace('MAP00000.raw')
    else:
        DeleteWorkspace('MAP'+str(runno[i])+'.raw')	




