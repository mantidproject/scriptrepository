"""
Script to perform absolute units data reduction for MARI
"""
from qtiGenie import *
from mantid.simpleapi import *
from mantid import config

import time

#instrument name:
inst=''
iliad_setup(inst)
ext='.raw'

# where to save resutls (usually specified in Mantid, data search directories)
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    save_dir = config.getString('defaultsave.directory')
    
print "Data will be saved into: ",save_dir
# map mask and cal file, again the values from Mantid, data search directories can be modified here
config.appendDataSearchDir('/usr/local/mprogs/InstrumentFiles/mari') 
# data (raw or nxs) run files -- values from data search directories can be modified here
config.appendDataSearchDir('/isisdatar55/NDXMARI/Instrument/data/cycle_05_1') 
config.appendDataSearchDir(r'd:/Data/MantidSystemTests/Data') 


maskfile='mar11015.msk' #'testMask2.msk'#hard mask out the edges of detectors, which tend to be noisy

#map file
mapfile='mari_res.map' # mapping file
#mapfile='/opt/Mantid/instrument/mapfiles/maps/parker_rings' #powder mapping file
mv_mapfile='4to1_mid_lowang'

# latest white beam vanadium file for bad detector diagnosis
wbvan=11060

#Run numbers can be specified as a list:
#runno=[17422,17423, etc]
runno=[11001] #[19399] #,[00004] 19402]

#Incident energy list e.g. ei=[20,30,40]
ei=12

#Sample info
sam_rmm = 435.96
#sam_mass= 2.465
sam_mass = 10
    
rebin_pars=[-11,0.05,11]
monovan=11015
argi = {};
argi['hardmaskPlus']=maskfile 
argi['sample_mass'] = sam_mass   
argi['sample_rmm']   =sam_rmm 



for i in range(len(runno)):
        #w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],rebin_pars,mapfile,mv_mapfile,**argi)
        # this does absolute units normalization as far as monovan is not None. Uses default map file, provided in MAPS_Parameters.xml file. Any changes from defaults should be provided here or above as parameters
        w1=dgreduce.arb_units(wbvan,runno[i],ei,rebin_pars,None,monovan,**argi)
        #w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,**argi)
        #w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,bkgd_range=[13000,19000],\
        #                     hardmaskPlus=maskfile,diag_remove_zero=False,save_format='none')

        
    #Alternative (abs units):
    #w1=iliad_abs(wbvan,runno[i],monovan[i],wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars[i]).strip('[]'),mapfile,mapfile,bkgd_range=[14000,19000],hardmaskPlus=maskfile,diag_remove_zero=False)
    save_file=inst+str(runno[i])+'_ei'+str(ei[i])  
    SaveNXSPE(w1,save_file+'Abs_DgrdOld.nxspe')
    #SaveNexus(w1,save_file+'newDgrd_NewQTG_NewDirectConv.nxs')	
    DeleteWorkspace(w1)




