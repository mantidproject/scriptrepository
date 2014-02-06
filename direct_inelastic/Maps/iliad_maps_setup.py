from qtiGenie import *
from mantid.simpleapi import *
from mantid import config

import time

def iliad_maps_setup():
	
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

def iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm):
	#instrument name:
	inst='map'
	iliad_setup(inst)
	ext='.raw'

	argi = {};
	if ( sam_rmm!=0 and sam_mass!=0 ) :
		abs_units=1
		argi['sample_mass'] = sam_mass   
		argi['sample_rmm']   =sam_rmm
	else:
		abs_units=0
		
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
	    
	
	argi['bkgd_range'] = [13000,19000]
	argi['hardmaskPlus']=maskfile 
	#argi['hardmaskOnly']=maskfile 
	argi['diag_remove_zero']=False
	argi['abs_units_van_range']=[-0.5*ei,0.7*ei]   
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

        if abs_units==1:
		w1=iliad_abs(wbvan,runno,monovan,wbvan,sam_rmm,sam_mass,ei,rebin_pars,mapfile,mv_mapfile,**argi)
	else:
		w1=iliad(wbvan,runno,ei,rebin_pars,mapfile,**argi)
		
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],rebin_pars,mapfile,mv_mapfile,**argi)
		# this does absolute units normalization as far as monovan is not None. Uses default map file, provided in MAPS_Parameters.xml file. Any changes from defaults should be provided here or above as parameters
		
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,**argi)
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,bkgd_range=[13000,19000],\
		#                     hardmaskPlus=maskfile,diag_remove_zero=False,save_format='none')
	
	save_file=inst+str(runno)+'_ei'+str(ei)  
	SaveNXSPE(w1,save_file+'.nxspe')
	SaveSPE(w1,save_file+'.spe')
	#SaveNexus(w1,save_file+'newDgrd_NewQTG_NewDirectConv.nxs')	
	DeleteWorkspace(w1)
	if runno==0:
		DeleteWorkspace('MAP00000.raw')
	else:
		DeleteWorkspace('MAP'+str(runno)+'.raw')	

def iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm):
	#instrument name:
	inst='map'
	iliad_setup(inst)
	ext='.raw'

	argi = {};
	if ( sam_rmm!=0 and sam_mass!=0 ) :
		abs_units=1
		argi['sample_mass'] = sam_mass   
		argi['sample_rmm']   =sam_rmm
	else:
		abs_units=0
		
	maskfile='4to1_022.msk' #'testMask2.msk'#hard mask out the edges of detectors, which tend to be noisy

	#map file
	mapfile='parker_rings' #powder mapping file
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

        if abs_units==1:
		w1=iliad_abs(wbvan,runno,monovan,wbvan,sam_rmm,sam_mass,ei,rebin_pars,mapfile,mv_mapfile,**argi)
	else:
		w1=iliad(wbvan,runno,ei,rebin_pars,mapfile,**argi)
		
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],rebin_pars,mapfile,mv_mapfile,**argi)
		# this does absolute units normalization as far as monovan is not None. Uses default map file, provided in MAPS_Parameters.xml file. Any changes from defaults should be provided here or above as parameters
		
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,**argi)
		#w1=iliad_abs(wbvan,runno[i],monovan,wbvan,sam_rmm,sam_mass,ei[i],str(rebin_pars).strip('[]'),mapfile,mv_mapfile,bkgd_range=[13000,19000],\
		#                     hardmaskPlus=maskfile,diag_remove_zero=False,save_format='none')
	
	save_file=inst+str(runno)+'_ei'+str(ei)  
	SaveNXSPE(w1,save_file+'_powder.nxspe')
	SaveSPE(w1,save_file+'_powder.spe')
	#SaveNexus(w1,save_file+'newDgrd_NewQTG_NewDirectConv.nxs')	
	#DeleteWorkspace(w1)
	if runno==0:
		DeleteWorkspace('MAP00000.raw')
	else:
		DeleteWorkspace('MAP'+str(runno)+'.raw')	



