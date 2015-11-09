import  MAPSReductionSample2015 as mpr
from mantid.simpleapi import *
from mantid import config

reload(mpr)
#rd = mpr.ReduceMAPS()
rd=mpr.MAPSReduction()


# set up advanced and main properties
rd.def_advanced_properties()
rd.def_main_properties()

#Ensure that default is not to sum runs
rd.reducer.prop_man.sum_runs=False

#Filename?
rd.set_custom_output_filename()

def iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range):

    rd.reducer.prop_man.map_file='4to1.map'
    rd.reducer.prop_man.hard_mask_file = "4to1_153.msk"
    #rd.reducer.prop_man.hard_mask_file = "4to1_143.msk"
    rd.reducer.prop_man.van_mass=30.1
    
    rd.reducer.prop_man.incident_energy=ei
    
    rd.reducer.prop_man.sample_run = runno
    rd.reducer.prop_man.wb_run=wbvan
    rd.reducer.prop_man.energy_bins=rebin_pars
    
    rd.reducer.prop_man.bkgd_range=bg_range
    
    if ( sam_rmm!=0 and sam_mass!=0 ) :
        abs_units=1
        rd.reducer.prop_man.sample_mass=sam_mass
        rd.reducer.prop_man.sample_rmm=sam_rmm
        rd.reducer.prop_man.monovan_run=monovan
    else:
        abs_units=0
        rd.reducer.prop_man.monovan_run='None'
        
     
    if isinstance(runno,(int,long)):
        runno=[runno]   
    
    if  len(runno)==1:
        runstring=str(runno)
        runstring=runstring.strip('[]')
        rd.reducer.prop_man.save_file_name='map'+runstring+'_ei'+str(int(round(ei)))
        rd.reducer.prop_man.sum_runs=False
    else:
        firstbit='map'
        lastbit='_ei'+str(int(round(ei)))
        middlebit=''
        for i in range(len(runno)):
            middlebit=middlebit+'_'+str(runno[i]).strip('[]')
        rd.reducer.prop_man.sum_runs=True
        rd.reducer.prop_man.save_file_name=firstbit+middlebit+lastbit
        
     #    
    rd.run_reduction()   
        
    
def iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range):
    
    #rd.reducer.prop_man.map_file='parker_rings.map'
    rd.reducer.prop_man.map_file='MAPS_rings.map'
    rd.reducer.prop_man.hard_mask_file = "4to1_153.msk"
    
    rd.reducer.prop_man.incident_energy=ei
    
    rd.reducer.prop_man.sample_run = runno
    rd.reducer.prop_man.wb_run=wbvan
    rd.reducer.prop_man.energy_bins=rebin_pars
    
    rd.reducer.prop_man.bkgd_range=bg_range
    rd.reducer.prop_man.save_format = 'nxspe' 
    
    if ( sam_rmm!=0 and sam_mass!=0 ) :
        abs_units=1
        rd.reducer.prop_man.sample_mass=sam_mass
        rd.reducer.prop_man.sample_rmm=sam_rmm
        rd.reducer.prop_man.monovan_run=monovan
    else:
        abs_units=0
        rd.reducer.prop_man.monovan_run='None'
    
    if isinstance(runno,(int,long)):
        runno=[runno]   
    
    if  len(runno)==1:
        runstring=str(runno)
        runstring=runstring.strip('[]')
        rd.reducer.prop_man.save_file_name='map'+runstring+'_ei'+str(int(round(ei)))+'_powder'
        rd.reducer.prop_man.sum_runs=False
    else:
        firstbit='map'
        lastbit='_ei'+str(int(round(ei)))+'_powder'
        middlebit=''
        for i in range(len(runno)):
            middlebit=middlebit+'_'+str(runno[i]).strip('[]')
        rd.reducer.prop_man.sum_runs=True
        rd.reducer.prop_man.save_file_name=firstbit+middlebit+lastbit
    #    
    rd.run_reduction()
    
