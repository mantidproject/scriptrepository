import  MERReductionSrRuO3 as mpr
from mantid.simpleapi import *
from mantid import config
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload

reload(mpr)
rd=mpr.MERLINReduction()

map_mask_dir = '/home/zmp58988/RB1510482'
config.appendDataSearchDir(map_mask_dir)

# set up advanced and main properties
rd.def_advanced_properties()
rd.def_main_properties()

#Filename?
rd.set_custom_output_filename()

def iliad_merlin_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range):

    rd.reducer.prop_man.map_file='one2one_143.map'
    #rd.reducer.prop_man.hardmaskOnly = "MER23698.msk"
    rd.reducer.prop_man.hardmaskPlus = "MER23698.msk"
    
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
        
     
    #rd.reducer.prop_man.save_file_name='mer'+str(runno)+'_ei'+str(int(round(ei)))
    rd.run_reduction()
    
def iliad_merlin_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range):
    
    rd.reducer.prop_man.map_file='rings_143.map'
    rd.reducer.prop_man.hard_mask_file = "MER23698.msk"
    
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

    #if  len(runno)==1:
    #    rd.reducer.prop_man.save_file_name='map'+str(runno)+'_ei'+str(int(round(ei)))+'_powder'
    #    rd.reducer.prop_man.sum_runs=False
    #else:
    #    firstbit='map'
    #    lastbit='_ei'+str(int(round(ei)))+'_powder'
    #    middlebit=''
    #    for i in range(len(runno)):
    #        middlebit=middlebit+'_'+str(runno[i])
    #    rd.reducer.prop_man.sum_runs=True
    #    rd.reducer.prop_man.save_file_name=firstbit+middlebit+lastbit
    #    
    rd.run_reduction()
    