import  MARI_Sample as mpr
from mantid.simpleapi import *
from mantid import config

reload(mpr)
rd = mpr.MARIReduction()

# set up advanced and main properties
rd.def_advanced_properties()
rd.def_main_properties()

#Filename?
rd.set_custom_output_filename()

def iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm):

    rd.reducer.prop_man.map_file="mari_res2013.map"
    rd.reducer.prop_man.hard_mask_file = "mari_mask2015.msk"
    
    rd.reducer.prop_man.incident_energy=ei
    
    rd.reducer.prop_man.sample_run = runno
    rd.reducer.prop_man.wb_run=wbvan
    rd.reducer.prop_man.energy_bins=[-ei,ei/100.,0.95*ei]
    
    if ( sam_rmm!=0 and sam_mass!=0 ) :
        abs_units=1
        rd.reducer.prop_man.sample_mass=sam_mass
        rd.reducer.prop_man.sample_rmm=sam_rmm
        rd.reducer.prop_man.monovan_run=monovan
    else:
        abs_units=0
        rd.reducer.prop_man.monovan_run=None
        
     
    rd.reducer.prop_man.save_file_name='mar'+str(runno)+'_ei'+str(int(round(ei)))
    rd.run_reduction()
    
