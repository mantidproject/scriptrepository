""" Sample MAPS reduction script """ 
from __future__ import (absolute_import, division, print_function, unicode_literals)
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import os
from mantid import *
from Direct.ReductionWrapper import *

# path and default directory
sys.path.insert(0,'/instrument/MAPS/RBNumber/USER_RB_FOLDER');
config['defaultsave.directory'] = '/instrument/MAPS/RBNumber/USER_RB_FOLDER'

class MAPSReduction(ReductionWrapper):
#------------------------------------------------------------------------------------#
   @MainProperties
   def def_main_properties(self):
       """ Define main properties used in reduction. These are the property 
           a user usually wants to change
       """ 
       prop = {}

       prop['sample_run']      = [37524,37525]
       prop['sum_runs']        = True        
       prop['wb_run']          = 37520
       
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       prop['incident_energy'] = [300]
       prop['energy_bins']     = [-0.2,0.002,0.9]
          
       # Absolute units reduction properties. Set prop['monovan_run']=None to do relative units
       prop['monovan_run']      = None        #  vanadium run in the same configuration as your sample 
       #prop['sample_mass']     = 41.104
       #prop['sample_rmm']      = 398.9439
       #prop['vanadium-mass']   = 32.94
       
       return prop
#------------------------------------------------------------------------------------#
   @AdvancedProperties
   def def_advanced_properties(self):
      """  Set up advanced properties, describing reduction.
           These are the properties, usually provided by an instrument 
           scientist. Separation between simple and advanced properties depends
           on scientist, experiment and user.   All are necessary for reduction 
           to work properly
      """
      prop = {}
      #prop['map_file']         = "4to1.map"
      prop['map_file']         = "MAPS_rings.map"
      prop['monovan_mapfile']  = "4to1_mid_lowang.map"
      prop['hardmaskOnly']     = "4to1_193_msk.xml" # disable diag, use only hard mask
      #prop['hard_mask_file']   = "4to1_193_msk.xml"
      prop['run_diagnostics']  = True
      prop['bkgd_range']       = [15000,19000]
      prop['normalise_method'] = 'current'
      prop['monovan_lo_frac']  = -0.5 # default is -0.6
      #prop['monovan_hi_frac']  = 0.7 # default is 0.7, no need to change
      #prop['abs_units_van_range']=[-40,40] # specify energy range directly, to
                                     #override relative default energy range
      prop['diag_remove_zero'] = False
      prop['wb_integr_range']  = [20,100] 
      prop['save_format']      = 'nxspe' # nxs or spe
      prop['data_file_ext']    ='.nxs' # if two input files with the same name and
                                       # different extension found, what to prefer.
      return prop
      #
#------------------------------------------------------------------------------------#
   @iliad
   def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file
          Overload only if custom reduction is needed or 
          special features are requested
      """
      results = ReductionWrapper.reduce(self,input_file,output_directory)
      return results
   #
   #
   def set_custom_output_filename(self):
      """ Define custom name of output files if standard one is not satisfactory
          In addition to that, example of accessing complex reduction properties
          Simple reduction properties can be accessed as e.g.: value= prop_man.sum_runs
      """
      def custom_name(prop_man):
            """ Sample function which builds filename from incident energy 
                and run number and adds some auxiliary information to it.
            """
            map_file = prop_man.map_file
            if 'rings' in map_file:
                ftype = '_powder'
            else:
                ftype = ''             

            ei      = PropertyManager.incident_energy.get_current()
            run_num = PropertyManager.sample_run.run_number()
            name    = "map{0}_ei{1:_<3.0f}meV{2}".format(run_num ,ei,ftype)
            return name
       
      # Uncomment this to use custom filename function
      return lambda : custom_name(self.reducer.prop_man)
      # Uncomment this to use standard file name generating function
      #return None
   #
   #
   def validation_file_place(self):
      """ Redefine this to the place, where validation file, used in conjunction with
         'validate_run' property, located. Here it defines the place to this script folder.
          but if this function is disabled, by default it looks for/places it 
          in a default save directory """
      return os.path.split(os.path.realpath(__file__))[0]
   
   def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'MAP',web_var)

#--------------------------------------------------------------------------------------------------------------------------- 
if __name__ == "__main__" or __name__ == "__builtin__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#

###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = MAPSReduction()
    rd.def_advanced_properties()
    rd.def_main_properties()

    run_dir = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(run_dir,'reduce_vars.py')
    rd.save_web_variables(file)
    
    rd.wait_for_file = 0  # waiting time interval in seconds
    
    rd.run_reduction()

