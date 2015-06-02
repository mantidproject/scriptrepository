""" Sample MAPS reduction script """ 
# Two rows necessary to run script outside of the mantid. You need also set up 
# appropriate python path-es
import os
#os.environ["PATH"] = r"c:/Mantid/Code/builds/br_master/bin/Release;"+\
#                     os.environ["PATH"]
#
from mantid import *
from Direct.ReductionWrapper import *

class MAPSReduction(ReductionWrapper):
#------------------------------------------------------------------------------------#
   @MainProperties
   def def_main_properties(self):
       """ Define main properties used in reduction. These are the property 
           a user usually wants to change
       """ 
       prop = {}
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       #
       prop['incident_energy'] = 450
       prop['energy_bins'] = [-50,2.5,425]

       # the range of files to reduce. This range ignored when deployed from autoreduction,
       # unless you going to sum these files. 
       # The range of numbers or run number is used when you run reduction from PC.
       prop['sample_run'] = [21384,21385]
       prop['wb_run'] = 21376
       #
       prop['sum_runs'] = False # set to true to sum everything provided to sample_run
       #                        # list
       # Absolute units reduction properties. Set prop['monovan_run']=None to do relative units
       prop['monovan_run'] = 21803  #  vanadium run in the same configuration as your sample 
       #prop['sample_mass'] = 41.104
       #prop['sample_rmm'] = 398.9439
       return prop
#------------------------------------------------------------------------------------#
   @AdvancedProperties
   def def_advanced_properties(self):
      """  Set up advanced properties, describing reduction.
           These are the properties, usually provided by an instrument 
           scientist
            
           separation between simple and advanced properties depends
           on scientist, experiment and user.   All are necessary for reduction 
           to work properly
      """
      prop = {}
      prop['map_file'] = "4to1.map"
      prop['monovan_mapfile'] = "4to1_mid_lowang.map"
      #prop['hardmaskOnly']=maskfile # disable diag, use only hard mask
      prop['hard_mask_file'] = "4to1_142.msk"
      prop['bkgd_range'] = [13000,19000]

      prop['monovan_lo_frac'] = -0.5 # default is -0.6
      #prop['monovan_hi_frac'] = 0.7 # default is 0.7, no need to change
      #prop['abs_units_van_range']=[-40,40] # specify energy range directly, to
                                     #override relative default energy range
      prop['diag_remove_zero'] = False
      prop['wb_integr_range'] = [20,100] 
      
      #prop['det_cal_file'] = "11060" what about calibration?
      prop['save_format'] = 'nxspe' # nxs or spe
      #prop['data_file_ext']='.nxs' # if two input files with the same name and
                                    #different extension found, what to prefer.
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
      #SaveNexus(ws,Filename = 'MARNewReduction.nxs')
      return results
   #
   #
   def set_custom_output_filename(self):
      """define custom name of output files if standard one is not satisfactory
        
          In addition to that, example of accessing complex reduction properties
          Simple reduction properties can be accessed as e.g.: value= prop_man.sum_runs
      """
      def custom_name(prop_man):
            """Sample function which builds filename from
              incident energy and run number and adds some auxiliary information
              to it.
            """
            # Note -- properties have the same names as the list of advanced and
            # main properties
            ei = PropertyManager.incident_energy.get_current()
            # sample run is more then just list of runs, so we use
            # the formalization below to access its methods
            run_num = PropertyManager.sample_run.run_number()
            name = "MAR{0}_Ei{1:<3.2f}meV_One2One".format(run_num ,ei)
            return name
       
      # Uncomment this to use custom filename function
      # Note: the properties are stored in prop_man class accessed as
        # below.
      return lambda : custom_name(self.reducer.prop_man)
      # use this method to use standard file name generating function
      #return None
   #
   #
   def validation_file_place(self):
      """Redefine this to the place, where validation file, used in conjunction with
         'validate_run' property, located. Here it defines the place to this script folder.
          but if this function is disabled, by default it looks for/places it 
          in a default save directory"""
      return os.path.split(os.path.realpath(__file__))[0]
   
   def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'MAP',web_var)

if __name__ == "__main__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    map_mask_dir = '/usr/local/mprogs/Libisis/InstrumentFiles/maps'
    # folder where input data can be found
    data_dir = '/home/maps/maps_data'
    # auxiliary folder with results
    #ref_data_dir = '/isisdatar55/ndxmaps/Instrument/data/cycle_09_05' 
    # Set input search path to values, specified above
    config.setDataSearchDirs('{0};{1}'.format(data_dir,map_mask_dir))
    # use appendDataSearch directory to add more locations to existing Mantid 
    # data search path
    #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
    # folder to save resulting spe/nxspe files.
    config['defaultsave.directory'] = '/home/maps/maps_users/Ewings/HoraceWorkshop/Fe_data/' #data_dir 

###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = MAPSReduction()
    # set up advanced and main properties
    rd.def_advanced_properties()
    rd.def_main_properties()

#### uncomment rows below to generate web variables and save then to transfer to   ###
    ## web services.
    run_dir = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(run_dir,'reduce_vars.py')
    rd.save_web_variables(file)

#### Set up time interval (sec) for reducer to check for input data file.         ####
    #  If this file is not present and this value is 0,reduction fails
    #  if this value >0 the reduction waits until file appears on the data
    #  search path checking after time specified below.
    rd.wait_for_file = 0  # waiting time interval in seconds

### Define a run number to validate reduction against future changes    #############
    # After reduction works well and all settings are done and verified, 
    # take a run number with good reduced results and build validation
    # for this result. 
    # Then place the validation run together with this reduction script.
    # Next time, the script will run reduction and compare the reduction results against
    # the results obtained earlier.
    #rd.validate_run_number = 21968  # Enabling this property disables normal reduction
    # and forces reduction to reduce run specified here and compares results against
    # validation file, processed earlier or calculate this file if run for the first time.
    #This would ensure that reduction script have not changed,
    #allow to identify the reason for changes if it was changed 
    # and would allow to recover the script,used to produce initial reduction
    #if changes are unacceptable.

####get reduction parameters from properties above, override what you want locally ###
   # and run reduction. Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g. 
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or 
   ## rd.reducer.prop_man.sum_runs = False
   # 
    
    rd.run_reduction()

