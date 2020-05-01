""" Sample LET reduction script """ 
import os,sys
#os.environ["PATH"] = r"c:/Mantid/Code/builds/br_master/bin/Release;"+os.environ["PATH"]
from mantid import *
from Direct.ReductionWrapper import *
# this exact name is necessary for web services to work until we implement factory
class LETReduction(ReductionWrapper):
#------------------------------------------------------------------------------------#++++++
   @MainProperties
   def def_main_properties(self):
       """ Define main properties used in reduction. These are the property 
           a user usually wants to change
       """ 
       prop = {}
#---------------------------------------------------------------------------#
#------------- To change for users ---------------------------------------#
#---------------------------------------------------------------------------#
      # multiple energies provided in the data file
       prop['incident_energy'] = [6.03,3.15,1.94] 
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       #
       prop['energy_bins'] = [-0.1,0.0025,0.8] #binning of the energy for the spe file. 
       #prop['energy_bins'] = [-2.5,0.01,4.9] #binning of the energy for the spe file. 
       #
       # the range of files to reduce.
       #
       # the range of files to reduce. This range ignored when deployed from autoreduction,
       # unless you going to sum these files. 
       # The range of numbers or run number is used when you run reduction from PC.
      # prop['sample_run'] = range(50359,52000) # 'LET18547.n001'
       prop['sample_run']=list(range(50776,50793)) #SET DATA RUN NUMBERS GIRI
       prop['wb_run'] = 45879   # monovan run number 
       #
       prop['sum_runs'] = True # set to true to sum everything provided to sample_run
       #   
       # Absolute units reduction properties. Set prop['monovan_run']=None to do arbitrary units
       prop['monovan_run'] = None # vanadium run in the same configuration as your sample
       #prop['sample_mass'] = 100 #  # mass of your sample
       #prop['sample_rmm'] = 200 #  molecular weight of scatterers in your sample
       prop['fixei']=True
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
      prop['map_file'] = 'LET_one2one_153.map'
      prop['det_cal_file'] = 'det_corrected_cycle153.dat'
      prop['norm_method']='current'
      prop['detector_van_range']=[4,20]
      prop['background_range'] = [92000,98000] # TOF range for the calculating flat background
      prop['hardmaskOnly']='hard_164.msk' # Use diag (hardmaskPlus option) to enhance hard masks
	  #prop['hardmaskOnly']='hard_164_Magnet_0-90.msk' # this is for magnet
      prop['check_background']=False
      prop['save_format'] = 'nxspe'
      # if two input files with the same name and  different extension found, what to prefer. 
      prop['data_file_ext']='.nxs' # for LET it may be choice between event and histo mode if 
      # raw file is written in histo, and nxs -- in event mode
      # Absolute units: map file to calculate monovan integrals
      prop['monovan_mapfile'] = 'LET_rings_153.map'
      prop['load_monitors_with_workspace']=False
      # change this to correct value and verify that motor_log_names refers correct and existing 
      # log name for crystal rotation to write correct psi value into nxspe files
      prop['motor_offset']=None
      # vanadium  mass valid from cycle 2013/5
      prop['vanadium-mass']=8.47
      #
      prop['monovan_lo_frac']=-0.4
      prop['monovan_hi_frac']= 0.4
      #RAE added
      
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
      #
      return results
   #
   def set_custom_output_filename(self):
      """ define custom name of output files if standard one is not satisfactory 
          In addition to that, example of accessing reduction properties 
          Changing them if necessary
      """ 
      def custom_name(prop_man):
          """ sample function which builds filename from 
              incident energy and run number and adds some auxiliary information 
              to it.
          """ 
          # Note -- properties have the same names  as the list of advanced and 
          # main properties
          ei = PropertyManager.incident_energy.get_current()
          # sample run is more then just list of runs, so we use 
          # the formalization below to access its methods
          run_num = PropertyManager.sample_run.run_number()
          name = "LET{0}_{1:<3.2f}meV_OneToOne".format(run_num ,ei)
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
       ReductionWrapper.__init__(self,'LET',web_var)
#
if __name__ == "__main__" or __name__ == "__builtin__" or __name__ == "mantidqt.widgets.codeeditor.execution":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
   
    map_mask_dir = '/usr/local/mprogs/InstrumentFiles/let'
    config.appendDataSearchDir(map_mask_dir)
    #LET_dat_dir='home/fep51756/Desktop/RB1810639/NXSPEfile'
    #config.appendDataSearchDir(LET_dat_dir)
    # folder where input data can be found
    #data_dir = 'd:/Data/Mantid_Testing/15_01_27/LET/data'
    # auxiliary folder with results
    #rez_dir = 'd:/Data/Mantid_Testing/15_01_27/LET'
    # Set input search path to values, specified above
    #config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,map_mask_dir,rez_dir))
    # use appendDataSearch directory to add more locations to existing Mantid 
    # data search path
    #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
    # folder to save resulting spe/nxspe files.
    #config['defaultsave.directory'] = rez_dir

###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = LETReduction()
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
    rd.wait_for_file = 120  # waiting time interval in seconds

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
   # 
    rd.run_reduction()

