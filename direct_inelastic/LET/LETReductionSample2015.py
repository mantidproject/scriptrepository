""" Sample LET reduction script """ 
# Two rows necessary to run script outside of the mantid. You need also set up 
# appropriate python path-es
import os,sys
os.environ["PATH"] = r"c:/Mantid/Code/builds/br_master/bin/Release;"+os.environ["PATH"]
sys.path.insert(0, "/opt/mantidnightly/bin") 
from mantid import *
from Direct.ReductionWrapper import *
# necessary for web services to work.  Will go with factory implemented
try:
    import reduce_vars as web_var
except:
    web_var = None

class ReduceLET_MultiRep2015(ReductionWrapper):
#------------------------------------------------------------------------------------#
   @MainProperties
   def def_main_properties(self):
       """ Define main properties used in reduction. These are the property 
           a user usually wants to change
       """ 
       prop = {}


       # multiple energies provided in the data file
       prop['incident_energy'] = [2.3] #,5.8]
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       #
       prop['energy_bins'] = [-0.25,0.005,0.9] #binning of the energy for the spe file. 
       #
       # the range of files to reduce. This range ignored when deployed from autoreduction,
       # unless you going to sum these files. 
       # The range of numbers or run number is used when you run reduction from PC.
       prop['sample_run'] = range(18417,18422)
       prop['wb_run'] = 18422
       #
       prop['sum_runs'] = False # set to true to sum everything provided to sample_run
       #                        # list
  
       # Absolute units reduction properties. Set prop['monovan_run']=None to do relative units
       prop['monovan_run'] = None # vanadium run in the same configuration as your sample
       #prop['sample_mass'] = 100 #  # mass of your sample
       #prop['sample_rmm'] = 200 #  molecular weight of scatterers in your sample
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
      prop['map_file'] = 'LET_one2one_cycle142.map'
      prop['det_cal_file'] = 'det_LET_cycle142.dat'
      prop['bleed'] = False
      prop['norm_method']='current'
      prop['detector_van_range']=[4.8,5.2]
      prop['background_range'] = [92000,98000] # TOF range for the calculating flat background
      prop['hardmaskOnly']='9Tmagnet_0to90_hard.msk' # diag does not work well on LET. At present only use a hard mask RIB has created

      prop['check_background']=False

      prop['save_format'] = 'nxspe'
       # if two input files with the same name and  different extension found, what to prefer. 
      prop['data_file_ext']='.nxs' # for LET it may be choice between event and histo mode if 
      # raw file is written in histo, and nxs -- in event mode
      # Absolute units: map file to calculate monovan integrals
      prop['monovan_mapfile'] = 'LET_rings_141.map'
      # change this to correct value and verify that motor_log_names refers correct and existing 
      # log name for crystal rotation to write correct psi value into nxspe files
      prop['motor_offset']=None
      return prop
      #
#------------------------------------------------------------------------------------#
   @iliad
   def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file

          Overload only if custom reduction is needed or 
          special features are requested
      """
      ws = ReductionWrapper.reduce(self,input_file,output_directory)
      #SaveNexus(ws,Filename = 'MARNewReduction.nxs')
      return ws
   #
   def validate_result(self,build_validation=False,Error=1.e-3,ToleranceRelErr=True):
      """ Change this method to verify different results     """
      # here we have: 
      #  18184                    run number with known reduction result
      # LET18184Ei2d30meV_Abs.nxs workspace for run above reduced earlier and we now 
      # validate against
      # build_validation -- if true, build and save new workspace rather
      #then validating the old one
      run_dir = os.path.dirname(os.path.realpath(__file__))
      # this row defines location of the validation file in this script folder
      validation_file = os.path.join(run_dir,"LET18184Ei2d30meV_Abs.nxs")
      rez,message = ReductionWrapper.build_or_validate_result(self,18184,
                                     validation_file,build_validation,
                                     Error,ToleranceRelErr)
      return rez,message
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
          ei = prop_man.incident_energy
          # sample run is more then just list of runs, so we use 
          # the formalization below to access its methods
          run_num = PropertyManager.sample_run.run_number()
          name = "RUN{0}atEi{1:<4.1f}meV_One2One".format(run_num ,ei)
          return name
       
      # Uncomment this to use custom filename function
      # Note: the properties are stored in prop_man class accessed as
      # below. 
      #return custom_name(self.reducer.prop_man)
      # use this method to use standard file name generating function
      return None
    #

   def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'LET',web_var)
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
def main(input_file=None,output_dir=None):
    """ This method is used to run code from web service
        and should not be touched unless you change the name of the
        particular ReductionWrapper class (e.g. ReduceLET_MultiRep2015 here)

        exception to change the output folder to save data to
    """
    # note web variables initialization
    rd = ReduceLET_MultiRep2015(web_var)
    rd.reduce(input_file,output_dir)
    
    # Define folder for web service to copy results to
    output_folder = ''
    return output_folder

if __name__ == "__main__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    map_mask_dir = 'c:/Users/wkc26243/Documents/work/Libisis/InstrumentFiles/let'
    # folder where input data can be found
    data_dir = 'd:/Data/Mantid_Testing/15_01_27/LET/data'
    # auxiliary folder with results
    rez_dir = 'd:/Data/Mantid_Testing/15_01_27/LET'
    # Set input search path to values, specified above
    config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,map_mask_dir,rez_dir))
    # use appendDataSearch directory to add more locations to existing Mantid 
    # data search path
    #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
    # folder to save resulting spe/nxspe files.
    config['defaultsave.directory'] = rez_dir

###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = ReduceLET_MultiRep2015()
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
    #  if this value >0 the reduction wait until file appears on the data 
    #  search path checking after time specified below.
    rd.wait_for_file = 0  # waiting time interval

####get reduction parameters from properties above, override what you want locally ###
   # and run reduction. Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g. 
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or 
   ## rd.reducer.prop_man.sum_runs = False
   # 
    rd.run_reduction()

#### Validate reduction result against known result, obtained earlier  ###
#    rez,mess=rd.validate_result()
#    if not rez:
#      raise RuntimeError("validation failed with error: {0}".format(mess))
#   else:
#     print "ALL Fine" 
