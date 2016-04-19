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
       prop['incident_energy'] = [2]
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       #
       prop['energy_bins'] = [-0.1,0.0025,0.6] #binning of the energy for the spe file. 
       #prop['energy_bins'] = [-2.5,0.01,4.9] #binning of the energy for the spe file. 
       #
       # the range of files to reduce.
       #
       # the range of files to reduce. This range ignored when deployed from autoreduction,
       # unless you going to sum these files. 
       # The range of numbers or run number is used when you run reduction from PC.
       prop['sample_run'] = 26771 # 'LET18547.n001'
       prop['wb_run'] = 25653   # monovan run number 
       #
       prop['sum_runs'] = False # set to true to sum everything provided to sample_run
       #   
       # Absolute units reduction properties. Set prop['monovan_run']=None to do arbitrary units
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
      prop['map_file'] = 'LET_one2one_153.map'
      prop['det_cal_file'] = 'det_corrected_cycle153.dat'
      prop['bleed'] = False
      prop['norm_method']='current'
      prop['detector_van_range']=[4.8,5.2]
      prop['background_range'] = [92000,98000] # TOF range for the calculating flat background
      prop['hardmaskOnly']='hard_2015_4_9T_0to90_test.msk' # Use diag (hardmaskPlus option) to enhance hard masks
      prop['check_background']=False
      prop['save_format'] = 'nxspe'
      # if two input files with the same name and  different extension found, what to prefer. 
      prop['data_file_ext']='.nxs' # for LET it may be choice between event and histo mode if 
      # raw file is written in histo, and nxs -- in event mode
      # Absolute units: map file to calculate monovan integrals
      prop['monovan_mapfile'] = 'LET_rings_153.map'
#
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
      prop['bleed_maxrate']=0.005
      return prop
      #
#------------------------------------------------------------------------------------#
   @iliad
   def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file

          Overload only if custom reduction is needed or 
          special features are requested
      """
      ws = PropertyManager.sample_run.get_workspace()
      wsf = filter_ts1_pulses(ws)
      self.reducer.prop_man.sample_run = wsf
      PropertyManager.sample_run.synchronize_ws()

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
       
def filter_ts1_pulses(ws,peak_period=20000,peak_res = 300):
    """ Function filters (removes) events occuring with periodicity cpecified by  peak_period value
        and within time interval peak_res (from each event period start)
        Time intervals are in units of uSec.
        
        Intended for filtering fast neutrons peaks generated by TS1 on LET
    
    """
    x_s = ws.readX(0)
    xmin = 0 #x_s[0]
    xmax = x_s[-1]
    #ws=Rebin(ws,OutputWorkspace=ws.name(),Params=[xmin,xmax,10],PreserveEvents=1)
    bin_ranges = range(int(xmin),int(xmax),peak_period)
    nb = len(bin_ranges)-1

    ws_list = []
    r_min = xmin
    print "Splitting input workspace {0} into {1} time intervals excluding specified periods".format(ws.name(),nb)
    for ind,x_min in enumerate(bin_ranges):
        r_max = x_min + peak_period-10
        if r_max > xmax:
            r_max = xmax
        if r_max<=r_min:
            continue
        print 'Step #{0}/{1}; r_min={2}, r_max={3}'.format(ind+1,nb,r_min,r_max)
        part  = FilterByXValue(ws,OutputWorkspace ='{0}_{1}'.format(ws.name(),ind+1),XMin=r_min,XMax=r_max)
        ws_list.append(part)
		
        r_min  = x_min + peak_period+peak_res
 
    targ_ws_name = '{0}_fltrd'.format(ws.name())
    MergeRuns(InputWorkspaces=ws_list,OutputWorkspace=targ_ws_name)
	# transfer monitor workspace
    wst = mtd[targ_ws_name]
    try:
       mon_ws = ws.getMonitorWorkspace()
    except RuntimeError: # May be old code or some problem connecting workspace with monitor workspace
      ws_name = ws.name()
      mon_ws_name = ws_name+'_monitors'
      if mon_ws_name in mtd:
         mon_ws = mtd[mon_ws_name]
      else:
         mon_ws = None
    if mon_ws:
       wst.setMonitorWorkspace(mon_ws) # connect workspace and monitors together  	
	   
    for ws in ws_list:
       DeleteWorkspace(ws)
	 
    return wst
       
#
if __name__ == "__main__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    map_mask_dir = '/usr/local/mprogs/InstrumentFiles/merlin'
    config.appendDataSearchDir(map_mask_dir)
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
   # 
    rd.run_reduction()

