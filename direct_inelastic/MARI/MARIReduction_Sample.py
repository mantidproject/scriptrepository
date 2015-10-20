#pylint: disable=invalid-name
""" Sample MARI reduction script """
import os,sys
#os.environ["PATH"] =\
#r"c:/Mantid/Code/builds/br_master/bin/Release;"+os.environ["PATH"]

""" Sample MARI reduction scrip used in testing ReductionWrapper """
from numpy import *
from mantid import *
from Direct.ReductionWrapper import *

class MARIReduction(ReductionWrapper):
    @MainProperties
    def def_main_properties(self):
        """Define main properties used in reduction. These are the property
           a user usually wants to change
		   
		   MARI Instrument scientist beware!!!!
		   -- the properties set up here may be owerridden in iliad_mari if you use it or 
		   in section name=='__main__' below if you do not use iliad_mari
        """
        prop = {}
        # if energy is specified as a list (even with single value e.g. ei=[81])
        # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
        # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
        #
        #prop['incident_energy'] = 50
        #prop['energy_bins'] = [-20,0.1,49]
        prop['incident_energy'] = 6.0
        prop['energy_bins'] = [-10,0.03,5.7]
        #prop['incident_energy'] = 10
        #prop['energy_bins'] = [-10,0.05,9]
        #
        # the range of files to reduce. This range ignored when deployed from autoreduction,
        # unless you going to sum these files. 
        # The range of numbers or run number is used when you run reduction from PC.

        prop['sample_run'] =[20181, 20182, 20183] #20176 #[20171, 20174]
        prop['wb_run'] = 19717

        #
        prop['sum_runs'] = False # set to true to sum everything provided to sample_run
        #                        # list
        # Absolute units reduction properties. Set prop['monovan_run']=None to do relative units
        prop['monovan_run'] = None
        #prop['sample_mass'] = 10
        #prop['sample_rmm'] = 10
        return prop

    @AdvancedProperties
    def def_advanced_properties(self):
        """Set up advanced properties, describing reduction.
           These are the properties, usually provided by an instrument
           scientist

           separation between simple and advanced properties depends
           on scientist, experiment and user.   All are necessary for reduction 
           to work properly
		   
		   MARI Instrument scientist beware!!!!
		   -- the properties set up here may be owerridden in iliad_mari if you use it or 
		   in section name=='__main__' below if you do not use iliad_mari		   
        """
        prop = {}
#       prop['sum_runs']=False #True
        prop['map_file'] = "mari_res2013.map"
        prop['monovan_mapfile'] = "mari_res2013.map"
        #prop['hardmaskOnly']=maskfile # disable diag, use only hard mask
        prop['hard_mask_file'] = "mari_mask2015.msk"
        prop['det_cal_file'] = 19717
        prop['save_format'] = 'nxspe'
        #
        #prop['wb_integr_range'] = [2,10]         
        #prop['data_file_ext']='.nxs' # if two input files with the same name and
                                    #different extension found, what to prefer.
        prop['load_monitors_with_workspace'] = True
        # change this to correct value and verify that motor_log_names refers correct and existing 
        # log name for crystal rotation to write correct psi value into nxspe files
        prop['motor_offset']=None
        return prop
      # 
    @iliad
    def reduce(self,input_file=None,output_directory=None):
        """Method executes reduction over single file

          Overload only if custom reduction is needed or 
          special features are requested
        """
        output = ReductionWrapper.reduce(self,input_file,output_directory)
        run_num = PropertyManager.sample_run.run_number()		
        NewName = 'MAR{0}Reduced'.format(run_num)
        RenameWorkspace(output,OutputWorkspace=NewName )		
        ei = PropertyManager.incident_energy.get_current()
        q_min = 0.04*sqrt(ei)		
        q_max = 1.3*sqrt(ei)		
        q_bins = str(q_min)+','+str(q_max/285.)+','+str(q_max)		
        SofQW3(InputWorkspace=NewName ,OutputWorkspace='MAR{0}Reduced'.format(run_num)+'_SQW',QAxisBinning=q_bins,Emode='Direct')
        Transpose(InputWorkspace='MAR{0}Reduced'.format(run_num)+'_SQW',OutputWorkspace='MAR{0}Reduced'.format(run_num)+'_SQW')
		#SaveNexus(outWS,Filename = 'MARNewReduction.nxs')
        return mtd[NewName]

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
            name = "MAR{0}_Ei{1:<3.2f}meV".format(run_num ,ei)
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
        ReductionWrapper.__init__(self,'MAR',web_var)
#------------------------------------------------------------------------------		
rd = MARIReduction()
# set up advanced and main properties
rd.def_advanced_properties()
rd.def_main_properties()

#Filename -- this generates dynamic name -- the method should be modified in MARIReduction_Sample
#rd.set_custom_output_filename()

def iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs=False):
    """Helper function, which allow to run MARIReduction in simple way
	inputs: 
	    runno -- one or list of run numbers to process
	    ei    -- incident energy
	    wbvan --  white beam vanadium run number or file name of the vanadium
        monovan-- monochromatic vanadium run number or file name
        sam_mass - mass of the sample under investigation
        sam_rmm  --rmm of sample under investigation
        sum_runs -- if true, all runs provided in runno list should be added together
	"""

    rd.reducer.prop_man.map_file="mari_res2013.map"
    rd.reducer.prop_man.hard_mask_file = "mari_mask2015.msk"
    
    rd.reducer.prop_man.incident_energy=ei
    rd.reducer.prop_man.sum_runs= sum_runs
    rd.reducer.prop_man.sample_run = runno
    rd.reducer.prop_man.wb_run=wbvan
    rd.reducer.prop_man.energy_bins=[-0.5*ei,ei/200.,0.97*ei]
    
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

if __name__ == "__main__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    #map_mask_dir = '/usr/local/mprogs/InstrumentFiles/maps'
    # folder where input data can be found
    #data_dir = r'\\isis\inst$\NDXMARI\Instrument\data\cycle_14_2'
    #config.appendDataSearchDir(map_mask_dir)
    #config.appendDataSearchDir(data_dir)

    root=os.path.dirname(os.path.realpath(__file__))
    #data_dir = os.path.join(root,r'data')

    #config.appendDataSearchDir(root)
    #config.appendDataSearchDir(data_dir)
    #config['defaultsave.directory']=root
###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = MARIReduction()
    rd.def_advanced_properties()
    rd.def_main_properties()

#### uncomment rows below to generate web variables and save then to transfer to ###
    ## web services.
    run_dir = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(run_dir,'reduce_vars.py')
    rd.save_web_variables(file)

#### Set up time interval (sec) for reducer to check for input data file.  ####
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
   # and run reduction.  Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g. 
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or 
   ## rd.reducer.prop_man.sum_runs = False
   # 
###### Run reduction over all run numbers or files assigned to ######
     # sample_run variable

    # return output workspace only if you are going to do
    # something with it here.  Running range of runs will return the array
    # of workspace pointers.
    #red_ws = rd.run_reduction()
    # usual way to go is to reduce workspace and save it internally
   
    rd.run_reduction()