
""" Sample LET reduction script """ 
import os,sys
#os.environ["PATH"] = r"c:/Mantid/Code/builds/br_master/bin/Release;"+os.environ["PATH"]
from mantid import *
from Direct.AbsorptionShapes import *
from Direct.ReductionWrapper import *
try:
    import mantidplot as mpl
except:
    mpl = None





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
        prop['sample_run']=range(50776,50793) #SET DATA RUN NUMBERS GIRI
        prop['wb_run'] = 45879   # monovan run number 
        #

        prop['sum_runs'] = True # set to true to sum everything provided to sample_run
        #
        # Absolute units reduction properties. Set prop['monovan_run']=None to do arbitrary units
        prop['monovan_run'] = None # vanadium run in the same configuration as your sample
        #prop['sample_mass'] = 100 #  # mass of your sample
        #prop['sample_rmm'] = 200 #  molecular weight of scatters in your sample
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
      
      # Uncomment properties below and set them up to proper background run numbers
      # to remove empty instrument background
      # the run, containing empty instrument background to remove from sample run
      #prop['empty_bg_run'] = None
      # the run, containing empty instrument background to remove from detector vanadium run
      #prop['empty_bg_run_for_wb'] = None
      # the run, containing empty instrument background to remove from monovan run
      # if not defined, the empty_bg_run will be used instead
      #prop['empty_bg_run_for_monovan'] = None
      # the run, containing empty instrument background to remove from WB for monovan run
      # if WB for monovan defined and this one not defined, the "empty_bg_run_for_wb"
      # will be used instead
      #prop['empty_bg_run_for_monoWb'] = None
      #
      #
      # Uncomment two following properties to correct for absorption
      # on sample or sample container during the experiment.
      # 1) Define the sample material and sample shape:
      #
      #prop['correct_absorption_on'] = Cylinder(['Fe'],{'Height':10,'Radius':2})
      #
      #  The shapes currently available are:
      #  Cylinder, FlatPlate, HollowCylinder and Sphere
      # Each class accepts two parameters, namely dictionary or list describing
      # sample material, as accepted by SetSampleMaterial algorithm
      # and Shape parameters as accepted by SetSample algorithm.

      # Go to iPython window, type from AbsorptionShapes import * and type
      # help(Cylinder),
      # help(Sphere) or help(anAbsorptionShape) to learn more about different
      # ways of setting the material and shape parameters:
      ##----------------------------
      # 2) Optionally provide additional parameters defining speed and accuracy of the
      #    absorption corrections algorithm in abs_corr_info dictionary.
      #
      #prop['abs_corr_info'] = {'is_mc':True}
      #
      #Two generic algorithms are currently available to
      #    correct for absorption:
      # a) MonteCarloAbsorption (invoked by 'is_mc':True key of abs_corr_info
      #    property and by default)
      # and
      # b) AbsorptionCorrections (invoked by 'is_fast':True key of
      #    abs_corr_info property). This algorithm has number of optimizations
      #    for the special shapes cases and need further scientific validations

      # All other properties, to be provided in the input dictionary of this
      # property are the non-sample related properties of the
      # correspondent correction algorithm.
      #
      # Test the speed and the accuracy of the selected algorithm
      # reducing all your data by running eval_absorption_corrections method below.
      ##-------------------------------
      return prop
      #
#------------------------------------------------------------------------------------#
   @iliad
   def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file
          Modify only if custom pre or post-processing is needed, namely:
      """
      """
        Define custom preprocessing procedure, to be applied to the whole
        run, obtained from the instrument before diagnostics is applied
        1) Early update the cycle variable, as the procedure is most likely invoked from
           the cycle over different file names:
        if (not input_file is None):
           self.reducer.sample_run = str(input_file)
        
        2) Get access to the pointer to the input workspace
           (the workspace will be loaded  and calibrated at this point if it has not been done yet)
           using sample_run property class:
        run_ws = PropertyManager.sample_run.get_workspace()
        
        3) Get access to any property defined in Instrument_Properties.xml file
            or redefined in the reduction script:
        properties = self.reducer.prop_man
            e.g:
        RunNumber = properties.sample_run
        ei = properties.incident_energy
        
        Perform custom preprocessing procedure (with the values you specified)
        preprocessed_ws = custom_preprocessing_procedure(run_ws,RunNumber,ei,...)
        
        4) Store preprocessed workspace in the sample_run property for further analysis from preprocessing
           for the case where the workspace name have changed in the Analysis Data Service during preprocessing
           (this situation is difficult to predict so better always do this operation)
        PropertyManager.sample_run.synchronize_ws(preprocessed_ws)
      """

      WS = ReductionWrapper.reduce(self,input_file,output_directory)

      """ Defined custom post-processing procedure, in the way, similar to
      the preprocessing procedure, using the WS pointer, returned by the reduction procedure
      above. If the run is reduced in multirep mode, the WS is the list of the reduced
      workspace pointers, so the procedure should be applied to each workspace.

      At this stage standard saving had already been done, so save your results if you need them
      in a future:
      for i in range(0:len(WS)):
        SaveNexus(WS[i],Filename = 'CustomFilenameN%d.nxs'.format(i))
      """
      return WS
      #

    def do_preprocessing(self,reducer,ws):
        """ Custom function, applied to each run or every workspace, the run is divided to
            in multirep mode
            Applied after diagnostics but before any further reduction is invoked.
            Inputs:
            self    -- initialized instance of the instrument reduction class
            reducer -- initialized instance of the reducer
                       (DirectEnergyConversion class initialized for specific reduction)
            ws         the workspace, describing the run or partial run in multirep mode
                       to preprocess

            By default, does nothing.
            Add code to do custom preprocessing.
            Must return pointer to the preprocessed workspace
        """
        return ws
      #
    def do_postprocessing(self,reducer,ws):
        """ Custom function, applied to each reduced run or every reduced workspace,
            the run is divided into, in multirep mode.
            Applied after reduction is completed but before saving the result.

            Inputs:
            self    -- initialized instance of the instrument reduction class
            reducer -- initialized instance of the reducer
                       (DirectEnergyConversion class initialized for specific reduction)
            ws         the workspace, describing the run or partial run in multirep mode
                       after reduction to postprocess
            By default, does nothing.
            Add code to do custom postprocessing.
            Must return pointer to the postprocessed workspace.
            The postprocessed workspace should be consistent with selected save method.
            (E.g. if you decide to convert workspace units to wavelength, you can not save result as nxspe)
        """
        return ws
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
    def eval_absorption_corrections(self,test_ws=None):
        """ The method to evaluate the speed and efficiency of the absorption corrections procedure,
            before applying your corrections to the whole workspace and all sample runs.

            The absorption correction procedure invoked with excessive accuracy can run for too
            long providing no real improvements in accuracy. This is why it is recommended to
            run this procedure evaluating absorption on selected detectors and
            deploy the corrections to the whole runs only after achieving satisfactory accuracy
            and execution time.
            The procedure evaluate and prints the expected time to run the absorption corrections
            on the whole run.

            Input:
            If provided, the pointer or the name of the workspace available in analysis data service.
            If it is not, the workspace is taken from PropertyManager.sample_run property

            Usage:
            Reduce single run and uncomment this method in the __main__ area to evaluate
            adsorption corrections.
            Change adsorption corrections parameters below to achieve best speed and
            acceptable accuracy
        """
        # Gain access to the property manager:
        propman =  rd.reducer.prop_man
        # Set up Sample as one of:
        # 1) Cylinder([Chem_formula],[Height,Radius])
        # 2) FlatPlate([Chem_formula],[Height,Width,Thick])
        # 3) HollowCylinder([Chem_formula],[Height,InnerRadius,OuterRadius])
        # 4) Sphere([[Chem_formula],Radius)
        # The units are in cm
        propman.correct_absorption_on = Cylinder('Fe',[10,2]) # Will be taken from def_advanced_properties
        #                                prop['correct_absorption_on'] =  if not defined here
        #
        # Use Monte-Carlo integration.  Take sparse energy points and a few integration attempts
        # to increase initial speed. Increase these numbers to achieve better accuracy.
        propman.abs_corr_info = {'EventsPerPoint':3000}#,'NumberOfWavelengthPoints':30}
        # See MonteCarloAbsorption for all possible properties description and possibility to define
        # a sparse instrument for speed.
        #
        # Gain access to the workspace. The workspace should contain Ei log, containing incident energy
        # (or be reduced)
        if test_ws is None:
            test_ws = PropertyManager.sample_run.get_workspace()
        # Define spectra list to test absorption on. Empty list will
        # define absorption on the whole workspace.
        check_spectra = [1,200]
        # Evaluate corrections on the selected spectra of the workspace and the time to obtain
        # the corrections on the whole workspace.
        corrections,time_to_correct_abs = self.evaluate_abs_corrections(test_ws,check_spectra)
        # When accuracy and speed of the corrections is satisfactory, copy chosen abs_corr_info
        # properties from above to the advanced_porperties area to run in reduction.
        if mpl is not None:
            n_spectra = len(check_spectra)
            if n_spectra == 0:
                n_specra = corrections.getNumberHistograms()
            mpl.plotSpectrum(corrections,range(0,n_spectra))
        #
        return corrections

#------------------------------------------------------------------------------------------------
    def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'LET',web_var)
       Mt = MethodType(self.do_preprocessing, self.reducer)
       DirectEnergyConversion.__setattr__(self.reducer,'do_preprocessing',Mt)
       Mt = MethodType(self.do_postprocessing, self.reducer)
       DirectEnergyConversion.__setattr__(self.reducer,'do_postprocessing',Mt)
#-------------------------------------------------------------------------------------------------#


def main(input_file=None,output_directory=None):
    """ This method is used to run code from web service
        and should not be touched unless you change the name of the
        particular ReductionWrapper class (e.g. ReduceMARI here)
        exception to change the output folder to save data to
    """
    # note web variables initialization
    rd = ReduceLET(web_var)
    rd.reduce(input_file,output_directory)
    # Define folder for web service to copy results to
    output_folder = ''
    return output_folder


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
####get reduction parameters from properties above, override what you want locally ###
   # and run reduction. Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g. 
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or 
   ## rd.reducer.prop_man.sum_runs = False
   # 
   # 
    rd.run_reduction()

###### Test absorption corrections to find optimal settings for corrections algorithm         ######
#     corr = rd.eval_absorption_corrections()




