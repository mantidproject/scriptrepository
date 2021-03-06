""" Sample HET reduction scrip  """
import os

from Direct.ReductionWrapper import *
try:
    import reduce_vars as web_var
except:
    web_var = None




class ReduceHET(ReductionWrapper):
#-------------------------------------------------------------------------------------------------#
    @MainProperties
    def def_main_properties(self):

       """ Define main properties used in reduction """
       prop = {}
       prop['sample_run'] = 5240
       prop['wb_run'] = 5223
       prop['incident_energy'] = 150
       prop['energy_bins'] = [-15,1.,145]

      # Absolute units reduction properties.
       prop['monovan_run'] = 5263
       prop['sample_mass'] = 10
       prop['sample_rmm'] = 79.55
       return prop
#-------------------------------------------------------------------------------------------------#
    @AdvancedProperties
    def def_advanced_properties(self):
      """  separation between simple and advanced properties depends
           on scientist, experiment and user.
           main properties override advanced properties.
      """
      prop = {}
      prop['map_file'] = 'HET_RingsNoPSD_2017.map' #
      prop['monovan_mapfile'] = "HET_WBVanGrouping_NoPSD_2017.map"
      prop['hard_mask_file'] = None #"mar11015.msk"
      prop['det_cal_file'] = 'HET_DETECTORS_CalFile.nxs'
      prop['save_format'] = '.nxspe'
      return prop
      #
#-------------------------------------------------------------------------------------------------#
    @iliad
    def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file
          Overload only if custom reduction is needed
      """
      ws = ReductionWrapper.reduce(self,input_file,output_directory)
      #SaveNexus(ws,Filename = 'MARNewReduction.nxs')
      return ws
    #
    def set_custom_output_filename(self):
        """ define custom name of output files if standard one is not satisfactory
          In addition to that, example of accessing reduction properties
          Changing them if necessary
        """
        def custom_name(prop_man):
            """Sample function which builds file-name from
              incident energy and run number and adds some auxiliary information
              to it.
            """
            # Note -- properties have the same names as the list of advanced and
            # main properties
            ei = prop_man.incident_energy
            # sample run is more then just list of runs, so we use
            # the formalization below to access its methods
            run_num = PropertyManager.sample_run.run_number()
            name = "RUN{0}atEi{1:<3.2f}meV_One2One".format(run_num ,ei)
            return name

        # Uncomment this to use custom file-name function
        # Note: the properties are stored in prop_man class accessed as
        # below.
        #return lambda : custom_name(self.reducer.prop_man)
        # use this method to use standard file name generating function
        return None
      

    def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'HET',web_var)
       

if __name__ == "__main__":
#-------------------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#-------------------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results             #####
    # It can be done here or from Mantid GUI
    # Folder where map files are located:
     #map_mask_dir = 'd:/Data/MantidSystemTests/Data'
    # folder where input data can be found
     #data_dir = 'd:/Data/Mantid_Testing/14_11_27'
     # auxiliary folder with results
     #ref_data_dir = 'd:/Data/MantidSystemTests/SystemTests/AnalysisTests/ReferenceResults'
     # Set input path to
     #config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,map_mask_dir,ref_data_dir))
     # use appendDataSearch directory to add to existing data search path
     #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
     # folder to save resulting spe/nxspe files.
     #config['defaultsave.directory'] = data_dir

###### Initialize reduction class above and set up reduction properties. Note no parameters  ######
     rd = ReduceHET()
     # set up advanced and main properties
     rd.def_advanced_properties()
     rd.def_main_properties()

     # uncomment rows below to save web variables to use in web services.
     #run_dir = os.path.dirname(os.path.realpath(__file__))
     #file = os.path.join(run_dir,'reduce_vars.py')
     #rd.save_web_variables(file)

     # Web-like reduction over sequence of files
     #rd.reduce()
###### Run reduction on all files provided as parameters ######
     red_ws = rd.run_reduction()
