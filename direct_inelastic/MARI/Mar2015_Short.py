import os
""" Sample MARI reduction script used  only locally)""" 
from Direct.ReductionWrapper import *

class ReduceMARIFromFile(ReductionWrapper):
    @MainProperties
    def def_main_properties(self):
        """ Define main properties used in reduction """ 
        prop = {}
        prop['sample_run'] = [19184,19185]
        prop['incident_energy'] = 15
        prop['energy_bins'] = [-10,0.03,15]

        #prop['sum_runs'] = False

        # Absolute units reduction properties.
        #prop['monovan_run'] = None
        #prop['sample_mass'] = 10
        #prop['sample_rmm'] = 435.96
        return prop

    @AdvancedProperties
    def def_advanced_properties(self):
        """  separation between simple and advanced properties depends
           on scientist, experiment and user.
           main properties override advanced properties.      
        """
        prop = {}
        prop['wb_run'] = 18622
        prop['map_file'] = "mari_res2013"
        prop['monovan_mapfile'] = "mari_res2013"
        prop['hard_mask_file'] = "mari_mask2014.msk"
        prop['det_cal_file'] = 'MAR18622.raw'
        prop['save_format'] = 'nxspe,spe'
        #Default for MARI
        #prop['data_file_ext']='.raw' # if two input files with the same name and
                                    #different extension found, what to prefer.
        return prop
      #
    @iliad
    def reduce(self,input_file=None,output_directory=None):
        """Method executes reduction over single file
         Overload only if custom reduction is needed
        """
        outWS = ReductionWrapper.reduce(self,input_file,output_directory)
        #SaveNexus(outWS,Filename = 'MARNewReduction.nxs')
        return outWS
 
    def __init__(self,web_var=None):
        """ sets properties defaults for the instrument with Name"""
        ReductionWrapper.__init__(self,'MAR',web_var)
#-------------------------------------------------------------------------------------------------#

if __name__ == "__main__":
    #data_dir = r'd:\Data\Mantid_Testing\15_03_01\data'
    #output_dir = r'd:\Data\Mantid_Testing\15_03_01'
    #config.setDataSearchDirs('{0};{1}'.format(data_dir,maps_dir,output_dir))
    maps_dir = r'/usr/local/mprogs/InstrumentFiles/mari'
    config.appendDataSearchDir(maps_dir)

    #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
    #config['defaultsave.directory'] = output_dir  # folder to save resulting spe/nxspe files.  Defaults are in

    # execute stuff from Mantid
    rd = ReduceMARIFromFile()
    #rd = ReduceMARIFromWorkspace()
    rd.def_advanced_properties()
    rd.def_main_properties()

#### Set up time interval (sec) for reducer to check for input data file.         ####
    #  If this file is not present and this value is 0,reduction fails 
    #  if this value >0 the reduction wait until file appears on the data 
    #  search path checking after time specified below.
    rd.wait_for_file = 0  # waiting time interval

###### Run reduction over all run numbers or files assigned to                   ######
    # sample_run  variable 
    rd.run_reduction()


