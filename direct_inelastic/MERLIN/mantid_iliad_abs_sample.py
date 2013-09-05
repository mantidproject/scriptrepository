"""
sample Direct inelastic reduction for MERLIN performed in absolute units
"""
from qtiGenie import *
#from PySlice2 import *

# Calculates 

inst='MER'
iliad_setup(inst)

#mapfile = 'rings_123.map'

# where to save resutls (usually specified in Mantid)
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
# map mask and cal files
config.appendDataSearchDir('../../../InstrumentFiles/merlin') 
# data (raw or nxs) run files
config.appendDataSearchDir('d:/Data/isis/Adroja') 
    
#load vanadium file    
whitebeamfile="13123"
LoadRaw(Filename=whitebeamfile,OutputWorkspace="wb_wksp",LoadLogFiles="0")
MonoVanWB="wb_wksp"
# set up single calibration file to use with all workspaces
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Mandatory positional parameters 
ei=25  # ei-guess
rebin_params=[-10,.1,23]
MonoVanRun=[13271] # mandatory for absolute units, for arbitrary should be None or absent in Iliad

#   Other positional  parameters
mapfile='one2one_123' # ring map file is used for powder.  if absend idf file value is used instead
# key-coded parameters
params={}
params['monovan_integr_range']=[-25,23]
params['norm_method']='current'
params['det_cal_file']='det_corr_123.dat'  #det_cal_file must be specified if the reduction sends out put to a workpsace
params['sample_mass']=32.62
params['sample_rmm']=50.94
params['monovan_mapfile']='rings_123.map' # good 
params['hardmaskPlus']='w1.msk'



#runs=[13151 13152 13153 13154 13155 13156 13157 13158 13159 13160 13161 13162 13163 13164 13165 13166 13167 13168 13169 13170 13171 13172 13173 13174 13175 13176 13177]
runs=[13271]
############## normal reduction####################

#save .nxspe file
for runfile in runs:
    save_file=inst+str(runfile)+'.nxspe'
    LoadRaw(Filename=str(runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")

    w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,MonoVanRun,**params)

    SaveNXSPE('w1',save_file)
print "All done"



