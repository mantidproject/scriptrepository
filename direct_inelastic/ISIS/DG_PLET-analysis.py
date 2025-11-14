#================================================================================
#   PLET reduction script v3
#   Uses PolCorr3 routines
#
#   Requires initial LET reduction using 1to1 mapping and .nxs files (not .nxspe)
#
#                                               JRS 12/11/25
#=================================================================================


import sys, os
sys.path.append(os.path.dirname(__file__))
from reduction_utils import Reduce
from mantid.simpleapi import AnalysisDataService as ADS

ADS.clear() # This clear all the workspaces (temporary fix for memory leak)

#======================================USER INPUT================================
rbno = 'USER_RB_FOLDER'    # RB number
NSF_first = False     # this flag should be true if you measured with the flipper off first (check JournalViewer)
pressure  = 0.75      # the pressure of the 3He gas
separate  = True      # this decides wether to output NSF/SF or COH/INC
he_mode   = 'fit'     # options: 'fit' or 'direct' - use 'fit' for set_helium_parameters
eis       = [8.61, 3.60, 1.97, 1.24]    # incident energies
PF        = [0.90, 0.92, 0.92, 0.92]    # Polarizer*Flipper efficiency
rings_map = 'RINGS_MAP_XML'             # Set to None for 1to1 output
#directory where reduced 1to1 .nxs files are stored
nxsdir    = f'/data/analysis/LET/RBNumber/RB{rbno}/'      

lbl='Quartz'
sample_runs = range(105952,105960)
cell = "3HeCal_105962-106023.txt"
   
#=================================END OF USER INPUT===================================
data = Reduce(sample_runs,
              eis,
              PF,
              he_pressure=pressure,
              he_mode=he_mode,
              NSF_first=NSF_first,
              separate=separate,
              label=lbl,
              rings_map=rings_map,
              nxsdir=nxsdir)
                       
#data.get_helium_parameters()
#data.set_helium_parameters(PHe0 = PHe, T1 = T1, T0run = T0run)
data.set_helium_parameters(cal=cell)
data.correct_data()
data.components()
data.output(type="nxspe")

#=====================================================================================
