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

ADS.clear()

# This script needs to be executed for each set of runs using one 3He cell 

#======================================USER INPUT================================
NSF_first = False   # this flag should be true if you measured with the flipper off first (check JournalViewer)
pressure = 0.75     # the pressure of the 3He gas
eis = [1.97]        # incident energy
PF  = [0.92]        # Polarizer*Flipper efficiency

sample_runs= range(105962,106024)   # cell #1, PHe0 = 0.529, T1 = 54.2 hours

#=================================END OF USER INPUT===================================
data = Reduce(sample_runs,
              eis,
              PF,
              he_pressure=pressure,
              he_mode='fit',
              NSF_first=NSF_first)          
data.get_helium_parameters(save=True)


