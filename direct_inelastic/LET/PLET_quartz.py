import PolCorr
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
reload(PolCorr)
from mantid import *

# This scripts reads reduced One2One .nxspe files of a sequence of NSF/SF 
# quartz calibration runs, for a sequence of LET energy reps, and uses the
# out-of-plane path length dependece of the polarization efficency of the He3 cell
# to calculate the Polarizer/Flipper efficency - and the "average" cell polarization
# Seems to work best when using quartz runs over a short time (not too much 
# variation in cell P and T  

# The calculated Polarizer/Flipper efficency PF should be used for 
# the data reduction in PLET_reduction_sample.py

sample_runs=list(range(57810,57820))
eis = [6.13,3.20,1.96] 
pressure = 0.9

quartz = PolCorr.Reduce(sample_runs,eis,[],he_pressure=pressure)

quartz.get_PF_from_quartz()


#print quartz.__dict__


