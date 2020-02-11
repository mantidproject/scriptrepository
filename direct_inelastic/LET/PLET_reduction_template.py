import PolCorr
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
reload(PolCorr)
from mantid import *

# This script reads reduced One2One .nxspe files and corrects / combines the data.
# Input files should be a sequence of NSF and SF runs using only one spin-filter
# data output as rings-mapped .nxspe files

NSF_first = True                    # set to False if flipper on is first run in sequence
pressure = 0.75                     # 3He pressure in banana cell
eis = [3.84,1.81,1.05]              # energy reps to reduce
PF = [0.95 - x*0.047 for x in eis]  # Polarizer/Flipper efficency for each ei.  
                                    # This can be found from a quartz run (see PLET_quartz.py)
# other possible arguments are:
#   he_path_length  - dafualt 0.06 m
#   polmon_distance - default 25.38 m (monitor in PLET coils in spectrum 7)
#   polmon_delay    - defualt 100 microsec
#   mask            - default 'PLET_184_msk.xml'

sample_runs=list(range(58310,58350)); lbl='PEO-D_375K'  # LET12

data = PolCorr.Reduce(sample_runs,eis,PF,he_pressure=pressure,NSF_first=NSF_first,label=lbl)

data.get_helium_parameters()    # Calculates the cell P_he(0) and T1 using the PLET monitor (uses first energy rep)

# data.get_helium_from_ROI(26700, 27000, roi="MaskWorkspace")   # uses Bragg peak defined by min and max TOF and ROI mask to calculate the cell paramters (uses first energy rep)

data.correct_data()             # corrects the data for finite polarization and energy dependent cell transmission

data.components()               # combines the NSF and SF components to produce Coherent and Incoherent XS

data.one2one_output()           # outputs the combined data (if components is run) or the NSF/SF data (if not)
data.rings_output()             # in rings or one2one format


#print test.__dict__

