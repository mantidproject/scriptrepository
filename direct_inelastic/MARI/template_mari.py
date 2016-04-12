from mantid import config
from MARIReduction_Sample import *
import time
import os
import datetime
import sys
try:
    #Note: due to the mantid-python implementation, one needs to run this 
    #script in Mantid  script window  TWICE!!!  to deploy the the changes made to MARIReduction_Sample.py file.
    reload(sys.modules['MARIReduction_Sample'])
except:
    print "*** WARNING can not reload MARIReduction_Sample file"
    pass

# Run number and Ei
runno=20942
sum_runs=False
ei=30

# White vanadium run number
wbvan=21334
# Default save directory
config['defaultsave.directory'] = '/instrument/MARI/RBNumber/USER_RB_FOLDER' #data_dir 

# Absolute normalisation parameters
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=0
sam_mass=0
sam_rmm=0

# If necessary, add any sequence of reduction paramerters defined in MARIParameters.xml file 
# to the end ot the illiad string using the form: property=value 
# (e.g.:  iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs,check_background=False)
iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs)

