from __future__ import print_function
from mantid import config
from MARIReduction_Sample import *
import time
import os
import datetime
import sys
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
try:
    #Note: due to the mantid-python implementation, one needs to run this 
    #script in Mantid  script window  TWICE!!!  to deploy the the changes made to MARIReduction_Sample.py file.
    sys.path.insert(0,'/instrument/MARI/RBNumber/USER_RB_FOLDER/')
    reload(sys.modules['MARIReduction_Sample'])
except:
    print("*** WARNING can not reload MARIReduction_Sample file")
    pass

# Run number and Ei
runno=26644
sum_runs=False
ei=[30, 11.8]

# White vanadium run number
wbvan=28041
# Default save directory
config['defaultsave.directory'] = '/instrument/MARI/RBNumber/USER_RB_FOLDER' #data_dir 

# Absolute normalisation parameters
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=0
sam_mass=0
sam_rmm=0

# Set to true to remove the constant ToF background from the data.
remove_bkg = True

# If necessary, add any sequence of reduction paramerters defined in MARIParameters.xml file 
# to the end ot the illiad string using the form: property=value 
# (e.g.:  iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs,check_background=False)
iliad_mari(runno, ei, wbvan, monovan, sam_mass, sam_rmm, sum_runs, check_background=remove_bkg,
    hard_mask_file='MASK_FILE_XML')

# To run reduction _and_ compute density of states together uncomment this and comment iliad_mari above
# bkgruns and runno can be lists, which means those runs will be summed, and the sum is reduced
#bkgruns = 20941
#iliad_dos(runno, wbvan, ei, monovan, sam_mass, sam_rmm, sum_runs, background=bkgrun, temperature=5)
