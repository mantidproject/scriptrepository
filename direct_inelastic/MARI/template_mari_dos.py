from __future__ import print_function
from mantid import config
from MARIReduction_Sample import *
import time
import os
import datetime
import sys
import numpy as np
if sys.version_info > (3,):
    if sys.version_info < (3,4):
        from imp import reload
    else:
        from importlib import reload
try:
    #Note: due to the mantid-python implementation, one needs to run this 
    #script in Mantid  script window  TWICE!!!  to deploy the the changes made to MARIReduction_2018_4.py file.
    sys.path.insert(0,'/instrument/MARI/RBNumber/USER_RB_FOLDER/')
    reload(sys.modules['MARIReduction_Sample'])
except:
    print("*** WARNING can not reload MARIReduction_Sample file")
    pass
    
# START EDITING HERE

# Incident energy(ies) of the run(s)
ei=[120,10]

# Energy bins and q-range for DOS calculation. Can be strings of numbers or lists. If using multiple Ei, can give
# list of strings or list of lists, with each element corresponding to each ei.
ebins = ['2, 1, 110', '0.5, 0.25, 9']
qrange = [['6,14'], ['2,4']]

# set the following to true to load the reduced data rather than run the reduction
load_reduce = True

# List of runs as a dictionary of dictionaries. The keys for the outer dictionary is the sample name.
# Keys for the inner dictionary are sample specific parameters (mass, molar mass) and the temperature index.
# In the following example, sample 1 was measured at 5K only and sample 2 at 5 and 300K.
# For each sample and temperature, if 'data' or 'background' is a list, all runs in the list will
# be summed together.
# The 'recalc' option: If True, forces reruns the reduction if load_reduce=False. (load_reduce takes precedent)
# E.g. if the workspace exists in the workspace list and 'recalc': False then the reduction will not be rerun and
# the script will use the workspace in memory.
# ssf is the self-shielding factor and msd is the mean-square displacement per sample and temperature
runs = {        
    'Sample_1': {'sam_mass':0, 'sam_rmm':0,
                           5: {'data': list(range(25478,25485)), 'background': list(range(25492,25500)), 'recalc':True, 'ssf':1., 'msd':0.}, 
                     },
    'Sample_2': {'sam_mass':0, 'sam_rmm':0,
                           5: {'data': list(range(25506,25516)), 'background': list(range(25492,25500)), 'recalc':True, 'ssf':1., 'msd':0.}, 
                           300: {'data': list(range(25526,25536)), 'background': list(range(25501,25510)), 'recalc':True, 'ssf':1., 'msd':0.}, 
                     },
}

# Set this to the monochromatic run number for absolute unit normalisation (or zero to ignore) 
monovan = 0

# set to 0 for no smoothing otherwise must be odd, and greater than or equal to 3.
nsmooth = 3

# set to True to put output files for each set of samples in a different folder
use_subdirs = True

# set to True to save text files of the calculated DOS
save_text = True

# END EDITING HERE


#White vanadium run number
wbvan=25779
#Default save directory
datadir = '/instrument/MARI/RBNumber/USER_RB_FOLDER/' #data_dir 

#Set to true to remove the constant ToF background from the data.
remove_bkg = False

# To get help
#print(help(iliad_dos))

iliad_dos(runs, wbvan, ei, monovan, check_background=remove_bkg, load_reduce=load_reduce,
         save_text=save_text, use_sub_directories=use_subdirs, nsmooth=nsmooth, save_folder=datadir)

# Alternative syntax without using the dictionaries, for a single sample at a single temperature
#iliad_dos([25478, 25479], ei=[120, 10], wbvan=25035, background=[25492, 25493], temperature=5)

