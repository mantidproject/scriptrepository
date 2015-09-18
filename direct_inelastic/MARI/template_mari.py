from mantid import config
from iliad_mari import *
import time
import os
import datetime

# Run number and Ei
runno=[19891]
sum_runs=False
ei=110

# White vanadium run number
wbvan=20465
# Default save directory
config['defaultsave.directory'] = 'home_folder'          #data_dir 

# Absolute normalisation parameters
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=0
sam_mass=0
sam_rmm=0


iliad_mar(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs)

