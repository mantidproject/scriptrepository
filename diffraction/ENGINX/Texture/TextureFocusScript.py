# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
from Engineering.texture.TextureUtils import find_all_files, run_focus_script

############### ENGINEERING DIFFRACTION INTERFACE FOCUS ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.

exp_name = "PostExp-ZrRingDiagScript"

# otherwise set root directory here:
root_dir = fr"C:\Users\kcd17618\Engineering_Mantid\User\{exp_name}"

# next, specify the folder with the files you would like to focus 
# (if you are using the standard scripts this might not need to change)
data_dir = fr"{root_dir}\AbsorptionCorrection"

# fill in the file paths for the vanadium and ceria runs (just run numbers might work if you are setup into the file system)
van_run = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day1\SteelDataset\DataFiles\ENGINX00361838.nxs"
ceria_run = "305738"

# set the path to the grouping file created by calibration
prm_path = None # fr"{root}\Calibration\ENGINX_305738_Texture30.prm"
grouping = "Texture30"

# Define some file paths
full_instr_calib = r"C:\Users\kcd17618\Documents\dev\mantid\mantid\scripts\Engineering\calib\ENGINX_full_instrument_calibration_193749.nxs"

######################### RUN SCRIPT ########################################

run_files = find_all_files(data_dir)

run_focus_script(wss = run_files, 
                 focus_dir = root_dir, 
                 van_run = van_run, 
                 ceria_run = ceria_run,
                 full_instr_calib = full_instr_calib, 
                 grouping = grouping,
                 prm_path = prm_path)

