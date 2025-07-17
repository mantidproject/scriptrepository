# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from os import path, makedirs, scandir
from Engineering.texture.TextureUtils import find_all_files, fit_all_peaks

############### ENGINEERING DIFFRACTION INTERFACE FITTING ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PostExp-ZrRingDiagScript"

# otherwise set root directory here:
root_dir = fr"C:\Users\kcd17618\Engineering_Mantid\User\{exp_name}"

# Next the folder contraining the workspaces you want to fit
file_folder = "Focus"
# These are likely within a sub-folder specified by the detector grouping
grouping = "Texture30"
prm_path = None
groupingfile_path = None

# You also need to specify a name for the folder the fit parameters will be saved in
fit_save_folder = "ScriptFitParameters"

# Finally, provide a list of peaks that you want to be fit within the spectra
#peaks = [2.03,1.44, 1.17, 0.91] # steel
peaks = [2.8, 2.575, 2.455, 1.89, 1.62, 1.46] # zr

######################### RUN SCRIPT ########################################

# create output directory
fit_save_dir = path.join(root_dir, fit_save_folder)
mk(fit_save_dir)

# find and load peaks

# get grouping directory name
calib_info = CalibrationInfo(group = GROUP(grouping))
if groupingfile_path:
    calib_info.set_grouping_file(groupingfile_path)
elif prm_path:
    calib_info.set_prm_filepath(prm_path) 
group_folder = calib_info.get_group_suffix()
focussed_data_dir = path.join(root_dir, file_folder, group_folder, "CombinedFiles")
focus_ws_paths = find_all_files(focussed_data_dir)
focus_wss = [path.splitext(path.basename(fp))[0] for fp in focus_ws_paths]
for iws, ws in enumerate(focus_wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_ws_paths[iws], OutputWorkspace= ws)

# execute the fitting                     
fit_all_peaks(focus_wss, peaks, 0.02, fit_save_dir)     
