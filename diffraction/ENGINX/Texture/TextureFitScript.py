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

# You also need to specify a name for the folder the fit parameters will be saved in
fit_save_folder = "ScriptFitParameters"

# Finally, provide a list of peaks that you want to be fit within the spectra
#peaks = [2.03,1.44, 1.17, 0.91] # steel
peaks = [2.8, 2.575, 2.455, 1.89, 1.62, 1.46] # zr

######################### RUN SCRIPT ########################################

# create output directory
save_dir = path.join(root_dir, fit_save_folder)
if not path.exists(save_dir):
    makedirs(save_dir)

# find and load peaks
focussed_data_dir = path.join(root_dir, file_folder, grouping, "CombinedFiles")
focus_wss = find_all_files(focussed_data_dir)
wss = [path.splitext(path.basename(fp))[0] for fp in focus_wss]
for iws, ws in enumerate(wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_wss[iws], OutputWorkspace= ws)

# execute the fitting                     
fit_all_peaks(wss, peaks, 0.02, save_dir)
