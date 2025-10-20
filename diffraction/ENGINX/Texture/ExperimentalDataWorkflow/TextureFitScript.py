# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from os import path, makedirs, scandir
from Engineering.texture.TextureUtils import find_all_files, fit_all_peaks, mk
from Engineering.common.calibration_info import CalibrationInfo
from Engineering.EnggUtils import GROUP

############### ENGINEERING DIFFRACTION INTERFACE FITTING ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PostExp-SteelCentre"

# otherwise set root directory here:
root_dir = fr"C:\Users\kcd17618\Engineering_Mantid\User\{exp_name}"

# Next the folder contraining the workspaces you want to fit
file_folder = "Focus"
# These are likely within a sub-folder specified by the detector grouping
grouping = "Texture30"
prm_path = None
groupingfile_path = None

# You also need to specify a name for the folder the fit parameters will be saved in
fit_save_folder = "ScriptFitParameters-FitTest"

# Provide a list of peaks that you want to be fit within the spectra
peaks = [2.03,1.44, 1.17, 0.91] # steel
#peaks = [2.8, 2.575, 2.455, 1.89, 1.62, 1.46] # zr

# The fitting has a couple of parameters that deal with when peaks are missing as a result of the texture
# The first parameter is 1_over_sigma_thresh - this determines the minimum value of I/sigma for a fit to be considered as for a valid peak
# any invalid peak will have parameters set to nan by default, but these nans can be overwritten by no_fit_value_dicts and nan_replacement
# no_fit_value_dict takes fitted parameter names and allows you to specify what the unfit value should be eg. {"I":0.0} - if you can't fit intensity
# set the value directly to 0.0
# nan_replacement then happens after this, if a nan_replacement method is given any parameters without an unfit_value provided will have nans replaced
# either with "zeros", or with the min/max/mean value of that parameter (Note: if all the values are nan, the value will remain nan)

i_over_sigma_thresh = 3.0
no_fit_value_dict = {"I": 0.0, "I_est": 0.0}
nan_replacement = "mean"

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
focus_ws_paths = find_all_files(focussed_data_dir)[:3]
focus_wss = [path.splitext(path.basename(fp))[0] for fp in focus_ws_paths]
for iws, ws in enumerate(focus_wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_ws_paths[iws], OutputWorkspace= ws)


# execute the fitting                     
fit_all_peaks(focus_wss, peaks, 0.02, fit_save_dir, i_over_sigma_thresh = i_over_sigma_thresh, nan_replacement = nan_replacement, no_fit_value_dict = no_fit_value_dict)     







