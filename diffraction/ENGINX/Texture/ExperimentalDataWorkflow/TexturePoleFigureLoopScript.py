# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from Engineering.texture.TextureUtils import find_all_files, create_pf_loop, get_xtal_structure
from Engineering.common.calibration_info import CalibrationInfo
from Engineering.EnggUtils import GROUP
import os

############### ENGINEERING DIFFRACTION INTERFACE POLE FIGURE ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################
# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PostExp-SteelCentre"

# otherwise set root directory here:
save_root = r"C:\Users\kcd17618\Engineering_Mantid"
root_dir = fr"{save_root}\User\{exp_name}"


ws_folder = "Focus"
fit_save_folder = "ScriptFitParameters-New"
# define the peaks of interest, NOTE these must correspond to sub folders in the fit directory
peaks = [2.03,1.44, 1.17, 0.91]
# define the columns you would like to create pole figures for
readout_columns = ["I", "X0"]
# you need to specify the detector grouping
grouping = "Texture30"
# and some grouping path if not using a standard
prm_path = None
groupingfile_path = None
# and the type of projection to plot
projection_method = "Azimuthal"

# you need to define the orientation of the intrinsic sample directions when the sample orientation matrix == I (no rotation)
# this should be the same as the reference state used in the absorption correction
#r2 = np.sqrt(2)/2
dir1 = np.array((1,0,0))
dir2 = np.array((0,1,0)) # projection axis
dir3 = np.array((0,0,1))
# you can also supply names for these three directions
dir_names = ["RD", "ND", "TD"]

# set whether you would like the plotted pole figure to be a scatter of experimental points or whether you would like to apply gaussian smoothing and
# plot a contour representation
scatter = "both"
# if contour, what should the kernel size of the gaussian be
kernel = 6.0

# do you want to include a scattering power correction
include_scatt_power = False
# if so what is the crystal structure, defined either by giving a cif file or supplying the lattice, space group and basis
xtal_input = None # "cif"/"array"/"string"
xtal_args = [] # for input "cif", require the cif filepath, for "array" array of lattice parameters, space group, basis
# for "string" lattice parameter string, space group and basis

# if you have set a crystal, you can also provide a set of hkls, the hkl_peaks dictionary is a useful way of assigning the peaks
hkl_peaks = {1.17: (1,1,2),1.43: (2,0,0),2.03: (1,1,0)} #Fe

chi2_thresh = 0.0   # max value of Chi^2 to be included as a point in the table
peak_thresh = 0.01   # max difference from either the HKL specified or the mean X0
scat_vol_pos = (0.0,0.0,0.0) # for now, can assume the gauge vol will be centred on origin

######################### RUN SCRIPT ########################################


# get grouping directory name
calib_info = CalibrationInfo(group = GROUP(grouping))
if groupingfile_path:
    calib_info.set_grouping_file(groupingfile_path)
elif prm_path:
    calib_info.set_prm_filepath(prm_path) 
group_folder = calib_info.get_group_suffix()
focussed_data_dir = os.path.join(root_dir, ws_folder, group_folder, "CombinedFiles")
focus_ws_paths = find_all_files(focussed_data_dir)
focus_wss = [os.path.splitext(os.path.basename(fp))[0] for fp in focus_ws_paths]
for iws, ws in enumerate(focus_wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_ws_paths[iws], OutputWorkspace= ws)

fit_load_dirs = [os.path.join(root_dir, fit_save_folder, group_folder, str(peak)) for peak in peaks]

hkls = [hkl_peaks[peak] if include_scatt_power else None for peak in peaks]

xtal = get_xtal_structure(xtal_input, *xtal_args) if xtal_input else None

fit_param_wss = []
for ifit, fit_folder in enumerate(fit_load_dirs):
    # get fit params
    fit_dir = os.path.join(root_dir, fit_folder)
    fit_wss = find_all_files(fit_dir)
    param_wss = [os.path.splitext(os.path.basename(fp))[0] for fp in fit_wss]
    fit_param_wss.append(param_wss)
    for iparam, param in enumerate(param_wss):
        if not ADS.doesExist(param):
            Load(Filename=fit_wss[iparam], OutputWorkspace=param)

create_pf_loop(wss = focus_wss,
               param_wss = fit_param_wss,
               include_scatt_power = include_scatt_power, 
               xtal = xtal,
               readout_columns = readout_columns, 
               hkls = hkls,
               dir1 = dir1, 
               dir2 = dir2, 
               dir3 = dir3, 
               dir_names = dir_names, 
               scatter = scatter,
               kernel = kernel, 
               scat_vol_pos = scat_vol_pos,
               chi2_thresh = chi2_thresh, 
               peak_thresh = peak_thresh, 
               save_root = save_root, 
               exp_name = exp_name, 
               projection_method = projection_method)

