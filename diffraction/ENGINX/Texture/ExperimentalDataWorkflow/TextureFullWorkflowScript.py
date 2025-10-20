from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from Engineering.texture.TextureUtils import find_all_files, mk, run_abs_corr, run_focus_script,  fit_all_peaks, create_pf_loop, get_xtal_structure
from Engineering.common.calibration_info import CalibrationInfo
from Engineering.EnggUtils import GROUP
import os

############### ENGINEERING DIFFRACTION INTERFACE ABSORPTION CORRECTION ANALOGUE #######################

# NOTE: for parameters that are not required, set as None
# In the example script some of these may already be set to none, with an example non-None parameter given as a comment

######################### EXPERIMENTAL INFORMATION ########################################

# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PR-TEST-ZR"

# otherwise set root directory here:
save_root = r"C:\Users\kcd17618\Engineering_Mantid"
root_dir = fr"{save_root}\User\{exp_name}"

# next, specify the folder with the files you would like to apply the absorption correction to 
corr_dir = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\ZrRing\DataFiles\Point2"

# For texture, it is expected that you have a single sample shape, that is reorientated between runs.
# this is handled by having a reference workspace with the shape in its neutral position 
# (position in the beamline when the goniometer is home)
# This reference workspace probably requires you to do some interacting and validating, so should be setup in the UI
# (Interfaces/Diffraction/Engineering Diffraction/Absorption Correction)

# if this is the case copy ref should be True and the ref_ws_path should be given
# otherwise, if set ref is true, it is assumed that the sample shapes are already present on the workspaces
copy_ref = True
ref_ws_path = os.path.join(root_dir, "ReferenceWorkspaces", f"{exp_name}_reference_workspace.nxs") 

# if using the reference you now need to reorientate the sample, this can be done using orientation files
# two standard types

# Euler Orientation (orient_file_is_euler = True)
# for this, euler_scheme and euler_axes_sense must be given to say which lab frame directions the goniometer axes are pointing along
# and where the rotations are counter-clockwise (1) or clockwise (-1)

# Matrix Orientation (orient_file_is_euler = False)
# for this the first 9 values in each row of the files are assumed to be flattened rotation matrix. 
# These are used to directly reorientate the samples 
orientation_file = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\ZrRing\Sscanss\Split\Zirc_ring_pose_matrices_mantid_point_1.txt"
orient_file_is_euler = False
euler_scheme = None # "YXY"
euler_axes_sense = None # "1,-1,1"

# Now you can specify information about the correction
include_abs_corr = True # whether to perform the correction based on absorption
monte_carlo_args = "SparseInstrument:True" # what arguments to pass to MonteCarloAbsorption alg
clear_ads_after = True # whether to remove the produced files from the ADS to free up RAM 
gauge_vol_preset = "4mmCube" # or "Custom" # the gauge volume being used
gauge_vol_shape_file = None # or "path/to/xml" # a custom gauge volume shape file

# There is also the option to output an attenuation table alongside correcting the data
# This will return a table of the attenuation coefficient at the point specified
include_atten_table = False
eval_point = "2.00"
eval_units = "dSpacing" #must be a valid argument for ConvertUnits

# Finally, you can add a divergence correction to the data, this is still a work in progress, so keep False for now
include_div_corr = False
div_hoz = None # 0.02
div_vert = None # 0.02
det_hoz = None # 0.02

############### ENGINEERING DIFFRACTION INTERFACE FOCUS ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# next, specify the folder with the files you would like to focus 
# (if you are using the standard scripts this might not need to change)
data_dir = fr"{root_dir}\AbsorptionCorrection"

# fill in the file paths for the vanadium and ceria runs (just run numbers might work if you are setup into the file system)
van_run = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day1\SteelDataset\DataFiles\ENGINX00361838.nxs"
ceria_run = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day1\SteelDataset\DataFiles\ENGINX00305738.nxs"

# set the path to the grouping file created by calibration
prm_path = None # r"C:\Users\kcd17618\Engineering_Mantid\Calibration\ENGINX_305738_Custom_block.prm"
grouping = "Texture30"
groupingfile_path = None

# Define some file paths
full_instr_calib = r"C:\Users\kcd17618\Documents\dev\mantid\mantid\scripts\Engineering\calib\ENGINX_full_instrument_calibration_193749.nxs"

############### ENGINEERING DIFFRACTION INTERFACE FITTING ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################


# Next the folder contraining the workspaces you want to fit
file_folder = "Focus"

# You also need to specify a name for the folder the fit parameters will be saved in
fit_save_folder = "ScriptFitParameters"

# Finally, provide a list of peaks that you want to be fit within the spectra
peaks = [2.80,2.59,2.46,1.89,1.62,1.46] # zr

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


############### ENGINEERING DIFFRACTION INTERFACE POLE FIGURE ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# define the columns you would like to create pole figures for
readout_columns = ["I", "X0"]
# and the type of projection to plot
projection_method = "Azimuthal"

# you need to define the orientation of the intrinsic sample directions when the sample orientation matrix == I (no rotation)
# this should be the same as the reference state used in the absorption correction
r2o2 = np.sqrt(2)/2
dir1 = np.array((r2o2,r2o2,0))
dir2 = np.array((0,0,1)) # projection axis
dir3 = np.array((r2o2,-r2o2,0))
# you can also supply names for these three directions
dir_names = ["RD", "AD", "HD"]

# set whether you would like the plotted pole figure to be a scatter of experimental points or whether you would like to apply gaussian smoothing and
# plot a contour representation
scatter = "both"
# if contour, what should the kernel size of the gaussian be
kernel = 6.0

# do you want to include a scattering power correction
include_scatt_power = False
# if so what is the crystal structure, defined either by giving a cif file or supplying the lattice, space group and basis
xtal_input = "cif" # "cif"/"array"/"string"
xtal_args = [r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\ZrRing\Zr.cif"] # for input "cif", require the cif filepath, for "array" array of lattice parameters, space group, basis
# for "string" lattice parameter string, space group and basis
# if you have set a crystal, you can also provide a set of hkls, the hkl_peaks dictionary is a useful way of assigning the peaks
hkl_peaks = dict(zip(peaks, [(1,0,0),(0,0,2),(1,0,1),(1,0,2),(2,-1,0),(1,0,3)]))

chi2_thresh = 0.0   # max value of Chi^2 to be included as a point in the table
peak_thresh = 0.01   # max difference from either the HKL specified or the mean X0
scat_vol_pos = (0.0,0.0,0.0) # for now, can assume the gauge vol will be centred on origin

######################### RUN SCRIPT ########################################

#~~~~~~~~~~~~~~~~~~~~~ ABSORPTION CORRECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load the ref workspace
ref_ws_str = os.path.splitext(os.path.basename(ref_ws_path))[0]
Load(Filename = ref_ws_path, OutputWorkspace = ref_ws_str)

# load data workspaces
corr_wss = find_all_files(corr_dir)
raw_wss = [os.path.splitext(os.path.basename(fp))[0] for fp in corr_wss]
for iws, ws in enumerate(raw_wss):
    if not ADS.doesExist(ws):
        Load(Filename = corr_wss[iws], OutputWorkspace= ws)

# run script
run_abs_corr(wss = raw_wss, 
             ref_ws = ref_ws_str, 
             orientation_file = orientation_file, 
             orient_file_is_euler = orient_file_is_euler, 
             euler_scheme = euler_scheme, 
             euler_axes_sense = euler_axes_sense, 
             copy_ref = copy_ref, 
             include_abs_corr = include_abs_corr, 
             monte_carlo_args = monte_carlo_args, 
             gauge_vol_preset = gauge_vol_preset, 
             gauge_vol_shape_file = gauge_vol_shape_file,
             include_atten_table = include_atten_table, 
             eval_point = eval_point, 
             eval_units = eval_units,
             exp_name = exp_name, 
             root_dir = root_dir, 
             include_div_corr = include_div_corr,
             div_hoz = div_hoz, 
             div_vert = div_vert, 
             det_hoz = det_hoz, 
             clear_ads_after = clear_ads_after)

#~~~~~~~~~~~~~~~~~~~~~~~~~~ FOCUSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
to_be_focussed_files = find_all_files(data_dir)

run_focus_script(wss = to_be_focussed_files, 
                 focus_dir = root_dir, 
                 van_run = van_run, 
                 ceria_run = ceria_run,
                 full_instr_calib = full_instr_calib, 
                 grouping = grouping,
                 prm_path = prm_path)        


#~~~~~~~~~~~~~~~~~~~~~ PEAK FITTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
# create output directory
fit_save_dir = os.path.join(root_dir, fit_save_folder)
mk(fit_save_dir)

# find and load peaks

# get grouping directory name
calib_info = CalibrationInfo(group = GROUP(grouping))
if groupingfile_path:
    calib_info.set_grouping_file(groupingfile_path)
elif prm_path:
    calib_info.set_prm_filepath(prm_path) 
group_folder = calib_info.get_group_suffix()
focussed_data_dir = os.path.join(root_dir, file_folder, group_folder, "CombinedFiles")
focus_ws_paths = find_all_files(focussed_data_dir)
focus_wss = [os.path.splitext(os.path.basename(fp))[0] for fp in focus_ws_paths]
for iws, ws in enumerate(focus_wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_ws_paths[iws], OutputWorkspace= ws)

# execute the fitting                     
fit_all_peaks(focus_wss, peaks, 0.02, fit_save_dir, i_over_sigma_thresh = i_over_sigma_thresh, nan_replacement = nan_replacement, no_fit_value_dict = no_fit_value_dict)     

#~~~~~~~~~~~~~~~~~~~~~ POLE FIGURE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_load_dirs = [os.path.join(root_dir, fit_save_folder, group_folder, str(peak)) for peak in peaks]

hkls = [hkl_peaks[peak] for peak in peaks]

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
               xtal_input = xtal_input,
               xtal_args = xtal_args,
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
