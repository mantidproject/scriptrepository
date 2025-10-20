# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from os import path, makedirs, scandir
from Engineering.texture.TextureUtils import find_all_files, run_abs_corr

############### ENGINEERING DIFFRACTION INTERFACE ABSORPTION CORRECTION ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################

# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PostExp-ZrRingDiagScript"

# otherwise set root directory here:
root_dir = fr"C:\Users\kcd17618\Engineering_Mantid\User\{exp_name}"

# next, specify the folder with the files you would like to apply the absorption correction to 
corr_dir = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\ZrRing\DataFiles\Point2"

# For texture, it is expected that you have a single ssmple shape, that is reorientated between runs.
# this is handled by having a reference workspace with the shape in its neutral position 
# (position in the beamline when the goniometer is home)
# This reference workspace probably requires you to do some interacting and validating, so should be setup in the UI
# (Interfaces/Diffraction/Engineering Diffraction/Absorption Correction)

# if this is the case copy ref should be True and the ref_ws_path should be given
# otherwise, if set ref is true, it is assumed that the sample shapes are already present on the workspaces
copy_ref = True
ref_ws_path = path.join(root_dir, "ReferenceWorkspaces", f"{exp_name}_reference_workspace.nxs") 

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
euler_scheme = "YXY"
euler_axes_sense = "1,-1,1"

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
div_hoz = 0.02
div_vert = 0.02
det_hoz = 0.02

######################### RUN SCRIPT ########################################

# load the ref workspace
ref_ws_str = path.splitext(path.basename(ref_ws_path))[0]
Load(Filename = ref_ws_path, OutputWorkspace = ref_ws_str)

# load data workspaces
corr_wss = find_all_files(corr_dir)
wss = [path.splitext(path.basename(fp))[0] for fp in corr_wss]
for iws, ws in enumerate(wss):
    if not ADS.doesExist(ws):
        Load(Filename = corr_wss[iws], OutputWorkspace= ws)

# run script
run_abs_corr(wss = wss, 
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
        
        
            
        
        
