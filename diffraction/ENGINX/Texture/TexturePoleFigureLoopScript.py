# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from Engineering.texture.TextureUtils import find_all_files, create_pf_loop

############### ENGINEERING DIFFRACTION INTERFACE POLE FIGURE ANALOGUE #######################

######################### EXPERIMENTAL INFORMATION ########################################
# First, you need to specify your file directories, If you are happy to use the same root, from experiment
# to experiment, you can just change this experiment name.
exp_name = "PostExp-ZrRingDiagScript"

# otherwise set root directory here:
save_root = r"C:\Users\kcd17618\Engineering_Mantid"
root_dir = fr"{save_root}\User\{exp_name}"


ws_folder = "Focus"
fit_folder = "ScriptFitParameters"
# define the peaks of interest, NOTE these must correspond to sub folders in the fit directory
peaks = [2.8, 2.575, 2.455, 1.89, 1.62, 1.46]
# define the columns you would like to create pole figures for
readout_columns = ["I", "I_est", "X0"]
# you need to specify the detector grouping
grouping = "Texture30"
# and the type of projection to plot
projection_method = "Azimuthal"

# you need to define the orientation of the intrinsic sample directions when the sample orientation matrix == I (no rotation)
# this should be the same as the reference state used in the absorption correction
r2 = np.sqrt(2)/2
dir1 = np.array((0,0,1))
dir2 = np.array((r2,r2,0)) # projection axis
dir3 = np.array((r2,-r2,0))
# you can also supply names for these three directions
dir_names = ["AD", "HD", "RD"]

# set whether you would like the plotted pole figure to be a scatter of experimental points or whether you would like to apply gaussian smoothing and
# plot a contour representation
scatter = True
# if contour, what should the kernel size of the gaussian be
kernel = 6.0

# do you want to include a scattering power correction
include_scatt_power = False
# if so what is the crystal structure, defined either by giving a cif file or supplying the lattice, space group and basis
cif = None
lattice = None #"2.8665  2.8665  2.8665"
space_group = None #"I m -3 m"
basis = None # "Fe 0 0 0 1.0 0.05; Fe 0.5 0.5 0.5 1.0 0.05"
# if you have set a crystal, you can also provide a set of hkls, the hkl_peaks dictionary is a useful way of assigning the peaks
hkl_peaks = {1.17: (1,1,2),1.43: (2,0,0),2.03: (1,1,0)} #Fe

chi2_thresh = 0.4   # max value of Chi^2 to be included as a point in the table
peak_thresh = 0.01   # max difference from either the HKL specified or the mean X0
scat_vol_pos = (0.0,0.0,0.0) # for now, can assume the gauge vol will be centred on origin

######################### RUN SCRIPT ########################################

hkls = [hkl_peaks[peak] if include_scatt_power else None for peak in peaks]
create_pf_loop(root_dir = root_dir, 
               ws_folder = ws_folder, 
               fit_folder = fit_folder, 
               peaks = peaks, 
               grouping = grouping, 
               include_scatt_power = include_scatt_power, 
               cif = cif, 
               lattice = lattice, 
               space_group = space_group, 
               basis = basis,
               hkls = hkls, 
               readout_columns = readout_columns, 
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
