# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.api import AnalysisDataService as ADS
from os import path, makedirs, scandir
from TexturePoleFigurePresenter import create_pf
from Engineering.texture.polefigure.polefigure_model import TextureProjection

from importlib import reload
import sys
reload(sys.modules["TexturePoleFigurePresenter"])

def find_all_files(directory):
    files = []
    with scandir(directory) as entries:
        for entry in entries:
            if entry.is_file():
                files.append(entry.path)
    return files


exp_name = "PostExp-ZrRingDiag"
save_root = r"C:\Users\kcd17618\Engineering_Mantid"
root_dir = fr"{save_root}\User\{exp_name}"
grouping = "Texture30"
ws_folder = "Focus"
fit_folder = "ScriptFitParameters"
peak = 2.575
readout_column = "I"
grouping = "Texture30"
projection_method = "Azimuthal"


r2 = np.sqrt(2)/2
dir1 = np.array((0,0,1))
dir2 = np.array((r2,r2,0)) # projection axis
dir3 = np.array((r2,-r2,0))
dir_names = ["AD", "HD", "RD"]

scatter = True
kernel = 6.0

include_scatt_power = False
cif = None
lattice = None #"2.8665  2.8665  2.8665"
space_group = None #"I m -3 m"
basis = None # "Fe 0 0 0 1.0 0.05; Fe 0.5 0.5 0.5 1.0 0.05"
# the value you set for peak when fit was done must be in this dict and correspond to the correct hkl
hkl_peaks = {1.17: (1,1,2),1.43: (2,0,0),2.03: (1,1,0)} #Fe
hkl = hkl_peaks[peak] if include_scatt_power else None

chi2_thresh = 0.0   # max value of Chi^2 to be included as a point in the table
peak_thresh = 0.0   # max difference from either the HKL specified or the mean X0
scat_vol_pos = (0.0,0.0,0.0) # for now, can assume the gauge vol will be centred on origin

# get ws paths
focus_dir = path.join(root_dir, ws_folder, grouping, "CombinedFiles")
focus_wss = find_all_files(focus_dir)
wss = [path.splitext(path.basename(fp))[0] for fp in focus_wss]

# get fit params
fit_dir = path.join(root_dir, fit_folder, grouping, str(peak))
fit_wss = find_all_files(fit_dir)
params = [path.splitext(path.basename(fp))[0] for fp in fit_wss]

for iws, ws in enumerate(wss):
    if not ADS.doesExist(ws):
        Load(Filename = focus_wss[iws], OutputWorkspace= ws)
    if not ADS.doesExist(params[iws]):
        Load(Filename = fit_wss[iws], OutputWorkspace= params[iws])


create_pf(wss, params, include_scatt_power, cif, lattice, space_group, basis,
                    hkl, readout_column, dir1, dir2, dir3, dir_names, scatter, kernel, scat_vol_pos, 
                    chi2_thresh, peak_thresh, save_root, exp_name, grouping, projection_method)
