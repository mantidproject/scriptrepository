# =======================================================
#
#      WHITE VAN REDUCTION SCRIPT
#      JRS 26/4/23
#
#      Reads white van run, subtracts empty, and
#      writes out white vanadium integrals to an
#      ascii file for use with DG_reduction.py
#
#      No mask is employed at this stage
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
from mantid.kernel.funcinspect import lhs_info
import numpy as np
import time

t = time.time()         #start the clock

#=======================User Inputs======================
whitevan       = 78715                  # white vanadium run
whitevan_bg    = None                   # background for the white vanadium
whitevan_trans = 1                      # transmission factor
#========================================================

#==================Local Contact Inputs==================
inst = 'INSTRUMENT_NAME'                # instrument
cycle = 'CYCLE_ID'                      # cycle
wv_lrange = [1,5]                       # wavelength integration limits for output of vanadium integrals
wv_detrange = [30000,60000]             # spectrum index range for average intensity calculation
idebug = False                          # keep itermediate workspaces for debugging
save_dir = f'/instrument/{inst}/RBNumber/USER_RB_FOLDER' # Set to None to avoid resetting
#========================================================

config['default.instrument'] = inst
if save_dir is not None:
    config['defaultsave.directory'] = save_dir
cycle_shortform = cycle[2:] if cycle.startswith('20') else cycle
data_dir = f'/archive/NDX{inst}/Instrument/data/cycle_{cycle_shortform}/'
config.appendDataSearchDir(data_dir)
if save_dir is not None:
    config.appendDataSearchDir(save_dir)

#========================================================
def try_load_no_mon(irun, ws_name=None):
    # LoadISISNexus needs LoadMonitors='Exclude'
    # LoadEventNexus needs LoadMonitors=False
    if ws_name is None:
        nret, ws_name = lhs_info('both')
        if nret > 1:
            raise RuntimeError('load: Only one output workspace name can be given')
        ws_name = ws_name[0]
    try:
        return Load(str(irun), LoadMonitors='Exclude', OutputWorkspace=ws_name)
    except ValueError:
        return Load(str(irun), LoadMonitors=False, OutputWorkspace=ws_name)

def load_sum(run_list):
    if isinstance(run_list, str) or not hasattr(run_list, '__iter__'):
        run_list = [run_list]
    for ii, irun in enumerate(run_list):
        try_load_no_mon(irun, 'ws')
        if ii == 0:
            w_buf = CloneWorkspace('ws')
            print(f'run #{irun} loaded')
        else:
            w_buf = Plus('w_buf', 'ws')
            print(f'run #{irun} added')
    nret, ws_name = lhs_info('both')
    if nret == 1:
        CloneWorkspace('w_buf', OutputWorkspace=ws_name[0])
#========================================================

print(f'\n======= {inst} white van reduction =======')

# =================load white van and background and subtract=====================
print(f'WHITE_VAN {inst}: Loading white vanadium run(s)')
wv = load_sum(whitevan)
wv = NormaliseByCurrent('wv')
wv_corrected = Scale('wv', 1/whitevan_trans, 'Multiply')
if (whitevan_bg is not None):
    print(f'... subtracting white vanadium background run(s)')
    wv_bg = load_sum(whitevan_bg)
    wv_bg = NormaliseByCurrent('wv_bg')
    wv_corrected = wv/whitevan_trans - wv_bg

# ==============create wv_norm workspace (white vanadium integrals)===============
print('... Calculating white vanadium integrals')
print(f'... Summing between {wv_lrange[0]:.1f} < \u03BB < {wv_lrange[1]:.1f} \u212B')
WV_normalised_integrals = ConvertUnits(wv_corrected, 'Wavelength')
WV_normalised_integrals = Rebin(WV_normalised_integrals, f'{wv_lrange[0]},100,{wv_lrange[1]}', PreserveEvents=False)
wv_normt = Transpose(WV_normalised_integrals)
if wv_detrange is not None:
    wv_normt = CropWorkspace(wv_normt, XMin=wv_detrange[0], Xmax=wv_detrange[1])
wv_scale = Integration(wv_normt)
scale_factor = len(wv_normt.dataY(0)) / wv_scale.dataY(0)[0]
WV_normalised_integrals = Scale(WV_normalised_integrals, scale_factor, 'Multiply')

# ===========================output integrals file================================
wv_name = whitevan if (isinstance(whitevan, str) or len(whitevan)==1) else whitevan[0]
ofile = f'WV_{wv_name}.txt'
print(f'WHITE_VAN {inst}: Saving integrals in {ofile}')
print(f'... average intensity = {1/scale_factor:.2f}')
SaveAscii(WV_normalised_integrals, ofile)

# ===================================cleanup======================================
ws_list = ADS.getObjectNames()
if (not idebug):
    nx_list = [ss for ss in ws_list if 'wv' in ss]
    for ss in nx_list:
        ADS.remove(ss)

# ==================================================================
print(f'\nWHITE_VAN {inst}: Reduction complete in {time.time() - t:.1f} seconds')

