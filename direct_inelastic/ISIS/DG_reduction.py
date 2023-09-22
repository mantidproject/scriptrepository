# =======================================================
#
#      EXCITATIONS INSTRUMENTS REDUCTION SCRIPT
#      JRS & MDL 3/7/23
#
#      Reads sample run(s), corrects for backgound,
#      normalises to a white vanadium run, and  the
#      absolute normalisation factor for each Ei (if MV
#      file is specified).
#      - Masks bad detectors with hard mask only.
#      - Converts time-of-flight to energy transfer
#      - performs Q-rebinning for QENS data ('_red.nxs' format)
#      - Outputs .nxspe, .nxs or _red.nxs files
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
import os, sys
import time
from importlib import reload

#=======================User Inputs=======================
powder         = True                        # powder or 1to1 map
sumruns        = False                       # set to True to sum sample runs
sample         = [93338]                     # sample runs (list)
sample_bg      = 93329                       # single background run
wv_file        = 'WV_91329.txt'              # white vanadium integral file (mandatory)
Ei_list        = [3.71,1.035,1.775]          # incident energies (Ei) - from PyChop
Erange         = [-0.8,0.0025,0.8]           # energy transfer range to output in fractions of Ei
trans          = [0.95,0.95,0.95]            # elastic line transmission factors for each Ei
mask           = 'MASK_FILE_XML'             # hard mask
# what to do if see run with same angle as previous run (and sumruns==False)
#  - 'replace': sum and replace first nxspe file created
#  - 'ignore': ignore and create an nxspe for each run
#  - 'accumulate': sum all runs with this angle when creating new nxspe
same_angle_action = 'ignore'
#========================================================

#==================Absolute Units Inputs=================
mv_file       = None        # pre-processed MV calibration
sample_mass   = 1           # mass of the sample
sample_fwt    = 50.9415     # formula weight of sample
monovan_mass  = None        # Mass of vanadium sample (ask local contact)
#========================================================

#==================Local Contact Inputs==================
inst = 'INSTRUMENT_NAME'
cycle = 'CYCLE_ID'              # cycle number
fixei = True                    # True for LET since no monitor 3
powdermap = 'RINGS_MAP_XML'     # rings mapping file - must be .xml format
file_wait = 30                  # wait for data file to appear (seconds)
keepworkspaces = True           # should be false for Horace scans
saveformat = '.nxspe'           # format of output, '.nxspe', '.nxs'
QENS = False                    # output Q-binned data for QENS data analysis '_red.nxs'
Qbins = 20                      # approximate number of Q-bins for QENS
theta_range = [5., 65.]         # useful theta range for Q-binning (QENS)
idebug = False                  # keep workspaces and check absolute units on elastic line
save_dir = f'/instrument/{inst}/RBNumber/USER_RB_FOLDER'  # Set to None to avoid reseting
psi_motor_name = 'rot'          # name of the rotation motor in the logs
angles_workspace = 'angles_ws'  # name of workspace to store previously seen angles
#========================================================

# ==============================setup directroties================================
config['default.instrument'] = inst
if save_dir is not None:
    config['defaultsave.directory'] = save_dir
if inst == 'MARI':
    source = 'Moderator'
    m2spec = 2                  # specID of monitor2 (pre-sample)
    m3spec = 3                  # specID of monitor3 (post-sample)
    #monovan_mass = 32.58       # mass of vanadium cylinder
    same_angle_action = 'ignore'
elif inst == 'MERLIN':
    source = 'undulator'
    m2spec = 69636              # specID of monitor2 (pre-sample)
    m3spec = 69640              # specID of monitor3 (post-sample)
    #monovan_mass = 32.62       # mass of vanadium cylinder
elif inst == 'MAPS':
    source = 'undulator'
    m2spec = 36867              # specID of monitor2 (pre-sample)
    m3spec = 36868              # specID of monitor3 (post-sample)
    #monovan_mass = 30.1        # mass of vanadium cylinder
elif inst == 'LET':
    source = 'undulator'
    m2spec = 98310              # specID of monitor2 (pre-sample)
    m3spec = None               # specID of monitor3 (post-sample)
    #monovan_mass = 3.56        # mass of vanadium cylinder
else:
    raise RuntimeError(f'Unrecognised instrument: {inst}')

cycle_shortform = cycle[2:] if cycle.startswith('20') else cycle
datadir = f'/archive/NDX{inst}/Instrument/data/cycle_{cycle_shortform}/'
instfiles  = '/usr/local/mprogs/InstrumentFiles/'
mapdir  = f'{instfiles}/{inst.swapcase()}/'
config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)
ws_list = ADS.getObjectNames()

# Try to load subsidiary utilities
sys.path.append(instfiles)
try:
    from reduction_utils import *
    utils_loaded = True
except ImportError:
    utils_loaded = False
    def rename_existing_ws(*a, **k): raise RuntimeError('Unable to load utils')
    def remove_extra_spectra_if_mari(*a, **k): raise RuntimeError('Unable to load utils')

#========================================================
# Helper functions
def tryload(irun):              # loops till data file appears
    if isinstance(irun, str) and irun in mtd:
        # If users give a string which matches an existing workspace, use that instead
        return rename_existing_ws(irun)
    while True:
        try:
            ws = Load(str(irun),LoadMonitors=True)
        except TypeError:
            ws = Load(str(irun),LoadMonitors='Separate')
        except ValueError:
            print(f'...waiting for run #{irun}')
            Pause(file_wait)
            continue
        break
    if inst == 'MARI':
        remove_extra_spectra_if_mari('ws')
    return mtd['ws']

def load_sum(run_list):
    for ii, irun in enumerate(run_list):
        tryload(irun)
        if ii == 0:
            w_buf = CloneWorkspace('ws')
            w_buf_monitors = CloneWorkspace('ws_monitors')
            print(f'run #{irun} loaded')
        else:
            w_buf = Plus('w_buf', 'ws')
            w_buf_monitors = Plus('w_buf_monitors', 'ws_monitors')
            print(f'run #{irun} added')
#========================================================


print(f'\n======= {inst} data reduction =======')
print(f'Working directory... {ConfigService.Instance().getString("defaultsave.directory")}\n')

# ============================create lists if necessary==========================
if isinstance(sample, str) or not hasattr(sample, '__iter__'):
    sample = [sample]
if sample_bg is not None and (isinstance(sample_bg, str) or not hasattr(sample_bg, '__iter__')):
    sample_bg = [sample_bg]

#==================================load hard mask================================
if mask is None:
    print(f'{inst}: WARNING - No hard mask!  Bad detectors wont be removed')
if mask not in ws_list and mask is not None:
    print(f'{inst}: Loading hard mask - {mask}')
    LoadMask(inst,mask,OutputWorkspace=mask)
else:
    print(f'{inst}: Using previously loaded hard mask - {mask}')

# ===============================load whitevan file==============================
if wv_file is None:
    raise RuntimeError(f'{inst}: ERROR - white vanadium calibration file missing')
if wv_file not in ws_list:
    print(f'{inst}: Loading white vanadium - {wv_file}')
    LoadAscii(Filename=wv_file, OutputWorkspace=wv_file)
    ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN', OutputWorkspace=wv_file)
else:
    print(f'{inst}: Using previously loaded white vanadium - {wv_file}')

# ===============================load monovan file===============================
mv_fac = []
if mv_file is not None and monovan_mass is not None:
    if mv_file not in ws_list:
        print(f'{inst}: Loading monovan calibration factors - {mv_file}')
        LoadAscii(mv_file,OutputWorkspace=mv_file)
    mv_ws = mtd[mv_file]
    mv_eis = mv_ws.readX(0)
    mv_cal = mv_ws.readY(0)
# check that monovan is compatible with Ei_list)
    Ei_diff = sum([x-y for x,y in zip(mv_eis,Ei_list)])
    if (abs(Ei_diff) > 0):
        print('----ERROR: Monovan file Eis not compatible with Ei_list')
    for Ei in Ei_list:
        ii = np.where(mv_eis == Ei)
        mvf = mv_cal[ii][0]
        van_fwt = 50.9415               # vanadium molar mass (g)
        van_xs = 5080. / 4. / np.pi     # vanadium differential cross-section (mbarns/st)
        mvf = (monovan_mass / van_fwt) / (sample_mass / sample_fwt) * van_xs / mvf
        mv_fac.append(mvf)
else:
    print(f'{inst}: Skipping absolute calibration')
    mv_fac = [x/x for x in Ei_list]     # monovan factors = 1 by default

# =======================load background runs and sum=========================
if sample_bg is not None:
    load_sum(sample_bg)
    ws_bg = CloneWorkspace('w_buf')
    ws_bg = NormaliseByCurrent('ws_bg')

# =======================sum sample runs if required=========================
if sumruns:
    load_sum(sample)
    ws = CloneWorkspace('w_buf')
    ws_monitors = CloneWorkspace('w_buf_monitors')
    sample = [sample[0]]

# =====================angles cache stuff====================================
if utils_loaded:
    build_angles_ws(sample, angles_workspace, psi_motor_name)
    runs_with_angles_already_seen = []
else:
    same_angle_action = 'ignore'

# =======================sample run loop=====================================
for irun in sample:
    t = time.time()         # start the clock

    # check loaded workspaces and remove old nxspe workspaces
    if not keepworkspaces:
        ws_list = ADS.getObjectNames()
        nx_list = [ss for ss in ws_list if saveformat in ss]
        for ss in nx_list:
            ADS.remove(ss)

    # if this run is at an angle already seen and we are doing 'replace', skip
    if same_angle_action.lower() != 'ignore' and irun in runs_with_angles_already_seen:
        continue

    print('============')
    if not sumruns:
        # Checks rotation angle
        if same_angle_action.lower() != 'ignore':
            runs_with_same_angles = get_angle(irun, angles_workspace, psi_motor_name, tryload)
            if len(runs_with_same_angles) > 1:
                load_sum(runs_with_same_angles)
                if same_angle_action.lower() == 'replace':
                    irun = runs_with_same_angles[0]
                    runs_with_angles_already_seen += runs_with_same_angles
        else:
            tryload(irun)
            print(f'Loading run# {irun}')
    ws = NormaliseByCurrent('ws')

# ============================= Ei loop =====================================
    for ienergy in range(len(Ei_list)):
        Ei  = Ei_list[ienergy]
        origEi = Ei
        tr  = trans[ienergy]
        mvf = mv_fac[ienergy]
        print(f'\n{inst}: Reducing data for Ei={Ei:.2f} meV')

        ws_corrected = Scale('ws', 1 / tr, 'Multiply')
        if sample_bg is not None:
            print(f'... subtracting background - transmission factor = {tr:.2f}')
            ws_corrected  = ws/tr - ws_bg

        # normalise to WB vanadium and apply fixed mask
        print('... normalising/masking data')
        ws_norm = Divide('ws_corrected', wv_file)       # white beam normalisation
        MaskDetectors(ws_norm,MaskedWorkspace=mask,ForceInstrumentMasking=True)

        # t2e section
        print('... t2e section')
        ws_monitors = mtd['ws_monitors']
        spectra = ws_monitors.getSpectrumNumbers()
        index = spectra.index(m2spec)
        m2pos = ws.detectorInfo().position(index)[2]

        if inst == 'MARI' and utils_loaded and origEi < 4.01:
            # Shifts data / monitors into second frame for MARI
            ws_norm, ws_monitors = shift_frame_for_mari_lowE(origEi, wsname='ws_norm', wsmon='ws_monitors')

        # this section shifts the time-of-flight such that the monitor2 peak
        # in the current monitor workspace (post monochromator) is at t=0 and L=0
        # note that the offest is not the predicted value due to energy dependence of the source position

        if m3spec is not None and not fixei:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor1Spec=m2spec,Monitor2Spec=m3spec,EnergyEstimate=Ei)
            print(f'... refined Ei={Ei:.2f} meV')
        else:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=fixei)

        print(f'... m2 tof={mon2_peak:.2f} mus, m2 pos={m2pos:.2f} m')

        ws_norm = ScaleX(ws_norm, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
        MoveInstrumentComponent(ws_norm, ComponentName=source, Z=m2pos, RelativePosition=False)

        ws_out = ConvertUnits(ws_norm, 'DeltaE', EMode='Direct', EFixed=Ei)
        ws_out = Rebin(ws_out, [x*Ei for x in Erange], PreserveEvents=False)
        ws_out = DetectorEfficiencyCor(ws_out, IncidentEnergy=Ei)
        ws_out = CorrectKiKf(ws_out, Efixed=Ei, EMode='Direct')

        # monovan scaling
        if mv_file is not None:
            print(f'... applying mono van calibration factor {mvf:.1f}')
        ws_out = Scale('ws_out', mvf, 'Multiply')

        # rings grouping if desired
        ofile_suffix='_1to1'
        if powder or inst == 'MARI' or QENS:
            ws_out=GroupDetectors(ws_out, MapFile=powdermap, Behaviour='Average')
            ofile_suffix = '_powder'
            if inst == 'MARI' or QENS:
                ofile_suffix = ''
            print(f'... powder grouping using {powdermap}')
        if sample_bg is not None and inst == 'MARI':
            ofile_suffix += '_sub'

        # output nxspe file
        ofile = f'{inst[:3]}{irun}_{origEi:g}meV{ofile_suffix}'

        # check elastic line (debug mode)
        if idebug:
            Rebin('ws_out', [-0.05*Ei, 100, 0.05*Ei], PreserveEvents=False, OutputWorkspace=ofile+'_elastic')
            Transpose(f'{ofile}_elastic', OutputWorkspace=f'{ofile}_elastic')

        # output appropriate formats
        ws_out = ConvertToDistribution(ws_out)
        if QENS:
            saveformat = '.nxs'
        print(f'{inst}: Writing {ofile}{saveformat}')
        if saveformat.lower() == '.nxspe':
            SaveNXSPE('ws_out', ofile+saveformat, Efixed=Ei, KiOverKfScaling=True)
            #if utils_loaded and not powder:
            #    copy_inst_info(ofile+saveformat, 'ws_out')   # Copies instrument info for Horace
        elif saveformat.lower() == '.nxs':
            rmlogs = {'events_log', 'frame_log', 'good_frame_log', 'period_log', 'proton_charge', 'raw_events_log'}
            RemoveLogs('ws_out', KeepLogs=','.join(set(mtd['ws_out'].run().keys()).difference(rmlogs)))
            SaveNexus('ws_out', ofile+saveformat)
        if QENS:
            print('... outputting QENS "_red" format')
            theta = np.array([theta_range[0], theta_range[1]])*np.pi/180.
            Q     = 1.39 * np.sqrt(Ei) * np.sin(theta)
            Q     = np.around(Q*10) / 10.
            Qbin  = int((Q[1] - Q[0])) / Qbins
            print(f'... Q-axis = [{Q[0]+Qbin:g},{Qbin:g},{Q[1]-Qbin:g}]')
            ws_out = SofQW3('ws_out', [Q[0]+Qbin, Qbin, Q[1]-Qbin], 'Direct', Efixed=Ei)
            # these lines are Anthony Lim's method of removing NaNs
            # NaN are changed to zeros, and then placed at the end of the energy range
            spectra = range(ws_out.getNumberHistograms())
            for ispec in spectra:
                x = ws_out.dataX(ispec)
                y = ws_out.dataY(ispec)
                e = ws_out.dataE(ispec)
                smidge = 1e-6
                for iq in range(len(x)-1):
                    if np.isnan(y[iq]):
                        x[iq] = np.max(x) + smidge
                        y[iq] = 0.0
                        e[iq] = 0.0
                        smidge += 1e-6
            SaveNexus('ws_out', ofile+"_red"+saveformat)

        CloneWorkspace('ws_out',OutputWorkspace=ofile+saveformat)

# ============================= End of Ei loop ================================

    print(f'\n{inst}: Reduction complete in {time.time() - t:.1f} seconds\n')

# ============================= End of run loop ================================

# cleanup
if not idebug:
    ws_list = ADS.getObjectNames()
    nx_list = [ss for ss in ws_list if 'w_buf' in ss or 'ws' in ss]
    for ss in nx_list:
        ADS.remove(ss)
