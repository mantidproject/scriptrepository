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
from mantid.kernel.funcinspect import lhs_info
import numpy as np
import os, sys
import time
from importlib import reload

#!begin_params
#=======================User Inputs=======================
powder         = True                        # powder or 1to1 map
sumruns        = False                       # set to True to sum sample runs
sample         = [93338]                     # sample runs (list)
sample_bg      = 93329                       # single background run
sample_cd      = None                        # Cadmium run for instrument background subtraction
wv_file        = 'WV_91329.txt'              # white vanadium integral file (mandatory)
Ei_list        = ['auto']                    # incident energies (Ei) - from PyChop
Erange         = [-0.8,0.0025,0.8]           # energy transfer range to output in fractions of Ei
trans          = [0.95,0.95,0.95]            # elastic line transmission factors for each Ei
mask           = 'MASK_FILE_XML'             # hard mask
# what to do if see run with same angle as previous run (and sumruns==False)
#  - 'replace': sum and replace first nxspe file created
#  - 'ignore': ignore and create an nxspe for each run
same_angle_action = 'ignore'
# to enable creating multiple reduced data files from a single
# "continuous scan" set the following variables to a valid log "block"
# name and the bin size (width) to something larger than zero.
# An optional unit can be given.
# If sumruns=True, the sample list is summed and then divided
# If sumruns=False and multiple sample runs are listed, an error is raised.
# This is also only compatible with same_angle_action='ignore'
cs_block = None
cs_block_unit = ''
cs_bin_size = 0
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
keepworkspaces = False          # should be false for Horace scans
saveformat = '.nxspe'           # format of output, '.nxspe', '.nxs'
QENS = False                    # output Q-binned data for QENS data analysis '_red.nxs'
Qbins = 20                      # approximate number of Q-bins for QENS
theta_range = [5., 65.]         # useful theta range for Q-binning (QENS)
idebug = False                  # keep workspaces and check absolute units on elastic line
save_dir = f'/data/analysis/{inst}/RBNumber/USER_RB_FOLDER'  # Set to None to avoid reseting
psi_motor_name = 'rot'          # name of the rotation motor in the logs
angles_workspace = 'angles_ws'  # name of workspace to store previously seen angles
sumruns_savemem = False         # Compresses event in summed ws to save memory
                                # (causes some loss of data so cannot use event filtering)
cs_smidge = 0.001               # Tolerance on continuous scan bins size
cs_conv_to_md = False           # Convert a continuous scan to an MDWorkspace
cs_conv_pars = {}               # A dict with 'lattice_pars', 'lattice_angles', 'u', 'v', 'psi0'
cs_single = False               # Outputs an EventWorkspace of a continuous scan (does not histogram)
#========================================================
#!end_params

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
    m2spec = 69634              # specID of monitor2 (pre-sample)
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
datadir2 = f'/data/instrument/{inst}/CYCLE20{cycle_shortform.replace("_","")}/USER_RB_FOLDER/'
instfiles  = '/usr/local/mprogs/InstrumentFiles/'
mapdir  = f'{instfiles}/{inst.swapcase()}/'
config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)
config.appendDataSearchDir(datadir2)
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
    def autoei(*a, **k):
        raise RuntimeError('Cannot use Auto-Ei: ' \
            'You need to copy the reduction_utils.py file somewhere on the Python path')

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
            try:
                ws = Load(str(irun),LoadMonitors='Separate')
            except ValueError: # Possibly when file is being copied over from NDX machine
                Pause(2)
                ws = Load(str(irun),LoadMonitors=True)
        except ValueError:
            print(f'...waiting for run #{irun}')
            Pause(file_wait)
            continue
        break
    if inst == 'MARI':
        remove_extra_spectra_if_mari('ws')
    return mtd['ws']

def load_sum(run_list, block_name=None):
    for ii, irun in enumerate(run_list):
        tryload(irun)
        if ii == 0:
            w_buf = CloneWorkspace('ws')
            w_buf_monitors = CloneWorkspace('ws_monitors')
            if block_name:
                bval = mtd['ws'].getRun().getLogData(block_name).filtered_value
            print(f'run #{irun} loaded')
        else:
            w_buf = Plus('w_buf', 'ws')
            w_buf_monitors = Plus('w_buf_monitors', 'ws_monitors')
            if block_name:
                bv2 = mtd['ws'].getRun().getLogData(block_name).filtered_value
                bval = np.concatenate((bval, bv2))
            print(f'run #{irun} added')
    ADS.remove('ws')
    ADS.remove('ws_monitors')
    wsout_name = lhs_info('names')[0]
    if wsout_name != 'w_buf':
        RenameWorkspace('w_buf_monitors', wsout_name+'_monitors')
        RenameWorkspace('w_buf', wsout_name)
    return (mtd[wsout_name], bval) if block_name else mtd[wsout_name]

#========================================================


print(f'\n======= {inst} data reduction =======')
print(f'Working directory... {ConfigService.Instance().getString("defaultsave.directory")}\n')

# ============================create lists if necessary==========================
if isinstance(sample, str) or not hasattr(sample, '__iter__'):
    sample = [sample]
if sample_bg is not None and (isinstance(sample_bg, str) or not hasattr(sample_bg, '__iter__')):
    sample_bg = [sample_bg]
if sample_cd is not None and (isinstance(sample_cd, str) or not hasattr(sample_cd, '__iter__')):
    sample_cd = [sample_cd]

#==================================load hard mask================================
if mask is None:
    print(f'{inst}: WARNING - No hard mask!  Bad detectors wont be removed')
if mask not in ws_list and mask is not None:
    assert mask.endswith('xml') or mask.endswith('msk'), 'Mask file should be .msk or .xml'
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

# =======================load background runs and sum=========================
if sample_cd is not None:
    ws_cd = load_sum(sample_cd)
    ws_cd = NormaliseByCurrent(ws_cd)
if sample_bg is not None:
    ws_bg = load_sum(sample_bg)
    ws_bg = NormaliseByCurrent(ws_bg)
    if sample_cd is not None:
        ws_bg = ws_bg - ws_cd

# ==========================continuous scan stuff=============================
if cs_block and cs_bin_size > 0:
    if len(sample) > 1 and not sumruns:
        raise RuntimeError('Continous scans for multiple (non-summed) files not supported')
    if same_angle_action.lower() != 'ignore':
        raise RuntimeError(f'Continous scans not compatible with same_angle_action="{same_angle_action}"')
    # Sets sumruns false so don't redo the sum below. Instead sum here to get filtered block values.
    # This is due to a bug in the Mantid "Plus" algorithm, which gives incorrect filtered values
    # https://github.com/mantidproject/mantid/issues/36194
    sumruns = False
    ws_full, bval = load_sum(sample, cs_block)
    ws = ws_full # So that auto-ei works (it assumes the workspace variable is "ws")
    ws_monitors = CloneWorkspace('ws_full_monitors')
    bval_range = max(bval) - min(bval)
    bval_nbins = int(bval_range / cs_bin_size)
    bval_remainder = bval_range - cs_bin_size * bval_nbins
    # Overrides the sample list in order to hijack the loop over run numbers below
    irun_orig = sample[0]
    sample = [x*cs_bin_size + min(bval) for x in range(bval_nbins)]
    unit = cs_block_unit
    if not unit:
        unit = ws_full.getRun().getLogData(cs_block).units
    print(f'{cs_block} = {min(bval):.1f} {unit} to {max(bval):.1f} {unit}')
    print(f'... filtering in {cs_bin_size:.1f} {unit} steps')
    print(f'... N={bval_nbins} bins with {bval_remainder:.2f} {unit} remainder')
if cs_conv_to_md:
    powder = False
    assert all([v in cs_conv_pars for v in ['lattice_pars', 'lattice_ang', 'u', 'v', 'psi0']]), \
        'Conversion to MDWorkspace needs parameters: "lattice_pars", "lattice_ang", "u", "v", "psi0"'

# =======================sum sample runs if required=========================
sumsuf = sumruns and len(sample) > 1
if sumruns:
    ws = load_sum(sample)
    sample = [sample[0]]

# =============================auto-Ei stuff=================================
is_auto = lambda x: isinstance(x, str) and 'auto' in x.lower()
if is_auto(Ei_list) or hasattr(Ei_list, '__iter__') and is_auto(Ei_list[0]):
    use_auto_ei = True
    try:
        Ei_list = autoei(ws)
    except (NameError, RuntimeError) as e:
        fn = str(sample[0]).split(inst[:3])[-1].replace('.raw', '.nxs').replace('.s','.n0')
        try:
            ws_tmp_mons = LoadNexusMonitors(fn, OutputWorkspace='ws_tmp_mons')
            Ei_list = autoei(ws_tmp_mons)
        except:
            pass
    if not Ei_list and mv_file is not None:
        raise RuntimeError(f'Invalid run(s) {sample}. Chopper(s) not running ' \
                            'or could not determine neutron energy')
    print(f"Automatically determined Ei's: {Ei_list}")
    if len(trans) < len(Ei_list):
        print(f'{inst}: WARNING - not enough transmision values for auto-Ei. ' \
              'Extending list with last (end) value')
        trans += [trans[-1]]*(len(Ei_list) - len(trans))
else:
    use_auto_ei = False

# ============================load monovan file==============================
mv_fac = []
if mv_file is not None and monovan_mass is not None:
    use_auto_ei = False
    if mv_file not in ws_list:
        print(f'{inst}: Loading monovan calibration factors - {mv_file}')
        LoadAscii(mv_file,OutputWorkspace=mv_file)
    mv_ws = mtd[mv_file]
    mv_eis = mv_ws.readX(0)
    mv_cal = mv_ws.readY(0)
# check that monovan is compatible with Ei_list)
    Ei_diff = sum([x-y for x,y in zip(mv_eis,Ei_list)])
    if (abs(Ei_diff) > 0):
        raise RuntimeError('----ERROR: Monovan file Eis not compatible with Ei_list')
    for Ei in Ei_list:
        ii = np.where(mv_eis == Ei)
        mvf = mv_cal[ii][0]
        van_fwt = 50.9415               # vanadium molar mass (g)
        van_xs = 5080. / 4. / np.pi     # vanadium differential cross-section (mbarns/st)
        mvf = (monovan_mass / van_fwt) / (sample_mass / sample_fwt) * van_xs / mvf
        mv_fac.append(mvf)
else:
    print(f'{inst}: Skipping absolute calibration')
    mv_fac = [1.0 for x in Ei_list]     # monovan factors = 1 by default

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

    # For continuous scans generate the workspace based on binned log values
    # NB. we replaced the sample (run) list with list of block min bin boundaries
    if cs_block and cs_bin_size > 0:
        val = round(irun + cs_bin_size/2,1)
        print(f"Filtering: {cs_block}={val:.1f} ± {round(cs_bin_size/2,2):.2f} {unit}")
        ws = FilterByLogValue(ws_full, cs_block, irun, irun + cs_bin_size - cs_smidge)

    # if this run is at an angle already seen and we are doing 'replace', skip
    if same_angle_action.lower() != 'ignore' and irun in runs_with_angles_already_seen:
        continue

    print('============')
    if not sumruns and not cs_block:
        # Checks rotation angle
        if same_angle_action.lower() != 'ignore':
            runs_with_same_angles = get_angle(irun, angles_workspace, psi_motor_name, tryload)
            if len(runs_with_same_angles) > 1:
                ws = load_sum(runs_with_same_angles)
                irun = runs_with_same_angles[0]
                runs_with_angles_already_seen += runs_with_same_angles
        else:
            tryload(irun)
            print(f'Loading run# {irun}')

    # Fixes the current if it's been corrupted by a bug in Mantid
    try:
        if mtd['ws'].getRun().getLogData('gd_prtn_chrg').value == 0:
            AddSampleLog('ws', 'gd_prtn_chrg', str(np.sum(mtd['ws'].getRun().getLogData('proton_charge').value)), 'Number')
        ws = NormaliseByCurrent('ws')
    except RuntimeError:
        print(f'{inst}: WARNING: Could not normalise run. No integrated current found or current is zero')
        ws = mtd['ws']
    if sumruns and sumruns_savemem:
        ws = CompressEvents(ws, Tolerance=1e-5)  # Tolerance in microseconds

    # instrument geometry to work out ToF ranges
    sampos = ws.getInstrument().getSample().getPos()
    l1 = (sampos - ws.getInstrument().getSource().getPos()).norm()
    l2 = (ws.getDetector(ws.getSpectrumNumbers()[1]).getPos() - sampos).norm()

    # Updates ei_loop if auto _and_ not using monovan
    if use_auto_ei:
        Ei_list = autoei(ws)
        if len(trans) < len(Ei_list): trans += [trans[-1]]*(len(Ei_list) - len(trans))
        if len(mv_fac) < len(Ei_list): mv_fac += [1]*(len(Ei_list) - len(mv_fac))

# ============================= Ei loop =====================================
    for ienergy in range(len(Ei_list)):
        Ei  = Ei_list[ienergy]
        origEi = Ei
        tr  = trans[ienergy]
        mvf = mv_fac[ienergy]
        print(f'\n{inst}: Reducing data for Ei={Ei:.2f} meV')

        # Low energy reps on MARI end up in the second frame
        t_shift = 20000 if inst == 'MARI' and origEi < 3.1 else 0

        tofs = ws.readX(0)
        if np.max(tofs) < 0.1 and isinstance(ws, mantid.dataobjects.EventWorkspace): # Probably live data don't crop
            ws_rep = CloneWorkspace(ws)
        else:
            tof_min = np.sqrt(l1**2 * 5.227e6 / Ei) - t_shift
            tof_max = tof_min + np.sqrt(l2**2 * 5.226e6 / (Ei*(1-Erange[-1])))
            ws_rep = CropWorkspace(ws, max(min(tofs), tof_min), min(max(tofs), tof_max))

        if sample_cd is not None:
            print(f'... subtracting Cd background')
            ws_rep = ws_rep - ws_cd

        if sample_bg is not None:
            print(f'... subtracting background - transmission factor = {tr:.2f}')
            ws_rep  = ws_rep/tr - ws_bg
        else:
            ws_rep = Scale('ws_rep', 1 / tr, 'Multiply')

        # normalise to WB vanadium and apply fixed mask
        print('... normalising/masking data')
        ws_rep = Divide('ws_rep', wv_file)       # white beam normalisation
        MaskDetectors(ws_rep, MaskedWorkspace=mask, ForceInstrumentMasking=True)

        # t2e section
        print('... t2e section')
        ws_monitors = mtd['ws_monitors']
        spectra = ws_monitors.getSpectrumNumbers()
        m2pos = ws.detectorInfo().position(spectra.index(m2spec))[2]

        if inst == 'MARI' and utils_loaded and origEi < 4.01:
            # Shifts data / monitors into second frame for MARI
            ws_rep, ws_monitors = shift_frame_for_mari_lowE(origEi, wsname='ws_rep', wsmon='ws_monitors')

        # this section shifts the time-of-flight such that the monitor2 peak
        # in the current monitor workspace (post monochromator) is at t=0 and L=0
        # note that the offest is not the predicted value due to energy dependence of the source position

        if m3spec is not None:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor1Spec=m2spec,Monitor2Spec=m3spec,EnergyEstimate=Ei,FixEi=fixei)
            print(f'... refined Ei={Ei:.2f} meV')
        else:
            (Ei,mon2_peak,_,_) = GetEi(ws_monitors,Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=fixei)

        print(f'... m2 tof={mon2_peak:.2f} mus, m2 pos={m2pos:.2f} m')

        ws_rep = ScaleX(ws_rep, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
        MoveInstrumentComponent(ws_rep, ComponentName=source, Z=m2pos, RelativePosition=False)

        ws_rep = ConvertUnits(ws_rep, 'DeltaE', EMode='Direct', EFixed=Ei)
        ws_out = Rebin(ws_rep, [x*origEi for x in Erange], PreserveEvents=cs_single)
        if not cs_single:
            ws_out = DetectorEfficiencyCor(ws_out, IncidentEnergy=Ei)
        ws_out = CorrectKiKf(ws_out, Efixed=Ei, EMode='Direct')
        ADS.remove('ws_rep')

        # monovan scaling
        if mv_file is not None:
            print(f'... applying mono van calibration factor {mvf:.1f}')
        ws_out = Scale('ws_out', mvf, 'Multiply')

        # Sets output file name
        if cs_block and cs_bin_size > 0:
            ofile_prefix = f'{inst[:3]}{irun_orig}'
            ofile_suffix = f"_{val:.1f}{unit}"
        else:
            ofile_prefix = f'{inst[:3]}{irun}'
            ofile_suffix = ''
        ofile_prefix = ofile_prefix[3:] if ofile_prefix[:3] == ofile_prefix[3:6] else ofile_prefix

        # rings grouping if desired
        if powder or inst == 'MARI' or QENS:
            ws_out=GroupDetectors(ws_out, MapFile=powdermap, Behaviour='Average')
            ofile_suffix += '_powder'
            if inst == 'MARI' or QENS:
                ofile_suffix = ''
            print(f'... powder grouping using {powdermap}')
        else:
            ofile_suffix += '_1to1'
        if inst == 'MARI':
            if sumsuf:
                ofile_suffix += 'sum'
            if sample_bg is not None:
                ofile_suffix += '_sub'

        # output nxspe file
        ofile = f'{ofile_prefix}_{origEi:g}meV{ofile_suffix}'

        # check elastic line (debug mode)
        if idebug:
            Rebin('ws_out', [-0.05*Ei, 100, 0.05*Ei], PreserveEvents=False, OutputWorkspace=ofile+'_elastic')
            Transpose(f'{ofile}_elastic', OutputWorkspace=f'{ofile}_elastic')

        # Convert to MD:
        if cs_conv_to_md:
            SetGoniometer('ws_out', Axis0=f'{cs_block},0,1,0,1', Axis1=f'{cs_conv_pars["psi0"]},0,1,0,-1')
            a, b, c = tuple(cs_conv_pars['lattice_pars'])
            alpha, beta, gamma = tuple(cs_conv_pars['lattice_ang'])
            SetUB('ws_out', a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, u=cs_conv_pars['u'], v=cs_conv_pars['v'])
            mn, mx = ConvertToMDMinMaxGlobal('ws_out', QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='HKL')
            ConvertToMD('ws_out', QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='HKL', QConversionScales='HKL',
                MinValues=mn, MaxValues=mx, PreprocDetectorsWS='-', OutputWorkspace=ofile+'_ang_md')
            continue

        if cs_single:
            CloneWorkspace('ws_out', OutputWorkspace='single_md')
            continue

        # output appropriate formats
        ws_out = ConvertToDistribution(ws_out)
        if QENS:
            saveformat = '.nxs'
        print(f'{inst}: Writing {ofile}{saveformat}')
        if saveformat.lower() == '.nxspe':
            SaveNXSPE('ws_out', ofile+saveformat, Efixed=Ei, KiOverKfScaling=True)
            if utils_loaded:
                copy_inst_info(ofile+saveformat, 'ws_out')   # Copies instrument info for Horace
        elif saveformat.lower() == '.nxs':
            rmlogs = {'events_log', 'frame_log', 'good_frame_log', 'period_log', 'proton_charge', 'raw_events_log'}
            RemoveLogs('ws_out', KeepLogs=','.join(set(mtd['ws_out'].run().keys()).difference(rmlogs)))
            SaveNexusProcessed('ws_out', ofile+saveformat, PreserveEvents=False, CompressNexus=True)
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
