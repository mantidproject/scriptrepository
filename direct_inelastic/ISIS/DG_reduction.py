# =======================================================
#
#      EXCITATIONS INSTRUMENTS REDUCTION SCRIPT
#      JRS 20/6/23
#
#      Reads sample run(s), corrects for backgound,
#      normalises to a white vanadium run, and  the
#      absolute normalisation factor for each Ei (if MV
#      file is specified). s
#      - Masks bad detectors with hard mask only.
#      - Converts time-of-flight to energy transfer
#      - performs Q-rebinning for QENS data only
#      - Outputs .nxspe, .nxs or _red.nxs fiile
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import matplotlib.pyplot as plt
import numpy as np
import re, os, h5py
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
#========================================================

#==================Absolute Units Inputs=================
mv_file       = None        # pre-processed MV calibration
sample_mass   = 1           # mass of the sample
sample_fwt    = 50.9415     # formula weight of sample
#========================================================

#==================Local Contact Inputs==================
inst = 'INSTRUMENT_NAME'
cycle = 'CYCLE_ID'              # cycle number
fixei = True                    # True for LET since no monitor 3
powdermap = 'RINGS_MAP_XML'     # rings mapping file - must be .xml format
file_wait = 30                  # wait for data file to appear (seconds)
keepworkspaces = True           # should be false for Horace scans
saveformat = '.nxspe'           # format of output, '.nxspe', '.nxs'
QENS = True                     # output Q-binned data for QENS data analysis '_red.nxs'
Qbins = 20                      # approximate number of Q-bins for QENS
idebug = False                  # keep workspaces and check absolute units on elastic line
save_dir = f'/instrument/{inst}/RBNumber/USER_RB_FOLDER'  # Set to None to avoid reseting
#========================================================

#========================================================
def tryload(irun):              # loops till data file appears
    while True:
        try:
            ws = Load(str(irun),LoadMonitors=True)
        except TypeError:
            ws = Load(str(irun),LoadMonitors='Separate')
        except ValueError:
            print(f'...waiting for run #{irun}')
            time.sleep(file_wait)
            continue
        break
#========================================================

#========================================================
# Functions to copy instrument info needed by HORACE
# Resolution convolution if it exists in the raw file
def get_raw_file_from_ws(ws):
    for alg in [h for h in ws.getHistory().getAlgorithmHistories() if 'Load' in h.name()]:
        for prp in [a for a in alg.getProperties() if 'Filename' in a.name()]:
            if re.search('[0-9]*.nxs', prp.value()) is not None:
                return prp.value()
    raise RuntimeError('Could not find raw NeXus file in workspace history')

def copy_inst_info(outfile, in_ws):
    print(get_raw_file_from_ws(mtd[in_ws]))
    if not os.path.exists(outfile):
        outfile = mantid.simpleapi.config['defaultsave.directory'] + '/' + os.path.basename(outfile)
    print(outfile)
    with h5py.File(get_raw_file_from_ws(mtd[in_ws]), 'r') as raw:
        with h5py.File(outfile, 'r+') as spe:
            exclude = ['dae', 'detector_1', 'name']
            to_copy = [k for k in raw['/raw_data_1/instrument'] if not any([x in k for x in exclude])]
            print(spe.keys())
            spe_root = list(spe.keys())[0]
            en0 = spe[f'{spe_root}/instrument/fermi/energy'][()]
            if 'fermi' in to_copy:
                del spe[f'{spe_root}/instrument/fermi']
            for grp in to_copy:
                print(grp)
                src = raw[f'/raw_data_1/instrument/{grp}']
                h5py.Group.copy(src, src, spe[f'{spe_root}/instrument/'])
            if 'fermi' in to_copy:
                spe[f'{spe_root}/instrument/fermi/energy'][()] = en0
            detroot = f'{spe_root}/instrument/detector_elements_1'
            spe.create_group(detroot)
            for df0, df1 in zip(['SPEC', 'UDET', 'DELT', 'LEN2', 'CODE', 'TTHE', 'UT01'], \
                ['spectrum_number', 'detector_number', 'delt', 'distance', 'detector_code', 'polar_angle', 'azimuthal_angle']):
                src = raw[f'/raw_data_1/isis_vms_compat/{df0}']
                h5py.Group.copy(src, src, spe[detroot], df1)
            for nn in range(raw['/raw_data_1/isis_vms_compat/NUSE'][0]):
                src = raw[f'/raw_data_1/isis_vms_compat/UT{nn+1:02d}']
                h5py.Group.copy(src, src, spe[detroot], f'user_table{nn+1:02d}')
#========================================================

# ==============================setup directroties================================
config['default.instrument'] = inst
if save_dir is not None:
    config['defaultsave.directory'] = save_dir
if inst == 'MARI':
    source = 'Moderator'
    m2spec = 2                  # specID of monitor2 (pre-sample)
    m3spec = 3                  # specID of monitor3 (post-sample)
    monovan_mass = 32.58        # mass of vanadium cylinder
elif inst == 'MERLIN':
    source = 'undulator'
    m2spec = 69636              # specID of monitor2 (pre-sample)
    m3spec = 69640              # specID of monitor3 (post-sample)
    monovan_mass = 32.62        # mass of vanadium cylinder
elif inst == 'MAPS':
    source = 'undulator'
    m2spec = 41475              # specID of monitor2 (pre-sample)
    m3spec = 41476              # specID of monitor3 (post-sample)
    monovan_mass = 30.1         # mass of vanadium cylinder
elif inst == 'LET':
    source = 'undulator'
    m2spec = 98310              # specID of monitor2 (pre-sample)
    m3spec = None               # specID of monitor3 (post-sample)
    monovan_mass = 3.56         # mass of vanadium cylinder
else:
    raise RuntimeError(f'Unrecognised instrument: {inst}')

cycle_shortform = cycle[2:] if cycle.startswith('20') else cycle
datadir = f'/archive/NDX{inst}/Instrument/data/cycle_{cycle_shortform}/'
mapdir  = f'/usr/local/mprogs/InstrumentFiles/{inst.swapcase()}/'
config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)
ws_list = ADS.getObjectNames()

print(f'\n======= {inst} data reduction =======')
print(f'Working directory... {ConfigService.Instance().getString("defaultsave.directory")}\n')

# ============================create lists if necessary==========================
if type(sample) is not list:
    sample = [sample]
if type(sample_bg) is not list and sample_bg is not None:
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
if wv_file not in ws_list:
    print(f'{inst}: Loading white vanadium - {wv_file}')
    LoadAscii(Filename=wv_file, OutputWorkspace=wv_file)
    ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN', OutputWorkspace=wv_file)
else:
    print(f'{inst}: Using previously loaded white vanadium - {wv_file}')

# ===============================load monovan file===============================
mv_fac = []
if mv_file is not None:
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
    for irun in sample_bg:
        ws_bg = Load(str(irun),LoadMonitors=False)
        if (irun == sample_bg[0]):
            print(f'background run #{irun} loaded')
            w_buf = CloneWorkspace('ws_bg')
        else:
            print(f'background run #{irun} added')
            w_buf = Plus('w_buf', 'ws_bg')
    ws_bg = CloneWorkspace('w_buf')
    ws_bg = NormaliseByCurrent('ws_bg')

# =======================sum sample runs if required=========================
if sumruns:
    for irun in sample:
        tryload(irun)
        if (irun == sample[0]):
            print(f'run #{irun} loaded')
            w_buf = CloneWorkspace('ws')
            w_buf_monitors = CloneWorkspace('ws_monitors')
        else:
            print(f'run #{irun} added')
            w_buf = Plus('w_buf', 'ws')
            w_buf_monitors = Plus('w_buf_monitors', 'ws_monitors')
    ws = CloneWorkspace('w_buf')
    ws_monitors = CloneWorkspace('w_buf_monitors')
    sample = [sample[0]]

# =======================sample run loop=====================================
for irun in sample:
    t = time.time()         # start the clock

    # check loaded workspaces and remove old nxspe workspaces
    if not keepworkspaces:
        ws_list = ADS.getObjectNames()
        nx_list = [ss for ss in ws_list if saveformat in ss]
        for ss in nx_list:
            ADS.remove(ss)

    print('============')
    if not sumruns:
        tryload(irun)
        print(f'Loading run# {irun}')
    if inst == 'MARI' and mtd['ws'].getNumberHistograms() > 918:
        ws = RemoveSpectra('ws',[0])
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
            if inst == 'MARI':
                ofile_suffix = ''
            if QENS:
                ofile_suffix = '_red'
            print(f'... powder grouping using {powdermap}')

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
            copy_inst_info(ofile, 'ws_out')
        elif saveformat.lower() == '.nxs' and not QENS:
            SaveNexus('ws_out', ofile+saveformat)
        if QENS:
            print('... outputting QENS "_red" format')
            theta = np.array([7.5,65.])*np.pi/180.
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
            SaveNexus('ws_out', ofile+saveformat)

        CloneWorkspace('ws_out',OutputWorkspace=ofile+saveformat)

# ============================= End of Ei loop ================================

    print(f'\n{inst}: Reduction complete in {time.time() - t:.1f} seconds\n')

# ============================= End of run loop ================================

# cleanup
if (not idebug):
    ws_list = ADS.getObjectNames()
    nx_list = [ss for ss in ws_list if 'w_buf' in ss or 'ws' in ss]
    for ss in nx_list:
        ADS.remove(ss)
