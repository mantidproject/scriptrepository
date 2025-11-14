from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
import re, os, sys, h5py
import json
import warnings
import importlib
import types
import scipy.optimize
from os.path import abspath, dirname
from mantid.kernel.funcinspect import lhs_info

#========================================================
# General utility functions

def rename_existing_ws(ws_name):
    # Gets an existing workspace and ensures we have its monitors for reduction
    ws = CloneWorkspace(ws_name, OutputWorkspace='ws')
    try:
        mon_ws = ws.getMonitorWorkspace()
    except RuntimeError:
        specInfo = ws.spectrumInfo()
        try:
            mon_list = [i for i in range(ws.getNumberHistograms()) if specInfo.isMonitor(i)]
        except RuntimeError:
            mon_list = []
        if len(mon_list) > 0:
            ExtractMonitors(ws_name, DetectorWorkspace='ws', MonitorWorkspace='ws_monitors')
            if np.max(mtd['ws_monitors'].readX(0)) < 0.1:
                _create_dummy_monitors(ws_name)
        else:
            _get_mon_from_history(ws_name)
    else:
        CloneWorkspace(mon_ws, OutputWorkspace='ws_monitors')
    return ws

def _get_mon_from_history(ws_name):
    # Tries to look in the history of a workspace for its original raw file
    # loads the monitors from that raw file.
    orig_file = None
    for hist in mtd[ws_name].getHistory().getAlgorithmHistories():
        if hist.name().startswith('StartLiveData'):
            return _create_dummy_monitors(ws_name)
        if hist.name().startswith('Load') and 'Filename' in [pp.name() for pp in hist.getProperties()]:
            orig_file = hist.getPropertyValue('Filename')
            break
    if orig_file is None:
        raise RuntimeError(f'Cannot find original file from workspace {ws_name} to load logs from')
    try:
        LoadEventNexus(orig_file, SpectrumMax=10, LoadMonitors=True, OutputWorkspace='tmp_mons')
    except (TypeError, ValueError) as err:
        LoadRaw(orig_file, SpectrumMax=10, LoadMonitors='Separate', OutputWorkspace='tmp_mons')
    ws_mon_name = f'{ws_name}_monitors'
    RenameWorkspace('tmp_mons_monitors', ws_mon_name)
    DeleteWorkspace('tmp_mons')
    mtd[ws_name].setMonitorWorkspace(mtd[ws_mon_name])
    CloneWorkspace(ws_mon_name, OutputWorkspace='ws_monitors')

MONDAT = {
    'MERLIN': {'ws':range(69632, 69641), 'l2':[3.258]+[1.504]*4+[4.247]*4, 'th':[180]*5+[0]*4},
    'MAPS': {'ws':range(36864, 36868), 'l2':[4.109,2.805,1.716,8.35], 'th':[180]*3+[0]},
    'LET': {'ws':range(98304,98312), 'l2':[17.758, 17.06, 16.558, 13.164, 9.255, 1.333, 1.088, 1.088], 'th':[180]*8}
}
def _create_dummy_monitors(ws_name):
    # Creates dummy monitors for a live data workspace
    inst = mtd[ws_name].getInstrument().getName()
    CreateSimulationWorkspace(inst, [100,100,19000], UnitX='TOF', OutputWorkspace='ws_monitors')
    ExtractSpectra('ws_monitors', WorkspaceIndexList=MONDAT[inst]['ws'], OutputWorkspace='ws_monitors')
    EditInstrumentGeometry('ws_monitors', L2=MONDAT[inst]['l2'], Polar=MONDAT[inst]['th'])

def get_angle(irun, angle_workspace='angle_ws', psi_motor_name='rot', tryload=None):
    # Checks if a workspace with previous angles exists and if we've seen this run before
    if angle_workspace not in mtd:
        CreateWorkspace(OutputWorkspace=angle_workspace, DataX=0, DataY=0)
    prev_angles = {}
    if mtd[angle_workspace].getRun().hasProperty('angles_seen'):
        prev_angles = json.loads(mtd[angle_workspace].getRun().getProperty('angles_seen').value)
    if irun in sum(prev_angles.values(), []):
        angle = float([k for k, v in prev_angles.items() if irun in v][0])
    else:
        if tryload is not None:
            ws = tryload(irun)
        else:
            try:
                ws = Load(orig_file, SpectrumMax=10, LoadMonitors=True, StoreInADS=False)
            except TypeError:
                ws = Load(orig_file, SpectrumMax=10, LoadMonitors='Separate', StoreInADS=False)
        angle = ws.getRun().getLogData(psi_motor_name).value[-1]
        # Reduce to 0.2 degree accuracy to check for equivalent angles
        angle = np.round(angle * 5.) / 5.
        print(f'Read logs from run {irun}, at rotation {angle} degrees')
    angles = str(angle)
    if angles not in prev_angles.keys():
        prev_angles[angles] = []
    if irun not in prev_angles[angles]:
        prev_angles[angles].append(irun)
    # Save previous angle information to workspace
    loginfo = json.dumps(prev_angles)
    mtd[angle_workspace].getRun().addProperty('angles_seen', loginfo, '', True)
    return prev_angles[angles]

def build_angles_ws(run_list, angle_workspace, psi_motor_name):
    for irun in run_list:
        # Just run through the list to build up list of angles in the angles_workspace
        try:
            get_angle(irun, angle_workspace=angle_workspace, psi_motor_name=psi_motor_name)
        except:
            pass

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
    print(f'Copying Instrument Info to file {outfile}')
    try:
        raw_file_name = get_raw_file_from_ws(mtd[in_ws])
    except RuntimeError:
        return
    if not os.path.exists(outfile):
        outfile = os.path.join(mantid.simpleapi.config['defaultsave.directory'], os.path.basename(outfile))
    en0 = mtd[in_ws].getEFixed(mtd[in_ws].getDetector(0).getID())
    with h5py.File(raw_file_name, 'r') as raw:
        exclude = ['dae', 'detector_1', 'name']
        try:
            to_copy = set([k for k in raw['/raw_data_1/instrument'] if not any([x in k for x in exclude])])
        except KeyError:  # Live data file
            return
        if 'aperture' not in to_copy and 'mono_chopper' not in to_copy:
            return
        reps = [k for k in to_copy if k.startswith('rep_')]
        to_copy = to_copy.difference(reps)
        has_fermi = 'fermi' in to_copy
        is_main_rep = False
        thisrep = None
        if has_fermi:
            if np.abs(en0 - float(raw['/raw_data_1/instrument/fermi/energy'][()])) < (en0/20):
                is_main_rep = True
            elif len(reps) == 0:
                return     # No rep information
            fermi_grp = raw['/raw_data_1/instrument/fermi/']
            to_copy = to_copy.difference(set(['fermi']))
        if not is_main_rep and len(reps) > 0:
            thisrep = reps[np.argsort([np.abs(en0 - float(k[4:])) for k in reps])[0]]
            if np.abs(en0 - float(thisrep[4:])) > (en0/20):
                return     # No rep information
            rep_copy = [k for k in raw[f'/raw_data_1/instrument/{thisrep}']]
            if 'fermi' in rep_copy:
                has_fermi = True
                fermi_grp = raw[f'/raw_data_1/instrument/{thisrep}/fermi']
                rep_copy = rep_copy.difference(['fermi'])
        n_spec = len(raw['/raw_data_1/instrument/detector_1/spectrum_index'])
        with h5py.File(outfile, 'r+') as spe:
            spe_root = list(spe.keys())[0]
            if has_fermi:
                fields = set([k for k in fermi_grp]).difference(['energy', 'rotation_speed'])
                for fd in fields:
                    src = fermi_grp[fd]
                    h5py.Group.copy(src, src, spe[f'{spe_root}/instrument/fermi'])
                if 'rotation_speed' in fermi_grp:
                    spe[f'{spe_root}/instrument/fermi'].create_dataset('rotation_speed', (), dtype='float64')
                    spe[f'{spe_root}/instrument/fermi/rotation_speed'][()] = fermi_grp['rotation_speed'][()]
                    for atr in fermi_grp['rotation_speed'].attrs:
                        spe[f'{spe_root}/instrument/fermi/rotation_speed'].attrs[atr] = fermi_grp['rotation_speed'].attrs[atr]
            if not is_main_rep and thisrep is not None:
                for grp in rep_copy:
                    src = raw[f'/raw_data_1/instrument/{thisrep}/{grp}']
                    h5py.Group.copy(src, src, spe[f'{spe_root}/instrument/'])
                    to_copy = to_copy.difference([grp])
            for grp in to_copy:
                src = raw[f'/raw_data_1/instrument/{grp}']
                h5py.Group.copy(src, src, spe[f'{spe_root}/instrument/'])
            detroot = f'{spe_root}/instrument/detector_elements_1'
            udet = np.array(raw['/raw_data_1/isis_vms_compat/UDET'])
            spe.create_group(detroot)
            spe[detroot].attrs['NX_class'] = np.array('NXdetector', dtype='S')
            for df0, df1 in zip(['UDET', 'DELT', 'LEN2', 'CODE', 'TTHE', 'UT01'], \
                ['detector_number', 'delt', 'distance', 'detector_code', 'polar_angle', 'azimuthal_angle']):
                src = raw[f'/raw_data_1/isis_vms_compat/{df0}']
                h5py.Group.copy(src, src, spe[detroot], df1)
            for nn in range(raw['/raw_data_1/isis_vms_compat/NUSE'][0]):
                src = raw[f'/raw_data_1/isis_vms_compat/UT{nn+1:02d}']
                h5py.Group.copy(src, src, spe[detroot], f'user_table{nn+1:02d}')
            ws = mtd[in_ws]
            s2, i2 = [], []
            for ii in range(ws.getNumberHistograms()):
                d_id = ws.getSpectrum(ii).getDetectorIDs()
                s2 += d_id
                i2 += [ii+1] * len(d_id)
            _, c1, c2 = np.intersect1d(udet, s2, return_indices=True)
            d2w = -np.ones(udet.shape, dtype=np.int32)
            d2w[c1] = np.array(i2)[c2]
            # Loads masks and process
            dI = ws.detectorInfo()
            iM = np.array([dI.isMasked(ii) for ii in range(len(dI))])
            spe[detroot].create_dataset('det2work', d2w.shape, dtype='i4', data=d2w)
            spe[detroot].create_dataset('det_mask', iM.shape, dtype=np.bool_, data=iM)

#========================================================
# MARI specific functions

def remove_extra_spectra_if_mari(wsname='ws'):
    ws = mtd[wsname]
    if ws.getNumberHistograms() == 919:
        RemoveSpectra(ws, [0], OutputWorkspace=wsname)
        try:
            ws_mon = ws.getMonitorWorkspace()
        except RuntimeError:
            try:
                ws_mon = mtd[f'{wsname}_monitors']
            except KeyError:
                ws_mon = None
        if ws_mon:
            if ws_mon.getNumberHistograms() > 3:
                RemoveSpectra(ws_mon, [3], OutputWorkspace=ws_mon.name())

def shift_frame_for_mari_lowE(origEi, wsname='ws_norm', wsmon='ws_monitors'):
    ws_norm, ws_monitors = (mtd[wsname], mtd[wsmon])
    ws_out, wsmon_out = lhs_info('names')
    if origEi < 4.01:
        # If Ei < 4, mon 3 is in 2nd frame so need to shift it by 20ms
        ws_monitors = ScaleX(wsmon, 20000, Operation='Add', IndexMin=2, IndexMax=2, OutputWorkspace=wsmon_out)
        if origEi < 3.1:
            # Additionally if Ei<3.1, data will also be in 2nd frame, shift all ToF by 20ms
            ws_norm = ScaleX(wsname, 20000, Operation='Add', IndexMin=0, IndexMax=ws_norm.getNumberHistograms()-1, OutputWorkspace=ws_out)
    return ws_norm, ws_monitors

def gen_ana_bkg(quietws='MAR29313', target_ws=None):
    # Generates an analytic background workspace from the quiet counts data
    # by fitting each spectra with a decaying exponential a*exp(-b*ToF)
    if quietws not in mtd:
        Load(Filename=f'{quietws}.nxs', OutputWorkspace=quietws)
        remove_extra_spectra_if_mari(quietws)
        NormaliseByCurrent(InputWorkspace=quietws, OutputWorkspace=quietws)
    bw, bx = (100, 1)
    ws = mtd[quietws]
    if target_ws:
        wx = RebinToWorkspace(ws, target_ws, PreserveEvents=False, OutputWorkspace='bkg_his')
    else:
        wx = Rebin(ws, f'1700,{bx},19000', PreserveEvents=False, OutputWorkspace='bkg_his')
    ws = Rebin(ws, f'1700,{bw},19000', OutputWorkspace='ws_ana_tmp')
    current = ws.run().getProtonCharge()
    xx = ws.extractX(); xx = (xx[:, 1:] + xx[:, :-1]) / 2.
    yy = ws.extractY()
    def fitfu(x, *pars):
        return pars[0] * np.exp(-x * pars[1])
    fdi = {}
    for spn in range(yy.shape[0]):
        if np.max(yy[spn,:]) > (0.002):
            popt, _ = scipy.optimize.curve_fit(fitfu, xx[spn,:], yy[spn,:], p0=[6e-4, 1./4000.])
            popt[0] /= bw
            fdi[spn] = popt
            x1 = wx.readX(spn)
            y1 = fitfu((x1[1:] + x1[:-1]) / 2., *popt)
            wx.setY(spn, y1)
            wx.setE(spn, np.sqrt(y1*bw*current) / (bw*current))
        else:
            wx.setY(spn, wx.readY(spn)*0)
            wx.setE(spn, wx.readY(spn)*0)
    if '28952' in quietws:
        for spn in range(697 - 4, 759 - 4):
            bb = spn - 259
            if bb in fdi.keys():
                x1 = wx.readX(bb)
                y1 = fitfu((x1[1:] + x1[:-1]) / 2., *fdi[bb])
                wx.setY(spn, y1)
                wx.setE(spn, np.sqrt(y1*bw*current) / (bw*current))
    bkg_ev = ConvertToEventWorkspace(wx)
    return bkg_ev, wx

def mari_remove_ana_bkg(wsname):
    if sub_ana is True and not sample_cd:
        if 'bkg_ev' not in mtd:
            gen_ana_bkg()
        current = mtd[wsname].run().getProtonCharge()
        if isinstance(mtd[wsname], mantid.dataobjects.EventWorkspace):
            ws = mtd[wsname] - mtd['bkg_ev'] * current
        else:
            try:
                ws = mtd[wsname] - mtd['bkg_his'] * current
            except ValueError:
                gen_ana_bkg(target_ws=mtd[wsname])
                ws = mtd[wsname] - mtd['bkg_his'] * current

#========================================================
# Auto-Ei routine


def mode(inp):
    return max(set(inp), key=list(inp).count)


def roundlog10(val):
    expn = int(min(0, np.floor(np.log10(val)) - 2))
    return np.round(val, decimals=-expn)


def autoei(ws):
    assert hasattr(ws, 'getInstrument') and hasattr(ws, 'getRun'), 'Input must be a Mantid workspace'
    inst = ws.getInstrument().getName()
    run = ws.getRun()

    def getLog(logname):
        logv = run.getProperty(logname)
        if 'Filtered' in str(type(logv)):
            return logv.filtered_value
        t0 = run.startTime().to_datetime64()
        t1 = run.endTime().to_datetime64()
        t = logv.times
        mask = np.logical_and(t0 <= t, t <= t1)
        return logv.value[mask]


    def getfracLog(logname, frac=0.25, expn=1.):
        val = getLog(logname)
        if expn == 1.:
            return mode(val[int(len(val) * frac):])
        else:
            return mode(np.round(val[int(len(val) * frac):] * expn)) / expn


    if inst == 'LET':
        try:
            c1_freq = getfracLog('Chopper1_Disk1_speed')
        except ValueError:
            return []
        c5_mode = run.getProperty('Chopper5_slits').value[-1]
        if 'open' in c5_mode.lower():    # Whitebeam mode
            return []
        f1, f2, c4_freq = (getfracLog(nm) for nm in ['Chopper5_Disk1_speed', 'Chopper3_speed', 'Chopper4_speed'])
        ph1a, ph1b, ph2, ph3, ph4, ph5a, ph5b = (getfracLog(nm) for nm in ['Chopper1_Disk1_phase', 'Chopper1_Disk2_phase',
            'Chopper2_phase', 'Chopper3_phase', 'Chopper4_phase', 'Chopper5_Disk1_phase', 'Chopper5_Disk2_phase'])
        # Moderator Chopper distances in m
        lm1a, lm1b, lm2, lm3, lm4, lm5a, lm5b = (7.839, 7.828, 8.4, 11.749, 15.664, 23.504, 23.496)
        # Chopper time offsets in us
        tc1a, tc1b, tc3, tc4 = ((slp / frq) + yc for slp, frq, yc in zip([1.3333e5, 0.7694e5, 186944, 115000], \
            [c1_freq, c1_freq, f2, c4_freq], [7.3, 2.4, -87.2, -82.607]))
        if 'resolution' in c5_mode.lower():
            tc5a = (1.5631e6 / f1) + 0.70
            tc5b = ( 1.3149e6 / f1) + 3.13
        elif 'flux' in c5_mode.lower():
            tc5a = (1.4520e6 / f1) + 2.56
            tc5b = (1.1204e6 / f1) + 3.33
        elif 'intermediate' in c5_mode.lower():
            tc5a = (1.0630e6 / f1) + 5.88
            tc5b = (1.8149e6 / f1) + 2.17
        else:
            warnings.warning(f'Unknown mode "{c5_mode}" using High Flux settings')
            tc5a = (1.4520e6 / f1) + 2.56
            tc5b = (1.1204e6 / f1) + 3.33
        # Determines the allowed reps through each set of choppers
        tfoc = 70.   # LET sets tfoc=80 for Ei<5 and 60 otherwise, we don't know which it is so use the average
        periods = [1.e6 / nslit / frq for nslit, frq in zip([6, 2, 6, 2], [c1_freq, f2, c4_freq, f1])]
        delays = [phase + tc - tfoc for phase, tc in zip([ph1a, ph3, ph4, ph5a], [tc1a, tc3, tc4, tc5a])]
        # Chopper 2 is the bandwidth chopper and its phase determines the min and max ToF.
        tf0 = 2500 if ph2 > 50000 else ph2   # Sets a minimum such that Ei~=50 is max Ei
        tf1 = (ph2 + 27000) % 100000         # Opening width of chopper 2 is ~26ms [TODO: tune this value!]
        def inrange(tf, l):
            lim = [2286.26 * l / sqrte for sqrte in [2286.26 * lm2 / t for t in [tf0, tf1]]]
            return tf >= lim[0] and tf < lim[1]
        # Determine the reps which make it through the other choppers 1, 3, 4 and 5
        eis, eisv = ([], [])
        for l, d, p, in zip([lm1a, lm3, lm4, lm5a], delays, periods):
            eisv.append(np.array([(((2286.26 * l) / tf)**2) for tf in [(d + s*p) for s in range(-30, 30)] if inrange(tf, l)]))
        # Calculates which reps which passes Chopper 5 also make it through all the others within 5% tolerance
        for ei in eisv[3]:
            istrans = [np.any((np.abs(eisn - ei) / eisn) < 0.05) for eisn in eisv[:3]]
            if all(istrans):
                eis.append(roundlog10(ei))
        return eis


    elif inst == 'MARI':
        try:
            freq = mode(getLog('Fermi_Speed'))
        except ValueError:
            return []
        delay = getfracLog('Fermi_delay')
        lmc = 10.04   # Moderator-Fermi distance
        period = 1.e6 / freq
        try:
            phase1, phase2 = (mode(getLog(nam)) for nam in ['Phase_Thick_1', 'Phase_Thick_2'])
        except RuntimeError:
            # Pre-upgrade MARI - single slot disk
            phase_disk = mode(getLog('Phase_Thick')) + 500 # Offset used in mari_utils.gcl
            eis_disk = ((2286.26*7.05) / (phase_disk))**2
            if mode(getLog('Freq_Thick')) > 0 or (run.hasProperty('nchannels') and run.getProperty('nchannels').value > 1000):
                return [eis_disk]
            else:
                delay_calc = ((2286.26 * lmc) / np.sqrt(eis_disk))
                eis = list({((2286.26*lmc) / (delay_calc + s*period))**2 for s in range(-10, 10)})
                return [roundlog10(ei) for ei in np.sort(eis)[::-1] if ei > 10]
        try:
            ei_nominal = mode(getLog('Ei_nominal'))
        except RuntimeError:  # Old file
            ei_nominal = ((2286.26 * lmc) / delay)**2
        sqrt_ei = np.sqrt(ei_nominal)
        delay_calc = ((2286.26 * lmc) / sqrt_ei)
        t_offset_ref = {'S':2033.3/freq-5.4, 'G':1339.9/freq-7.3, 'A':-4790.69/freq+17.7}
        t_offset = delay - (delay_calc % period)
        chopper_type = min(t_offset_ref.keys(), key=lambda x:np.abs(t_offset - t_offset_ref[x]))
        nom_disk1, nom_disk2 = (((2286.26 * l) / sqrt_ei) - c for l, c in zip([7.861, 7.904], [5879., 6041.]))
        delt_disk1, delt_disk2 = (ph - nom for ph, nom in zip([phase1, phase2], [nom_disk1, nom_disk2]))
        disk_delta = delt_disk2 - delt_disk1
        slots_delta = np.round(disk_delta / 202.11) / 10
        assert slots_delta % 1.0 < 0.2, 'Bad slots calculation'
        slots = {0:[0,1,2,4], 1:[0,1], 2:[0,2], 3:[0], 4:[0]}[abs(int(round(slots_delta)))]
        disk_ref = 6 - (np.round(delt_disk1 / 202.11) / 10)
        assert disk_ref % 1.0 < 0.2, f'Bad disk calculation'
        disk = {0:disk_ref, 1:disk_ref-1, 2:1 if disk_ref==2 else 0, 3:0, 4:0}[abs(int(round(slots_delta)))]
        reps = [d-disk for d in slots]
        eis_disk = {((2286.26*lmc) / (delay_calc + s*2500.))**2 for s in reps}
        period = period / 2. if 'G' in chopper_type.upper() else period
        eis = {((2286.26*lmc) / (delay_calc + s*period))**2 for s in range(-10, 10) if (delay_calc+s*period) > 0}
        # If disk is off, assume open and let all reps through
        if abs(mode(getLog('Freq_Thick_1'))) < 1:
            eis_disk = [ei for ei in eis if ei > (2.9 if 'G' in chopper_type.upper() else 40)]
        inrange = lambda x: (x > 1 and x < 2.9) or (x > 4 and x < 1000)
        return [roundlog10(ei) for ei in np.sort(list(eis.intersection(eis_disk)))[::-1] if inrange(ei)]


    elif inst == 'MAPS':
        try:
            freq = mode(getLog('Fermi_Speed'))
        except ValueError:
            return []
        period = 1.e6 / freq
        delay_angle = getfracLog('Fermi_Delay', expn=10.)
        disk_delay = mode(getLog('Disc_Delay'))
        lmc = 10.143   # Mod-Fermi distance
        ldc = 8.831    # Mod-Disk distance
        pickup_position = 180.
        t_offsets = [20000 * (slit - pickup_position) / 360. for slit in [180., -39.1, 0., 39.1]]
        disk_eis = [((2286.372 * ldc) / tof)**2 for tof in [disk_delay + t_offset for t_offset in t_offsets] if tof > 0]
        fermi_eis = []
        for pickup_pos, t_off in zip([195.1, 14.429], [26.2, 24.53]):  # 'S', and 'A' chopper respectively
            delay = period * (delay_angle - pickup_pos) / 360.
            fermi_eis += [((2286.372 * lmc) / tf)**2 for tf in [(delay + s*period) for s in range(-10, 10)] if tf > 100 and tf < 12000]
        fermi_eis = np.array(fermi_eis)
        eis = []
        for ei in disk_eis:
            ed = np.abs(fermi_eis - ei) / fermi_eis
            if np.any(ed < 0.1):
                eis.append(roundlog10(ei))
        return eis


    elif inst == 'MERLIN':
        run = ws.getRun()
        try:
            freq = mode(getLog('Chopper_Speed'))
        except ValueError:
            return []
        delay = getfracLog('Chopper_Delay')
        disk_delay = mode(getLog('Disc_Delay'))
        rrm_mode = np.abs(disk_delay - 13700) < 10 or np.abs(disk_delay - 12400) < 10
        lmc = 9.995    # Mod-Fermi distance
        ldc = 9.2176   # Mod-Disk distance
        if not rrm_mode:
            return [roundlog10(((2286.26*ldc) / disk_delay)**2)]
        else:          # Assume using Gd chopper (t_offset = 2000/freq-5, for 'S' it is 300/freq-8.1)
            # The opto on the Fermi was replaced on 28/11/24 - need to use different constants
            m, c = (2000, -5) if run.endTime().to_datetime64() < np.datetime64('2024-11-27') else (6483, -5.5)
            period = 20000 * (25. / freq)
            tof = (2286.26 * lmc) / np.sqrt(roundlog10(((2286.26 * lmc) / (delay - m/freq - c))**2))
            tfmx = 7500 if np.abs(disk_delay - 12400) < 10 else 8500
            return [roundlog10(((2286.26*lmc) / tf)**2) for tf in [(tof + s*period) for s in range(-10, 10)] if tf > 1500 and tf < tfmx]


    else:
        raise RuntimeError(f'Instrument {inst} not supported')


#========================================================
# Continuous rotation routines
def controt_fill_in_log(ws_full, cs_block):
    print('## Continuous rotation log interval larger than 1s!')
    print('## Reconstructing log by interpolation, this can take up to a minute.')
    onesec = np.timedelta64(1, 's')
    logval = ws_full.getRun().getLogData(cs_block)
    start = ws_full.getRun().startTime().to_datetime64()
    logs = [[logval.times[ii], logval.value[ii]] for ii in range(len(logval.value)) if logval.times[ii] > start]
    AddTimeSeriesLog(ws_full, cs_block, str(logs[0][0]), logs[0][1], DeleteExisting=True)
    for ii in range(len(logs)-1):
        tdif, vdif = (logs[ii+1][0] - logs[ii][0], logs[ii+1][1] - logs[ii][1])
        if tdif > onesec and abs(vdif) > 0:
            newstep = int(tdif / onesec)
            v0 = logs[ii][1]
            vdif = vdif / newstep
            for jj in range(1, newstep):
                tim = logs[ii][0] + onesec*jj
                AddTimeSeriesLog(ws_full, cs_block, str(tim), v0 + vdif*jj)
        else:
            AddTimeSeriesLog(ws_full, cs_block, str(logs[ii][0]), logs[ii][1])


#========================================================
# Iliad driver routines

_DGRED = None

class DG_reduction_wrapper:

    subsdict = {'\nconfig':'\n#config',
                'save_dir = ':'save_dir = None #',
                'INSTRUMENT_NAME':'MARI',
                'MASK_FILE_XML':'mari_mask2023_1.xml',
                'RINGS_MAP_XML':'mari_res2013.map',
                'whitevan\s*=\s*[0-9]*':'whitevan = 28580',
                'sample\s*=\s*\\[*[\\]0-9,]+':'sample = [28727, 28728]',
                'sample_bg\s*=\s*\\[*[\\]0-9,]+':'sample_bg = None',
                'wv_file\s*=\s*[\\\'A-z0-9\\.]*':'wv_file = \'WV_28580.txt\'',
                'wv_detrange\s*=\s*[\\[\\]0-9,]*':'wv_detrange = None',
                'Ei_list\s*=\s*[\\[\\]\\.0-9,]+.*':'Ei_list = [1.84, 1.1]'}

    def __init__(self):
        self.curdir = abspath(dirname(__file__))
        self.reduction = self.load_code(os.path.join(self.curdir, 'DG_reduction.py'))
        self.whitevan = self.load_code(os.path.join(self.curdir, 'DG_whitevan.py'))
        self.monovan = self.load_code(os.path.join(self.curdir, 'DG_monovan.py'))

    def load_code(self, filename):
        loader = importlib.machinery.SourceFileLoader('<dgredwrapper>', filename)
        # Loads the main reduction script and compiles it to bytecode
        src = loader.get_data(filename).decode()
        if 'INSTRUMENT_NAME' in src:
            for ky, val in self.subsdict.items():
                src = re.sub(ky, val, src)
        x0, x1 = (src.find('#!begin_params'), src.find('#!end_params'))
        src_param = src[x0:x1]
        src_body = (src[:x0] + re.sub('\n', '\n#', src_param) + src[x1:]).encode()
        code = loader.source_to_code(src_body, filename)
        # Parses the variables section
        tmp_mod = types.ModuleType('dgredwrapper_param')
        exec(src_param, tmp_mod.__dict__)
        params = {k:getattr(tmp_mod, k) for k in dir(tmp_mod) if not k.startswith('_')}
        return code, params

    def __call__(self, mod='reduction', **kwargs):
        code, params = getattr(self, mod)
        env = params
        tmp_mod = types.ModuleType('DG_red_exec')
        tmp_mod.__file__ = os.path.join(self.curdir, f'DG_{mod}.py')
        env.update(tmp_mod.__dict__)
        env.update(kwargs)
        exec(code, env)

def run_reduction(mod='reduction', **kwargs):
    global _DGRED
    if _DGRED is None:
        _DGRED = DG_reduction_wrapper()
    if 'inst' in kwargs:
        config['default.instrument'] = kwargs['inst']
    _DGRED(mod=mod, **kwargs)

def run_whitevan(**kwargs):
    run_reduction(mod='whitevan', **kwargs)

def run_monovan(**kwargs):
    run_reduction(mod='monovan', **kwargs)

def _tryload(runno):
    if isinstance(runno, str) and os.path.exists(runno):
        print(f'{runno} file exists - pre-loading.')
        outname = os.path.basename(runno).split('.')[0]
        try:
            Load(runno, OutputWorkspace=outname, LoadMonitors=True)
        except TypeError:
            Load(runno, OutputWorkspace=outname, LoadMonitors='Separate')
        return outname
    return runno

def iliad(runno, ei, wbvan, monovan=None, sam_mass=None, sam_rmm=None, sum_runs=False, **kwargs):
    wv_name = wbvan if (isinstance(wbvan, (str, int, float)) or len(wbvan)==1) else wbvan[0]
    wv_file = f'WV_{wv_name}.txt'
    wv_args = {}
    if 'inst' in kwargs:
        kwargs['inst'] = kwargs['inst'].upper()
        config['default.instrument'] = kwargs['inst']
        wv_args = {'inst':kwargs['inst']}
    try:
        LoadAscii(wv_file, OutputWorkspace=wv_file)
        ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN', OutputWorkspace=wv_file)
    except ValueError:
        run_whitevan(whitevan=wbvan, **wv_args)
        ws = mtd['WV_normalised_integrals']
        if ws.getNumberHistograms() == 919:
            RemoveSpectra(ws, [0], OutputWorkspace=ws.name())
            SaveAscii(ws, wv_file)
    Ei_list = ei if hasattr(ei, '__iter__') else [ei]
    if 'hard_mask_file' in kwargs:
        kwargs['mask'] = kwargs.pop('hard_mask_file')
    if 'sumruns' not in kwargs and sum_runs is True:
        kwargs['sumruns'] = True
    if monovan is not None and isinstance(monovan, (int, float)) and monovan > 0:
        mv_file = f'MV_{monovan}.txt'
        mvkw = {}
        for pp in [v for v in ['mask', 'inst'] if v in kwargs]:
            mvkw[pp] = kwargs[pp]
        try:
            LoadAscii(mv_file, OutputWorkspace=mv_file)
            ReplaceSpecialValues(mv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN', OutputWorkspace=mv_file)
        except ValueError:
            run_monovan(monovan=monovan, Ei_list=Ei_list, wv_file=wv_file, **mvkw)
        kwargs['sample_mass'] = sam_mass
        kwargs['sample_fwt'] = sam_rmm
        kwargs['mv_file'] = mv_file
    runno = [_tryload(r) for r in runno] if isinstance(runno, list) else _tryload(runno)
    run_reduction(sample=runno, Ei_list=ei if hasattr(ei, '__iter__') else [ei], wv_file=wv_file, **kwargs)


#=================================================================================
#
#   PolCorr3
#
#   Routines used in DG_PLET-calibration and DG_PLET-analysis  scripts
#   
#                                              JRS, Gino Cassella and Gøran Nilsen
#                                              14/11/25
#=================================================================================

from mantid import *
from mantid.simpleapi import *
from scipy.optimize import curve_fit
import numpy as np
import os

class Reduce():
    def __init__ (self, sample, energies, PF,
                  he_pressure = 0.75,
                  he_path_length = 0.06,
                  he_mode = 'direct',
                  polmon_distance = 25.39,
                  polmon_delay = 100,
                  polmon_spectrum = 98311,
                  p0 = np.array([0.85,0.90,0.90,0.5]),
                  name_format = 'LET{0}_{1:g}meV_1to1.nxs',
                  nxsdir = '',
                  label = '',
                  mask = None,
                  rings_map = None,
                  NSF_first = False,
                  separate = False,
                  sum_runs = True,
                  event_mode = True
                  ):

        self.energies = energies
        self.sample_runs = sample
        self.PF = PF
        self.he_pressure = he_pressure          # in bar
        self.he_path_length = he_path_length    # in metres
        self.he_mode = he_mode                  # either 'direct' or 'fit' - 'direct' takes 3He polarization straight from monitor, 'fit' fits to an exponential decay - 'direct' only recommended for short runs
        self.polmon_distance = polmon_distance  # in metres
        self.polmon_delay = polmon_delay        # in microsecs
        self.polmon_spectrum = polmon_spectrum
        self.polmon_wsindex = polmon_spectrum-98305
        self.p0 = p0
        self.NSF_first = NSF_first
        self.separate = separate
        self.sum_runs = sum_runs
        self.event_mode = event_mode
        self.label = label
        self.mask = mask
        self.rings_map = rings_map
        self.name_format = name_format
        self.nxsdir = nxsdir

#=================================================================================
#---------------------------------------------------------------------------

    def generate_dummy(self, run, ei):
        """Clears out a workspace with the correct dimensions to be used
           as an empty workspace that can be cloned and populated"""
        print("      ... Generating dummy workspace...")
        Dummy = LoadNexus(self.nxsdir + '/' + self.name_format.format(run, ei))
        DummyWorkspace = mtd['Dummy']*0.0
        self.DummyWorkspace = mtd['DummyWorkspace']
        DeleteWorkspace(Dummy)

#=================================================================================
#---------------------------------------------------------------------------

    def get_PF_from_quartz(self): 
        """Calculates PF (product of polariser and flipper efficiencies) using the out of plane
           angular dependence of the flipping ratio of quartz through the 3He analyser. Several
           Ei are required to generate a reliable set of PF"""

        FAP_tofit = np.array([])
        FAP_tofite = np.array([])
        gamma = np.array([])
        wl = np.array([])
        
        for ei in self.energies:
            wl = np.append(wl,9.045 / np.sqrt(ei))
            self.generate_dummy(self.sample_runs[0],ei)

            NSF_quartz_total = CloneWorkspace(self.DummyWorkspace)
            SF_quartz_total  = CloneWorkspace(self.DummyWorkspace)

# array of out-of-plane angles
            gamma = np.append(gamma,np.deg2rad(np.linspace(-30,30,num=256)))
            
            print("****************************************")        

# read in NSF and SF quartz and totalize
            for run in self.sample_runs[::2]:
                if self.NSF_first:
                    NSF_run = run
                    SF_run = run + 1
                else:
                    NSF_run = run + 1
                    SF_run = run

                print("Loading quartz runs {0} (NSF) and {1} (SF) at {2:<3.2f}meV".format(NSF_run,SF_run, ei))
                NSF_quartz = LoadNXSPE(self.name_format.format(NSF_run, ei))
                SF_quartz  = LoadNXSPE(self.name_format.format(SF_run, ei))
                
                NSF_quartz_total += NSF_quartz
                SF_quartz_total  += SF_quartz 

            if (not "Masking" in mtd) and (self.mask is not none):
                LoadMask('let',InputFile=self.mask,RefWorkspace=NSF_quartz_total,OutputWorkspace='Masking')

# integrate over the elastic line and mask
            NSF_quartz_total = Integration(NSF_quartz_total,RangeLower=-ei*0.03,RangeUpper=ei*0.03)
            SF_quartz_total  = Integration(SF_quartz_total, RangeLower=-ei*0.03,RangeUpper=ei*0.03)
            if "Masking" in mtd:
                MaskDetectors(NSF_quartz_total, MaskedWorkspace='Masking')
                MaskDetectors(SF_quartz_total,  MaskedWorkspace='Masking')

# LET_gamma_grouping sums over 2-theta to leave workspaces as function of gamma
            NSF_quartz_total = GroupDetectors('NSF_quartz_total', MapFile='LET_gamma_grouping.xml')
            SF_quartz_total  = GroupDetectors('SF_quartz_total',  MapFile='LET_gamma_grouping.xml')

            FAP = (NSF_quartz_total - SF_quartz_total) / (NSF_quartz_total + SF_quartz_total)
            FAP = Transpose(FAP)
            FAP_tofit  = np.append(FAP_tofit,FAP.extractY()[0])
            FAP_tofite = np.append(FAP_tofite,FAP.extractE()[0])

# neutron polarization vs. gamma            
        
        def FAP_gamma(gamma,wl,PF,PHe):
            return  PF*np.tanh(7.33*wl*self.he_path_length*self.he_pressure*PHe*(1.0/np.cos(gamma)))
            
        def FAP_fit(gamma, *params):
            y_fit = np.array([])
            PHe = params[-1]
            gg = int(len(gamma)/len(self.energies))
            
            for i in range(len(self.energies)):
                PF = params[i]
                extract = gamma[i*gg:(i+1)*gg]
                y_fit = np.append(y_fit,FAP_gamma(extract,wl[i],PF,PHe))
            
            return y_fit

# Only fit to values where our secant approximations holds well and with physically meaningful numbers
# of counts
        indices = np.logical_and((np.abs(gamma) > 0.01),(np.abs(gamma) < 0.3))
        indices = np.logical_and((FAP_tofit > 0.), indices)

        try:
            import matplotlib.pyplot as pyplot
        except ModuleNotFoundError:
            pass
        else:
            fig,ax = pyplot.subplots()
            ax.plot(gamma[indices],FAP_tofit[indices])
            popt, pcov = curve_fit(FAP_fit, gamma[indices], FAP_tofit[indices],p0=self.p0, sigma=FAP_tofite[indices],bounds=(0.9*self.p0,1.1*self.p0))
            fit_y = FAP_fit(gamma[indices],popt[0],popt[1],popt[2],popt[3])
            ax.plot(gamma[indices],fit_y)
            pyplot.show()

        PF = np.zeros_like(self.energies)
        PFe = PF
        PHe = popt[-1]
        PHee = np.sqrt(np.diag(pcov))[-1]

        print("****************************************")
        for i in range(len(PF)):
            PF[i] = popt[i]
            PFe[i] = np.sqrt(np.diag(pcov))[i]
            print("PF   = {0:1.3f} +/- {1:1.3f}".format(popt[i], PFe[i]))
        print("P_He = {0:1.3f} +/- {1:1.3f}".format(PHe, PHee))
        print("****************************************") 

#=================================================================================
#-------------------------------------------------------------------------------------
# Finds the initial PHe0 and T1 of a range of runs for a single 3He cell.  If "save" is 
# True then a calibration file is used by set_helium_parameters.  Otherwise the helium
# parameters are set and may used by following routines


    def get_helium_parameters(self, save=False):
        """Finds the initial PHe0 and T1 of a range of runs for a single 3He cell.  If "save" is 
           True then a calibration file is used by set_helium_parameters.  Otherwise the helium
           parameters are set and may used by following routines"""

        ei = self.energies[0]
        PF = self.PF[0]
        wl = 9.045 / np.sqrt(ei)
        TOF = 251.9 * wl * self.polmon_distance + self.polmon_delay

# Now we have PF we can use the cell transmission monitor in a sample run to fit the lifetime
# and polarisation of the 3He cell, calculate the flipping ratio for each run, and apply the pol
# corrections. It is also possible to directly calculate the 3He polarisation for short runs.
        first = True
        t0    = 0
        PHes  = []
        PHese = []
        times = []

        print("\nStarting get_helium_parameters...")
        print("Using {0} meV rep with PF={1}".format(ei,PF))
        print("****************************************************")        
        for run in self.sample_runs[::2]:
            if self.NSF_first:
                NSF_run = run
                SF_run = run + 1
            else:
                NSF_run = run + 1
                SF_run = run

            print("Calculating P_He for runs NSF:{0} and SF:{1} at {2:<3.2f} meV at t={3:.0f} \u03BCs".format(NSF_run,SF_run,ei,TOF))

#load monitors and calculate flipping ratios
            if self.event_mode:
                NSF_monitors = LoadNexusMonitors("LET000{0}.nxs".format(NSF_run))
                SF_monitors  = LoadNexusMonitors("LET000{0}.nxs".format(SF_run))
            else:
                print('Warning: histogram mode!')
                NSF_monitors = LoadNexus("LET000{0}.nxs".format(NSF_run),self.polmon_spectrum,self.polmon_spectrum)
                SF_monitors  = LoadNexus("LET000{0}.nxs".format(SF_run),self.polmon_spectrum,self.polmon_spectrum)
            
            NSF_monitors = NormaliseByCurrent(NSF_monitors,RecalculatePCharge=True)
            SF_monitors  = NormaliseByCurrent(SF_monitors,RecalculatePCharge=True)
            
            start_time = NSF_monitors.getSampleDetails().startTime().to_datetime64()
            end_time   =  SF_monitors.getSampleDetails().endTime().to_datetime64()
            
            time = start_time + (end_time - start_time)/2.0
            
            if first:
                t0 = time
                first = False

# Integrate the monitors over the appropriate TOF     
            if self.event_mode:
                NSF_int = Integration(NSF_monitors,RangeLower=TOF-100,RangeUpper=TOF+100,StartWorkspaceIndex=self.polmon_wsindex)                                      
                SF_int  = Integration(SF_monitors, RangeLower=TOF-100,RangeUpper=TOF+100,StartWorkspaceIndex=self.polmon_wsindex)
            else:
                NSF_int = Integration(NSF_monitors,RangeLower=TOF-1000,RangeUpper=TOF+1000)
                SF_int = Integration(SF_monitors,  RangeLower=TOF-1000,RangeUpper=TOF+1000)
               
            FAP = (NSF_int - SF_int) / (NSF_int + SF_int)
            A   = FAP.readY(0)[0] / PF
            Ae_fractional = FAP.readE(0)[0] / FAP.readY(0)[0]
            PHe = np.abs(np.arctanh(A) / (7.33*wl*self.he_path_length*self.he_pressure))
            
            if not np.isfinite(PHe):
                print("ERROR: Bad P_He value...")
                return
            
            PHes.append(PHe)
            PHese.append(Ae_fractional*PHe)
            times.append(float((time - t0))*1e-9/60./60.)

        if self.he_mode == 'direct':
            self.PHe0 = PHes
            self.PHe0e = PHese
            self.times = times
            np.set_printoptions(precision=3)
            print("\nDirect cell polarizations P={}+/-{}".format(self.PHe0,self.PHe0e))
            print("************************************************************************") 
            
        elif self.he_mode == 'fit':
            print('\nFitting to extract P0 and T1')
# fit the helium polarization vs time
            def exp_decay(t,T1,P0):
                return P0 * np.exp(t/T1)
            
            popt, pcov = curve_fit(exp_decay, times, PHes, p0=[-20, 0.5], sigma=PHese)
            self.T1    = popt[0]
            self.T1e   = np.sqrt(np.diag(pcov))[0]
            self.PHe0  = popt[1]
            self.PHe0e = np.sqrt(np.diag(pcov))[1]
            self.t0    = t0
            self.times = times
            self.Phes  = PHes
            he_fit     = exp_decay(times, self.T1, self.PHe0)

# Plot PHe vs t and the fit
            try:
                import matplotlib.pyplot as pyplot
            except ModuleNotFoundError:
                pass
            else:
                pyplot.errorbar(np.array(times),np.array(PHes),yerr=np.array(np.abs(PHese)),fmt='o')
                pyplot.plot(np.array(times), np.array(he_fit))
                pyplot.title(f"T1 from cell runs: {self.sample_runs[0]} to {self.sample_runs[-1]}")
                pyplot.xlabel("Time (hours since installation)")
                pyplot.ylabel("3He polarization")
                pyplot.figtext(0.48,0.78,
                    f"PHe0 = {self.PHe0:.3f} ± {self.PHe0e:.4f} \nT1     = {-self.T1:.2f} ±  {self.T1e:.2f} hours",
                    fontsize=12)
                pyplot.show()

# Save out the calibration file
            if save:
                os.chdir(config["defaultsave.directory"])
                cal_out = f"3HeCal_{self.sample_runs[0]}-{self.sample_runs[-1]}.txt"
                file = open(cal_out, "w")
                print(self.PHe0,-self.T1,self.sample_runs[0], file=file)
                file.close

            print("\nInitial Cell polarization P0={0:.3f} ± {1:.4f} with lifetime T1={2:.2f} ± {3:.2f} hours".format(self.PHe0,self.PHe0e,-self.T1,self.T1e))
            print(f"\nCell calibration written to {cal_out}")
            print("************************************************************************") 
            
#=================================================================================
#--------------------------------------------------------------------------
# Sets the 3He cell parameters according to a specified calibration file (written by
# get_helium_parameters above) or according to explictly given parameters.

    def set_helium_parameters(self, cal=None, PHe0=0, T1=0, T0run = 0):

        """Sets the 3He cell parameters according to a specified calibration file (written by
           get_helium_parameters above) or according to explictly given parameters."""

        times = []
        print("\nStarting set_helium_parameters...")
        print("****************************************************")

        if cal is not None:
            os.chdir(config["defaultsave.directory"])
            file = open(cal, "r")
            s = file.readline()
            res = [float(x) for x in s.split()]
            PHe0  = res[0]
            T1    = res[1]
            T0run = int(res[2])
            
        if (PHe0 == 0 or T1 == 0 or T0run == 0):
            print("Must specify 3He parameters directly or in  calibration file")

#load first cell run and extract T0
        if self.event_mode:
            NSF_monitors = LoadNexusMonitors("LET000{0}.nxs".format(T0run+1))
            SF_monitors  = LoadNexusMonitors("LET000{0}.nxs".format(T0run))
        else:
            print('Warning: histogram mode!')
            NSF_monitors = LoadNexus("LET000{0}.nxs".format(T0run+1),self.polmon_spectrum,self.polmon_spectrum)
            SF_monitors  = LoadNexus("LET000{0}.nxs".format(T0run),self.polmon_spectrum,self.polmon_spectrum)
            
        start_time = NSF_monitors.getSampleDetails().startTime().to_datetime64()
        end_time   =  SF_monitors.getSampleDetails().endTime().to_datetime64()
        time = start_time + (end_time - start_time)/2.0

        t0 = time

        for run in self.sample_runs[::2]:
            if self.NSF_first:
                NSF_run = run
                SF_run = run + 1
            else:
                NSF_run = run + 1
                SF_run = run

#load monitors and calculate flipping ratios
            if self.event_mode:
                NSF_monitors = LoadNexusMonitors("LET000{0}.nxs".format(NSF_run))
                SF_monitors  = LoadNexusMonitors("LET000{0}.nxs".format(SF_run))
            else:
                print('Warning: histogram mode!')
                NSF_monitors = LoadNexus("LET000{0}.nxs".format(NSF_run),self.polmon_spectrum,self.polmon_spectrum)
                SF_monitors  = LoadNexus("LET000{0}.nxs".format(SF_run),self.polmon_spectrum,self.polmon_spectrum)
            
            start_time = NSF_monitors.getSampleDetails().startTime().to_datetime64()
            end_time   =  SF_monitors.getSampleDetails().endTime().to_datetime64()
            time = start_time + (end_time - start_time)/2.0
            time = float((time - t0))*1e-9/60./60.

            times.append(time)
            print(f"Setting run times for NSF:{NSF_run} and SF:{SF_run} to {time:.2f} hours")
        
        print(f"Setting PHe_0 = {PHe0:3f} and T1 = {T1:.2f} hours")

        self.PHe0 = PHe0
        self.T1 = -T1
        self.times = times

#=================================================================================
#--------------------------------------------------------------------------
# Corrects the data for the time-dependent cell polarization (spin leakage 
# corrections according to Scharpf/Williams) and time-dependent 3He cell transmission


    def correct_data(self):
        """Using PF, PHe0, and T1, calculate the flipping ratio for all (gamma, t) and correct runs,
           summing and averaging at the end for powder data."""

        PF_iter = iter(self.PF)

        for ei in self.energies:

            PF = next(PF_iter)
            print("****************************************")
            print("\ncorrect_data: Correcting {0:<3.2f}meV rep with PF={1} \n".format(ei,PF))
            NSF_out = "PLET_{0}_{1:<3.2f}meV_NSF".format(self.label,ei)
            SF_out  = "PLET_{0}_{1:<3.2f}meV_SF".format(self.label,ei)   

            self.generate_dummy(self.sample_runs[0],ei)
            NSF_total       = CloneWorkspace(self.DummyWorkspace)
            SF_total        = CloneWorkspace(self.DummyWorkspace)
            Scharpf_ws      = CloneWorkspace(self.DummyWorkspace)
            transmission_ws = CloneWorkspace(self.DummyWorkspace)

# Make wavelength and gamma arrays, same shape as the DummyWorkspace
            y_shape = NSF_total.extractY().shape   
            delta_e_binning = NSF_total.extractX()
            final_energy_binning = np.abs(delta_e_binning[:y_shape[0],:y_shape[1]] - ei)
            wl = 9.045 / np.sqrt(final_energy_binning)
            dummy_gammaspace      = np.zeros(y_shape).T
            dummy_gammaspace[:,:] = np.tile(np.deg2rad(np.linspace(-30,30,num=256)), 384)
            gammaspace            = dummy_gammaspace.T

# Calculate cell opacity
            opacity = 7.33 * wl * self.he_path_length * self.he_pressure
            itime = 0

# Cycle through runs - uses same "times" array as get_helium_parrameters
            for run in self.sample_runs[::2]:
                if self.NSF_first:
                    NSF_run = run
                    SF_run = run + 1
                else:
                    NSF_run = run + 1
                    SF_run = run

                print("Correcting runs NSF:{0} and SF:{1} at {2:<3.2f}meV".format(NSF_run,SF_run,ei))

                NSF_One2One = LoadNexus(self.nxsdir + self.name_format.format(NSF_run, ei))
                SF_One2One  = LoadNexus(self.nxsdir + self.name_format.format(SF_run, ei))   

# Calculate FAP, cell transmission and Scharpf correction factor for these runs

                if self.he_mode == 'direct':
                    PHe = self.PHe0[itime]
                    time = self.times[itime]
                    itime += 1
                elif self.he_mode == 'fit':
                    time = self.times[itime]
                    itime += 1
                    PHe = self.PHe0 * np.exp(time / self.T1)

                FAP = PF * np.tanh(opacity*(1.0/np.cos(gammaspace)*PHe))
                print("After {0:1.2f} hours, cell polarization is {1:1.3f}".format(time,PHe))
                flipping_ratio = (1.0 + FAP) / (1.0 - FAP)
                transmission = np.exp(-opacity) * np.cosh(opacity * PHe)
                Scharpf = (1.0 / (flipping_ratio - 1.0))

# Populate Scharpf and transmission workspaces
                for i in range(y_shape[0]):
                    Scharpf_ws.setY(i, Scharpf[i,:y_shape[1]])
                    Scharpf_ws.setX(i, delta_e_binning[i,:])
                for i in range(y_shape[0]):
                    transmission_ws.setY(i, transmission[i,:y_shape[1]])
                    transmission_ws.setX(i, delta_e_binning[i,:])
                 
# Apply corrections, totalise and rename workspace to save it
                Diff     = (NSF_One2One - SF_One2One) * Scharpf_ws
                NSF_corr = (NSF_One2One + Diff) / transmission_ws
                SF_corr  = (SF_One2One  - Diff) / transmission_ws
                
                if self.sum_runs == True:
                    NSF_total.setYUnit('')
                    SF_total.setYUnit('')
                    NSF_total = NSF_total + NSF_corr
                    SF_total  = SF_total  + SF_corr
                else:
                    RenameWorkspace(NSF_corr,OutputWorkspace=NSF_out)
                    RenameWorkspace(SF_corr,OutputWorkspace=SF_out)
                    self.one2one_output(NSF_run,SF_run,ei)
                    
            if self.sum_runs == True:
                NSF_total /= len(self.sample_runs[::2])             
                SF_total  /= len(self.sample_runs[::2])
                RenameWorkspace(NSF_total,OutputWorkspace=NSF_out)
                RenameWorkspace(SF_total,OutputWorkspace=SF_out)                
            
            print("****************************************")

#=================================================================================
#----------------------------------------------------------------------------
# combines the SF and NSF scattering to obtain the Coherent and Incoherent scattering

    def components(self):
        
        """combines the SF and NSF scattering to obtain the Coherent and Incoherent scattering"""

        self.separate = True
        print("Combining components...")
        for ei in self.energies:

            NSF_in  = "PLET_{0}_{1:<3.2f}meV_NSF".format(self.label,ei)
            SF_in   = "PLET_{0}_{1:<3.2f}meV_SF".format(self.label,ei)            
            coh_out = "PLET_{0}_{1:<3.2f}meV_coh".format(self.label,ei)
            inc_out = "PLET_{0}_{1:<3.2f}meV_inc".format(self.label,ei)

            NSF = mtd[NSF_in]
            SF  = mtd[SF_in]

            coh = NSF - 0.5 * SF 
            inc = 1.5 * SF

            RenameWorkspace(coh,OutputWorkspace=coh_out)
            RenameWorkspace(inc,OutputWorkspace=inc_out)

#=================================================================================
#---------------------------------------------------------------------------
# Normalises the data for each energy to the level of the incoherent scattering
#                    **DOESN'T WORK, DO NOT USE**

    def norm_inc(self):
        
        if not self.separate:
            print("Must separate into components before running norm_inc")
            
        print("Normalising to the incoherent scattering level...")
        for ei in self.energies:

            inc_in  = "PLET_{0}_{1:<3.2f}meV_inc".format(self.label,ei)
            coh_in  = "PLET_{0}_{1:<3.2f}meV_coh".format(self.label,ei) 
            
            inc = mtd[inc_in]
            coh = mtd[coh_in]
            
            int_inc = Integration(inc, RangeLower = -0.03 * ei, RangeUpper = 0.03 * ei, StartWorkspaceIndex=45000, EndWorkspaceIndex=47000)
            norm_inc_factor = SumSpectra(int_inc, WeightedSum=True)
            
            factor = norm_inc_factor.readY(0)[0]

            inc = Scale(inc, Factor = 2000/factor)
            coh = Scale(coh, Factor = 2000/factor)
            
            RenameWorkspace(coh,OutputWorkspace=coh_in)
            RenameWorkspace(inc,OutputWorkspace=inc_in)
            
#=================================================================================
#---------------------------------------------------------------------------
# Outputs the data in nxspe or nxs format.  If no rings map is specified, then output
# is 1to1 format.

    def output(self, type='nxspe'):
        """Outputs the data in nxspe or nxs format.  If no rings map is specified, then output
           is 1to1 format.
        """

        if self.rings_map is not None:
            format = "rings"
        else:
            format = "1to1"
            
        print(f"Writing {format} {type} files...")
        
        if self.separate:
            ext1 = "coh"
            ext2 = "inc"
        else:
            ext1 = "NSF"
            ext2 = "SF"

        for ei in self.energies:
            NSF_in  = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext1)
            SF_in   = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext2)
            NSF_out  = "PLET_{0}_{1:<3.2f}meV_{2}_{3}.{4}".format(self.label,ei,ext1,format,type)
            SF_out   = "PLET_{0}_{1:<3.2f}meV_{2}_{3}.{4}".format(self.label,ei,ext2,format,type)
            NSF = mtd[NSF_in]
            SF  = mtd[SF_in]
            
            if format == "rings":
                grouped_NSF = GroupDetectors(NSF,MapFile=self.rings_map,PreserveEvents=False,Behaviour='Average')
                grouped_SF  = GroupDetectors(SF, MapFile=self.rings_map,PreserveEvents=False,Behaviour='Average')
                if type == "nxspe":
                    SaveNXSPE(grouped_NSF,Filename=NSF_out,Efixed=ei, KiOverKfScaling=False)
                    SaveNXSPE(grouped_SF, Filename=SF_out, Efixed=ei, KiOverKfScaling=False)
                else:
                    SaveNexus(grouped_NSF, Filename=NSF_out)
                    SaveNexus(grouped_SF,  Filename=SF_out)
                    
            else:
                if type == "nxspe":
                    SaveNXSPE(NSF, Filename=NSF_out, Efixed=ei, KiOverKfScaling=False)
                    SaveNXSPE(SF,  Filename=SF_out,  Efixed=ei, KiOverKfScaling=False)
                else:
                    SaveNexus(NSF, Filename=NSF_out)
                    SaveNexus(SF,  Filename=SF_out)

#=================================================================================
#---------------------------------------------------------------------------

