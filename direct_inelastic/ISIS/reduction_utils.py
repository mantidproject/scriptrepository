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
        if hist.name().startswith('Load') and 'Filename' in [pp.name() for pp in hist.getProperties()]:
            orig_file = hist.getPropertyValue('Filename')
            break
    if orig_file is None:
        raise RuntimeError(f'Cannot find original file from workspace {ws_name} to load logs from')
    try:
        Load(orig_file, SpectrumMax=10, LoadMonitors=True, OutputWorkspace='tmp_mons')
    except TypeError:
        Load(orig_file, SpectrumMax=10, LoadMonitors='Separate', OutputWorkspace='tmp_mons')
    ws_mon_name = f'{ws_name}_monitors'
    RenameWorkspace('tmp_mons_monitors', ws_mon_name)
    DeleteWorkspace('tmp_mons')
    mtd[ws_name].setMonitorWorkspace(mtd[ws_mon_name])
    CloneWorkspace(ws_mon_name, OutputWorkspace='ws_monitors')

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
        to_copy = set([k for k in raw['/raw_data_1/instrument'] if not any([x in k for x in exclude])])
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
        t_offset_ref = {'S':2033.3/freq-5.4, 'G':1339.9/freq-7.3}
        t_offset = delay - (delay_calc % period)
        chopper_type = min(t_offset_ref.keys(), key=lambda x:np.abs(t_offset - t_offset_ref[x]))
        nom_disk1, nom_disk2 = (((2286.26 * l) / sqrt_ei) - c for l, c in zip([7.861, 7.904], [5879., 6041.]))
        delt_disk1, delt_disk2 = (ph - nom for ph, nom in zip([phase1, phase2], [nom_disk1, nom_disk2]))
        disk_delta = delt_disk2 - delt_disk1
        slots_delta = np.round(disk_delta / 202.11) / 10
        assert slots_delta % 1.0 < 0.2, 'Bad slots calculation'
        slots = {0:[0,1,2,4], 1:[0,1], 2:[0,2], 3:[0], 4:[0]}[abs(int(slots_delta))]
        disk_ref = 6 - (np.round(delt_disk1 / 202.11) / 10)
        assert disk_ref % 1.0 < 0.2, f'Bad disk calculation'
        disk = {0:disk_ref, 1:disk_ref-1, 2:1 if disk_ref==2 else 0, 3:0, 4:0}[abs(int(slots_delta))]
        reps = [d-disk for d in slots]
        eis_disk = {((2286.26*lmc) / (delay_calc + s*2500.))**2 for s in reps}
        period = period / 2. if 'G' in chopper_type.upper() else period
        eis = {((2286.26*lmc) / (delay_calc + s*period))**2 for s in range(-10, 10)}
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

def iliad(runno, ei, wbvan, monovan=None, sam_mass=None, sam_rmm=None, sum_runs=False, **kwargs):
    wv_name = wbvan if (isinstance(wbvan, (str, int, float)) or len(wbvan)==1) else wbvan[0]
    wv_file = f'WV_{wv_name}.txt'
    wv_args = {}
    if 'inst' in kwargs:
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
    run_reduction(sample=runno, Ei_list=ei if hasattr(ei, '__iter__') else [ei], wv_file=wv_file, **kwargs)
