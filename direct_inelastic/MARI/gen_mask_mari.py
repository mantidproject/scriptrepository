from mantid.simpleapi import mtd, Load, CreateWorkspace, LoadInstrument
import numpy as np
from six import string_types

def get_det_ids(ws):
    detids = []
    specnum = []
    offset = -1
    for idx in range(ws.getNumberHistograms()):
        try:
            #detids.append(ws.getDetector(idx).getID())
            spec = ws.getSpectrum(idx)
            detids.append(spec.getDetectorIDs()[0])
            specnum.append(spec.getSpectrumNo())
        except IndexError:
            detids.append(-1)
            specnum.append(-1)
        if offset == -1 and detids[-1] > 1000:
            offset = idx + 1
    return np.array(detids), np.array(specnum), offset

def loadspec(run_id):
    if hasattr(run_id, 'getSpectrum'):
        ws = run_id
    else:
        ws_name = run_id if isinstance(run_id, string_types) else 'MAR{}'.format(run_id)
        if ws_name in mtd:
            ws = mtd[ws_name]
        else:
            ws = Load(ws_name, OutputWorkspace=ws_name)
    return ws

def map_mask_specnum(spec_list, white_run, quiet_run=None):
    ws_white = loadspec(white_run)
    ws_quiet = ws_white if quiet_run is None else loadspec(quiet_run)
    ws_xml = CreateWorkspace([0]*921, [1]*921, NSpec=921, OutputWorkspace='ws_xml')
    LoadInstrument(Workspace=ws_xml, Filename='/home/vqq25957/.mantid/instrument/MARI_Definition.xml', MonitorList='1-3', RewriteSpectraMap='True')
    
    detids_white, specnum_white, off_white = get_det_ids(ws_white)
    detids_quiet, specnum_quiet, off_quiet = get_det_ids(ws_quiet)
    detids_xml, specnum_xml, off_xml = get_det_ids(ws_xml)

    out_list = []
    for snum in spec_list:
        id_spec = np.where(specnum_quiet == snum)[0]
        d_quiet = detids_quiet[id_spec]
        id_quiet = np.where(detids_white == d_quiet)[0]
        d_xml = detids_xml[id_quiet]
        id_remaped = np.where(detids_white == d_xml)[0]
        # Check mapped id is ok
        d_white = detids_white[id_remaped - off_white]
        id_xml = np.where(detids_xml == d_white)[0]
        d_remaped = detids_white[id_xml + off_white]
        print(d_quiet)
        if d_remaped == d_quiet:
            out_list.append(id_remaped[0])
        else:
            raise RuntimeError('Wrong spectra')
    return out_list

print(map_mask_specnum([535, 647, 709], 25779, 26190))
