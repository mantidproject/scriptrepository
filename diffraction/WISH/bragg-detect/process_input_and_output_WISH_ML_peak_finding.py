from mantid.simpleapi import * 
import matplotlib.pyplot as plt
import numpy as np
import h5py
from os import path
from bragg_detect import detect_bragg_peaks
from bragg_utils import make_3d_array, createPeaksWorkspaceFromIndices

# functions to process data to use as ML input

def load_raw_run(runno):
    """
    Load raw data and pre-process
    :param runno: WISH run number to load
    :return: MatrixWorkspace
    """
    LoadRaw(Filename='WISH000'+str(runno)+'.raw', OutputWorkspace=str(runno), LoadLogFiles=False)
    CropWorkspace(InputWorkspace=str(runno), OutputWorkspace=str(runno), XMin=6000, XMax=99000)
    NormaliseByCurrent(InputWorkspace=str(runno), OutputWorkspace=str(runno))
    ConvertUnits(InputWorkspace=str(runno), OutputWorkspace=str(runno), Target='Wavelength', 
        ConvertFromPointData=False)
    ws = ConvertToPointData(InputWorkspace=str(runno), OutputWorkspace=str(runno))
    return ws

    
def save_data(data, filepath):
    """
    Save 3d array of WISH data in hdf5 file
    :param data: 3d array of WISH data
    :param filepath: filepath to save
    """
    _, fname = path.split(filepath)
    with h5py.File(filepath, "w") as f:
        dset = f.create_dataset(fname, data=data.astype('float32'), compression="gzip", compression_opts=9)
    return

def read_data(filepath):
    """
    Read 3d array of WISH data in hdf5 file
    :param filepath: filepath to read
    :return: 3d array of WISH data
    """
    with h5py.File(filepath, 'r') as f:
        dset = list(f.keys())[0]
        data = f[dset][:]
    return data

def save_found_peaks(filepath, indices):
    np.savetxt(filepath, indices, fmt='%d')
    
    
if __name__ == "__main__" or "mantidqt.widgets.codeeditor.execution":
    # 1. reduce data and save 3d numpy array from ML algorithm (bragg-detect)
    runno = 49505
    do_load_hdf5 = False
    do_save_hdf5 = False
    hdf5_path = ""  # path to save WISH data array

    if not do_load_hdf5:
        # produce 3d array from raw data
        ws = load_raw_run(runno)
        data = make_3d_array(ws)
        if do_save_hdf5 and hdf5_path:
            save_data(data, hdf5_path)
    elif hdf5_path:
        # read 3d array from file
        data = read_data(hdf5_path)
        # create instrument ws for PeakWorkspace creation
        ws = LoadEmptyInstrument(InstrumentName='WISH', OutputWorkspace='WISH')
        axis = ws.getAxis(0)
        axis.setUnit("Wavelength")

    # 2. Run ML algorithm here (or load results)
    do_load_peaks = True
    do_save_peaks = False
    peaks_path = r"C:\Users\xhg73778\Documents\bragg-detect-main\bragg-detect-main\example\result\peak_locations.txt"

    if do_load_peaks and peaks_path:
        # load from existing
        indices = np.loadtxt(peaks_path, dtype=int)
    else:
        indices = detect_bragg_peaks(data, large_peak_size=[10, 10, 50], threshold=.2, workers=1)
        if do_save_peaks and peaks_path:
            save_found_peaks(indices)
                           
    # 3. Process output of ML algorithm and create PeaksWorkspace
    peaks = createPeaksWorkspaceFromIndices(ws, 'ML_peaks', indices, data)
    
