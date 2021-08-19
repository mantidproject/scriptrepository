from mantid.simpleapi import * 
import matplotlib.pyplot as plt
import numpy as np
import h5py
from os import path
from bragg_detect import detect_bragg_peaks


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

def make_3d_array(ws, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=5):
    """
    Extract Ydata from WISH MatrixWorkspace and shape into 3d array
    :param ws: MatrixWorkspace of WISH data with xunit wavelength
    :param ntubes: number of tubes in instrument (WISH)
    :param npix_per_tube: number of detector pixels in each tube
    :param ntubes_per_bank: number of tubes per bank
    :param nmonitors: number of monitor spectra (assumed first spectra in ws)
    :return: 3d numpy array ntubes x npix per tube x nbins
    """
    y = ws.extractY()[nmonitors:,:]  # exclude monitors - alternatively load with monitors separate?
    nbins = ws.blocksize()  # 4451
    y = np.reshape(y, (ntubes, npix_per_tube, nbins))[:,::-1,:]  # ntubes x npix x nbins (note flipped pix along tube)
    # reverse order of tubes in each bank
    nbanks = ntubes//ntubes_per_bank
    for ibank in range(0, nbanks):
        istart = ibank*ntubes_per_bank
        iend = (ibank+1)*ntubes_per_bank
        y[istart:iend,:,:] = y[istart:iend,:,:][::-1,:,:]
    return y
    
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
    
# functions to process output of ML peak finding

def findSpectrumIndex(indices, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=5):
    """
    :param indices: indices of found peaks in 3d array
    :param ntubes: number of tubes in instrument (WISH)
    :param npix_per_tube: number of detector pixels in each tube
    :param ntubes_per_bank: number of tubes per bank
    :param nmonitors: number of monitor spectra (assumed first spectra in ws)
    :return: list of spectrum indindices of workspace corresponding to indices of found peaks in 3d array
    """
    ibank = np.floor(indices[:,0]/(ntubes_per_bank-1))
    itube = ibank*ntubes_per_bank + ((ntubes_per_bank-1) - indices[:,0] % ntubes_per_bank)
    ipix = (npix_per_tube-1) - indices[:,1]
    specIndex = np.ravel_multi_index((itube.astype(int), ipix), 
        dims=(ntubes, npix_per_tube), order='C') + nmonitors
    return specIndex.tolist()  # so elements are of type int not numpy.int32

def createPeaksWorkspaceFromIndices(ws, peak_wsname, indices, data):
    """
    Create a PeaksWorkspace using indices of peaks found in 3d array. WISH has 4 detectors so put peak in central pixel.
    Could add peak using avg QLab for peaks at that lambda in all detectors but peak placed in nearest detector anyway
    for more details see https://github.com/mantidproject/mantid/issues/31944
    :param ws: MatrixWorkspace of WISH data with xunit wavelength
    :param peak_wsname: Output name of peaks workspace created
    :param indices: indices of peaks found in 3d array
    :param data: 3d array of data
    :return: PeaksWorkspace
    """
    ispec = findSpectrumIndex(indices, *data.shape[0:2])
    peaks = CreatePeaksWorkspace(InstrumentWorkspace=ws, NumberOfPeaks=0, 
                OutputWorkspace=peak_wsname)
    for ipk in range(len(ispec)):
        wavelength = ws.readX(ispec[ipk])[indices[ipk,2]]
        # four detectors per spectrum so use one of the central ones
        detIDs = ws.getSpectrum(ispec[ipk]).getDetectorIDs()
        idet = (len(detIDs)-1)//2  # pick central pixel
        AddPeak(PeaksWorkspace=peaks, RunWorkspace=ws, TOF=wavelength, DetectorID=detIDs[idet],
            BinCount=data[indices[ipk,0], indices[ipk,1], indices[ipk,2]])
    return peaks
    
if __name__ == "__main__" or "mantidqt.widgets.codeeditor.execution":
    # 1. reduce data and save 3d numpy array from ML algorithm (bragg-detect)
    runno = 49505
    do_load_hdf5 = False
    do_save_hdf5 = False
    hdf5_path = ""  # path to save WISH data array

    if not do_load_hdf5:
        # produce 3d array from raw data
        ws = load_raw_run(49505)
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
    so_save_peaks = False
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

