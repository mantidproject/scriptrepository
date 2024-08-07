from mantid.simpleapi import CreatePeaksWorkspace, AddPeak
import numpy as np

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
                OutputWorkspace=peak_wsname, EnableLogging=False)
    for ipk in range(len(ispec)):
        wavelength = ws.readX(ispec[ipk])[indices[ipk,2]]
        # four detectors per spectrum so use one of the central ones
        detIDs = ws.getSpectrum(ispec[ipk]).getDetectorIDs()
        idet = (len(detIDs)-1)//2  # pick central pixel
        AddPeak(PeaksWorkspace=peaks, RunWorkspace=ws, TOF=wavelength, DetectorID=detIDs[idet],
            BinCount=data[indices[ipk,0], indices[ipk,1], indices[ipk,2]], EnableLogging=False)
    return peaks


def findSpectrumIndex(indices, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=5):
    """
    :param indices: indices of found peaks in 3d array
    :param ntubes: number of tubes in instrument (WISH)
    :param npix_per_tube: number of detector pixels in each tube
    :param ntubes_per_bank: number of tubes per bank
    :param nmonitors: number of monitor spectra (assumed first spectra in ws)
    :return: list of spectrum indindices of workspace corresponding to indices of found peaks in 3d array
    """
    # find bank and then reverse order of tubes in bank
    ibank = np.floor(indices[:,0]/ntubes_per_bank)
    itube = ibank*ntubes_per_bank + ((ntubes_per_bank-1) - indices[:,0] % ntubes_per_bank)
    # flip tube 
    ipix = (npix_per_tube-1) - indices[:,1]
    # get spectrum index
    specIndex = np.ravel_multi_index((itube.astype(int), ipix), 
        dims=(ntubes, npix_per_tube), order='C') + nmonitors
    return specIndex.tolist()  # so elements are of type int not numpy.int32