# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from IntegratePeaksSkew import InstrumentArrayConverter, PeakData
from mantid.api import FunctionFactory

def make_3d_array(ws, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=5):
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

# powder taken with MR-SF settings (same typically used for single-crystal)
ws = Load(Filename='WISH00056081.raw', OutputWorkspace='WISH00056081')
# reduce size of data and improve stats - increase TOF bin-width by factor 3 (typically fine for peak finding on WISH)
ws = Rebunch(InputWorkspace=ws, NBunch=3, OutputWorkspace=ws.name())
ws = ConvertUnits(InputWorkspace=ws, OutputWorkspace=ws.name(), Target='dSpacing')

# set a UB matrix (taken from PdCrO2 - could use any!)
ub = np.array([[ 0.08749108,  0.37636055, -0.00169743],
               [-0.01788245,  0.00531973,  0.05517389],
               [ 0.383845  ,  0.11677283,  0.00295732]])
SetUB(ws, UB=ub)
# set goniometer (roatets UB matrix by 5 degrees around vertical axis - can add more generate a series of datasets!)

axis0_rotations = [f"{angle},0,1,0,-1" for angle in range(0, 360, 5)]

for rot_i, rot in enumerate(axis0_rotations):
    cloned_ws = ws.clone()
    SetGoniometer(cloned_ws, Axis0=rot)
    # predict peaks
    peaks = PredictPeaks(InputWorkspace=cloned_ws, WavelengthMin=0.8, WavelengthMax=10, MinDSpacing=0.8, MaxDSpacing=15, 
                         ReflectionCondition='Primitive', OutputWorkspace='peaks')
    array_converter = InstrumentArrayConverter(cloned_ws)
    minI, maxI = 1,50
    nfwhm=6
    sig_rows = 1
    sig_cols = 2
    np.random.seed(1)
    # each row in labels
    peak_labels = np.zeros((peaks.getNumberPeaks(), 4)) # each row contains [ispec, TOF index of peak cen, TOF index of peak min, TOF index of peak max]
    frac_cutoff = 0.02  # cutoff in integrated intensity used to find peak extents for labelling
    for ipk, pk in enumerate(peaks):
        # get data array in window around peak region
        peak_data = array_converter.get_peak_data(pk, pk.getDetectorID(), peaks.column("BankName")[ipk], 15, 15, nrows_edge=1, ncols_edge=1)
        ispecs = np.array(cloned_ws.getIndicesFromDetectorIDs([int(d) for d in peak_data.detids.flatten()])).reshape(peak_data.detids.shape)
        # get arrays representing distance of each pixel from peak position on the detector/instrument view
        rows, cols = np.meshgrid(np.arange(ispecs.shape[1])-peak_data.icol, np.arange(ispecs.shape[0])-peak_data.irow)  # cols, rows
        gaussian = (1/(2*np.pi*sig_rows*sig_cols))*np.exp(-0.5*(rows/sig_rows)**2 -0.5*(cols/sig_cols)**2 )
        # estimate parameters
        func = FunctionFactory.Instance().createPeakFunction("BackToBackExponential")
        func.setParameter("X0", pk.getDSpacing())  # set centre
        func.setMatrixWorkspace(cloned_ws, int(ispecs[peak_data.irow, peak_data.icol]), 0.0, 0.0)
        func_callable = FunctionWrapper(func)
        # loop over all spectra and simulate peak in given d-spacing range
        fwhm = func.fwhm()
        dmin = pk.getDSpacing() - nfwhm*fwhm
        dmax = pk.getDSpacing() + nfwhm*fwhm
        intensity = minI + np.random.rand()*(maxI-minI)
        for irow in range(ispecs.shape[0]):
            for icol in range(ispecs.shape[1]):
                ispec = int(ispecs[irow, icol])
                # get d-spacing range to evaluate peak
                istart = cloned_ws.yIndexOfX(dmin, ispec) if dmin > cloned_ws.readX(ispec)[0] else 0
                iend = cloned_ws.yIndexOfX(dmax, ispec) if dmax < cloned_ws.readX(ispec)[-1] else cloned_ws.blocksize()
                # get d-spacing values
                dspac = cloned_ws.readX(ispec)
                dspac = 0.5*(dspac[istart:iend] + dspac[istart+1:iend+1])  # convert to bin centers from edgess
                # simulate data and add to workspace
                func_callable.setParameter('I', gaussian[irow,icol]*intensity)
                ypk = func_callable(dspac)
                if irow == peak_data.irow and icol == peak_data.icol:
                    # spectrum corresponding to predicted peak position
                    # find peak extents for labeling - use cutoff in integrated intensity
                    icen = cloned_ws.yIndexOfX(pk.getDSpacing(), ispec)
                    y_cdf = np.cumsum(ypk)  # cumulative distribtuion function
                    y_cdf /= y_cdf[-1]
                    ilo = np.clip(istart + np.argmin(abs(y_cdf - frac_cutoff)), a_min=0, a_max=icen)
                    ihi = np.clip(istart + np.argmin(abs(y_cdf - (1-frac_cutoff))), a_min=icen, a_max=cloned_ws.blocksize())
                    peak_labels[ipk] = [ispec, icen, ilo, ihi]
                cloned_ws.dataY(ispec)[istart:iend] += np.random.poisson(lam=ypk)
    # convert back to TOF
    ConvertUnits(InputWorkspace=cloned_ws, OutputWorkspace=cloned_ws.name(), Target='TOF')
    labels_path = r"/home/wj1132075/Desktop/CNN_Data/Labels/labels_{}.txt".format(rot_i)
    np_ws_path = r"/home/wj1132075/Desktop/CNN_Data/Workspaces/workspace_{}.npz".format(rot_i)
    np.savetxt(labels_path, peak_labels, fmt='%d', delimiter=",")
    npws = make_3d_array(cloned_ws)
    np.savez_compressed(np_ws_path, npws)
    DeleteWorkspace(cloned_ws)