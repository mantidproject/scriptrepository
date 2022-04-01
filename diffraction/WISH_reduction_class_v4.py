from typing import Sequence, Optional
import numpy as np
from enum import Enum
from mantid.simpleapi import (Load, LoadParameterFile, CropWorkspace, NormaliseByCurrent, AnalysisDataService as ADS,
                              Rebunch, SmoothData, Minus, RebinToWorkspace, Divide, SetSample, SphericalAbsorption,
                              MonteCarloAbsorption, SaveReflections, IntegratePeaksMD, ConvertUnits, SetGoniometer,
                              CreatePeaksWorkspace, MaskDetectors, CombinePeaksWorkspaces, CloneWorkspace,
                              DeleteWorkspaces, FindSXPeaks, MaskBTP, ClearMaskFlag, LoadIsawUB, IndexPeaks,
                              FindUBUsingLatticeParameters, FindGlobalBMatrix, OptimizeCrystalPlacement,
                              CalculateUMatrix, ConvertToDiffractionMDWorkspace, PredictPeaks, 
                              PredictFractionalPeaks, SaveNexus, FilterPeaks,MoveInstrumentComponent, 
                              RotateInstrumentComponent, LoadParameterFile, SCDCalibratePanels,
                              ClearInstrumentParameters, RenameWorkspace, RenameWorkspace, DeleteWorkspace,
                              TransformHKL, CopySample, DeleteTableRows, ReplaceSpecialValues, NormaliseToMonitor,
                              ExtractMonitors, SmoothNeighbours, SaveIsawUB)
from FindGoniometerFromUB import getSignMaxAbsValInCol
from mantid.geometry import UnitCell, CrystalStructure
from os import path
from scipy.signal import  convolve2d
from scipy.ndimage import label
from scipy.stats import moment

##############
# Helper funcs for peak finding
#############

def findArrayIndex(ws, peak_ws, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=5):
    """
    :return: list of indices in 3d array
    """
    # convert det ids to spec indices
    ispec = np.array(ADS.retrieve(ws).getIndicesFromDetectorIDs(
        ADS.retrieve(peak_ws).column('DetID')))
    indices = np.zeros((len(ispec),3), dtype=int)
    # get tube and pixel coords
    itube, ipix = np.unravel_index(ispec - nmonitors, 
        shape=(ntubes, npix_per_tube), order='C')
    indices[:,1] = (npix_per_tube-1) - ipix  # flip tube
    # find bank and then reverse order of tubes in bank
    ibank = np.floor(itube/ntubes_per_bank)
    indices[:,0] = ibank*ntubes_per_bank + ((ntubes_per_bank-1) - itube % ntubes_per_bank)
    # find the index in the xdimension (assuming TOF)
    for irow, iw in enumerate(ispec):
        xbins = ADS.retrieve(ws).readX(int(iw))
        indices[irow,2] = np.argmin( abs(xbins - 
            ADS.retrieve(peak_ws).getPeak(irow).getTOF()) )
    return indices

def find_bg_pts_seek_skew(I, ibg):
    # sort and grow seed
    isort = np.argsort(-I[ibg]) # descending order
    iend = len(I)
    prev_skew =  moment(I[ibg[isort[:iend]]], 3)
    for istart in range(1, iend):
        this_skew  =  moment(I[ibg[isort[istart:iend]]], 3)
        if this_skew <= 0 or this_skew >= prev_skew:
            break
        else:
            prev_skew = this_skew
    return ibg[isort[istart:iend]], ibg[isort[:istart]] # bg, peak

def first_deriv_filter(std):
    x = np.linspace(-int(3.5*std), int(3.5*std), 7*std)
    return -x*np.exp(-0.5*((x/std)**2))/(np.sqrt(2*np.pi)*std**3)
    
def smooth_filter(std):
    x = np.linspace(-int(3.5*std), int(3.5*std), 7*std)
    return np.exp(-0.5*((x/std)**2))/(np.sqrt(2*np.pi)*std)

def find_peak_bounds(yy, ee_sq, istart, iend, ntol = 5, scale=1.15):                          
    # get initial seed from nbg_pts around max in half window
    yderiv = np.convolve(yy[istart:iend], first_deriv_filter(3), 'same')
    iend = np.argmin(yderiv) + istart
    istart = np.argmax(yderiv) + istart
    prev_IoverSig = np.sum(yy[istart:iend])/np.sqrt(np.sum(ee_sq[istart:iend]))
    nbad = 0
    for ishift in range(1,istart):
        this_IoverSig = np.sum(yy[istart-ishift:iend])/(
            np.sqrt(np.sum(ee_sq[istart-ishift:iend])))
        if this_IoverSig > prev_IoverSig:
            prev_IoverSig = this_IoverSig
            nbad = 0
        elif nbad < ntol:
            nbad = nbad + 1
        else:
            break
    istart = istart - (ishift - nbad)
    prev_IoverSig = np.sum(yy[istart:iend])/np.sqrt(np.sum(ee_sq[istart:iend]))
    nbad = 0
    for ishift in range(1, yy.size - iend):
        this_IoverSig = np.sum(yy[istart:iend+ishift])/(
            np.sqrt(np.sum(ee_sq[istart:iend+ishift])))
        if this_IoverSig > prev_IoverSig:
            prev_IoverSig = this_IoverSig
            nbad = 0
        elif nbad < ntol:
            nbad = nbad + 1
        else:
            break
    iend = iend + (ishift - nbad)
    # according to https://journals.iucr.org/a/issues/1974/04/00/a10706/a10706.pdf
    # I/sig always underestimates peak limits (for Gauss?) by ~10-15%
    pad = int(scale*(iend-istart)/2)
    return max(istart-int(0.5*pad), 0), min(iend+pad, len(yy)-1)

def is_peak_mask_valid(peak_mask, npk_min=3, density_min=0.35, 
        ntubes_max = 8, npix_max = 15, nholes_max = 2, ncavity_max=0):
    if peak_mask.sum() < npk_min:
        return False
    if np.sum(peak_mask.sum(axis=1)>0) > ntubes_max:
        return False
    if np.sum(peak_mask.sum(axis=0)>0) > npix_max:
        return False
    density = peak_mask.sum()/( np.sum(peak_mask.sum(axis=0)>0)*np.sum(peak_mask.sum(axis=1)>0) )
    if density < density_min:
        return False
    # check for holes
    kernel=np.array([[0,1,0],[1,0,1],[0,1,0]])/4
    hole_mask = convolve2d(peak_mask, kernel, mode='same')
    norm = convolve2d(np.ones(hole_mask.shape), kernel, mode='same')
    nholes = np.sum(np.logical_and((hole_mask/norm) > 0.9, ~peak_mask))
    if nholes > nholes_max:
        return False
    # check for cavities (2 adjacent bg points)
    kernel=np.array([[0,1,1,0],[1,0,0,1],[0,1,1,0]])/6
    hole_mask = convolve2d(peak_mask, kernel, mode='same')
    norm = convolve2d(np.ones(hole_mask.shape), kernel, mode='same')
    ncavity = np.sum(np.logical_and((hole_mask/norm) > 0.9, ~peak_mask))
    if ncavity > ncavity_max:
        return False
    kernel=np.array([[0,1,0],[1,0,1],[1,0,1],[0,1,0]])/6
    hole_mask = convolve2d(peak_mask, kernel, mode='same')
    norm = convolve2d(np.ones(hole_mask.shape), kernel, mode='same')
    ncavity = np.sum(np.logical_and((hole_mask/norm) > 0.9, ~peak_mask))
    if ncavity > ncavity_max:
        return False 
    else:
        return True
        
def make_3d_array(ws, ntubes=1520, npix_per_tube=128, ntubes_per_bank=152, nmonitors=0): # 5 monitors if not excluded
    x = ws.extractX()[nmonitors:,:]  # exclude monitors - alternatively load with monitors separate?
    y = ws.extractY()[nmonitors:,:]  # exclude monitors - alternatively load with monitors separate?
    e = ws.extractE()[nmonitors:,:]  # exclude monitors - alternatively load with monitors separate?
    nbins = ws.blocksize()  # 4451
    x = np.reshape(x, (ntubes, npix_per_tube, nbins+1))[:,::-1,:]  # ntubes x npix x nbins (note flipped pix along tube)
    y = np.reshape(y, (ntubes, npix_per_tube, nbins))[:,::-1,:]  # ntubes x npix x nbins (note flipped pix along tube)
    e = np.reshape(e, (ntubes, npix_per_tube, nbins))[:,::-1,:]  # ntubes x npix x nbins (note flipped pix along tube)
    # reverse order of tubes in each bank
    nbanks = ntubes//ntubes_per_bank
    for ibank in range(0, nbanks):
        istart = ibank*ntubes_per_bank
        iend = (ibank+1)*ntubes_per_bank
        x[istart:iend,:,:] = x[istart:iend,:,:][::-1,:,:]
        y[istart:iend,:,:] = y[istart:iend,:,:][::-1,:,:]
        e[istart:iend,:,:] = e[istart:iend,:,:][::-1,:,:]
    return x, y, e
      
def integrate_peaks_table_skew(ws_run, pk_ws, xbank, ybank, ebank,
    dpixel = 8, integrate_on_edge = False, dt0_over_t0 = 0.04, dth = 0.015, 
    ipks = None, do_plot = False, out_pk_wsname = 'peaks_int'):
    # Empty table workspace
    pk_ws_int = CreatePeaksWorkspace(InstrumentWorkspace=ws_run, NumberOfPeaks=0, 
        OutputWorkspace=out_pk_wsname)
    CopySample(InputWorkspace=pk_ws, OutputWorkspace=pk_ws_int, CopyName=False, 
        CopyMaterial=False, CopyEnvironment=False, CopyShape=False)  # copy UB
    # extract data
    indices = findArrayIndex(ws_run.name(), pk_ws.name(), nmonitors=0)
    # loop over each peak
    ipks = range(pk_ws.getNumberPeaks()) if ipks is None else ipks
    for ipk in ipks:
        pk = pk_ws.getPeak(ipk)
        # get TOF and TOF window
        pkTOF = pk.getTOF()
        dTOF = pkTOF*np.sqrt(dt0_over_t0**2 + (dth/np.tan(pk.getScattering()/2))**2)
        # find TOF indices over which to itnegrate
        ii, jj  = indices[ipk,0], indices[ipk,1] # peak center indices in data array
        xpk = 0.5*(xbank[ii,jj,:-1] + xbank[ii,jj,1:]) 
        frac_dTOF = 1 # so as not to pick up sattelites with similar d-spacing)
        istart = np.argmin(abs(xpk - (pkTOF - 0.4*frac_dTOF*dTOF)))
        iend = np.argmin(abs(xpk - (pkTOF + 0.6*frac_dTOF*dTOF)))
        # view window on detector integrated over TOF
        jj_min, jj_max = max(jj-dpixel, 0), jj+dpixel+1
        bank_limit = int(ybank.shape[0]/2)
        if ii < bank_limit:
            ii_min, ii_max = max(ii-dpixel,0), min(ii+dpixel+1, bank_limit + 1)
        else:
            ii_min, ii_max = max(ii-dpixel,bank_limit), ii+dpixel+1
        # find peaks on instrument vied
        ypeak2D = ybank[ii_min:ii_max,jj_min:jj_max, istart:iend].sum(axis=2)
        ibg2D, ipeak2D = find_bg_pts_seek_skew(ypeak2D.flatten(), np.arange(ypeak2D.size))
        ipeak2D_x, ipeak2D_y = np.unravel_index(ipeak2D, ypeak2D.shape)
        # find contigous region containing peak
        mask = np.zeros(ypeak2D.shape, dtype=bool)
        mask[ipeak2D_x, ipeak2D_y] = True
        labeled_array, num_features = label(mask)
        peak_mask = labeled_array == labeled_array[ii-ii_min, jj-jj_min]  # label at peak pos
        # check if mask touches edge of detector
        on_edge = False
        if ii < bank_limit:
            if ii_min == 0:
                on_edge = peak_mask.sum(axis=1)[0] > 0
            if ii_max > bank_limit:
                on_edge = peak_mask.sum(axis=1)[-1] > 0
        else:
            if ii_min == bank_limit:
                on_edge = peak_mask.sum(axis=1)[0] > 0
            if ii_max > ybank.shape[0]:
                on_edge = peak_mask.sum(axis=1)[-1] > 0
        if jj_min == 0:
            on_edge = peak_mask.sum(axis=0)[0] > 0
        if jj_max > ybank.shape[1]:
            on_edge = peak_mask.sum(axis=0)[-1] > 0
        if (on_edge and not integrate_on_edge) or not is_peak_mask_valid(peak_mask) or not mask[ii-ii_min, jj-jj_min]:
            intens, sig = 0.0, 0.0
        else:
            # increase mask to include adjacent pixels (increase stats and ensure got all peak)
            kernel= np.ones((3,3))
            bg_mask = convolve2d(peak_mask, kernel, mode='same')
            norm = convolve2d(np.ones(bg_mask.shape), kernel, mode='same')
            bg_mask = (bg_mask/norm) > 0
            bg_mask = np.logical_and(bg_mask, ~mask)  # mask includes all non bg points - not just peak label
            # proceed to find peak in 1D
            ypk = ybank[ii_min:ii_max,jj_min:jj_max, :][peak_mask].sum(axis=0)
            epk = np.sqrt((ebank[ii_min:ii_max,jj_min:jj_max, :][peak_mask]**2).sum(axis=0))
            # subtract background from bgshell normalised by ratio of npixels
            ybg = ybank[ii_min:ii_max,jj_min:jj_max, :][bg_mask].sum(axis=0)*peak_mask.sum()/bg_mask.sum()
            ebg_sq = np.sum(ebank[ii_min:ii_max,jj_min:jj_max, :][bg_mask]**2, axis=0)*((peak_mask.sum()/bg_mask.sum())**2)
            # replace with avg. error in same quarter of spectrum
            ix_pk = np.argmin(abs(xpk - pkTOF))
            width = len(xpk)/4
            iquarter = ix_pk//width
            ebg_quarter = np.sqrt(ebg_sq[int(iquarter*width):int((iquarter+1)*width)])
            epk[epk==0] = np.mean(ebg_quarter[ebg_quarter>0])  # look in quarter spectrum 
            # smooth bg and subtract
            ybg = np.convolve(ybg, smooth_filter(1), 'same')
            ebg_sq = np.convolve(ebg_sq, smooth_filter(1), 'same')
            ypk = ypk - ybg  # removes background and powder lines (roughly)
            epk = np.sqrt(epk**2 + ebg_sq)
            # normalise by bin width
            dx = np.diff(xbank[ii,jj,:])
            ypk = ypk*dx
            epk = epk*dx
            # find bg and pk bins in focused spectrum by maximising I/sig
            frac_dTOF = 1.2 # so as not to pick up sattelites with similar d-spacing)
            istart_seed = np.argmin(abs(xpk - (pkTOF - 0.4*frac_dTOF*dTOF)))
            iend_seed = np.argmin(abs(xpk - (pkTOF + 0.6*frac_dTOF*dTOF)))
            istart_pk, iend_pk = find_peak_bounds(ypk, epk**2, istart_seed, iend_seed, ntol = 10) 
            intens = np.sum(ypk[istart_pk:iend_pk])
            sig = np.sqrt(np.sum(epk[istart_pk:iend_pk]**2))
            if do_plot:
                fig, ax = plt.subplots(1,2,subplot_kw={'projection': 'mantid'})
                vmax = 0.5*(ypeak2D[peak_mask].min() + ypeak2D[peak_mask].mean())
                img = ax[0].imshow(ypeak2D, vmax=vmax)  # data.sum(axis=2).T[::-1,::-1] to look like cylindrical Y
                ax[0].plot(*np.where(bg_mask.T), 'xw')
                ax[0].plot(*np.where(np.logical_and(peak_mask, ~bg_mask).T), 'xr')
                ax[0].plot(jj-jj_min, ii-ii_min, 'or')
                ax[0].set_title(str(ipk) + ' ' + 
                    str(pk.getIntHKL()).replace('[','(').replace(']',')') + 
                    ' lambda = ' +  str(np.round(pk.getWavelength(),2)) + ' Ang') 
                ax[1].errorbar(xpk, ypk, yerr=epk, marker='o', markersize=3, capsize=2, ls='', color='k', alpha=0.5)
                ax[1].axvline(xpk[istart], ls='--', color='r')
                ax[1].axvline(xpk[iend], ls='--', color='r')
                ax[1].axvline(xpk[istart_pk], ls='-', color='g')
                ax[1].axvline(xpk[iend_pk], ls='-', color='g')
                ax[1].axhline(0, ls=':', color='r')
                # set axes limits
                ipad = int((iend-istart)/5)
                imin = max(min(istart_pk,istart)-ipad, 0)
                imax = min(max(iend, iend_pk)+ipad, len(xpk)-1)
                ax[1].set_xlim(xpk[imin], xpk[imax])
                ymin = np.min(ypk[imin:imax]-2*epk[imin:imax])
                ymax = np.max(ypk[imin:imax]+2*epk[imin:imax])
                ax[1].set_ylim(ymin, ymax)
                ax[1].set_title('I/sig = '+  str(np.round(intens/sig,2)))
                fig.colorbar(img, orientation = 'horizontal', ax=ax[0])
                fig.tight_layout()
                fig.show()
        # update peak (incl. Lorz factor)
        th = pk.getScattering()/2
        wl = pk.getWavelength()
        L = (np.sin(th)**2)/(wl**4)
        pk.setIntensity(intens*L)
        pk.setSigmaIntensity(sig*L)
        pk_ws_int.addPeak(pk)
    return pk_ws_int 
            
##############
# Reduction class
#############        

class PEAK_TYPE(Enum):
    FOUND = 0
    PREDICTED = 1
    PREDICTED_SAT = 2
    
class WISHReduction:
    def __init__(self,
                 vanadium_runno: str):
        self.runs = dict()  # runno: {MD:, peaks_found}
        self.van_runno = vanadium_runno
        self.van_ws = None
        self.sphere_shape = '''<sphere id="sphere">
                               <centre x="0.0"  y="0.0" z="0.0" />
                               <radius val="0.0025"/>
                               </sphere>'''  # sphere radius 2.5mm
        self.xtal_structure = None
        
    def set_van_ws(self, van_ws):
        self.van_ws = van_ws
        
    def set_run_data(self, run, data_dict):
        self.runs[str(run)] = data_dict
        
    def set_xtal_structure(*args, **kwargs):
        self.xtal_structure = CrystalStructure(*args, **kwargs)
        
    def delete_run_data(self, run):
        run = str(run)
        for field in ['ws', 'MD']:
            if self.runs[run][field]:
                DeleteWorkspace(self.runs[run][field])
                self.runs[run][field] = None
            
    def process_data(self,
                     runs: Sequence[str],
                     find_peaks = True,
                     bg = 0.5,
                     sample_dict: Optional[dict] = None,
                     nEventsMC = 1200,
                     gonio_motor = "wccr",
                     make_MD = True):

        if not self._is_vanadium_processed():
            return

        for run in runs:
            wsname = self.load_run(run)
            # set goniometer
            if ADS.retrieve(wsname).run().hasProperty(gonio_motor):
                wccr = ADS.retrieve(wsname).run().getPropertyAsSingleValueWithTimeAveragedMean(gonio_motor)
                SetGoniometer(Workspace=wsname, Axis0=str(wccr) +',0,1,0,1')  # vertical
            # find peaks before normalised
            peaks = self.find_sx_peaks(wsname, bg=bg, out_pk_wsname=wsname+'_peaks') if find_peaks else None
            # correct for empty counts and normalise by vanadium
            ConvertUnits(InputWorkspace=wsname, OutputWorkspace=wsname, Target='Wavelength')
            self._divide_workspaces(wsname, self.van_ws)
            # set sample geometry and material
            if sample_dict:
                SetSample(wsname, **sample_dict)
                # correct for attenuation
                if "<sphere" in ADS.retrieve(wsname).sample().getShape().getShapeXML():
                    transmission = SphericalAbsorption(InputWorkspace=wsname, 
                        OutputWorkspace='transmission')
                else:
                    transmission = MonteCarloAbsorption(InputWorkspace=wsname,  
                        OutputWorkspace='transmission', EventsPerPoint=nEventsMC)
                self._divide_workspaces(wsname, transmission)
                DeleteWorkspace(transmission)
            # convert to md
            if make_MD:
                md_name = wsname + '_MD'
                self._normalise_by_bin_width(wsname) # normalise by bin-width in K = 2pi/lambda
                ConvertToDiffractionMDWorkspace(InputWorkspace=wsname, OutputWorkspace=md_name,
                                                LorentzCorrection=True, OneEventPerBin=False)
                # undo normalisation by bvin width
                self._normalise_by_bin_width(wsname, undo=True)
            else:
                md_name = None
            ConvertUnits(InputWorkspace=wsname, OutputWorkspace=wsname, Target='TOF')
            
            # save results in dictionary
            self.runs[str(run)] = {'ws': wsname, 'MD': wsname + '_MD', 'found_pks': peaks}
        
    def load_isaw_ub(self, isaw_files: Sequence[str], runs: Optional[Sequence[str]] = None):
        if runs is None:
            runs = self.runs.keys()
        if len(isaw_files) == 1 and runs is None:
            isaw_files = len(runs)*isaw_files  # i.e. use same file for all runs
        if len(isaw_files) == len(runs):
            for run, isaw_file in zip(runs, isaw_files):
                try:
                    LoadIsawUB(InputWorkspace=self.runs[run]['ws'], Filename=isaw_file)
                    if self.runs[run]['found_pks'] is not None:
                        LoadIsawUB(InputWorkspace=self.runs[run]['found_pks'], Filename=isaw_file)
                except:
                    print('LoadIsawUB failed for run ' + run)
            
    def find_ub_using_lattice_params(self, global_B, tol,  *args, **kwargs):
        if global_B:
            FindGlobalBMatrix(PeakWorkspaces=[data['found_pks'] for data in self.runs.values()], *args, **kwargs)
        else:
            for run, data in self.runs.items():
                FindUBUsingLatticeParameters(PeaksWorkspace=data['found_pks'], Tolerance=tol, 
                    *args, **kwargs)
                IndexPeaks(PeaksWorkspace=data['found_pks'], Tolerance=tol, RoundHKLs=True)
                
    def calc_U_matrix(self, *args, **kwargs):
        for run, data in self.runs.items():
            CalculateUMatrix(PeaksWorkspace=data['found_pks'], *args, **kwargs)
    
    def calibrate_sample_pos(self, tol=0.1):
        for run, data in self.runs.items():
            IndexPeaks(PeaksWorkspace=data['found_pks'], Tolerance=tol, RoundHKLs=True)
            OptimizeCrystalPlacement(PeaksWorkspace=data['found_pks'],
                                     ModifiedPeaksWorkspace=data['found_pks'],
                                     FitInfoTable=run + '_sample_pos_fit_info',
                                     AdjustSampleOffsets=True,
                                     OptimizeGoniometerTilt=True,
                                     MaxSamplePositionChangeMeters=0.025,
                                     MaxIndexingError=tol)
            IndexPeaks(PeaksWorkspace=data['found_pks'], Tolerance=tol, RoundHKLs=True)
                                         
    def predict_peaks(self, *args, **kwargs):
        for run, data in self.runs.items():
            ws = None
            if data['found_pks'] is not None and ADS.retrieve(data['found_pks']).sample().hasOrientedLattice():
                ws = data['found_pks']
            elif data['ws'] is not None:
                ws = data['ws']
            if ws is None:
                continue  # skip
            pred_pks = PredictPeaks(InputWorkspace=ws, 
                OutputWorkspace=ADS.retrieve(ws).name() + '_predicted', *args, **kwargs)
            self.runs[run]['predicted_pks'] = pred_pks.name()
            
    def predict_satellite_peaks(self, *args, **kwargs):
        for run, data in self.runs.items():
            ws = None
            if data['found_pks'] is not None and ADS.retrieve(data['found_pks']).sample().hasOrientedLattice():
                ws = data['found_pks']
            elif data['ws'] is not None:
                ws = data['ws']
            if ws is None:
                continue  # skip
            pred_pks = PredictFractionalPeaks(InputWorkspace=ws, 
                OutputWorkspace=ADS.retrieve(ws).name() + '_predicted_satellite', *args, **kwargs)
            self.runs[run]['predicted_satellite_pks'] = pred_pks.name()
            
    def integrate_data(self, peak_type, tol, *args, **kwargs):
        for run, data in self.runs.items():
            if data['MD'] is not None:
                if peak_type == PEAK_TYPE.FOUND and tol > 0:
                    IndexPeaks(PeaksWorkspace=data['found_pks'], Tolerance=tol, RoundHKLs=True)
                    FilterPeaks(InputWorkspace=data['found_pks'], OutputWorkspace=data['found_pks'], 
                        FilterVariable='h^2+k^2+l^2', FilterValue=0, Operator='>')
                    pk_table =  data['found_pks']
                elif peak_type == PEAK_TYPE.PREDICTED:
                    pk_table = data['predicted_pks']
                else:
                    pk_table = data['predicted_satellite_pks']
                WISHReduction.remove_peaks_on_edge(pk_table)
                # only mask tubes on panels with edge facing beam in/out
                MaskBTP(Workspace=peaks, Bank='5-6', Tube='152')
                MaskBTP(Workspace=peaks, Bank='1,10', Tube='1')
                # integrate
                peaks_int = IntegratePeaksMD(InputWorkspace=data['MD'], 
                    PeaksWorkspace=pk_table, OutputWorkspace=ADS.retrieve(pk_table).name() + '_int',
                    MaskEdgeTubes=False, IntegrateIfOnEdge=False, UseOnePercentBackgroundCorrection=False,
                    *args, ** kwargs)
                self.runs[run]['integrated_pks'] = peaks_int.name()
            
    def integrate_data_skew(self, peak_type, tol, *args, **kwargs):
        for run, data in self.runs.items():
            if data['ws'] is not None:
                if peak_type == PEAK_TYPE.FOUND and tol > 0:
                    IndexPeaks(PeaksWorkspace=data['found_pks'], Tolerance=tol, RoundHKLs=True)
                    FilterPeaks(InputWorkspace=data['found_pks'], OutputWorkspace=data['found_pks'], 
                        FilterVariable='h^2+k^2+l^2', FilterValue=0, Operator='>')
                    pk_table =  data['found_pks']
                elif peak_type == PEAK_TYPE.PREDICTED:
                    pk_table = data['predicted_pks']
                else:
                    pk_table = data['predicted_satellite_pks']
                WISHReduction.remove_peaks_on_edge(pk_table)
                # extract data
                xbank, ybank, ebank =  make_3d_array(ADS.retrieve(data['ws']))
                peaks_int = integrate_peaks_table_skew(ADS.retrieve(data['ws']), 
                    ADS.retrieve(pk_table), xbank, ybank, ebank, 
                    out_pk_wsname = ADS.retrieve(pk_table).name() + '_int_skew', 
                    *args, **kwargs)                
                self.runs[run]['integrated_pks'] = peaks_int.name() 

    def save_integrated_peaks(self, save_dir: str, save_format: str,
                              split_files: Optional[bool] = False):
        all_integrated_peaks = CreatePeaksWorkspace(InstrumentWorkspace=self.van_ws, NumberOfPeaks=0)
        # get first run - use as reference UB and copy sample
        ws_ref = next(iter(self.runs.keys()))
        for run, data in self.runs.items():
            # ensure indexing is consistent
            # self.make_UB_consistent(self.runs[ws_ref]['integrated_pks'], 
            #     data['integrated_pks'])
            # save reflections
            filepath = path.join(save_dir, '_'.join([data['integrated_pks'], save_format])) + '.int'
            SaveReflections(InputWorkspace=data['integrated_pks'], Filename=filepath,
                Format=save_format, SplitFiles=split_files)
            SaveNexus(InputWorkspace=data['integrated_pks'], Filename=filepath[:-3] + 'nxs')
            # append to combined table
            all_integrated_peaks = CombinePeaksWorkspaces(LHSWorkspace=all_integrated_peaks, 
                RHSWorkspace=data['integrated_pks'])
        if len(self.runs) > 1:
            # copy lattice and UB from last run
            CopySample(InputWorkspace=data['integrated_pks'], OutputWorkspace=all_integrated_peaks , 
                CopyName=False, CopyMaterial=False, CopyEnvironment=False, CopyShape=False)
            # save with run range in filename
            min_ws = min(self.runs.keys(), key = lambda k: int("".join(filter(str.isdigit, k))))
            max_ws = max(self.runs.keys(), key = lambda k: int("".join(filter(str.isdigit, k))))
            filename = 'WISH000' + '-'.join([min_ws,max_ws]) + data['integrated_pks'][len(min_ws)+7:]
            filepath = path.join(save_dir, '_'.join([filename, save_format]) + '.int')
            SaveReflections(InputWorkspace=all_integrated_peaks, Filename=filepath,
                    Format=save_format, SplitFiles=split_files)
            SaveNexus(InputWorkspace=all_integrated_peaks, Filename=filepath[:-3] + 'nxs')
            RenameWorkspace(InputWorkspace=all_integrated_peaks, OutputWorkspace=filename)
        else:
            DeleteWorkspace(all_integrated_peaks)
            
    def remove_forbidden_peaks(self, peak_type = PEAK_TYPE.FOUND):
        if self.xtal_structure is None:
            return
        for run, data in self.runs.items():
            if peak_type == PEAK_TYPE.FOUND:
                pk_table =  data['found_pks']
            elif peak_type == PEAK_TYPE.PREDICTED:
                pk_table = data['predicted_pks']
            else:
                pk_table = data['predicted_satellite_pks']
            iforbidden = []
            for ipk, pk in enumerate(pk_table):
                if not self.xtal_structure.getSpaceGroup().isAllowedReflection(pk.getIntHKL()):
                    iforbidden.append(ipk)
            DeleteTableRows(TableWorkspace=pk_table, Rows=iforbidden)
            
    def process_vanadium(self, npoints=301):
        # vanadium
        self.van_ws = self.load_run(self.van_runno)
        ConvertUnits(InputWorkspace=self.van_ws, OutputWorkspace=self.van_ws, Target='Wavelength')
        SmoothNeighbours(InputWorkspace=self.van_ws, OutputWorkspace=self.van_ws,Radius=3)
        SmoothData(InputWorkspace=self.van_ws, OutputWorkspace=self.van_ws, NPoints=npoints)
        # SmoothData(InputWorkspace=self.van_ws, OutputWorkspace=self.van_ws, NPoints=int(npoints/5))
        # correct vanadium for absorption
        SetSample(self.van_ws, Geometry={'Shape': 'CSG', 'Value': self.sphere_shape},
                  Material={'ChemicalFormula': 'V0.95-Nb0.05', 'SampleNumberDensity': 0.0722})
        vanadium_transmission = SphericalAbsorption(InputWorkspace=self.van_ws, OutputWorkspace='vanadium_transmission')
        self._divide_workspaces(self.van_ws, vanadium_transmission)
        DeleteWorkspace(vanadium_transmission)
        
    def load_run(self, runno):
        wsname = 'WISH000' + str(runno)
        ws = Load(Filename=wsname+ '.raw', OutputWorkspace=wsname)
        # replace empty bin errors with 1 count
        si = ws.spectrumInfo()
        for ispec in range(ws.getNumberHistograms()):
            if si.hasDetectors(ispec) and not si.isMonitor(ispec):
                y = ws.readY(ispec)
                e = ws.dataE(ispec)
                e[y<1e-15] = 1
                ws.setE(ispec, e)
        CropWorkspace(InputWorkspace=wsname, OutputWorkspace=wsname, XMin=6000, XMax=99000)
        NormaliseToMonitor(InputWorkspace=wsname, OutputWorkspace=wsname, MonitorID=4)
        ReplaceSpecialValues(InputWorkspace=wsname, OutputWorkspace=wsname, 
            NaNValue=0, InfinityValue=0)
        NormaliseByCurrent(InputWorkspace=wsname, OutputWorkspace=wsname)
        ExtractMonitors(InputWorkspace=wsname, DetectorWorkspace=wsname, MonitorWorkspace=wsname + '_mon')
        return wsname

    def _is_vanadium_processed(self):
        return self.van_ws is not None and ADS.doesExist(self.van_ws)

    @staticmethod
    def _normalise_by_bin_width(wsname, undo=False):
        ws = ConvertUnits(InputWorkspace=wsname, OutputWorkspace=wsname, Target='Momentum')
        si = ws.spectrumInfo()
        for ispec in range(ws.getNumberHistograms()):
            if si.hasDetectors(ispec) and not si.isMonitor(ispec):
                # divide y and e by bin-width
                xedges = ws.readX(ispec)
                dx = np.diff(xedges)
                scale = dx if not undo else 1/dx
                ws.setY(ispec, ws.readY(ispec) / scale)
                ws.setE(ispec, ws.readE(ispec) / scale)

    @staticmethod
    def _minus_workspaces(ws_lhs, ws_rhs):
        RebinToWorkspace(WorkspaceToRebin=ws_rhs, WorkspaceToMatch=ws_lhs,
                         OutputWorkspace=ws_rhs, PreserveEvents=False)
        Minus(LHSWorkspace=ws_lhs, RHSWorkspace=ws_rhs, OutputWorkspace=ws_lhs)
        ReplaceSpecialValues(InputWorkspace=ws_lhs, OutputWorkspace=ws_lhs, 
            NaNValue=0, InfinityValue=0, BigNumberThreshold=1e15, SmallNumberThreshold=1e-15)
            
    @staticmethod
    def _divide_workspaces(ws_lhs, ws_rhs):
        RebinToWorkspace(WorkspaceToRebin=ws_rhs, WorkspaceToMatch=ws_lhs,
                         OutputWorkspace=ws_rhs, PreserveEvents=False)
        Divide(LHSWorkspace=ws_lhs, RHSWorkspace=ws_rhs, OutputWorkspace=ws_lhs)
        ReplaceSpecialValues(InputWorkspace=ws_lhs, OutputWorkspace=ws_lhs, 
            NaNValue=0, InfinityValue=0, BigNumberThreshold=1e15, SmallNumberThreshold=1e-15)        

    @staticmethod
    def make_UB_consistent(ws_ref, ws):
        # compare U matrix to perform TransformHKL to preserve indexing
        U_ref = ADS.retrieve(ws_ref).sample().getOrientedLattice().getU()
        U = ADS.retrieve(ws).sample().getOrientedLattice().getU()
        # find transform required  ( U_ref = U T^-1) - see TransformHKL docs for details
        transform = np.linalg.inv(getSignMaxAbsValInCol(np.matmul(np.linalg.inv(U), U_ref)))
        TransformHKL(PeaksWorkspace=ws, HKLTransform=transform, FindError=False)

    @staticmethod
    def find_sx_peaks(wsname, bg, out_pk_wsname='peaks'):
        ws = ADS.retrieve(wsname)
        # get unit to convert back to after peaks found
        xunit = None
        if not ws.getXDimension().name == 'Time-of-flight':
            xunit = ws.getXDimension().name.replace('-',"")
            ConvertUnits(InputWorkspace=wsname, OutputWorkspace=wsname, Target='TOF')  # FindSXPeaks requries TOF
        # extract y data (to use to determine to detemrine threshold)
        FindSXPeaks(InputWorkspace=wsname, PeakFindingStrategy='AllPeaks', 
            AbsoluteBackground=bg, ResolutionStrategy='AbsoluteResolution', 
            XResolution=400, PhiResolution=3, TwoThetaResolution=2,
            OutputWorkspace=out_pk_wsname)  # StartWorkspaceIndex=spec_min, EndWorkspaceIndex=spec_max,
        WISHReduction.remove_peaks_on_edge(out_pk_wsname)
        if xunit:
            ConvertUnits(InputWorkspace=wsname, OutputWorkspace=wsname, Target=xunit)
        return out_pk_wsname
        
    @staticmethod
    def remove_peaks_on_edge(pk_wsname, nedge_tube = 3, nedge_pix = 20):
        peaks = ADS.retrieve(pk_wsname)
        # filter peaks on edge of tubes
        row = np.array(peaks.column("Row"))
        iremove = np.where(np.logical_or(row < nedge_pix, row > 512 - nedge_pix))[0]
        DeleteTableRows(TableWorkspace=peaks, Rows=iremove)
        # filter out peaks on tubes at edge of detector banks 
        col = np.array(peaks.column("Col"))
        bank = np.array([int(name[-2:]) for name in peaks.column("BankName")])
        iremove = np.where(np.logical_and(col<nedge_tube, bank==1))[0]
        DeleteTableRows(TableWorkspace=peaks, Rows=iremove)
        iremove = np.where(np.logical_and(col>152 - nedge_tube, bank==5))[0]
        DeleteTableRows(TableWorkspace=peaks, Rows=iremove)
        iremove = np.where(np.logical_and(col<nedge_tube, bank==10))[0]
        DeleteTableRows(TableWorkspace=peaks, Rows=iremove)
        iremove = np.where(np.logical_and(col>152 - nedge_tube, bank==6))[0]
        DeleteTableRows(TableWorkspace=peaks, Rows=iremove)
       

sphere_GGG = '''<sphere id="sphere">
                <centre x="0.0"  y="0.0" z="0.0" />
                <radius val="0.0009"/>
                </sphere>'''  # sphere radius 0.9mm
sample_dict = {'Geometry': {'Shape': 'CSG', 'Value': sphere_GGG}, 
    'Material': {'ChemicalFormula': 'Ca3-Ga2-Ge3-O12', 'SampleNumberDensity': 0.00433}}
    
wish = WISHReduction(vanadium_runno="43526")
wish.process_vanadium()
for run in range(42728, 42735):
    print('RUN = ', run)
    wish.process_data([run], sample_dict=sample_dict, find_peaks=False)
    wish.load_isaw_ub(['/babylon/Public/RWaite/WISH000' + str(run) + '_RW.mat'],
        runs = [str(run)])
    wish.predict_peaks(MinDSpacing=0.5, WavelengthMin=0.75, 
        ReflectionCondition='Body centred')
    wish.integrate_data_skew(PEAK_TYPE.PREDICTED, tol=0)
    wish.delete_run_data(run)
wish.save_integrated_peaks(save_dir='/babylon/Public/RWaite/', save_format='Jana')


# intPeaksMDArgs = {'peakRadius': 0.15, 'backgroundInnerRadius': 0.155, 'backgroundOuterRadius': 0.195,
#                   'ellipsoid': True, 'fixQAxis': True, 'fixMajorAxisLength': True, 'useCentroid': True}
# # intPeaksMDArgs = {'peakRadius': 0.3, 'backgroundInnerRadius': 0, 'backgroundOuterRadius': 0,
# #                   'ellipsoid': False} # no bg shell most like SXD2001 shoebox
# # wish.integrate_data(PEAK_TYPE.PREDICTED, tol=0.1, **intPeaksMDArgs)
# rho = 	1.424 # g/cm^3  [1,2]
# Mf = 89.09 # g/mol
# NA = 6.022E23 # atoms/mol
# n = rho*NA/Mf # g/cm^3
# # 1 Ang = 1E-10m = 1E-8cm -> 1E-8 cm/Ang
# n = n*((1E-8)**3) # 0.0096255 atoms/Ang^3
#  
# vol = 5.8033*5.9749*12.3259
# # 4 formula units per unit cell [2]
# n = 4/vol # 0.0093591 atoms/Ang^3
# 
# # [1] https://pubs.acs.org/doi/pdf/10.1021/ja01853a020
# # [2] https://pubchem.ncbi.nlm.nih.gov/compound/Alanine

print(8/(12.27**3))