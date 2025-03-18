# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from plugins.algorithms.peakdata_utils import InstrumentArrayConverter, PeakData
from scipy.optimize import minimize, curve_fit
from mantid.api import FunctionFactory

def back_to_back_exp(x, I, A, B, X0, S, bg):
    """IFunc = BackToBackExponential()"""
    if np.isnan([I, A, B, X0, S, bg]).any():
        return np.ones(x.shape)*1e-4
    func = BackToBackExponential(I=I, A=A, B=B, X0=X0, S=S)
    out = func(x)
    out[~np.isfinite(out)] = 1e-4 # don't return NaN (often far from peak)
    out = out + bg
    return out

   
def estimate_intensity_and_background_of_spectra(x, y):
    ibg, ipeak = PeakData.find_bg_pts_seed_skew(np.squeeze(y))
    bg = y[ibg].mean()
    intens = np.sum(y[ipeak]-bg)*np.diff(x).mean()
    return intens, bg

runs = range(32863, 32895)
ws_list = []
pks_list = []
for run in runs:
    ws_list.append(f"SXD{run}")
    Load(Filename=f'{ws_list[-1]}.raw', OutputWorkspace=ws_list[-1])
    FindSXPeaks(InputWorkspace=ws_list[-1], PeakFindingStrategy='AllPeaks', 
                AbsoluteBackground=25, ResolutionStrategy='AbsoluteResolution', 
                XResolution=250, PhiResolution=2, TwoThetaResolution=2, 
                OutputWorkspace=ws_list[-1] + '_peaks')
    ConvertUnits(InputWorkspace=ws_list[-1], OutputWorkspace=ws_list[-1], Target='dSpacing')

nrows, ncols = 15, 15
dtof_over_tof = 0.12
bank_groups = ([1,6,10], [2,5,11,7,9], [3,4], [8])

wl_banks = {str(banks): [] for banks in bank_groups}
tth_banks = {str(banks): [] for banks in bank_groups}
ds_banks = {str(banks): [] for banks in bank_groups}
tofs_banks = {str(banks): [] for banks in bank_groups}
sigmas_banks = {str(banks): [] for banks in bank_groups}
sigmas_err_banks = {str(banks): [] for banks in bank_groups}
betas_banks = {str(banks): [] for banks in bank_groups}
betas_err_banks = {str(banks): [] for banks in bank_groups}

for wsname in ws_list:
    ws = mtd[wsname]
    array_converter = InstrumentArrayConverter(ws)
    peaks = mtd[wsname + '_peaks']
    for banks in bank_groups:
        # filter by bank
        peaks_bank = FilterPeaks(InputWorkspace=peaks, 
                                 Criterion="=", BankName=f"bank{banks[0]}")
        for bank in banks[1:]:
            peaks_bank_tmp = FilterPeaks(InputWorkspace=peaks, 
                                         Criterion="=", BankName=f"bank{bank}")
            peaks_bank = CombinePeaksWorkspaces(peaks_bank, peaks_bank_tmp, CombineMatchingPeaks=False,
                                                OutputWorkspace=peaks_bank.name())
            DeleteWorkspace(peaks_bank_tmp)
        # fit each peak indep.
        for ipk in range(peaks_bank.getNumberPeaks()):
            pk = peaks_bank.getPeak(ipk)
            bank_name = peaks_bank.column('BankName')[ipk]
            peak_data = array_converter.get_peak_data(pk, pk.getDetectorID(), bank_name, nrows, ncols, 
                                                      nrows_edge=1, ncols_edge=1)
            peak_data.peak_mask = ~peak_data.det_edges
            # print(peak_data._focus_detids(peak_data.detids.flatten()))
            if not np.any(peak_data.det_edges):
                # focus peak
                peak_data.peak_mask = np.ones(peak_data.detids.shape, dtype=bool)
                peak_data.peak_mask[0,:]=False
                peak_data.peak_mask[-1,:]=False
                peak_data.peak_mask[:,0]=False
                peak_data.peak_mask[:,-1]=False
                peak_data.non_bg_mask = peak_data.peak_mask
                peak_data.focus_data_in_detector_mask()
                # extract data, crop and convert to TOF
                tof = pk.getTOF()
                dcen = pk.getDSpacing()
                dx = dcen*dtof_over_tof
                x = peak_data.xpk
                ikeep = np.logical_and(x > dcen - dx/2, x < dcen + dx)
                x = x[ikeep] * (pk.getTOF()/dcen)
                y = peak_data.ypk[ikeep]
                e = np.sqrt(peak_data.epk_sq[ikeep])
                if y.sum()/np.sqrt(np.sum(e**2)) > 3:
                    # get initial params
                    intens, _ = estimate_intensity_and_background_of_spectra(x, y)
                    pinit = [max(intens, np.sum(y)*(x[1]-x[0])), 2, 0.03, tof, 10, 0]
                    # fit peak
                    inonzero = e > 0
                    #              I                    A        B         X0              S          bg
                    pbounds = [(0.0, np.inf), (2-1e-6,2+1e-6), (0, 2), (x[0],x[-1]), (0.0, 200), (-np.inf, np.inf)]
                    popt, pcov = curve_fit(back_to_back_exp, x[inonzero], y[inonzero], pinit, 
                                            sigma=e[inonzero], bounds = list(zip(*pbounds)),
                                            absolute_sigma=True, maxfev=5000)
                    perr = np.sqrt(np.diag(pcov))
                    # add chi-squared criteria
                    chisq = np.sum(((y - back_to_back_exp(x, *popt))/e)**2)
                    chisq_red = chisq/(len(y)-len(popt))
                    if popt[0]/perr[0] > 20 and abs(popt[2]/perr[2]) > 1 and abs(popt[4]/perr[4]) > 1 and chisq_red < 15:
                        tofs_banks[str(banks)].append(tof)
                        ds_banks[str(banks)].append(dcen)
                        sigmas_banks[str(banks)].append(popt[4])
                        sigmas_err_banks[str(banks)].append(perr[4])
                        betas_banks[str(banks)].append(popt[2])
                        betas_err_banks[str(banks)].append(perr[2])
                        wl_banks[str(banks)].append(pk.getWavelength())
                        tth_banks[str(banks)].append(pk.getScattering())
                        # # plot
                        # fig, ax = plt.subplots()
                        # ax.errorbar(x, y, yerr=e, capsize=2, marker="o", color="k", label='data')
                        # ax.plot(x, back_to_back_exp(x, *pinit), '--b', label='initial guess')
                        # ax.plot(x, back_to_back_exp(x, *popt), '-r', label='fit')
                        # ax.set_title(f"ipk = {ipk}\n I/sig = {np.round(popt[0]/perrs[0],1)}\n{wsname}, {bank_name}, chisq_red = {np.round(chisq_red,2)}")
                        # ax.axvline(pk.getTOF(), color=3*[0.5], ls='--')
                        # ax.set_xlabel('TOF (mus)')
                        # ax.set_ylabel('Intensity (a.u.)')
                        # ax.legend()
                        # fig.show()

####
# sig & beta for each bank (group)
####

linear_func = lambda x, m, c: m*x + c

for ibank in range(len(bank_groups)):
    ds = np.array(ds_banks[str(bank_groups[ibank])])
    tofs = np.array(tofs_banks[str(bank_groups[ibank])])
    sigmas = np.array(sigmas_banks[str(bank_groups[ibank])])
    sigmas_err = np.array(sigmas_err_banks[str(bank_groups[ibank])])
    betas = np.array(betas_banks[str(bank_groups[ibank])])
    betas_err = np.array(betas_err_banks[str(bank_groups[ibank])])
    # fit back to back exponential coeficients
    
    # sigma
    s_coef, s_cov = curve_fit(linear_func, tofs, sigmas, [np.mean(sigmas/tofs), 0], bounds = ([0, -1e-10], [np.inf, 1e-10]),
                              sigma=sigmas_err, absolute_sigma=False)
    # B
    b_coef, b_cov = curve_fit(linear_func, 1/(ds**4), betas, [np.mean((betas-min(betas))*(ds**4)), min(betas)], 
                              sigma=betas_err, absolute_sigma=False)
    # plot

    fig, ax = plt.subplots(1,2)
    ax[0].errorbar(tofs, sigmas, yerr=sigmas_err, ls='', marker='o', color='k', capsize=2, alpha=0.5)
    ax[0].plot(tofs, linear_func(tofs, *s_coef), '-r')
    ax[0].set_xlabel('TOF (mus)')
    ax[0].set_ylabel('sigma (mus)')
    ax[0].set_title(f"sig_1 = {np.round(s_coef[0],4)}")
    ax[1].errorbar(1/(ds**4), betas, yerr=betas_err, ls='', marker='o', color='k', capsize=2, alpha=0.5)
    ax[1].plot(1/(ds**4), linear_func(1/(ds**4), *b_coef), '-r')
    ax[1].set_title(f"beta_0 = {np.round(b_coef[1],4)}\nbeta_1 = {np.round(b_coef[0],4)}")
    ax[1].set_xlabel('1/d^4 (Ang^-4)')
    ax[1].set_ylabel('beta (1/mus)')
    fig.suptitle(f"banks: {str(bank_groups[ibank])}")
    fig.tight_layout()
    fig.show()         
            
####
# all aigmas
####

tth_all = []
tofs_all = []
sigmas_all = []

fig, ax = plt.subplots(1,2)           
for ibank in range(len(bank_groups)):
    tth = np.array(tth_banks[str(bank_groups[ibank])])
    wl = np.array(wl_banks[str(bank_groups[ibank])])
    tofs = np.array(tofs_banks[str(bank_groups[ibank])])
    sigmas = np.array(sigmas_banks[str(bank_groups[ibank])])
    ax[0].scatter(np.degrees(tth/2), sigmas/tofs, c=wl, vmin=0.5, vmax=3)
    line = ax[1].scatter((1/np.tan(tth/2))**2, (sigmas/tofs)**2, c=wl, vmin=0.5, vmax=3)
    tth_all.extend(tth.tolist())
    tofs_all.extend(tofs.tolist())
    sigmas_all.extend(sigmas.tolist())
# fit resolution to ,linear
tth_all = np.array(tth_all)
tofs_all = np.array(tofs_all)
sigmas_all = np.array(sigmas_all)
xvals = (1/np.tan(tth_all/2))**2
yvals = (sigmas_all/tofs_all)**2
res_coef, res_cov = curve_fit(linear_func, xvals, yvals, 
                              [np.mean(yvals/xvals), 0])
fobj = lambda x, m, c: np.sqrt(m * (1 / np.tan(x)**2) + c)
res_coef, res_cov = curve_fit(fobj, tth_all/2, sigmas_all/tofs_all, res_coef)
xfits = np.radians(np.linspace(11,80))
ax[1].plot((1/np.tan(xfits))**2, linear_func((1/np.tan(xfits))**2, *res_coef), '-r')               
ax[0].plot(np.degrees(xfits), np.sqrt(res_coef[0] * (1 / np.tan(xfits)**2) + res_coef[1]), '-r')                             
ax[0].set_ylabel('sigma/TOF')
ax[0].set_xlabel('tth/2')
ax[1].set_ylabel('(sigma/TOF)^2')
ax[1].set_xlabel('cot(tth/2)^2')
fig.suptitle(f"dth = {np.round(np.sqrt(res_coef[0]),4)}\n dT0/T0 = {np.round(np.sqrt(res_coef[1]),4)}")
fig.tight_layout()
fig.subplots_adjust(bottom=0.35)
cbar_ax = fig.add_axes([0.2,0.15,0.6,0.05])
fig.colorbar(line, orientation="horizontal", label='Wavelength', cax=cbar_ax)
fig.show()
   
print(np.min(sigmas_all/tofs_all))   
####
# all betas
####

wl_all = []
betas_all = []
betas_err_all = []

fig, ax = plt.subplots()           
for ibank in range(len(bank_groups)):
    wl = np.array(wl_banks[str(bank_groups[ibank])])
    tth = np.array(tth_banks[str(bank_groups[ibank])])
    betas = np.array(betas_banks[str(bank_groups[ibank])])
    betas_err = np.array(betas_banks[str(bank_groups[ibank])])
    line = ax.scatter(1/(wl**4), betas, c=np.degrees(tth/2), vmin=10, vmax=80)
    wl_all.extend(wl)
    betas_all.extend(betas)
    betas_err_all.extend(betas_err)
wl_all = np.array(wl_all)
betas_all = np.array(betas_all)
betas_err_all = np.array(betas_err_all)
res_coef, res_cov = curve_fit(linear_func, 1/(wl_all**4), betas_all, 
                              [np.mean(betas_all*wl_all**4), 0], sigma=betas_err_all)
xs = np.array([0, ax.get_xlim()[-1]])
ax.plot(xs, linear_func(xs, *res_coef), '-r')   
ax.set_xlabel('1/lambda^4')
ax.set_ylabel('beta')
fig.colorbar(line, label='tth/2')
fig.suptitle(f"beta_0 = {np.round(res_coef[1],4)}\n beta_1 = {np.round(res_coef[0],4)}")
fig.show()

#####
# Full extent of peak
#####
frac = 0.01

tth_all = []
tofs_all = []
dtof_all = []

fig, ax = plt.subplots(1,2)           
for ibank in range(len(bank_groups)):
    tth = np.array(tth_banks[str(bank_groups[ibank])])
    wl = np.array(wl_banks[str(bank_groups[ibank])])
    tofs = np.array(tofs_banks[str(bank_groups[ibank])])
    sigmas = np.array(sigmas_banks[str(bank_groups[ibank])])
    betas = np.array(betas_banks[str(bank_groups[ibank])])
    # find extent of peak
    dtof = np.zeros(tofs.shape)
    for ipk in range(len(dtof)):
        # evaluate peak function
        p = [1, 2, betas[ipk], tofs[ipk], sigmas[ipk], 0]
        dx = 1.2*tofs[ipk]*dtof_over_tof
        nx = 501
        x = np.linspace(tofs[ipk] - dx, tofs[ipk]+dx,nx)
        yfit = back_to_back_exp(x, *p)
        # find xlo and xhi based on frac
        ymax = yfit.max()
        ixlo = np.argmin(abs(yfit[0:nx//2] - frac*ymax))
        ixhi = nx//2 + np.argmin(abs(yfit[nx//2:] - frac*ymax))
        dtof[ipk] = x[ixhi] - x[ixlo]
        # fig,ax = plt.subplots()
        # ax.plot(x, yfit, 'ok')
        # ax.axhline(frac*ymax)
        # ax.axvline(x[ixlo])
        # ax.axvline(x[ixhi])
        # fig.show()
    ax[0].scatter(np.degrees(tth/2), dtof/tofs, c=wl, vmin=0.5, vmax=3)
    line = ax[1].scatter((1/wl**1)*(1/np.tan(tth/2))**2, (dtof/tofs)**2, c=wl, vmin=0.5, vmax=3)
    tth_all.extend(tth.tolist())
    tofs_all.extend(tofs.tolist())
    dtof_all.extend(dtof.tolist())
# fit resolution to ,linear
tth_all = np.array(tth_all)
tofs_all = np.array(tofs_all)
dtof_all = np.array(sigmas_all)
xvals = (1/np.tan(tth_all/2))**2
yvals = (dtof_all/tofs_all)**2
res_coef, res_cov = curve_fit(linear_func, xvals, yvals, 
                              [np.mean(yvals/xvals), 0])
fobj = lambda x, m, c: np.sqrt(m * (1 / np.tan(x)**2) + c)
res_coef, res_cov = curve_fit(fobj, tth_all/2, dtof_all/tofs_all, res_coef)
xfits = np.radians(np.linspace(11,80))
ax[1].plot((1/np.tan(xfits))**2, linear_func((1/np.tan(xfits))**2, *res_coef), '-r')               
ax[0].plot(np.degrees(xfits), np.sqrt(res_coef[0] * (1 / np.tan(xfits)**2) + res_coef[1]), '-r')                             
ax[0].set_ylabel('dTOF/TOF')
ax[0].set_xlabel('tth/2')
ax[1].set_ylabel('(sigma/TOF)^2')
ax[1].set_xlabel('cot(tth/2)^2')
fig.suptitle(f"dth = {np.round(np.sqrt(res_coef[0]),4)}\n dT0/T0 = {np.round(np.sqrt(res_coef[1]),4)}")
fig.tight_layout()
fig.subplots_adjust(bottom=0.35)
cbar_ax = fig.add_axes([0.2,0.15,0.6,0.05])
fig.colorbar(line, orientation="horizontal", label='Wavelength', cax=cbar_ax)
fig.show()


# print out nominal sig_0 for all banks
dth_sq = 1.34892856e-05
dt0_over_t0_sq = 2.43489688e-06
fsig = lambda th: np.sqrt(dth_sq * (1 / np.tan(th)**2) + dt0_over_t0_sq)
 
print(res_coef)
for ipk, pk in enumerate(mtd['SingleCrystalPeakTable']):
    th = pk.getScattering()/2
    print(mtd['SingleCrystalPeakTable'].column('BankName')[ipk], np.degrees(th), np.round(fsig(th),4))
        