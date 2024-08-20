# =======================================================
#
#      MONO VANADIUM REDUCTION SCRIPT
#      JRS 21/6/23
#
#      Reads mono vanadium run, corrects for backgound,
#      normalises to a White vanadium run, and finds the
#      absolute normalisation factor for each Ei.  These
#      are saved to an assci file, 'MV_<runno>.txt'
#
#      Optionally takes account of the vanadium DW factor.
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import numpy as np
import time

t = time.time()         #start the clock

#!begin_params
#============================User Input===========================
# Assumes that the thickness of the vanadium and the sample are the same
monovan       = 59138                       #mono van run
monovan_bg    = None                        #background for mono van
Ei_list       = [55]                        #incident energies - from PyChop
monovan_trans = [1]                         #transmission for each Ei
monovan_temp  = 300                         #vanadium temperature (K) for Debye-Waller correction
mask          = 'MASK_FILE_XML'             #mask
wv_file       = 'WV_57083.txt'              #whitevan integrals file
#=================================================================

#======================Local Contact Input========================
# MARI monitors:    m2spec=2,     m3spec=3     - fixei = False
# MERLIN monitors:  m2spec=69636, m3spec=69640 - fixei = False
# MAPS monitors:    m2spec=41475, m3spec=41476 - fixei = False
# LET monitors:     m2spec=98310, m3spec= None - fixei = True

inst = 'INSTRUMENT_NAME'
cycle = 'CYCLE_ID'
fixei   = False                 #true for LET since no m3
powdermap = 'RINGS_MAP_XML'     #rings mapping file
monovan_range = 0.1             #integration range (+/- fraction of Ei)
min_theta = 5.                  #minimum theta to avoid low Q detectors (SANS)
idebug  = False                 #keep workspaces
#=================================================================
#!end_params

config['default.instrument'] = inst
if inst == 'MARI':
    source = 'Moderator'
    m2spec = 2                  # specID of monitor2 (pre-sample)
    m3spec = 3                  # specID of monitor3 (post-sample)
elif inst == 'MERLIN':
    source = 'undulator'
    m2spec = 69636              # specID of monitor2 (pre-sample)
    m3spec = 69640              # specID of monitor3 (post-sample)
elif inst == 'MAPS':
    source = 'undulator'
    m2spec = 36867              # specID of monitor2 (pre-sample)
    m3spec = 36868              # specID of monitor3 (post-sample)
elif inst == 'LET':
    source = 'undulator'
    m2spec = 98310              # specID of monitor2 (pre-sample)
    m3spec = None               # specID of monitor3 (post-sample)
else:
    raise RuntimeError(f'Unrecognised instrument: {inst}')

if cycle.startswith('20'):
    cycle = cycle[2:]
datadir = '/archive/NDX'+inst+'/instrument/data/cycle_'+cycle+'/'
mapdir  = '/usr/local/mprogs/InstrumentFiles/'+inst.swapcase()+'/'
config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)
ws_list = ADS.getObjectNames()

print("\n======= "+inst+" MonoVan reduction =======")

#=============================load hard mask=============================================
if mask is None:
    print(inst+": WARNING - No hard mask!  Bad detectors wont be removed")
if mask not in ws_list:
    print(inst+": Loading hard mask - %s" % mask)
    LoadMask(inst,mask,OutputWorkspace=mask)
else:
    print(inst+": Using previously loaded hard mask - %s" % mask)

# load whitevan file
if wv_file not in ws_list:
    print(inst+": Loading white vanadium - %s" % wv_file)
    LoadAscii(Filename=wv_file,OutputWorkspace=wv_file)
    ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN',OutputWorkspace=wv_file)
else:
    print(inst+": Using previously loaded white vanadium - %s" % wv_file)

# =====load the mono vanadium and find absolute normalisation factor for each Ei  =======
print(inst + " MONO_VAN: Loading mono vanadium run #%i" % monovan)
mv = Load(str(monovan),LoadMonitors=True)
mv = NormaliseByCurrent('mv')
if monovan_bg is not None:
    if iprint:
        print("... subtracting mono vanadium background run #%i" % monovan_bg)
    mv_bg = Load(str(monovan_bg),LoadMonitors=False)
    mv_bg = NormaliseByCurrent('mv_bg')

mv_fac=[]; mv_error=[]
for ienergy in range(len(Ei_list)):
    tr = monovan_trans[ienergy]
    Ei = Ei_list[ienergy]
    print("\n"+inst+" MONO_VAN: calibration factor for Ei=%g meV" % Ei)
    print("... summing  |\u0394E| < %.2f meV" %  (monovan_range*Ei))

    mv_corrected = Scale('mv',1/tr,"Multiply")
    if (monovan_bg is not None):
        if iprint:
            print("... transmission factor = %.2f" % tr)
        mv_corrected = mv/tr - mv_bg
    mv_norm = Divide('mv_corrected',wv_file)
    MaskDetectors(mv_norm,MaskedWorkspace=mask)

    mv_monitors = mtd['mv_monitors']
    spectra = mv_monitors.getSpectrumNumbers()
    index = spectra.index(m2spec)
    m2pos = mv.detectorInfo().position(index)[2]

# this section shifts the time-of-flight such that the monitor2 peak
# in the current monitor workspace (post monochromator) is at t=0 and L=0
# note that the offest is not the predicted value due to energy dependence of the source position

    Ei_orig = Ei
    if m3spec is not None and not fixei:
        (Ei,mon2_peak,_,_) = GetEi(mv_monitors,Monitor1Spec=m2spec,Monitor2Spec=m3spec,EnergyEstimate=Ei,FixEi=fixei)
        print("... refined Ei=%.2f meV" % Ei)
    else:
        (Ei,mon2_peak,_,_) = GetEi(mv_monitors,Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=fixei)
    if np.isnan(Ei):
        Ei = Ei_orig
    print("... m2 tof=%.2f mus, m2 pos=%.2f m" % (mon2_peak,m2pos))

    mv_norm = ScaleX(mv_norm, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
    MoveInstrumentComponent(mv_norm, ComponentName=source, Z=m2pos, RelativePosition=False)

# convert to energy transfer and correct for energy dependent detector efficiency
    mv_out = ConvertUnits(mv_norm,'DeltaE',EMode='Direct',EFixed=Ei)
    mv_out = Rebin(mv_out, [-monovan_range*Ei,Ei/100.,monovan_range*Ei], PreserveEvents=False)
    mv_out = DetectorEfficiencyCor(mv_out, IncidentEnergy=Ei)

#rings group and Q-axis
    mv_out = GroupDetectors(mv_out, MapFile=powdermap, Behaviour='Average')
    mv_out = RemoveMaskedSpectra(mv_out)
    mv_out = ConvertSpectrumAxis(mv_out,'ElasticQ',EFixed=Ei)
#integrate under the elastic line between the specified limits
    mv_out = Integration(mv_out, -monovan_range*Ei, monovan_range*Ei)
    mv_out = Transpose(mv_out)
    Q = mv_out.readX(0)
    if idebug:
        print(Q)

#=============================Debye-Waller factor correction==================================

#  See F. Rieutord INX Manual, ILL Technical Report No. ILL90RI17T, 1990.
#  2W=2*[hbar^2/2M]*[2k/(kT_Debye_vana)^2]*[T*300/300]*Q^2=alpha*Q^2
#  with alpha = 0.0067 * T / 300
#  <=> alpha = 2* [hbar^2/2M]*[2k/(kT_Debye_vana)^2]
#  T_Debye_vana = 355 K
#  M = 50.9415*1.67265E-27 kg

    if monovan_temp is not None:
        DW_fac = np.exp(0.0067*monovan_temp*(Q**2)/300.)
        DW = CreateWorkspace(DataX=Q,DataY=DW_fac,UnitX='MomentumTransfer')
        mv_out = Multiply(mv_out,DW)

    Qmax = np.min([Q[-3],2.5])                               # max Q is 2.6Å^-1 or maxQ (avoiding edge) whichever is least
    theta_min = min_theta*np.pi/180.                         # minimum theta angle to avoid low angle detectors
    Qmin = 4*np.pi*np.sin(theta_min)*np.sqrt(Ei)/9.045       # don't use low Q detectors
    print("... average over %.2f < |Q| < %.2f \u212B^-1" % (Qmin, Qmax))
    mv_out = CropWorkspace(mv_out,Qmin,Qmax)                 # extract data not contaminated by low Q, or Al Bragg peaks
    CloneWorkspace(mv_out,OutputWorkspace='elastic_monovan_'+str(Ei))

# =================================take averages of data=====================================
    signal = mv_out.dataY(0)
    mv_fac_unweighted = np.sum(signal) / len(signal)
    errors = mv_out.dataE(0)
    weights = 1/(errors**2)
    mv_fac_weighted = np.sum(np.multiply(weights,signal)) / np.sum(weights)
    print("... unweighted average: %.4f" % mv_fac_unweighted)
    print("... weighted average: %.4f" % mv_fac_weighted)

    mv_fac.append(mv_fac_weighted)

#get standard errors
    mv_error.append(np.sqrt(np.sum(errors**2))/len(signal))

# ================================save out calibration factors==============================
monovan_factors = CreateWorkspace(DataX=Ei_list,DataY=mv_fac,DataE=mv_error,UnitX='Energy',YUnitLabel='MonoVan factor')
ofile = 'MV_'+str(monovan)+'.txt'
print("\n" + inst + " MONO_VAN: saving calibration factors to %s" % ofile)
SaveAscii(monovan_factors,ofile)

# ==================================================================

print(inst + " MONO_VAN: Reduction complete in %.1f seconds\n" % (time.time() - t))

# ==================================================================

# cleanup
if (not idebug):
    DeleteWorkspaces(['DW'])
    ws_list = ADS.getObjectNames()
    nx_list = [ss for ss in ws_list if 'mv' in ss]
    for ss in nx_list:
        ADS.remove(ss)
