# =======================================================
#
#      MONO VANADIUM REDUCTION SCRIPT
#      JRS 19/2/22
#
#      Reads mono vanadium run, corrects for backgound,
#      normalises to a White vanadium run, and finds the
#      absolute normalisation factor for each Ei.  These
#      are save to an assci file, 'MV_<runno>.txt'
#
#      Optionally takes account of the vanadium DW factor.
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import time

t = time.time()         #start the clock

#============================User Input===========================
# Assumes that the thickness of the vanadium and the sample are the same
monovan       = 72316           #mono van, 240/80
monovan_bg    = None
Ei_list       = [12.14,3.7,1.77] #240/80 Eis
monovan_trans = [1,1,1]          #transmission for each Ei
monovan_range = [-0.1,0.1]       #integration range (fraction of Ei)
monovan_temp  = 300              #vanadium temperature (K) for DW correction
mask          = 'LET_mask_212_base.xml'
wv_file       = 'WV_76865.txt'
#=================================================================

#======================Local Contact Input========================
datadir = '/archive/cycle_21_1/NDXLET/'
mapdir  = '/usr/local/mprogs/InstrumentFiles/let/'
maskdir = './'
m2spec  = 98310                 #specID of monitor2 (post monochromator)
m2pos   = 1.333                 #position of monitor2 wrt the sample position
iprint  = True                  #turn on output
idebug  = False                  #keep workspaces
powdermap = 'LET_rings_153.map'
#=================================================================

config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)

#load hard mask and create norm workspace
if iprint:
    print("MONO_VAN: Loading hard mask %s" % mask)
LoadMask('LET',maskdir+mask,OutputWorkspace=mask)

# load white van file
if iprint:
    print("MONO_VAN: Loading white vanadium integrals %s" % wv_file)
wv_norm = LoadAscii(Filename=wv_file)
wv_norm = ReplaceSpecialValues(wv_norm, SmallNumberThreshold=1e-20, SmallNumberValue='NaN')

# load the mono vanadium and find absolute normalisation factor for each Ei  
if iprint:
    print("MONO_VAN: Loading mono vanadium run #%i" % monovan)
mv = Load(str(monovan),LoadMonitors=True)
mv = NormaliseByCurrent('mv')
if (monovan_bg is not None):
    if iprint:
        print("... subtracting mono vanadium background run #%i" % monovan_bg)
    mv_bg = Load(str(monovan_bg),LoadMonitors=False)
    mv_bg = NormaliseByCurrent('mv_bg')

mv_fac=[]
for ienergy in range(len(Ei_list)):
    tr = monovan_trans[ienergy]
    Ei = Ei_list[ienergy]
    
    if iprint:
        print("MONO_VAN: calibration factor for Ei=%.2f" % Ei)
        print("... summing  |deltaE| < %.2f meV" %  (monovan_range[1]*Ei))
    
    mv_corrected = Scale('mv',1/tr,"Multiply")
    if (monovan_bg is not None):
        print("... transmission factor = %.2f" % tr)
        mv_corrected = mv/tr - mv_bg
    mv_norm = Divide('mv_corrected','wv_norm')
    MaskDetectors(mv_norm,MaskedWorkspace=mask)
        
# this section shifts the time-of-flight such that the monitor2 peak 
# in the current monitor workspace (post monochromator) is at t=0 and L=0
# note that the offest is not the predicted value due to energy dependence of the source position        
    (_,mon2_peak,_,_) = GetEi('mv_monitors',Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=True)
    mv_norm = ScaleX(mv_norm, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
    MoveInstrumentComponent(mv_norm, ComponentName='undulator', Z=-m2pos, RelativePosition=False)

# convert to energy transfer and integrate under the elastic line
    mv_out = ConvertUnits(mv_norm,'DeltaE',EMode='Direct',EFixed=Ei)
    mv_out = Rebin(mv_out, [monovan_range[0]*Ei,100,monovan_range[1]*Ei], PreserveEvents=False)
    mv_out = DetectorEfficiencyCor(mv_out, IncidentEnergy=Ei)
    MaskDetectors(mv_out, MaskedWorkspace=mask)
#rings group and Q-axis
    mv_out = GroupDetectors(mv_out, MapFile=powdermap, Behaviour='Average')
    mv_out = ConvertSpectrumAxis(mv_out,'ElasticQ',EFixed=Ei)
    mv_out = Transpose(mv_out)
    Q = mv_out.readX(0)

#Debye-Waller factor correction

#  See F. Rieutord INX Manual, ILL Technical Report No. ILL90RI17T, 1990.
#  2W=2*[hbar^2/2M]*[2k/(kT_Debye_vana)^2]*[T*300/300]*Q^2=alpha*Q^2
#  with alpha = 0.0067 * T / 300
#  <=> alpha = 2* [hbar^2/2M]*[2k/(kT_Debye_vana)^2]
#  T_Debye_vana = 355 K
#  M = 50.9415*1.67265E-27 kg

    if (monovan_temp is not None):
        DW_fac = np.exp(0.0067*monovan_temp*(Q**2)/300.)
        DW = CreateWorkspace(DataX=Q,DataY=DW_fac,UnitX='MomentumTransfer')
        mv_out = Multiply(mv_out,DW)

    Qmax = np.min([Q[-1]-0.3,2.5])   # max Q is 2.5Å^-1 or  [maxQ - 0.3Å^-1] whichever is least
    Qmin = Q[0]+0.3                  # don't use low Q detectors
    mv_out = CropWorkspace(mv_out,Qmin,Qmax)      # extract data not contaminated by low Q, or Al Bragg peaks
    CloneWorkspace(mv_out,OutputWorkspace='mv_out_'+str(Ei))
    
# take averages of data
    signal = mv_out.dataY(0)
    mv_fac_unweighted = np.sum(signal) / len(signal)
    errors = mv_out.dataE(0)
    weights = 1/(errors**2)
    mv_fac_weighted = np.sum(np.multiply(weights,signal)) / np.sum(weights)
    
    print("... unweighted average: %.4f" % mv_fac_unweighted)
    print("... weighted average: %.4f" % mv_fac_weighted)
 
    mv_fac.append(mv_fac_weighted)

# save out calibration factors
mv_ws = CreateWorkspace(DataX=Ei_list,DataY=mv_fac)
ofile = 'MV_'+str(monovan)+'.txt'
if iprint:
    print("MONO_VAN: saving calibration factors to %s" % ofile)
SaveAscii(mv_ws,ofile)
if (not idebug):
    DeleteWorkspaces(['mv_ws','mv_out','DW','mv','mv_corrected','mv_monitors','mv_norm'])

# ==================================================================
if iprint:
    print("MONO_VAN: Reduction complete in %.1f seconds\n" % (time.time() - t))