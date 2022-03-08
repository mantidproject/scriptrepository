# =======================================================
#
#      SAMPLE REDUCTION SCRIPT
#      JRS 19/2/22
#
#      Reads sample run, corrects for backgound,
#      normalises to a White vanadium run, and  the
#      absolute normalisation factor for each Ei (if MV
#      file is specified).  Outputs .nxspe file
#
#========================================================

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
import matplotlib.pyplot as plt
import numpy as np
import time
from importlib import reload

t = time.time()         #start the clock

#=======================User Input=======================
powder         = True
sumruns        = True
sample         = [76889,76892]                 
sample_bg      = 76894
wv_file        = 'WV_76865.txt'          #white vanadium integral file
Ei_list        = [22.78,7.52,3.7,2.2]        #240/80 Eis
Erange         = [-0.2,0.002,0.8]
trans          = [1,1,1,1]           #elastic line transmission factors for each Ei
mask           = 'LET_mask_212_base.xml'
#========================================================

#==================Absolute Units Input==================
mv_file       = None        # pre-processed MV calibration
monovan_mass  = 6.62
sample_mass   = 6.62
sample_fwt    = 50.9415
#========================================================

#==================Local Contact Input===================
datadir = '/archive/cycle_21_1/NDXLET/'
mapdir  = '/usr/local/mprogs/InstrumentFiles/let/'
maskdir = './'
m2spec  = 98310                 #specID of monitor2 (post monochromator)
m2pos   = 1.333                 #position of monitor2 wrt the sample position
iprint  = True                  #turn on output
idebug  = True                 #keep workspaces
powdermap = 'LET_rings_212_map.xml'
#========================================================

config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)

#check loaded workspaces and remove old nxspe workspaces
ws_list = ADS.getObjectNames()
nx_list = [ss for ss in ws_list if '.nxspe' in ss]
for ss in nx_list:
    ADS.remove(ss)

#load hard mask and create norm workspace
if mask not in ws_list:
    if iprint:
        print("LET: Loading hard mask %s" % mask)
    LoadMask('LET',maskdir+mask,OutputWorkspace=mask)
else:
    if iprint:
        print("LET: Using hard mask %s" % mask)
        
# load white van files
if wv_file not in ws_list:
    if iprint:
        print("LET: Loading white vanadium integrals %s" % wv_file)
    LoadAscii(Filename=wv_file,OutputWorkspace=wv_file)
    ReplaceSpecialValues(wv_file, SmallNumberThreshold=1e-20, SmallNumberValue='NaN',OutputWorkspace=wv_file)
else:
    if iprint:
        print("LET: using white vanadium integrals from %s" % wv_file)

# load monovan file
mv_fac = []
if (mv_file is not None):
    if mv_file not in ws_list:
        if iprint:
            print("LET Loading mono van calibration %s" % mv_file)
        LoadAscii(mv_file,OutputWorkspace=mv_file)
    mv_ws = mtd[mv_file]
    mv_eis = mv_ws.readX(0)
    mv_cal = mv_ws.readY(0)
# check that monovan is compatible with Ei_list)
    Ei_diff = sum([x-y for x,y in zip(mv_eis,Ei_list)])
    if (abs(Ei_diff) > 0):
        print("ERROR: Monovan file Eis not compatible with Ei_list")
    for Ei in Ei_list:
        ii = np.where(mv_eis == Ei)
        mvf = mv_cal[ii][0]
        mvf = (monovan_mass / 50.9415) / (sample_mass / sample_fwt) * (5080. / 4. / np.pi) / mvf
        mv_fac.append(mvf)
else:
    if iprint:
        print("LET: Skipping absolute calibration")
    mv_fac = [x/x for x in Ei_list]     # monovan factors = 1 by default

# =======================sum sample runs=====================================

if sumruns:
    for irun in sample:
        w_buf = Load(str(irun),LoadMonitors=True)
        if (irun == sample[0]):
            print("LET: Loading run #%i" % irun)
            ws = CloneWorkspace('w_buf'); ws_monitors = CloneWorkspace('w_buf_monitors')
        else:
            print("... Adding run #%i" % irun)
            ws = Plus('w_buf', 'ws'); w_sam_monitors = Plus('w_buf_monitors', 'ws_monitors')
    sample = [sample[0]]

# =======================sample run loop=====================================

for irun in sample:
    if not sumruns:
        print("Processing run #%i" % irun)
        ws = Load(str(irun),LoadMonitors=True)
    ws = NormaliseByCurrent('ws')

    if (sample_bg is not None):
        if iprint:
            print("... using sample background run #%i" % sample_bg)
        ws_bg = Load(str(sample_bg),LoadMonitors=False)
        ws_bg = NormaliseByCurrent('ws_bg')
        
# ============================= Ei loop =====================================  
    for ienergy in range(len(Ei_list)):
        Ei  = Ei_list[ienergy]
        tr  = trans[ienergy]
        mvf = mv_fac[ienergy]
        if iprint:
            print("\nLET:  Reducing data for Ei=%.2f meV" % Ei)

        ws_corrected = Scale('ws',1/tr,"Multiply")
        if (sample_bg is not None):
            print("... sample transmission factor = %.2f " % tr)
            ws_corrected  = ws/tr - ws_bg

# normalise to WB vanadium and apply fixed mask
        if iprint:
            print("... normalising/masking data")
        ws_norm = Divide('ws_corrected',wv_file)          # white beam normalisation
        MaskDetectors(ws_norm,MaskedWorkspace=mask)
        
# t2e section
        if iprint:
            print("... t2e section")
# this section shifts the time-of-flight such that the monitor2 peak 
# in the current monitor workspace (post monochromator) is at t=0 and L=0
# note that the offest is not the predicted value due to energy dependence of the source position        
        (_,mon2_peak,_,_) = GetEi('ws_monitors',Monitor2Spec=m2spec,EnergyEstimate=Ei,FixEi=True)
        ws_norm = ScaleX(ws_norm, Factor=-mon2_peak, Operation='Add', InstrumentParameter='DelayTime', Combine=True)
        MoveInstrumentComponent(ws_norm, ComponentName='undulator', Z=-m2pos, RelativePosition=False)

        ws_out = ConvertUnits(ws_norm,'DeltaE',EMode='Direct',EFixed=Ei)
        ws_out = Rebin(ws_out, [x*Ei for x in Erange], PreserveEvents=False)
        ws_out = DetectorEfficiencyCor(ws_out, IncidentEnergy=Ei)
        ws_out = CorrectKiKf(ws_out, Efixed=Ei, EMode='Direct')

# monovan scaling
        if iprint and mv_file is not None:
            print("... applying mono van calibration factor %.1f " % mvf)
        ws_out = Scale('ws_out',mvf,'Multiply') 

# rings grouping if desired
        ofile_suffix=''
        if powder:
            if iprint:
                print("... powder map grouping")
            ws_out=GroupDetectors(ws_out, MapFile=powdermap, Behaviour='Average')
            ofile_suffix = '_powder'

# output nxspe file
        ofile = 'RS_LET'+str(irun)+'_'+str(Ei)+'meV'+ofile_suffix+'.nxspe'
        if iprint:
            print("Writing %s" % ofile)
        
# check elastic line
        if idebug:
            Rebin('ws_out',[-0.05*Ei,100,0.05*Ei], PreserveEvents=False, OutputWorkspace=ofile+'_elastic')
            Transpose(ofile+'_elastic',OutputWorkspace=ofile+'_elastic')

# output nxspe
        ws_out = ConvertToDistribution(ws_out)
        SaveNXSPE('ws_out',ofile,Efixed=Ei,KiOverKfScaling=True)
        CloneWorkspace('ws_out',OutputWorkspace=ofile)

# cleanup
if (not idebug):
    DeleteWorkspaces(['ws_out','ws_norm','ws','ws_bg','ws_corrected','ws_monitors'])

# ============================= End of Ei loop ================================
if iprint:
    print("\nReduction complete in %.1f seconds\n" % (time.time() - t))

