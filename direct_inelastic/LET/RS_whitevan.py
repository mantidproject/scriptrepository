# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import time

t = time.time()         #start the clock

#=======================User Input=======================
whitevan       = 76865
whitevan_bg    = 76868
whitevan_trans = 0.8
mask           = 'LET_mask_212_base.xml'
#========================================================

#==================Local Contact Input===================
datadir = '/archive/cycle_21_2/NDXLET/'
mapdir  = '/usr/local/mprogs/InstrumentFiles/let/'
maskdir = './'
wv_lrange = [1,5]     #wavelength range for the white van integrals
iprint = True           #turn on output
#========================================================

config.appendDataSearchDir(mapdir)
config.appendDataSearchDir(datadir)

#load hard mask and create norm workspace
if iprint:
    print("WHITE_VAN: Loading hard mask %s" % mask)
LoadMask('LET',maskdir+mask,OutputWorkspace=mask)

# load white van and background and subtract
if iprint:
    print("WHITE_VAN: Loading white vanadium run #%i" % whitevan)
wv    = Load(str(whitevan),LoadMonitors='Exclude')
wv    = CropWorkspace(wv,XMax=99720)
wv    = NormaliseByCurrent('wv')
wv_corrected = Scale('wv',1/whitevan_trans,"Multiply")
if (whitevan_bg is not None):
    if iprint:
        print("... subtracting white vanadium background run #%i" % whitevan_bg)
    wv_bg = Load(str(whitevan_bg),LoadMonitors='Exclude')
    wv_bg = CropWorkspace(wv_bg,XMax=99720)
    wv_bg = NormaliseByCurrent('wv_bg')
    wv_corrected = wv/whitevan_trans - wv_bg
    
# create wv_norm workspace (white vanadium integrals)    
if iprint:
    print("WHITE_VAN: Calculating white vanadium integrals")
wv_norm = ConvertUnits(wv_corrected,'Wavelength')
wv_norm = Rebin(wv_norm,str(wv_lrange[0])+',100,'+str(wv_lrange[1]))
MaskDetectors(wv_norm, MaskedWorkspace=mask)
normt = Transpose(wv_norm); scale = Integration(normt); scale_factor = len(normt.dataY(0))/scale.dataY(0)[0]
wv_norm = Scale(wv_norm, scale_factor, 'Multiply')
DeleteWorkspaces(['scale','normt'])
ofile = 'WV_'+str(whitevan)+'.txt'
if iprint:
    print("WHITE_VAN: Saving integrals in %s\n" % ofile)
SaveAscii(wv_norm,ofile)

# ==================================================================
if iprint:
    print("WHITE_VAN: Reduction complete in %.1f seconds" % (time.time() - t))

