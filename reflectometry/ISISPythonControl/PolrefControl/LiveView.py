import os, sys
sys.path.append(r'C:\MantidInstall')
sys.path.append(r'C:\MantidInstall\bin')
sys.path.append(r'u:\mantid\TestingPointDetectorReduction')
sys.path.append(r'u:\mantid')
from PlotScan import PlotScan, ChopperChecks
from mantid.simpleapi import *
from WrappedReduction import *

r=Runs()

r2l=refl(r[0],twotheta=4.0);
r2l=Rebin(r2l, '2.5,0.05,11')