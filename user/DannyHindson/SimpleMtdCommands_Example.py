#
# Example: Basic Mantid commands
#
from __future__ import print_function
# import mantid algorithms
from mantid.simpleapi import *


# Print the available algorithms
#mtdHelp()

# Print information on a specific algorithm
#mtdHelp(LoadRaw)

# Perform some algorithms
LoadRaw("HET15869.raw", OutputWorkspace="test")
ConvertUnits("test","dSpacing", OutputWorkspace="converted")
Rebin("converted","0.1,0.001,5", OutputWorkspace="rebinned")

# clear up intermediate workspaces
DeleteWorkspace("test")
DeleteWorkspace("converted")

# extract the one we want
wksp = mtd['rebinned']

print("Rebinned workspace has " + str(wksp.getNumberHistograms()) + " histograms")
print("Spectrum 450's X data size = " + str(len(wksp.readX(450))) + " bin boundaries")
