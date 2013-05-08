"""
Merge all MD nexus files in a folder
"""
import os,glob

inputdirectory="/SNS/CNCS/IPTS-7236/shared/12mev_1p5K_step1_noV/"
outputfile=inputdirectory+"all.nxs"
os.chdir(inputdirectory)
Filenames=glob.glob("*MD.nxs")
s=", "+inputdirectory
MergeMDFiles(Filenames=inputdirectory+s.join(Filenames),OutputFilename=outputfile,OutputWorkspace="md1")
