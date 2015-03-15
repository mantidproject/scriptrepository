from qtiGenie import *
from PySlice2 import *

inst='mar'
iliad_setup(inst)
ext='.raw'
mapfile='mari_res2013'
#det_cal_file must be specified if the reduction sends out put to a workpsace
cal_file='MAR18622.raw'
#hard mask file
mask_file='mari_mask2014.msk'
#load vanadium file
whitebeamfile="18622"
LoadRaw(Filename=whitebeamfile,OutputWorkspace="wb_wksp",LoadLogFiles="0")
#---------------------------------------------------------------------------------------------------------------


#load run
############## normal reduction####################
runs=[19184]; ei=15; rebin_params='-10,0.03,15'; sum=False 


if sum == True:
    nxsp_file=inst+str(runs[0])+'sum.nxspe'
    for i in range(len(runs)):
        LoadRaw(Filename=str(runs[i]),OutputWorkspace="w"+str(i+1),LoadLogFiles="0")
        if i == 0: run_wksp = mtd['w1']
        else: run_wksp = run_wksp+mtd['w'+str(i+1)]
    w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,fixei=False,det_cal_file=cal_file,norm_method='current')#,hardmaskOnly=mask_file)
    RenameWorkspace(InputWorkspace='w1',OutputWorkspace=inst+str(runs[0])+'sum')
    SaveNXSPE(InputWorkspace='mar'+str(runs[0])+'sum',Filename=nxspe_file, KiOverKfScaling=True)
else:
    for runfile in runs:
        nxspe_file=inst+str(runfile)+'.nxspe'
        LoadRaw(Filename=str(runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")
        w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,fixei=False,det_cal_file=cal_file,norm_method='current')#,hardmaskOnly=mask_file)
        RenameWorkspace(InputWorkspace='w1',OutputWorkspace=inst+str(runfile))
        SaveNXSPE(InputWorkspace='mar'+str(runfile),Filename=nxspe_file, KiOverKfScaling=True)



