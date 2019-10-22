from mantid import config
from MAPSReduction_Sample import *
import time
import os
import datetime
import sys
try:
    #Note: due to the mantid-python implementation, one needs to run this 
    #script in Mantid  script window  TWICE!!!  to deploy  the the changes made to MAPSReduction_Sample.py file.
    sys.path.insert(0,'/instrument/MAPS/RBNumber/USER_RB_FOLDER');    #<- template parameter, modified on isiscompute to correct user
    reload(sys.modules['MAPSReduction_Sample'])
except:
    print "*** WARNING can not reload MAPSReduction_Sample file"
    pass

config['defaultsave.directory'] = '/instrument/MAPS/RBNumber/USER_RB_FOLDER' #data_dir 

#=========================================================================================================================================

#run number
runno=[37395]
sum_files = False

#incident neutron energy
ei=[200]

#white beam vanadium run number
wbvan=36501

#Energy bins in processed data [lo,step,hi]
rebin_pars=[-0.5,0.002,0.95]

#Monochromatic vanadium run number (for absolute units normalisation) - put [] if not using
#Sample mass (g), sample relative molar mass (both zero if not using)
monovan=[]
sam_mass=0
sam_rmm=0

# hard mask file for the cycle
hard_mask_file = "4to1_192_msk.xml"

#Process data as single crystal (4-to-1 detector mapping)
wait_for_file=0
iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,sum_files,wait_for_file=wait_for_file,map_file="4to1.map",hard_mask_file=hard_mask_file)

#Process data as powder (rings mapping)
#iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,sum_files,hard_mask_file=hard_mask_file)


