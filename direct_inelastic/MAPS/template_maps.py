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

#run number
runno=21385

#incident neutron energy
ei=450

#white beam vanadium run number
wbvan=21376

#Energy bins in processed data [lo,step,hi]
rebin_pars=[-50,2.5,425]

#Monochromatic vanadium run number (for absolute units normalisation) - put [] if not using
#Sample mass (g), sample relative molar mass (both zero if not using)
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
sum_files = False

monovan=[]
sam_mass=0
sam_rmm=0

#Set background time-of-flight range (suggested values [13000,19000]  - only change if you know what you are doing!)
bg_range=[13000,19000]

#Process data as single crystal (4-to-1 detector mapping)
# set          wait_for_file=time, to wait for files to appear on the data search path
iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,sum_files,bkgd_range=bg_range)

#Process data as powder (rings mapping)
iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,sum_files,bkgd_range=bg_range)


