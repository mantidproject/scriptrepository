from mantid import config
from iliad_maps import *
import time
import os
import datetime


config['defaultsave.directory'] = '/home/maps/maps_users/Hutchings/March2015/SPE'#data save director

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

monovan=[]
sam_mass=0
sam_rmm=0

#Set background time-of-flight range (suggested values [13000,19000]  - only change if you know what you are doing!)
bg_range=[13000,19000]

#Process data as single crystal (4-to-1 detector mapping)
iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)

#Process data as powder (rings mapping)
iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
