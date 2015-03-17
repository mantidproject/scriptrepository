from mantid import config
from iliad_maps import *
import time
import os
import datetime


config['defaultsave.directory'] = '/home/maps/maps_users/Hayden/July2014/SPE/Alex/'#data_dir 

runno=21385
ei=450
wbvan=21376
rebin_pars=[-50,2.5,425]
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=[]
sam_mass=0
sam_rmm=0

iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm)

iliad_maps_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm)
