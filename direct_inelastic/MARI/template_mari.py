from mantid import config
from iliad_mari import *
import time
import os
import datetime

config['defaultsave.directory'] = '/home/mari/Users/MARI_team'          #data_dir 

runno=19841
ei=750
wbvan=19717
rebin_pars=[-100,2,700]
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=0
sam_mass=0
sam_rmm=0

iliad_mari(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm)

