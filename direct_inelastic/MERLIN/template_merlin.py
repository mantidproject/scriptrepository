from mantid import config
from iliad_merlin import *
import time
import os
import datetime

#reload(iliad_merlin)

config['defaultsave.directory'] = '/home/zmp58988'#data save directory

#run number
runno=24434

#incident neutron energy
#ei=[180,68,35]


#white beam vanadium run number
wbvan=23684

#Energy bins in processed data [lo,step,hi]
rebin_pars=[-0.25,0.005,0.85]

#Monochromatic vanadium run number (for absolute units normalisation) - put [] if not using
#Sample mass (g), sample relative molar mass (both zero if not using)
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439

monovan=[]
sam_mass=0
sam_rmm=0

#Set background time-of-flight range (suggested values [13000,19000]  - only change if you know what you are doing!)
bg_range=[18000,19000]


############
#Process data as single crystal (4-to-1 detector mapping)
#iliad_merlin_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
#iliad_merlin_crystal(24516,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
############


############
#Process data as powder (rings mapping)
#iliad_merlin_powder(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
############

############
#Process multiple runs, waiting for data to arrive, with a specified end time (can also be done using the master script
#but it is easier to see what is going on here)

#Run number is range start, end+1 (due to annoying python syntax)
runno=range(24742,25515)
ei=[120,43,22]

for i in range(len(runno)):
	latest_file='/archiver55/cycle_14_3/NDXMERLIN/MER'+str(runno[i])+'.raw'
	latest_run=runno[i]
	#Check we have not reached some specified time
	timenow=datetime.datetime.now()
	#Must specify the end time in this year,month,date,hr,min,sec format
	endtime=datetime.datetime(2015,4,14,22,10,0)
	if timenow>=endtime:
		break
		
	while latest_run==runno[i]:
		if os.path.isfile(latest_file): 
			try:
				Pause(1)
				iliad_merlin_crystal(runno[i],ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm,bg_range)
				latest_run=latest_run+1
			except:
				print 'Skipped run numr '+str(latest_run)
				latest_run=latest_run+1
		else:
			Pause(600)

sam_mass=0



