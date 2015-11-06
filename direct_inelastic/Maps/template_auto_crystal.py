from iliad_maps_setup import *
import time
import os
import datetime

##COPY THIS FILE TO YOUR HOME DIRECTORY BEFORE EDITING - THIS IS THE MASTER COPY!!

iliad_maps_setup

##White beam vanadium run number
wbvan=19588

##Run number range you wish to process in this loop. Note that this should be sensible, as breaking out of loops is tricky!
runno=range(19739,19741)

##Incident energy - can be a single number, or a list, e.g. ei=60, or ei=[60,120,300], etc
ei=300

##Energy binning - [lo,step,hi]. If multiple ei specified above, then can have this as an array
##e.g. rebin_pars=[lo,step,hi], or rebin_pars=[[lo1,step1,hi1],[lo2,step2,hi2],...]
rebin_pars=[-30,3,285]

##Sample info - set these values to 0 if no monochromatic vanadium run performed yet
#sam_rmm=350.653
#sam_mass=29.8
sam_rmm=0
sam_mass=0


##monochromatic vanadium run number - set to [] if not yet done
monovan=[]

##Run the analysis routine, where files will be processed in order if they exist. If not, the routine waits 30 minutes and tries again
for i in range(len(runno)):
	latest_file='/home/maps/maps_data/MAP'+str(runno[i])+'.raw'
	latest_run=runno[i]
	#Check we have not reached some specified time
	timenow=datetime.datetime.now()
	#Must specify the end time in this year,month,date,hr,min,sec format
	endtime=datetime.datetime(2013,11,18,10,0,0)
	if timenow>=endtime:
		break
		
	while latest_run==runno[i]:
		if os.path.isfile(latest_file): 
			try:
				time.sleep(60)
				iliad_maps_crystal(runno[i],ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm)
				latest_run=latest_run+1
			except:
				print 'Skipped run number '+str(latest_run)
				latest_run=latest_run+1
		else:
			time.sleep(1800)

sam_mass=0
##You can also specify things like the rebinning parameters and ei inside the for-loop if you have lots of runs with the same configuration
##And there are plenty of other ways of doing this...
#for i in range(len(runno)):
#	myei=60
#	myrebin=[-6,0.6,54]
#	iliad_maps_crystal(runno[i],ei[i],wbvan,rebin_pars[i],monovan,sam_mass,sam_rmm)