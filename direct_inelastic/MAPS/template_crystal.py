from iliad_maps_setup import *

##COPY THIS FILE TO YOUR HOME DIRECTORY BEFORE EDITING - THIS IS THE MASTER COPY!!

iliad_maps_setup

##White beam vanadium run number
wbvan=19327

##Run numbers. Can be either a single number, e.g. runno=12345, or an array, e.g. runno=[12345,67890,...]
runno=[19399,19402]

##Incident energy - can be a single number, or a list, e.g. ei=60, or ei=[60,120,300], etc
ei=[60,60]

##Energy binning - [lo,step,hi]. If multiple ei specified above, then can have this as an array
##e.g. rebin_pars=[lo,step,hi], or rebin_pars=[[lo1,step1,hi1],[lo2,step2,hi2],...]
rebin_pars=[[-6,0.6,54],[-6,0.6,54]]

##Sample info - set these values to 0 if no monochromatic vanadium run performed yet
#sam_rmm=350.653
#sam_mass=2.465
sam_rmm=0
sam_mass=0


##monochromatic vanadium run number - set to [] if not yet done
#monovan=19403
monovan=[]

##Run the analysis routine for a single run, ei and set of rebin parameters
#iliad_maps_crystal(runno,ei,wbvan,rebin_pars,monovan,sam_mass,sam_rmm)

##If a list of ei is given, then can run a for-loop as below (ensure the array of ei, run numbers and rebin parameters are all the same length)
for i in range(len(runno)):
    iliad_maps_crystal(runno[i],ei[i],wbvan,rebin_pars[i],monovan,sam_mass,sam_rmm)
	
##You can also specify things like the rebinning parameters and ei inside the for-loop if you have lots of runs with the same configuration
##And there are plenty of other ways of doing this...
#for i in range(len(runno)):
#	myei=60
#	myrebin=[-6,0.6,54]
#	iliad_maps_crystal(runno[i],ei[i],wbvan,rebin_pars[i],monovan,sam_mass,sam_rmm)