"""
Print incident energy and t0. Useful to check statistics of monitors in individual runs.
"""

runs=range(10634,10690)							      #list if runs that have the same energy
Eguess=450									              #initial energy guess
datadir="/SNS/SEQ/IPTS-4081/data/"			  #Data directory	


def GetEiT0(filename,EiGuess):
	LoadNexusMonitors(Filename=filename,OutputWorkspace="monitor_ws")			#Load the monitors
	alg=GetEi(InputWorkspace="monitor_ws")			                          #Run GetEi algorithm
	[Ei,Tzero]=[float(alg[0],float(alg[3])]					                      #Extract incident energy and T0
	return [Ei,Tzero]	

Esum = 0.0
T0sum=0.0
idx = 0

for run in runs:
	filename=datadir+"SEQ_"+str(run)+"_event.nxs"
	[Efixed,T0]=GetEiT0(filename,Eguess)	     #Get Ei and -T0 using the function defined before
	idx+=1
	Esum += Efixed
	T0sum += T0
	print idx,run,Efixed,T0,Esum,T0sum
	
print "Average values ",Esum/idx,T0sum/idx
