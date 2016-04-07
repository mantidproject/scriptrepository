#Example script for using PyChop

#All examples are for MAPS, but the same routines can be used for MARI on MERLIN

#Import modules for inst desired:
from PyChopMAPS_CalibWB import *

#Also import numpy, as needed for some plotting
import numpy as np

#Run the calcuation for Ei=50meV, freq=250Hz, sloppy chopper:
Ei=500
freq=600
[res_elastic,res_inelastic,flux,tmod,tchp,ttot]=calc_chop(Ei,freq,'a')

#NB the tmod, tchp and v_van terms are internal debug, giving the moderator and chopper time widths (in seconds), and the total time width (in seconds) respectively
print(flux)
#Gives the flux in n / cm^2 / s / 160uAhrs
print(res_elastic)
#Gives the elastic resolution in meV

#################################################
#################################################

#Plot the inelastic resolution. It is calculated by PyChop for energy transfers in the range 0.05 to 0.95 Ei, in steps of 0.05Ei
plot(np.arange(0.05*Ei,Ei,0.05*Ei),np.flipud(res_inelastic),linewidth=2,color='blue')
title('MAPS Ei = '+str(Ei)+' meV, freq = '+str(freq)+' A-chopper')
xlabel('Energy transfer (meV)')
ylabel('Resolution (meV)')
grid('on')

#################################################
#################################################

#Plot the flux and elastic resolution vs incident energy (more complex example)
Ei=np.arange(80,610,10)
freq=600

res_elastic_scan=np.zeros(len(Ei))
flux_scan=np.zeros(len(Ei))
for i in range(len(Ei)):
    [res_elastic_scan[i],res_inelastic,flux_scan[i],tmod,tchp,ttot]=calc_chop(Ei[i],freq,'a')

#Flux vs Ei plot
plot(Ei,flux_scan,linewidth=2,color='red')
title('MAPS flux vs Ei, freq = '+str(freq)+' S-chopper')
xlabel('Incident energy (meV)')
ylabel('Flux (n / cm^2 / s / 160uAhrs)')
grid('on')

#Resolution vs Ei plot
plot(Ei,res_elastic_scan,linewidth=2,color='black')
title('MAPS resolution vs Ei, freq = '+str(freq)+' S-chopper')
xlabel('Incident energy (meV)')
ylabel('Elastic resolution (meV)')
grid('on')
