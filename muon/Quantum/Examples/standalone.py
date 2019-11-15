# try Quantum without Mantid
from __future__ import print_function
from quantumtabletools import RunModelledSystem
pars={
	"spins":("Mu","e"),
	"a(Mu)":(4463,),
	"measure":("integral",),
	"loop0par":"bmag",
	"loop0range":(0.0,1.0,5),
	}
q=RunModelledSystem(pars)

x=q[0][0]
y=q[1]
e=q[2]
for i in range(len(x)):
	print(x[i],y[0,i])
