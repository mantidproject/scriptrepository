from __future__ import print_function
from numpy import *
import time
units=1   #1 for hours, 60 for minutes, 3600 for secs
pl=6.1         #enter pressure length of cell

nl=0.0733*pl

##8893-8974

##Load Comparision Run in units lamda
Load(Filename=r'\\britannic\3He\LET Data\0213\LET00008922.raw',OutputWorkspace='FileI',SpectrumList='40966')
NormaliseByCurrent(InputWorkspace='FileI',OutputWorkspace='FileI')
ConvertUnits(InputWorkspace='FileI',OutputWorkspace='FileI',Target='Wavelength')

##Extract Initial Run Time
CreateLogPropertyTable(InputWorkspaces='FileI',OutputWorkspace='LogI',LogPropertyNames='run_start,run_end',GroupPolicy='All')
stI=mtd['LogI'].cell(0,0)
stI=time.strptime(stI, '%Y-%m-%dT%H:%M:%S')
stI=float(stI[7])*24+float(stI[3])+float(stI[4])/60+float(stI[3])/3600
endI=mtd['LogI'].cell(0,1)
endI=time.strptime(endI, '%Y-%m-%dT%H:%M:%S')
endI=float(endI[7])*24+float(endI[3])+float(endI[4])/60+float(endI[3])/3600

##Open file to write polarisations into
A_NPol= open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NPol', 'w')

##Open files,normalise to current
for i in range(8893,8922+1):
	print(i)
	Load(Filename=r'\\britannic\3He\LET Data\0213\LET0000'+str(i)+'.raw',OutputWorkspace='Filei',SpectrumList='40966')
	try:
		NormaliseByCurrent(InputWorkspace='Filei',OutputWorkspace='Filei')
	
## Extract Run Times (relative to beginning of first run)
		CreateLogPropertyTable(InputWorkspaces='filei',OutputWorkspace='Logi',LogPropertyNames='run_start,run_end',GroupPolicy='All')
		Endi = mtd['Logi'].cell(0,1)
		Endi = time.strptime(Endi, '%Y-%m-%dT%H:%M:%S')
		Endi=units*(float(Endi[7])*24+float(Endi[3])+float(Endi[4])/60+float(Endi[3])/3600-stI)
	
##Calculate He3 Polarisation
		ConvertUnits(InputWorkspace='Filei',OutputWorkspace='Filei',Target='Wavelength')
		Divide(LHSWorkspace='FileI',RHSWorkspace='Filei',OutputWorkspace='TnI/Tni',AllowDifferentNumberSpectra='1')
	
##Fit and write to file
		Fit(Function='name=UserFunction,Formula=(exp(2*'+str(nl)+'*x*PHeI)+1)/(exp('+str(nl)+'*x*(PHeI+PHei))+exp('+str(nl)+'*x*(PHeI-PHei))),PHeI=0.01,PHei=0.5,constraints=(1.00>=PHeI>=0.00, 1.00>=PHei>=0.00)',InputWorkspace='TnI/Tni'  ,Output='TnI/Tni',StartX='1.5',EndX='12.5')
		table = mtd['TnI/Tni_Parameters']
		A_NPol.write(str(Endi)+'	'+ str(abs(table.cell(1,1)))+'	' + str(table.cell(1,2))+'\n')
	except RuntimeError:
		pass
	
A_NPol.close()

LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NPol',OutputWorkspace='NPol',Unit='Time')