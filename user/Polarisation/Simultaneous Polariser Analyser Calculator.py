from __future__ import print_function
from numpy import *
import time
#Opacity value needs to be changed in function for different cells
#The mask applies to direct transmission on LET
#error reduced by having a greater difference in polarisation to the calibration file-
beammonitorindex=40965
Initial_File=8918
Final_File=8974
units=1   #1 for hours, 60 for minutes, 3600 for secs
#T=15   #Time between runs in mins
plpol=6.1   #pressure length of polariser
plana=8.64   #pressure length of analyser
#could enter expected polarisation
######################################################################################
UP='name=UserFunction,Formula=A*(1-exp( -(x+T0)/Tup)),T0=680,A=0.000286216,Tup=10.00,constraints=(0<=Tup)'
DOWN='name=UserFunction,Formula=A*(exp( -(x+T0)/T1)),T0=680,A=0.000286216,T1=31.11,constraints=(0<=T1)'
nl1=0.0733*plpol
nl2=0.0733*plana

##Create Mask
LoadRaw(Filename=r'\\isis\inst$\ndxlet\instrument\data\cycle_12_5\LET0000'+str(Initial_File)+'.raw',OutputWorkspace='Initial File')
MaskDetectors(Workspace='Initial File',WorkspaceIndexList='40965,5244,5244,5245,5245,5245,5245,5246,5246,5246,5246,5247,5247,5247,5247,5248,5248,5248,5248,5249,5249,5249,5249,5250,5250,5250,5250,5251,5251,5251,5251,5252,5252,5501,5501,5501,5502,5502,5502,5502,5503,5503,5503,5503,5504,5504,5504,5504,5505,5505,5505,5505,5506,5506,5506,5506,5507,5507,5507,5507,5508,5508,5508,5508,5509,5509,5757,5757,5757,5758,5758,5758,5758,5759,5759,5759,5759,5760,5760,5760,5760,5761,5761,5761,5761,5762,5762,5762,5762,5763,5763,5763,5763,5764,5764,5764,5764,5765')
ExtractMask(InputWorkspace='Initial File',OutputWorkspace='MaskWorkspace')
InvertMask(InputWorkspace='MaskWorkspace',OutputWorkspace='MaskWorkspace')

##Load calibration file and apply mask
LoadRaw(Filename=r'\\isis\inst$\ndxlet\instrument\data\cycle_12_5\LET0000'+str(Initial_File)+'.raw',OutputWorkspace='Initial File')
MaskDetectors(Workspace='Initial File',MaskedWorkspace='MaskWorkspace')
NormaliseByCurrent(InputWorkspace='Initial File',OutputWorkspace='Initial File')
ConvertUnits(InputWorkspace='Initial File',OutputWorkspace='Initial File_lam',Target='Wavelength')

##Extract Initial Run Time
CreateLogPropertyTable(InputWorkspaces='Initial File',OutputWorkspace='LogI',LogPropertyNames='run_start,run_end',GroupPolicy='All')
stI=mtd['LogI'].cell(0,0)
stI=time.strptime(stI, '%Y-%m-%dT%H:%M:%S')
stI=float(stI[7])*24+float(stI[3])+float(stI[4])/60+float(stI[3])/3600
endI=mtd['LogI'].cell(0,1)
endI=time.strptime(endI, '%Y-%m-%dT%H:%M:%S')
endI=float(endI[7])*24+float(endI[3])+float(endI[4])/60+float(endI[3])/3600

##Rebin to sum spectra and divide
Rebin(InputWorkspace='Initial File_lam',OutputWorkspace='Initial File_lam',Params='1.5,0.1,12.5')

##Separate out beam monitor
ExtractSingleSpectrum(InputWorkspace='Initial File_lam',OutputWorkspace='BeamMonI',WorkspaceIndex='40965')

SumSpectra(InputWorkspace='Initial File_lam',OutputWorkspace='Initial File_sum',IncludeMonitors=False)

##Open file to write polarisations into
A_NPol= open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NPol', 'w')
A_Ana= open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_Ana', 'w')

diagPol=open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\diagPol.csv', 'w')
diagAna=open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\diagAna.csv', 'w')
##Load file and apply mask
for i in range(Initial_File,Final_File+1):
	print(i)
	LoadRaw(Filename=r'\\isis\inst$\ndxlet\instrument\data\cycle_12_5\LET0000'+str(i)+'.raw',OutputWorkspace='filei')
	MaskDetectors(Workspace='filei',MaskedWorkspace='MaskWorkspace')
	try:
		NormaliseByCurrent(InputWorkspace='filei',OutputWorkspace='filei')
		ConvertUnits(InputWorkspace='filei',OutputWorkspace='filei_lam',Target='Wavelength')
## Extract Run Times (relative to beginning of first run)
		CreateLogPropertyTable(InputWorkspaces='filei',OutputWorkspace='Logi',LogPropertyNames='run_start,run_end',GroupPolicy='All')
		Endi = mtd['Logi'].cell(0,1)
		Endi = time.strptime(Endi, '%Y-%m-%dT%H:%M:%S')
		Endi=units*(float(Endi[7])*24+float(Endi[3])+float(Endi[4])/60+float(Endi[3])/3600-stI)
		
##Rebin to sum spectra and divide
		Rebin(InputWorkspace='filei_lam',OutputWorkspace='filei_lam',Params='1.5,0.1,12.5')

##Separate out beam monitor and calculate Polariser Phe
		ExtractSingleSpectrum(InputWorkspace='filei_lam',OutputWorkspace='BeamMoni',WorkspaceIndex='40965')
		Divide(LHSWorkspace='BeamMonI',RHSWorkspace='BeamMoni',OutputWorkspace='T1nI/T1ni',AllowDifferentNumberSpectra='1')
		Fit(Function='name=UserFunction,Formula=(exp(2*'+str(nl1)+'*x*PHeI)+1)/(exp('+str(nl1)+'*x*(PHeI+PHei))+exp('+str(nl1)+'*x*(PHeI-PHei))),PHeI=0.1,PHei=0.5,constraints=(1.00>=PHeI>=0.00, 1.00>=PHei>=0.00)',InputWorkspace='T1nI/T1ni'  ,Output='T1nI/T1ni',StartX='1.5',EndX='6.5')
		table = mtd['T1nI/T1ni_Parameters']
		A_NPol.write(str(Endi)+'	'+ str(abs(table.cell(1,1)))+'	' + str(table.cell(1,2))+'\n')
		if i == Initial_File:
			diagPol.write('run end'+','+ 'P1HeI'+',' + 'EP1HeI'+','+ 'P1Hei'+',' + 'EP1Hei'+'\n')
		diagPol.write(str(Endi)+','+ str(abs(table.cell(0,1)))+',' + str(table.cell(0,2))+','+ str(abs(table.cell(1,1)))+',' + str(table.cell(1,2))+'\n')
##Sum spectra and claculate RELATIVE transmission(=Earlier/Later)
		SumSpectra(InputWorkspace='filei_lam',OutputWorkspace='filei_sum',IncludeMonitors=False)
		Divide(LHSWorkspace='Initial File_sum',RHSWorkspace='filei_sum',OutputWorkspace='T2nI/T2ni',AllowDifferentNumberSpectra='1')

##Extract Polariser variables
		P1HeI=str(table.cell(0,1))
		P1Hei=str(table.cell(1,1))
		
##Setup Expected Polarisation
		

##Fit to find Analyser PHe3
		Fit(Function='name=UserFunction,Formula=(exp(2*x*('+str(nl1)+'*'+P1HeI+'+'+str(nl2)+'*P2HeI))+1)/(exp(x*('+str(nl1)+'*('+P1HeI+'+'+P1Hei+')+'+str(nl2)+'*(P2HeI+P2Hei)))+exp(x*('+str(nl1)+'*('+P1HeI+'-'+P1Hei+')+'+str(nl2)+'*(P2HeI-P2Hei)))),P2HeI=0.5,P2Hei=0.5,constraints=(1.00>=P2HeI>=0.00, 1.00>=P2Hei>=0.00)',InputWorkspace='T2nI/T2ni'  ,Output='T2nI/T2ni',StartX='1.5',EndX='6.5')
		table = mtd['T2nI/T2ni_Parameters']
		A_Ana.write(str(Endi)+'	'+ str(abs(table.cell(1,1)))+'	' + str(table.cell(1,2))+'\n')
		if i == Initial_File:
			diagAna.write('run end'+','+'P2HeI'+',' + 'EP2HeI'+','+ 'P2Hei'+',' + 'EP2Hei'+'\n')
		diagAna.write(str(Endi)+','+ str(abs(table.cell(0,1)))+',' + str(table.cell(0,2))+','+ str(abs(table.cell(1,1)))+',' + str(table.cell(1,2))+'\n')
	except RuntimeError:
		pass
		
A_NPol.close()
A_Ana.close()

diagPol.close()
diagAna.close()

LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NPol',OutputWorkspace='HePol',Unit='Time')
Fit(UP,InputWorkspace='HePol',Output='HePol',StartX='0',EndX=str(Endi))

LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_Ana',OutputWorkspace='HeAna',Unit='Time')
Fit(DOWN,InputWorkspace='HeAna',Output='HeAna',StartX='0',EndX=str(Endi))