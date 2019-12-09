from __future__ import print_function
import os
# Quantum/Mantid regression tests!
# set of source Tables, all called xxx_Tab.nxs
# run QuantumTableDrivenSimulation() on them
# set of outputs that should be matched, corresponding xxx_Result.nxs
# 
# Option to create reference files if they don't exist
createMissing=False
forceCreation=False # even if file exists and may be different!

directory=os.path.dirname(__file__)
print(directory)
files=list(map(str.lower,os.listdir(directory)))
nsims=0
nfits=0
ncreations=0
for f in files:
	if(f[-8:]=="_tab.nxs"):
		fw=f[:-8]+"_tab"
		rw=f[:-8]+"_reference"
		r=f[:-8]+"_reference.nxs"
		sw=f[:-8]+"_simulation"
		if((r in files) or createMissing):
			print("testing simulation ",fw)
			LoadNexus(Filename=os.path.join(directory,f),OutputWorkspace=fw)
			if(r in files):
				LoadNexus(Filename=os.path.join(directory,r),OutputWorkspace=rw)
			QuantumTableDrivenSimulation(ModelTable=fw,Results=sw) # might throw an error, let it be printed
			if(r in files):
				cm=CompareWorkspaces(rw,sw,Tolerance=1.E-7,CheckSpectraMap=False)
				if(cm[0]):
					print("simulation output as expected")
					nsims=nsims+1
				else:
					for cmr in cm[1]:
						print(cmr)
					if(not forceCreation):
						raise Exception(fw+" ran, but output wrong")
			if((forceCreation or not(r in files)) and createMissing):
				SaveNexus(InputWorkspace=sw,Filename=os.path.join(directory,r))
				print("created output for ",fw)
				ncreations=ncreations+1
			# Plot it (intelligent)
			if(mtd[sw].getNumberHistograms()>5 and mtd[sw].getAxis(1).isNumeric()):
				plot2D(sw)
			else:
				plotSpectrum(sw,list(range(mtd[sw].getNumberHistograms())))

# files for use as Fit Function
# table xxx_Fn.nxs, fit real data xxx_Data.nxs (spectrum 0)
# need hint of starting params xxx_InitialPars.nxs
# could be too dependent on performance of Fit()
# so evaluate function, but with params xxx_InitialPars.nxs and x values from xxx_Data.nxs and get xxx_Data(data,calc,diff)
# go on to fit properly.
	if(f[-7:]=="_fn.nxs"):
		ipw=f[:-7]+"_initialpars" # the xxx_Parameters table, perhaps edited
		ip=f[:-7]+"_initialpars.nxs" # the xxx_Parameters table, perhaps edited
		dw=f[:-7]+"_data" # also reference guess, if num of spectra > 1, ie an old xxx_Workspace
		d=f[:-7]+"_data.nxs" # also reference guess, if num of spectra > 1, ie an old xxx_Workspace
		gw=f[:-7]+"_guess" # for output
		bn=f[:-7] # base for Fit()'s _Parameters, etc
		ww=f[:-7]+"_Workspace"
		if((d in files) and (ip in files)):
			LoadNexus(Filename=os.path.join(directory,f),OutputWorkspace="Tab")
			LoadNexus(Filename=os.path.join(directory,ip),OutputWorkspace=ipw)
			LoadNexus(Filename=os.path.join(directory,d),OutputWorkspace=dw)
			# extract initial pars from table
			Pars=[]
			for row in mtd[ipw]:
				Pars.append(row['Name'].split(".")[-1]+"="+str(row['Value']))
			FnString="name=QuantumTableDrivenFunction"+str(len(Pars)-3)+","+",".join(Pars[:-1]) # note there's Scale and Baseline, and omit Cost Function
			print("testing fit function from ",f," with ",len(Pars)-3," parameters")
			Fit(Function=FnString,InputWorkspace=dw,WorkspaceIndex=0,Output=bn,MaxIterations=0,CreateOutput=True)
			RenameWorkspace(bn+"_Workspace",gw)
			if(mtd[dw].getNumberHistograms() > 1):
				cm=CompareWorkspaces(dw,gw,Tolerance=1.E-7,CheckSpectraMap=False)
				if(cm[0]):
					print("Fit ran, initial guess curve as expected")
					nfits=nfits+1
				else:
					for cmr in cm[1]:
						print(cmr)
					if(not forceCreation):
						raise Exception(f+" ran, but output wrong")
			if((forceCreation or not(mtd[dw].getNumberHistograms() > 1)) and createMissing):
				SaveNexus(InputWorkspace=gw,Filename=os.path.join(directory,d))
				print("created output for ",f)
				ncreations=ncreations+1
			# and now fit properly
			Fit(Function=FnString,InputWorkspace=dw,WorkspaceIndex=0,Output=bn,CreateOutput=True)
			# plot function and fit
			plotSpectrum(ww,[0,1])

print(nsims," simulation and ",nfits," fit function tests successfully completed!")
if(ncreations>0):
	print(ncreations," tests run for the firat time and results saved")
