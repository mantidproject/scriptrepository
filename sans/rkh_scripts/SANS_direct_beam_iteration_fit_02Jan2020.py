# extra stuff needed for new monoGaussCoil  function
from __future__ import (absolute_import, division, print_function)
import numpy as np
import math as math
#
from mantid.simpleapi import *
from shutil import copy as shcopy
from gauss_coil_fit1 import *

#import nxs as nxs
import os
# there are a lot of useful controls over data reduction hiiding in mantidinstall\scripts\SANS\ISISCommandInerface.py
from ISISCommandInterface import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace
from reducer_singleton import ReductionSingleton
#
# RKH's Direct Beam file iterator that refines the shape of the direct beam files using a strongly scattering sample with
# not too much incoherent scattering and not too much multiple scatter, such as a d/h polystyrene polymer standard
# 19 Aug 2019, changed "maskfile" to "userfile" for consistency with current usage.
# 21 Aug 2019 All relevant parameters are now passed through the call (however the instrument name appears at lines ~ 77 & ~179
#                 though no doubt there might be some way for python to deal with this?),  lowlam & highlam now computed automatically, 
#                 db file range now checked, zeros in ratios now avoided, plot is now switched to log-log
#
# 11/12/19 (with help from David F),  this avoids having to keep regenerating mask files for each layer.
def MaskMASKFILE(details):
    ReductionSingleton().settings["MaskFiles"] = details
    print('in MaskMASKFILE details=', details )

# 16/12/19 (wth help from Matt A) now set the DB file name, no more need to keep multiple edits of user files.
# BUT see also SetCorrectionFile() which RKH wrote in March 2015 !
def SetMonDirect(details):
    ReductionSingleton().instrument.cur_detector().correction_file = details
    ReductionSingleton().instrument.other_detector().correction_file = details
    print('in SetMonDirect - Correction file is {}'.format(ReductionSingleton().instrument.cur_detector().correction_file))

#
def iterateDBFile(rsam,rcan,tsam,tdb,tcan,Now,Next,niterations,filepath,userfilename,DBfilename,q1longW,q2longW,q1shortW,q2shortW,expectI0,fitparams,tieparams,flatbkg=0.0,fitmodel="polyGaussCoil",wlist=[1.75,1.95,2.25,2.75,2.95,3.35,3.75,4.15,4.55,4.95,5.35,5.75,6.15,6.95,7.75,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.5]):
	#
	# PLEASE also read the detailed comments below from line ~89
	#
	# Assume here that if Now=1 we are starting afresh, so can copy DBfilebase.txt to DBfilebase_v1.txt
	# if Now=2, 3, 4 etc than assume that the DBfilebase_v2.txt already exists.
	# now, next are variables here, but passed as integer constants, so play safe:
	now=Now
	next=Next
	# flat background to be subtracted from everything, IF flat background is independent of wavelength this should
	# not make any difference to results, but does make the log(I(Q) plot look better. Other versions of this code have
	# attempted to use a background that varies with wavelength, see example for "shift" below.
	#flatbkg=0.25
	#nIterations=3
	nIterations=niterations
	# DBfilename DOES NOT have _v1 etc at end
	# if starting afresh generate the first DBfile
	DBfilebase=DBfilename.rsplit(".",1)[0]
	DBfileext  ="."+DBfilename.rsplit(".",1)[1]
	if (now==1):shcopy(filepath+DBfilename,filepath+DBfilebase +"_v" +str(now) +DBfileext )
	#
	# initial read direct beam file so can check its wavelength range and set uo wext1 etc outisde the iterations loop
	W1=wlist[0]
	W2=wlist[-1]
	directbeam=LoadRKH(filepath+ DBfilebase+'_v'+str(now)+DBfileext)
	w1db=directbeam.readX(0)[0]
	w2db=directbeam.readX(0)[-1]
	if(wlist[0]<w1db): 
		print( "========= ERROR ====== direct beam file needs to go to lower wavelength ==========")
		exit
	if(wlist[-1]>w2db):
		print( "========= ERROR ====== direct beam file needs to go to longer wavelength ==========")
		exit
	wext1=w1db
	wext2=0.5*(wext1 + W1) 
	wext4=w2db
	wext3=0.5*(wext4 + W2)
	print( 'wext1,2,3,4 = ',wext1,wext2,wext3,wext4 )
	# lowlam is just above wext2, highlam just below wext3
	lowlam = 0.9*wext2 + 0.1*W1
	highlam = 0.1*W2   + 0.9*wext3
	# bit messy here as saveRKH takes mid points of histogram bins, 
	step = (w2db - w1db)/len(directbeam.readY(0) -1)
	# set up parameters for InterpolatingRebin used below
	interpolatestring=str(lowlam)+','+str(step)+','+str(highlam)
	print( "lowlam, step, highlam = ", interpolatestring)
	#
	i=1
	for ii in range(nIterations):
		SANS2DTUBES()
		#Set reduction to 1D (note that if this is left out, 1D is the default)
		Set1D()
		# Read a user file.   NOTE direct beam file is then over-ridden
		# pull in the direct beam that we need to modify
		directbeam=LoadRKH(filepath+ DBfilebase+'_v'+str(now)+DBfileext)
		newdirectfile=filepath+DBfilebase+'_v'+str(next)+DBfileext
		#
		# First we do a reduction over the full wavelength range, then over the ranges in wlist[] array
		# NOTE if the userfile masks out time ranges for say Bragg peaks then some steps in wlist[] may
		# need to be adjusted if empty I(Q) appear. The steps must be small enough to cover well for Bragg dips
		# in beam windows etc, but wide enough to give reasonable statistics.on I(Q)
		# A reduction for each wavelength band is divided by the full reduction, results appearing in "ratio"
		# Each ratio is averaged over the Q range QA to QB, but QA and QB move by geometric progression 
		# (i.e. equal steps in the log(Q) plot) between q1shortW and q2shortW at short wavelength and q1longW 
		# and q2longW at long wavelength. Thus the range used moves to smaller Q towards longer wavelengths.
		#
		# I(Q) and ratio for each wavelength band are plotted for you ( the plot is now automatically switched to log-log axes)
		# Ranges q1shortW, q2shortW, q1longW & q2longW will have to be adjusted to suit the sample, the counting 
		# statistics and Q binning chosen.
		# e.g. if q1shortW is too low, there could be zeros in the lower part ofthe "ratio" workspaces at short wavelengths
		# which will upset the averaging! Since this is not always obvious, especially on a log plot, this version of the code
		# now removes the zeros.
		# IF the I(Q) from all wavelengths overlap nicely,thus the iterated direct beam file is converging, then the 
		# ratios will all be flat and tending to 1.0, try to avoid noisy data and regions dominated by background 
		#(which likely varies with wavelength ).
		# Once the direct beam file is "converged", edit a  normal user file to make use of it, but either it or the ovwerall 
		# scale factor will likely need a rescale to get the correct absolut intensities!
		# The averaged ratios are collected as a function of wavelength in workspace "scale".
		# The direct beam file itself should have quite small (and equal due to LoadRKH) wavelength bins, and MUST 
		# extend beyond the usual wavelength range as extrapolation of "scale" is required in order to 
		# keep the cubic spline interpolation function well behaved. The call to InterpolatingRebin() near end of routine
		# needs to have an extra two points below wlist[0] at wext1 & wext2, and two points past the last value in 
		# wlist (which is wlist[-1] ) at wext3 & wext4. The interpolation can run from just above wext2 to just below wext3.
		#
		# Extrapolations to "scale" previously have some hard coded wavelength values, but the version here 
		# tries to automate this, see wext1 & wext2 at short wavelength, wext3 & wext4 at long wavelength.
		#
		# Workspace "scale2" contains the final correction function to multiply the previous direct beam, plot
		# this on top of "scale" to check that all looks well.
		# Workspaces direct_beam_hist and direct_beam_new_hist should be self explanatory.
		#
		# We ought to try to fit the transmissions (at least for ESS), this might perhaps remove the need for the smoothing carried out by
		# InterpolatingRebin, but risks hiding any genuine Bragg dips, which are sharp on a short pulse source.
		#
		# NOTE if your direct beam curve gets zeros in it they will stay there forever as the adjustments
		# here are multiplicative, so if say  you extend the wavlength range you may have to manually edit in some
		# values or otherwise generate a new curve.
		#
		# (NOTE - there is also much horrible messing about between Mantid distributions and histograms.)
		#
		MaskFile(userfilename)
		# can add this here   MaskMASKFILE('Mask_tubeends_deadtubes.xml,Mask_beamstop.xml,Mask_extra.xml')
		# Change default limits if we want to (in principle the whole reduction can be set up by python calls rather
		# than having a user file)
		#LimitsR(50.,170.)
		#   Here the arguments are min, max, step, step type
		#LimitsWav(1.75,16.5,0.1, 'LIN')
		#LimitsQ(0.004, 1.0, 0.002, 'LIN')
		# command to switch banks in from Python is Detector('bank-name'), i.e. Detector('front-detector').
		#Detector('front-detector')
		# Gravity(True) or Gravity(False) that will switch accounting for gravity on and off 
		# in the next call to WavRangeReduction.
		#
		AssignSample(rsam+'.nxs')
		AssignCan(rcan+'.nxs')
		TransFit('Off',lambdamin=W1,lambdamax=W2)
		TransmissionSample(tsam+'.nxs', tdb+'.nxs')
		TransmissionCan(tcan+'.nxs', tdb+'.nxs')
		nloops=len(wlist)-1
		# tried set DB file immediately after the user file, but something else reinitialised it to empty string (needs investigation!), works here right before WavRangeReduction
		SetMonDirect(filepath+ DBfilebase+'_v'+str(now)+'.txt')
		#MaskMASKFILE(view_layer_base+'.xml')
		print('Correction file is {}'.format(ReductionSingleton().instrument.cur_detector().correction_file))
		print('MASKFILE= {}'.format(ReductionSingleton().settings["MaskFiles"]))
		reducedAll = WavRangeReduction(wlist[0], wlist[nloops], False)
		
		# RKH 16/10/19 try fiiting the full range here, then rescale the reduced data AND the initial direct beam ==========================================
		fitstring="name="+fitmodel+","+fitparams
		Fit(fitstring,reducedAll,Ties=tieparams,StartX=q1longW,EndX=q2shortW,CreateOutput=True)

		aa=mtd[reducedAll+"_Parameters"].row(0).values()
		#print aa
		I0_fit=aa[1]
		fit_by_expect = I0_fit / expectI0
		print( "I0_fit = ",I0_fit," +- ",aa[2],"    expectI0 = ",expectI0," fit_by_expect = ",fit_by_expect)
		Scale(InputWorkspace=reducedAll,OutputWorkspace=reducedAll,Factor=1.0/fit_by_expect,Operation='Multiply')
		Scale(InputWorkspace=directbeam,OutputWorkspace=directbeam,Factor=fit_by_expect,Operation='Multiply')
		
		# ==========================================================================================================================================
		SaveRKH(reducedAll,reducedAll,Append=False)
		shift=-flatbkg
		Scale(InputWorkspace=reducedAll,OutputWorkspace=reducedAll,Factor=shift,Operation='Add')
		plt = plotSpectrum(reducedAll,0)
		data = mtd[reducedAll].readY(0)
		plotmin = min(data)/5.0
		plotmin = max(plotmin,1e-04)
		plotmax = max(data)*2.0
		activelayer = plt.activeLayer()
		# force log axes, but annonyingly this explicitly needs the ranges
		activelayer.setAxisScale(Layer.Left, plotmin,plotmax, Layer.Log10)
		activelayer.setAxisScale(Layer.Bottom,mtd[reducedAll].readX(0)[0]*0.95 ,mtd[reducedAll].readX(0)[-1]*1.05, Layer.Log10)
		#
		outlist=[]
		scaleX=[]
		scaleY=[]
		scaleE=[]
		#  must go to unity at extremes of ranges, these get modified again below
		scaleX.append(wext1)
		scaleY.append(1.0)
		scaleE.append(0.0)
		scaleX.append(wext2)
		scaleY.append(1.0)
		scaleE.append(0.0)
		#
		for i in range(nloops):
			SANS2DTUBES()
			#Set reduction to 1D (note that if this is left out, 1D is the default)
			Set1D()
			MaskFile(userfilename)
			# can over-ride here e.g. MaskMASKFILE('Mask_tubeends_deadtubes.xml,Mask_beamstop.xml,Mask_extra.xml')
			# Assign run numbers (.nxs for nexus)
			AssignSample(rsam+'.nxs')
			AssignCan(rcan+'.nxs')
			wav1 = wlist[i]
			wav2 = wlist[i+1]
			# must make sure that if fit transmissions it is done over same wavelength range as for full reduction
			# so that same numerical values are used at same wavelengths
			TransFit('Off',lambdamin=W1,lambdamax=W2)
			TransmissionSample(tsam+'.nxs', tdb+'.nxs')
			TransmissionCan(tcan+'.nxs', tdb+'.nxs')
			SetMonDirect(filepath+ DBfilebase+'_v'+str(now)+'.txt')
			#    MaskMASKFILE(view_layer_base+'_'+str(layer)+'.xml')
			print('Correction file is {}'.format(ReductionSingleton().instrument.cur_detector().correction_file))
			#    print('MASKFILE= {}'.format(ReductionSingleton().settings["MaskFiles"]))
			reduced = WavRangeReduction(wav1, wav2, False)
			SaveRKH(reduced,reduced,Append=False)
			# BKG larger at long wavelength, shift be change in bkg
			#wav3=(wav1+wav2)*0.5
			# cubic fit to  batch mode flat bkg fits using fixed scale & Rg in SASVIEW, where I(Q=0)=51.922	
			#shift=-0.42121+wav3*0.012356-wav3*wav3*0.0027909 +flatbkg
			#
			# RKH 16/10/19, though we scaled the directbeam, the non-scaled one is still being used for the reduction!
			Scale(InputWorkspace=reduced,OutputWorkspace=reduced,Factor=1.0/fit_by_expect,Operation='Multiply')
			shift=-flatbkg
			Scale(InputWorkspace=reduced,OutputWorkspace=reduced,Factor=shift,Operation='Add')
			SaveRKH(reduced,reduced+"_shifted",Append=False)
			#    
			outlist.append(reduced)
			mergePlots(plt,plotSpectrum(reduced,0))
			# set up varying Q range to integrate over
			QA = q1longW +(1./wav2-1./W2)*(q1shortW-q1longW)/(1./W1-1./W2)
			QB = q2longW  +(1./wav1-1./W2) *(q2shortW-q2longW)/(1./W1-1./W2)
			QAA =mtd[reduced].readX(0)[0]
			# avoid getting zeros in ratio which then bias the average!
			if(QAA > QA):
				print( "=======  increased QA from ",QA," to ",QAA)
				QA=QAA
			QBB =  mtd[reduced].readX(0)[-1]
			if(QBB < QB):
				print( "=======  decreased QB from ",QB," to ",QBB)
				QB=QBB
			# ratio full wavelength data over chosen Q range
			print( "wav1, wav2 = ", wav1,wav2,"  QA, QB = ", QA,QB)
			ratio='ratio_'+str(wav1)+'_'+str(wav2)
			# If the division falls over here it is usually beacuse you have empty workspace due to Bragg time masks
			CropWorkspace(reducedAll,OutputWorkspace=ratio,Xmin=QA,Xmax=QB)
			#numerator=CropWorkspace(reduced,Xmin=QA,Xmax=QB)
			numerator=RebinToWorkspace(reduced,ratio)
			Divide("numerator",ratio,OutputWorkspace=ratio)
			mergePlots(plt,plotSpectrum(ratio,0))
			#  now integrate the ratio to get the average, NOTE it rounds in QA and QB, so to get correct average need to read back the actual range
			# "ratio" is a histogram, so to get the mean value we have to divide by the Q range
			Integration(ratio,Outputworkspace="integral",RangeLower=QA,RangeUpper=QB)
			QBB= mtd["integral"].readX(0)[1]
			QAA= mtd["integral"].readX(0)[0]
			norm = mtd["integral"].readY(0)[0]/(QBB-QAA)
			normE = mtd["integral"].readE(0)[0]/(QBB-QAA)
			print( QAA,QBB,norm,normE)
			# "scale" will be a histogram, so have to multiply by the wavelength range
			scaleX.append(wav1)
			int=mtd["integral"].readY(0)[0]
			scaleY.append(norm*(wav2-wav1))
			scaleE.append(normE*(wav2-wav1))

		scaleX.append(wav2)
		print( "wav2 ",wav2)
		#  add two dummy points of same value (in distribution) at end of list
		scaleX.append(wext3)
		scaleY.append(norm*(wext3-wav2))
		scaleE.append(normE*(wext3-wav2))
		scaleX.append(wext4)
		scaleY.append(norm*(wext4-wext3))
		scaleE.append(normE*(wext4-wext3))
		# now change the two dummy points at the start 
		scaleY[0]=scaleY[2]*(wext2-wext1)/(wlist[1]-wlist[0])
		scaleY[1]=scaleY[2]*(wlist[0]-wext2)/(wlist[1]-wlist[0])
		# ======================================================================================
		#
		CreateWorkspace(scaleX,scaleY,scaleE,UnitX="Wavelength",OutputWorkspace="scale")
		# small problem, LoadRKH brings in a distribution,, InterpolatingRebin needs a histogram
		# ConvertToHistogram leaves Y values unchanged, just splits the x values, adding one extra x 
		directbeam_hist=ConvertToHistogram("directbeam")
		#
		# interpolating rebin expects a proper histogram with "counts and time", ought to write own function here,
		# but it does some clever stuff with cubic splines,
		# it seems to need a big over-run at the ends of the range - at least 2 spare points each end
		#  BEWARE - hard coded wavelength bin size BUT this does not really matter, as RebinToWorkspace below shifts results around again below!
		scale2=InterpolatingRebin("scale",interpolatestring)
		# ConvertToDistribution only divides Y values by bin width, leaves original histogram type X values !
		ConvertToDistribution("scale")
		ConvertToDistribution("scale2")
		#CropWorkspace("directbeam_hist","directbeam_hist_crop",XMin=1.0,XMax=15.0)
		#
		# Suspect that this is putting in ZEROS when outside the original wavelength range in scale2
		scale2=RebinToWorkspace(scale2,"directbeam_hist")
		directbeam_new_hist=Multiply("directbeam_hist","scale2")
		SaveRKH("directbeam_new_hist",newdirectfile,Append=False)
		print( "Wrote new direct beam file:",newdirectfile)
		# now get it back as point mode (there are some small cumulative rounding errors introduced here on the X values)
		new=LoadRKH(newdirectfile)
		RenameWorkspace('new','directbeam'+'_v'+str(next))
		RenameWorkspace('scale','scale'+'_v'+str(now))
		RenameWorkspace('scale2','scale2'+'_v'+str(now))
		now=next
		next=next+1
		i=i+1
	print( "NOTE final direct beam is only properly scaled to I0 from fit if it is converged!")

# ============================== main code runs from here  =====================================================================================
#
# this routine cannot be run blindly, it will need careful adjustments to the parameters - see comments in code above.
# NOTE you must have a userfilebase.txt  or .dat depending on which beam line you are on ....  where the directbeam does NOT have a version number, i.e. is the same as DBfilebase.dat
#
#  OOPS some DB files are .txt, some are .dat so pass the full name with extension
#  OOPS expected I(Q=0) is wrong here as I didn't know which polymer this run is for !
#iterateDBFile(         rsam,            rcan         tsam,        tdb,        tcan,           now,next,niterations,filepath,userfilename,DBfilename,
iterateDBFile("57200","57202","57185","57187","57187",1,12,2,"c:/mantidinstall_OLD/my_output/","USER_SANS2D_184F_2p4_4m_M3_Andre_12mm_Changer.txt","DIRECT_SANS2D_REAR_57200-add_4m_12mm_20x20_24Feb19.txt",q1longW=0.007,q2longW=0.05,q1shortW=0.05,q2shortW=0.175,expectI0=50.0,fitparams="I0=50.0,Mw_by_Mn=1.02,Background=0.05",tieparams="Mw_by_Mn=1.02",flatbkg=0.05,fitmodel="polyGaussCoil",wlist=[1.75,2.0,13.0,16.5])
# RTi polymer used on Larmor is  I0 55.77, Rg 60, Pd 1.02
#CloneWorkspace("41596rear_1D_1.0_12.9",OutputWorkspace="41596rear_1D_1.0_12.9"+"_dist")
#ConvertToDistribution("41596rear_1D_1.0_12.9"+"_dist")
# the simple inline formula here does not give correct results and fails to converge (note 500 iterations limit reached) ,, suspect rounding errors????
# so need to create own function, based on one from sasview
#Fit("name=UserFunction, Formula=bkg+I0*(exp(-x*Rg)+x*Rg -1.0)/(0.5*x*x*Rg*Rg),I0=60,Rg=50,bkg=0.25","41596rear_1D_1.0_12.9",Constraints="40.0<I0<100.0,40.0<Rg<80.0,0.1<bkg<0.7",StartX=0.001,EndX=0.4,CreateOutput=True)
'''
# some test code for the fitting:
LoadRKH("c:/MantidInstall/my_output/41596rear_1D_1.0_12.9",OutputWorkspace="test")
#
# test new monoGaussCoil - same result as sasview on one example so far ....
Fit("name=monoGaussCoil,I0=60,Rg=50,Background=0.25","test",StartX=0.001,EndX=0.4,CreateOutput=True)
#
# test new polyGaussCoil - same result as sasview to ~ 3 significant figures  on one example so far ....
Fit("name=polyGaussCoil,I0=40,Rg=50,Mw_by_Mn=1.01,Background=0.25","test",Ties="Mw_by_Mn=1.02",StartX=0.001,EndX=0.4,CreateOutput=True)
Fit("name=polyGaussCoil,I0=40,Rg=50,Mw_by_Mn=1.02,Background=0.25","test",StartX=0.001,EndX=0.4,CreateOutput=True)
'''