from __future__ import print_function
from mantid.simpleapi import *
from shutil import copyfile
from math import *
import nxs as nxs
import os
from ISISCommandInterface import *
from mantid.api import WorkspaceGroup
from mantid.api import IEventWorkspace

# Taken from the offspec scripts. It outputs strings but seems to work with addruns
'''
parse a text string of the format "1-6:2+8+9,10+11+12+13-19:3,20-24"
to return a structure containing the separated lists [1, 3, 5, 8, 9], 
[10, 11, 12, 13, 16, 19] and [20, 21, 22, 23, 24]
as integer lists that addRuns can handle.
'''
def parseRunList(istring):
	if len(istring) >0: 
		s1=istring.split(',')
		rlist1=[]
		for i in range(len(s1)):
			tstr=s1[i].strip()
			if len(tstr) > 0:
				rlist1.append(tstr)
		rlist=[]
		for i in range(len(rlist1)):
			rlist2=[]
			if rlist1[i].find('+') >= 0:
				tstr=rlist1[i].split('+')
				for j in range(len(tstr)):
					if tstr[j].find(':') >=0 and tstr[j].find('-') >=0:
						tstr[j].strip()
						tstr2=tstr[j].split('-')
						tstr3=tstr2[1].split(':')
						r1=list(range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1])))
						for k in r1:
							rlist2.append(str(k))
					elif tstr[j].find('-') >=0:
						tstr[j].strip()
						tstr2=tstr[j].split('-')
						r1=list(range(int(tstr2[0]),int(tstr2[1])+1))
						for k in r1:
							rlist2.append(str(k))
					else:
						rlist2.append(tstr[j])
			else:
				if rlist1[i].find(':') >=0 and rlist1[i].find('-')>=0:
					rlist1[i].strip()
					tstr2=rlist1[i].split('-')
					tstr3=tstr2[1].split(':')
					r1=list(range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1])))
					for k in r1:
						rlist2.append(str(k))
				elif rlist1[i].find('-')>=0:
					rlist1[i].strip()
					tstr2=rlist1[i].split('-')
					r1=list(range(int(tstr2[0]),int(tstr2[1])+1))
					for k in r1:
						rlist2.append(str(k))
				else:
					rlist2.append(rlist1[i])
			rlist.append(rlist2)
	return rlist

'''
SANS add runs from RKH
'''

def add_runs(runlist,pathout,instrument='LARMOR',keepwksp=0,savewksp=1,eventbinning=0):
	if(os.path.isdir(pathout)=='False'):
		return
	pfix=instrument
	runlist=parseRunList(runlist)
	print("runlist=",runlist)
	runlist=runlist[0]
	Load(runlist[0],LoadMonitors='1',OutputWorkspace='added')
	#print "          uampHr = ", added.getRun().getProtonCharge()
	if isinstance(mtd['added'],IEventWorkspace):
	  if(eventbinning==1):
	    #Rebin('added','5.0,20.0,100000.0',PreserveEvents=False,OutputWorkspace='addedreb')
	    #Rebin('added_monitors','5.0,20.0,100000.0',OutputWorkspace='addedmonreb')
	    #ConjoinWorkspaces('addedmonreb','addedreb',CheckOverlapping=False)
	    #RenameWorkspace('addedmonreb',OutputWorkspace='added')
	    #DeleteWorkspace('added_monitors')
	    pass        
	  else:
	    Rebin('added','5.0,100.0,100000.0',PreserveEvents=False,OutputWorkspace='addedreb')
	    Rebin('added_monitors','5.0,100.0,100000.0',OutputWorkspace='addedmonreb')
	    ConjoinWorkspaces('addedmonreb','addedreb',CheckOverlapping=False)
	    RenameWorkspace('addedmonreb',OutputWorkspace='added')
	    DeleteWorkspace('added_monitors')
	else:
		pass
		ConjoinWorkspaces('added_monitors','added',CheckOverlapping=False)
		RenameWorkspace('added_monitors',OutputWorkspace='added')

	
	if(len(runlist)>1):
		for i in range(1,len(runlist)):
			#print "adding run ",str(i)        
			Load(runlist[i],LoadMonitors='1',OutputWorkspace='wtemp')
			if isinstance(mtd['wtemp'],IEventWorkspace):
			  if(eventbinning==1):
			    #Rebin('wtemp','5.0,20.0,100000.0',PreserveEvents=False,OutputWorkspace='wtempreb')
			    #Rebin('wtemp_monitors','5.0,20.0,100000.0',OutputWorkspace='wtempmonreb')
			    #ConjoinWorkspaces('wtempmonreb','wtempreb',CheckOverlapping=False)
			    #RenameWorkspace('wtempmonreb',OutputWorkspace='wtemp')
			    #DeleteWorkspace('wtemp_monitors')
			    pass
			  else:
			    Rebin('wtemp','5.0,100.0,100000.0',PreserveEvents=False,OutputWorkspace='wtempreb')
			    Rebin('wtemp_monitors','5.0,100.0,100000.0',OutputWorkspace='wtempmonreb')
			    ConjoinWorkspaces('wtempmonreb','wtempreb',CheckOverlapping=False)
			    RenameWorkspace('wtempmonreb',OutputWorkspace='wtemp')
			    DeleteWorkspace('wtemp_monitors')
			else:
				#pass
				ConjoinWorkspaces('wtemp_monitors','wtemp',CheckOverlapping=False)
				RenameWorkspace('wtemp_monitors',OutputWorkspace='wtemp')
			#print "          uampHr = ", wtemp.getRun().getProtonCharge()
			Plus('added','wtemp',OutputWorkspace='added')
			#Plus('added_monitors','wtemp_monitors',OutputWorkspace='added_monitors')
		DeleteWorkspace("wtemp")
	
	if savewksp==1:
		nzeros=8-len(runlist[0])
		fpad=""
		for ii in range(nzeros):
			fpad+="0"
		print("writing file:   "+pathout+pfix+fpad+runlist[0]+"-add.nxs")
		SaveNexusProcessed("added",pathout+pfix+fpad+runlist[0]+"-add.nxs")
		added=mtd['added']
		#print "    total uampHr = ", added.getRun().getProtonCharge()
    
	retwksp=''
	if keepwksp==0:
	    try:
	        DeleteWorkspace("added")
	        DeleteWorkspace("added_monitors")
	    except:
	        pass
	else:
	    try:
	        RenameWorkspace("added",OutputWorkspace=runlist[0]+'-add')
	        retwksp=[runlist[0]+'-add']
	        RenameWorkspace("added_monitors",OutputWorkspace=runlist[0]+'_monitors-add')
	        retwksp=[runlist[0]+'-add',runlist[0]+'_monitors-add']
	    except:
	        pass
	return retwksp

def larmor1D(rsample,rcan,tsample,tdb,tcan,tdbcan="",maskfile="",wkspname="",lmin=0.9,lmax=12.5,setthickness=0,thickness=1.0,setwidth=0,width=6.0,setheight=0,height=8.0,diagnostics=0,periods=[-1,-1,-1,-1,-1,-1],dirname='c:/Data/Processed/',saveFile=1):
    '''
    periods array is defined in order sample sans, can sans, sample trans, sample db, can trans, can db
    if tdbcan is not defined then tdb is used
    '''
    LARMOR()
    #Set reduction to 1D (note that if this is left out, 1D is the default)
    Set1D()
    MaskFile(maskfile)
    # Assign run numbers (.nxs for nexus)
    if(periods[0]>0):
        AssignSample(rsample+'.nxs',period=periods[0])
    else:
        AssignSample(rsample+'.nxs')
    a1=rsample.split('-')
    sroot=a1[0]
    if(setthickness!=0):
        ReductionSingleton().get_sample().geometry.thickness=thickness
    if(setwidth!=0):
        ReductionSingleton().get_sample().geometry.width=width
    if(setheight!=0):
        ReductionSingleton().get_sample().geometry.height=height
    if(len(rcan)>0):
        if(periods[1]>0):
            AssignCan(rcan+'.nxs',period=periods[1])
        else:
            AssignCan(rcan+'.nxs')
        a1=rcan.split('-')
        croot=a1[0]
    TransFit('On',lambdamin=lmin,lambdamax=lmax)
    if(len(tsample)>0):
        TransmissionSample(tsample+'.nxs', tdb+'.nxs',period_t=periods[2],period_d=periods[3])
        a1=tsample.split('-')
        tsroot=a1[0]
        a1=tdb.split('-')
        tdbroot=a1[0]
    if(len(tcan)>0):
        if(len(tdbcan)>0):
            TransmissionCan(tcan+'.nxs', tdbcan+'.nxs',period_t=periods[4],period_d=periods[5])
        else:
            TransmissionCan(tcan+'.nxs', tdb+'.nxs',period_t=periods[4],period_d=periods[5])
        a1=tdb.split('-')
        tdbroot=a1[0]
        a1=tcan.split('-')
        tcroot=a1[0]
    reduced = WavRangeReduction(lmin, lmax, False)
    if(periods[0]>0):
        sroot=sroot+'p'+str(periods[0])
        croot=croot+'p'+str(periods[0])
        tsroot=tsroot+'p'+str(periods[0])
        tcroot=tsroot+'p'+str(periods[0])
        tdbroot=tsroot+'p'+str(periods[0])
    if(len(wkspname) > 0):
        try:
            RenameWorkspace(sroot+'rear_1D_'+str(lmin)+'_'+str(lmax),OutputWorkspace=wkspname)
        except:
            pass
    if(diagnostics==0):
        try:
            DeleteWorkspace(sroot+'rear_1D_'+str(lmin)+'_'+str(lmax)+'_incident_monitor')
        except:
            pass
        try:
            DeleteWorkspace(sroot+'_sans_nxs')
        except:
            pass
        if(len(tsample)>0):
            try:
                DeleteWorkspace(tsroot+'_trans_nxs')
            except:
                pass
            try:
                DeleteWorkspace(tsroot+'_trans_sample_'+str(lmin)+'_'+str(lmax)+'_unfitted')
            except:
                pass
        if(len(rcan)>0):
            try:
                DeleteWorkspace(croot+'_sans_nxs')
            except:
                pass
            try:
                DeleteWorkspace(sroot+'rear_1D_'+str(lmin)+'_'+str(lmax)+'_can_tmp_incident_monitor')
            except:
                pass
        if(len(tcan)>0):
            try:
                DeleteWorkspace(tcroot+'_trans_nxs')
            except:
                pass
            try:
                DeleteWorkspace(tsroot+'_trans_can_'+str(lmin)+'_'+str(lmax)+'_unfitted')
            except:
                pass
        if(len(tsample)>0 or len(tcan)>0):
            try:
                DeleteWorkspace(tdbroot+'_trans_nxs')
            except:
                pass
    elif(diagnostics==1):
        try:
            DeleteWorkspace(sroot+'rear_1D_'+str(lmin)+'_'+str(lmax)+'_incident_monitor')
        except:
            pass
        try:
            DeleteWorkspace(sroot+'_sans_nxs')
        except:
            pass
        if(len(tsample)>0):
            try:
                DeleteWorkspace(tsroot+'_trans_nxs')
            except:
                pass
        if(len(rcan)>0):
            try:
                DeleteWorkspace(croot+'_sans_nxs')
            except:
                pass
            try:    
                DeleteWorkspace(sroot+'rear_1D_'+str(lmin)+'_'+str(lmax)+'_can_tmp_incident_monitor')
            except:
                pass
        if(len(tcan)>0):
            DeleteWorkspace(tcroot+'_trans_nxs')
        if(len(tsample)>0 or len(tcan)>0):
            try:
                DeleteWorkspace(tdbroot+'_trans_nxs')
            except:
                pass
    if saveFile==1:
        # save to file but make sure we don't have any : characters in the filename as windows doesn't like it
        wkspname1=string.replace(wkspname,':','_')
        SaveRKH(wkspname,dirname+'r'+rsample+'_'+wkspname1+'.txt',Append=0)
        SaveCanSAS1D(wkspname,dirname+'r'+rsample+'_'+wkspname1+'.xml',RadiationSource='Spallation Neutron Source',Append=0)
    

def larmor2D(rsample,rcan,tsample,tdb,tcan,tdbcan="",maskfile="",wkspname="",lmin=0.9,lmax=12.5,setthickness=0,thickness=1.0,setwidth=0,width=6.0,setheight=0,height=8.0,diagnostics=0,periods=[-1,-1,-1,-1,-1,-1],dirname='c:/Data/Processed/',saveFile=1):
    '''
    periods array is defined in order sample sans, can sans, sample trans, sample db, can trans, can db
    if tdbcan is not defined then tdb is used
    '''
    LARMOR()
    #Set reduction to 1D (note that if this is left out, 1D is the default)
    Set2D()
    MaskFile(maskfile)
    # Assign run numbers (.nxs for nexus)
    if(periods[0]>0):
        AssignSample(rsample+'.nxs',period=periods[0])
    else:
        AssignSample(rsample+'.nxs')
    a1=rsample.split('-')
    sroot=a1[0]
    if(setthickness!=0):
        ReductionSingleton().get_sample().geometry.thickness=thickness
    if(setwidth!=0):
        ReductionSingleton().get_sample().geometry.width=width
    if(setheight!=0):
        ReductionSingleton().get_sample().geometry.height=height
    if(len(rcan)>0):
        if(periods[1]>0):
            AssignCan(rcan+'.nxs',period=periods[1])
        else:
            AssignCan(rcan+'.nxs')
        a1=rcan.split('-')
        croot=a1[0]
    TransFit('On',lambdamin=lmin,lambdamax=lmax)
    if(len(tsample)>0):
        TransmissionSample(tsample+'.nxs', tdb+'.nxs',period_t=periods[2],period_d=periods[3])
        a1=tsample.split('-')
        tsroot=a1[0]
        a1=tdb.split('-')
        tdbroot=a1[0]
    if(len(tcan)>0):
        if(len(tdbcan)>0):
            TransmissionCan(tcan+'.nxs', tdbcan+'.nxs',period_t=periods[4],period_d=periods[5])
        else:
            TransmissionCan(tcan+'.nxs', tdb+'.nxs',period_t=periods[4],period_d=periods[5])
        a1=tdb.split('-')
        tdbroot=a1[0]
        a1=tcan.split('-')
        tcroot=a1[0]    
    reduced = WavRangeReduction(lmin, lmax, False)

def calibrateTubes(wkspName,calibrationfile='8tubeCalibration_25-03-2015_r2284-2296.nxs'):
   #
   # pixel by pixel efficiency correction for the linear detector
   #
   
   flood_wksp = "8tubeCalibration_25-03-2015_r2284-2296"
   if  flood_wksp not in mtd.getObjectNames():
		#Load("Vanadium_tube_calib_1to1_May14.nxs",OutputWorkspace="Vanadium_tube_calib_1to1_May14")
		#Load("Vanadium_tube_calib_1to1_25062014.nxs",OutputWorkspace="Vanadium_tube_calib_1to1_25062014")
		#Load(calibrationfile,OutputWorkspace="Vanadium_tube_calib_1to1_30062014")
		Load(calibrationfile,OutputWorkspace="8tubeCalibration_25-03-2015_r2284-2296")

   CopyInstrumentParameters("8tubeCalibration_25-03-2015_r2284-2296",wkspName)

def virtualslitxml(ypos):
	# Active detector length is ~640mm
	pixelsize=0.640/500.0
	shapexml='<cuboid id="shape">'
	shapexml+='<left-front-bottom-point x="-0.4" y="'+str(ypos-0.5*pixelsize)+'" z="29.6"/>'
	shapexml+='<left-front-top-point x="-0.4" y="'+str(ypos+0.5*pixelsize)+'" z="29.6"/>'
	shapexml+='<left-back-bottom-point x="-0.4" y="'+str(ypos-0.5*pixelsize)+'" z="29.8"/>'
	shapexml+='<right-front-bottom-point x="0.4" y="'+str(ypos-0.5*pixelsize)+'" z="29.6"/>'
	shapexml+='</cuboid>'
	return shapexml

def createMapFile(rnum,fname):
	'''
	Generate the mapfile required for create1DPSD
	e.g. createMapFile(1203,'w:\\Users\Masks\Integrate_tubes_map_file.xml')
	'''
	w1=Load(str(rnum))
	calibrateTubes(w1,"")
	det=CropWorkspace(w1,StartWorkspaceIndex=10)
	f=open(fname,'w')
	s='<?xml version="1.0" encoding="UTF-8" ?>\n'
	f.write(s)
	s='<detector-grouping>\n'
	f.write(s)
	for i in range(1,11):
		s='<group name="m'+str(i)+'"><detids val="'+str(i)+'"/></group>\n'
		f.write(s)
	pixelsize=0.640/500.0
	for i in range(500):
		s='<group name="'+str(i)+'"> <detids val="'
		shapexml=virtualslitxml(-0.32+i*pixelsize)
		detlist=FindDetectorsInShape('det',shapexml)
		for j in range(len(detlist)-1):
			s+=str(detlist[j])+','
		s+=str(detlist[len(detlist)-1])+'"/> </group>\n'
		print(s)
		f.write(s)
	s='</detector-grouping>\n'
	f.write(s)
	f.close()

def create1DPSD(rnum,lmin=0.8,lmax=13.0,Mapfile='w:\\Users\Masks\Integrate_tubes_map_file.xml'):
	'''
	Sum the detector in horizontal stripes according to the setup within a mapfile
	'''
	w1=Load(str(rnum))
	calibrateTubes(w1,"")
	w1=ConvertUnits(w1,'Wavelength',AlignBins=1)
	w1=CropWorkspace(w1,lmin,lmax)
	w2=CropWorkspace(w1,StartWorkspaceIndex=0,EndWorkspaceIndex=0)
	detgrp=GroupDetectors(w1,MapFile=Mapfile)
	detgrp=CropWorkspace(detgrp,StartWorkspaceIndex=10)
	w3=detgrp/w2
	w4=Integration(w3,RangeLower=lmin,RangeUpper=lmax)
	x=[]
	y=[]
	e=[]
	pixelsize=0.64/500.0
	for i in range(500):
		x.append(-0.32-0.5*pixelsize+i*pixelsize)
		y.append(w4.readY(i)[0])
		e.append(w4.readE(i)[0])
	x.append(-0.32-0.5*pixelsize+500*pixelsize)
	CreateWorkspace(x,y,e,1,OutputWorkspace=str(rnum)+'_d_h')
	RenameWorkspace(w3,OutputWorkspace=str(rnum)+'_1DPSD')
	return [x,y,e]

def createDBFile(rnumsans,filename,path='W:/Users/Masks/',monitor=0,diagnostics=0):
    # direct beam file is simply the ratio of an upstream monitor to the main detector
    sanswksp=add_runs(rnumsans,'',instrument='LARMOR',keepwksp=1,savewksp=0)
    
    if (len(sanswksp) > 1):
        EBSANS=mtd[sanswksp[0]]
        EBSANS=ConvertUnits(EBSANS,'Wavelength',AlignBins=1)
        EBSANS=Rebin(EBSANS,'0.0,0.04,16.0')
        EBSANSsum=SumSpectra(EBSANS)
        EBSANS=mtd[sanswksp[1]]
        EBSANS=ConvertUnits(EBSANS,'Wavelength',AlignBins=1)
        EBSANS=Rebin(EBSANS,'0.0,0.04,16.0')
        EBSANSm1=SumSpectra(EBSANS,StartWorkspaceIndex=monitor-1,EndWorkspaceIndex=monitor-1)
    else:
        EBSANS=mtd[sanswksp]
        EBSANS=ConvertUnits(EBSANS,'Wavelength',AlignBins=1)
        EBSANS=Rebin(EBSANS,'0.0,0.04,16.0')
        EBSANSsum=SumSpectra(EBSANS,StartWorkspaceIndex=11)
        EBSANSm1=SumSpectra(EBSANS,StartWorkspaceIndex=monitor-1,EndWorkspaceIndex=monitor-1)

    EBSANSnorm=0.1*EBSANSsum/EBSANSm1

    EBSANSnorm_cut=Rebin(EBSANSnorm,'0.9,0.082,12.9')
	
    SaveRKH(EBSANSnorm_cut,path+filename,Append=0)
    RenameWorkspace(EBSANSnorm_cut,OutputWorkspace=filename)
    if diagnostics == 0:
		DeleteWorkspace(EBSANS)
		DeleteWorkspace(EBTRANS)
		DeleteWorkspace(EBSANSsum)
		DeleteWorkspace(EBSANSm1)
		DeleteWorkspace(EBSANSnorm)

#
# Direct Beam file iterator that refines the shape of the direct beam files using a strongly scattering sample with
# not to much incoherent scattering such as a polymer standard
#
def calculateDBFile(rsam,rcan,tsam,tdb,tcan,niterations,filepath,maskfilebase,DBfilebase,wlist=[1.0,1.4,1.8,2.2,3.0,4.0,5.0,7.0,9.0,11.0,12.9 ]):
	now=1
	next=2
	#nIterations=3
	nIterations=niterations
	# generate the extra maskfiles and extra DB files required
	#filepath='W:/Users/Masks/'
	#maskfilebase='MASKLARMOR_142A_LarmorTeam_8mm_AfterNewHome_BR0deg'
	#DBfilebase='DIRECT_787_13August14_1761_corrected'
	for i in range(nIterations):
		shl.copy(filepath+maskfilebase+'.txt',filepath+maskfilebase+str(i+1)+'.txt')
		shl.copy(filepath+DBfilebase+'.dat',filepath+DBfilebase+'_v'+str(i+1)+'.dat')
		for line in finp.input(filepath+maskfilebase+str(i+1)+'.txt',inplace=True):
			print(line.replace('MON/DIRECT='+DBfilebase+'.dat','MON/DIRECT='+DBfilebase+'_v'+str(i+1)+'.dat'))

	#rsam="1761-add"
	#rcan="1757-add"
	#tsam="1762-add"
	#tdb="1758-add"
	#tcan="1758-add"

	i=1
	for ii in range(nIterations):
		LARMOR()
		#Set reduction to 1D (note that if this is left out, 1D is the default)
		Set1D()
		# Read a mask file.   NOTE same direct beam file needs to be mentioned below !
		maskfilename=maskfilebase+str(now)+'.txt'
		# pull in the direct beam that we need to modify MAKE SURE IT IS THE SAME AS IN THE MASK FILE !
		directbeam=LoadRKH(filepath+ DBfilebase+'_v'+str(now)+'.dat')
		newdirectfile=filepath+DBfilebase+'_v'+str(next)+'.dat'
		fitallbkg=0.254
		#
		# 09/03/12     2.375 to 2.625  is empty  direct beam file goes from 1.75 in steps of 0.2 ang
		wlist=[1.0,1.4,1.8,2.2,3.0,4.0,5.0,7.0,9.0,11.0,12.9 ]
		#wlist=[1.0,1.4,2.0,3.5,5.0,7.0,10.0,12.9 ]
		MaskFile(maskfilename)
		#LimitsQ(0.002,1.0,0.08,'LOG')
		W1=1.0
		Q1A=.06
		Q1B=.52
		W2=12.9
		# was .004,to 0.02 but gc drooping near b/s due to multiple scatter??l Try 0.01
		Q2A=0.008
		Q2B=0.03
		# Change default limits if we want to
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
		# this loads the workspaces but does not do the calc, so you have to do one data reduction to get it calculated, then
		# modify it and do the data reduction again
		TransFit('Off',lambdamin=0.9,lambdamax=12.9)
		TransmissionSample(tsam+'.nxs', tdb+'.nxs')
		TransmissionCan(tcan+'.nxs', tdb+'.nxs')
		nloops=len(wlist)-1
		reducedAll = WavRangeReduction(wlist[0], wlist[nloops], False)
		SaveRKH(reducedAll,reducedAll,Append=False)
		shift=-fitallbkg
		Scale(InputWorkspace=reducedAll,OutputWorkspace=reducedAll,Factor=shift,Operation='Add')
		plt = plotSpectrum(reducedAll,0)
		#
		# reduce data
		outlist=[]
		scaleX=[]
		scaleY=[]
		scaleE=[]
		#  must go to unity at extremes of ranges, these get modified again below
		scaleX.append(0.5)
		scaleY.append(1.0)
		scaleE.append(0.0)
		scaleX.append(1.125)
		scaleY.append(1.0)
		scaleE.append(0.0)
		#
		for i in range(nloops):
			LARMOR()
			#Set reduction to 1D (note that if this is left out, 1D is the default)
			Set1D()
			MaskFile(maskfilename)
			# Assign run numbers (.nxs for nexus)
			AssignSample(rsam+'.nxs')
			AssignCan(rcan+'.nxs')
			wav1 = wlist[i]
			wav2 = wlist[i+1]
			TransFit('Off',lambdamin=0.9,lambdamax=12.9)
			TransmissionSample(tsam+'.nxs', tdb+'.nxs')
			TransmissionCan(tcan+'.nxs', tdb+'.nxs')
			#  reduced = WavRangeReduction(wav1, wav2, DefaultTrans, no_clean=True)
			reduced = WavRangeReduction(wav1, wav2, False)
			SaveRKH(reduced,reduced,Append=False)
			# BKG larger at long wavelength, shift be change in bkg
			wav3=(wav1+wav2)*0.5
			# cubic fit to  batch mode flat bkg fits using fixed scale & Rg in SASVIEW, where I(Q=0)=51.922	
			#shift=-0.42121+wav3*0.012356-wav3*wav3*0.0027909 +fitallbkg
			shift=-fitallbkg
			Scale(InputWorkspace=reduced,OutputWorkspace=reduced,Factor=shift,Operation='Add')
			SaveRKH(reduced,reduced+"_shifted",Append=False)
			#    
			outlist.append(reduced)
			mergePlots(plt,plotSpectrum(reduced,0))
			# set up varying Q range to integrate over
			QA=Q2A+(1./wav2-1./W2)*(Q1A-Q2A)/(1./W1-1./W2)
			QB=Q2B+(1./wav1-1./W2)*(Q1B-Q2B)/(1./W1-1./W2)
			# ratio full wavelength data over chosen Q range
			print(wav1,wav2,QA,QB)
			ratio='ratio_'+str(wav1)+'_'+str(wav2)
			# NEW x3 lines, if the division falls over here it is usually beacuse you have empty workspace due to Bragg time masks
			CropWorkspace(reducedAll,OutputWorkspace=ratio,Xmin=QA,Xmax=QB)
			numerator=CropWorkspace(reduced,Xmin=QA,Xmax=QB)
			Divide("numerator",ratio,OutputWorkspace=ratio)
			mergePlots(plt,plotSpectrum(ratio,0))
			#  now integrate the ratio, NOTE it rounds in QA and QB, so to get correct average need to read back the actual range
			#NEW
			Integration(ratio,Outputworkspace="integral",RangeLower=QA,RangeUpper=QB)
			QBB= mtd["integral"].readX(0)[1]
			QAA= mtd["integral"].readX(0)[0]
			norm = mtd["integral"].readY(0)[0]/(QBB-QAA)
			normE = mtd["integral"].readE(0)[0]/(QBB-QAA)
			print(QAA,QBB,norm,normE)
			scaleX.append(wav1)
			int=mtd["integral"].readY(0)[0]
			scaleY.append(norm*(wav2-wav1))
			scaleE.append(normE*(wav2-wav1))

		scaleX.append(wav2)
		#  must go to unity or same as last value at extremes of ranges 
		scaleX.append(19.0)
		scaleY.append(norm*(19.0-wav2))
		scaleE.append(normE*(19.0-wav2))
		scaleX.append(20.0)
		scaleY.append(norm*(20.0-19.0))
		scaleE.append(normE*(20.0-19.0))
		# now change the two dummy points at the start 
		scaleY[0]=scaleY[2]*(1.125-0.5)/(wlist[1]-wlist[0])
		scaleY[1]=scaleY[2]*(wlist[0]-1.125)/(wlist[1]-wlist[0])
		# ======================================================================================
		# need sort some errors for this
		# NEW
		CreateWorkspace(scaleX,scaleY,scaleE,UnitX="Wavelength",OutputWorkspace="scale")
		# small problem, LoadRKH brings in a distribution,, InterpolatingRebin needs a histogram
		# ConvertToHistogram leaves Y values unchanged, just splits the x values, adding one extra x 
		# NEW
		directbeam_hist=ConvertToHistogram("directbeam")
		# interpolating rebin expects a proper histogram wtih counts and time,
		# it seems to need a big over-run at the ends of the range!
		# NEW
		scale2=InterpolatingRebin("scale","0.9,0.082,12.9")
		# ConvertToDistribution only divides by bin width, leave original X values !
		ConvertToDistribution("scale")
		ConvertToDistribution("scale2")
		#CropWorkspace("directbeam_hist","directbeam_hist_crop",XMin=1.0,XMax=15.0)
		#NEW
		scale2=RebinToWorkspace(scale2,"directbeam_hist")
		directbeam_new_hist=Multiply("directbeam_hist","scale2")
		SaveRKH("directbeam_new_hist",newdirectfile,Append=False)
		now=next
		next=next+1
		i=i+1
