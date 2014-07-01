from mantid.simpleapi import *
from shutil import copyfile
from math import *
import nxs as nxs
import os
from ISISCommandInterface import *

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
						r1=range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1]))
						for k in r1:
							rlist2.append(str(k))
					elif tstr[j].find('-') >=0:
						tstr[j].strip()
						tstr2=tstr[j].split('-')
						r1=range(int(tstr2[0]),int(tstr2[1])+1)
						for k in r1:
							rlist2.append(str(k))
					else:
						rlist2.append(tstr[j])
			else:
				if rlist1[i].find(':') >=0 and rlist1[i].find('-')>=0:
					rlist1[i].strip()
					tstr2=rlist1[i].split('-')
					tstr3=tstr2[1].split(':')
					r1=range(int(tstr2[0]),int(tstr3[0])+1,int(tstr3[1]))
					for k in r1:
						rlist2.append(str(k))
				elif rlist1[i].find('-')>=0:
					rlist1[i].strip()
					tstr2=rlist1[i].split('-')
					r1=range(int(tstr2[0]),int(tstr2[1])+1)
					for k in r1:
						rlist2.append(str(k))
				else:
					rlist2.append(rlist1[i])
			rlist.append(rlist2)
	return rlist

'''
SANS add runs from RKH
'''

def add_runs(runlist,pathout,instrument='LARMOR',keepwksp=0,savewksp=1):
	if(os.path.isdir(pathout)=='False'):
		return
	pfix=instrument
	runlist=parseRunList(runlist)
	runlist=runlist[0]
	added=Load(runlist[0])
	print "          uampHr = ", added.getRun().getProtonCharge()
	
	if(len(runlist)>1):
		for i in range(1,len(runlist)):
			wtemp=Load(runlist[i])
			print "          uampHr = ", wtemp.getRun().getProtonCharge()
			added=added+wtemp
		DeleteWorkspace("wtemp")
	
	if savewksp==1:
		nzeros=8-len(runlist[0])
		fpad=""
		for ii in range(nzeros):
			fpad+="0"
		print "writing file:   "+pathout+pfix+fpad+runlist[0]+"-add.nxs"
		SaveNexusProcessed("added",pathout+pfix+fpad+runlist[0]+"-add.nxs")  
		print "    total uampHr = ", added.getRun().getProtonCharge()
	if keepwksp==0:
		DeleteWorkspace("added")
	else:
		RenameWorkspace("added",OutputWorkspace=runlist[0]+'-add')
	return mtd[runlist[0]+'-add']

def larmor1D(rsample,rcan,tsample,tdb,tcan,maskfile,wkspname,lmin=0.9,lmax=12.5):
	LARMOR()
	#Set reduction to 1D (note that if this is left out, 1D is the default)
	Set1D()
	MaskFile(maskfile)
	# Assign run numbers (.nxs for nexus)
	AssignSample(rsample+'.nxs')
	AssignCan(rcan+'.nxs')
	TransFit('Off',lambdamin=lmin,lambdamax=lmax)
	TransmissionSample(tsample+'.nxs', tdb+'.nxs')
	TransmissionCan(tcan+'.nxs', tdb+'.nxs')
	reduced = WavRangeReduction(lmin, lmax, False)
	if(len(wkspname) > 0):
		a1=rsample.split('-')
		wksproot=a1[0]
		RenameWorkspace(wksproot+'rear_1D_'+str(lmin)+'_'+str(lmax),OutputWorkspace=wkspname)

def calibrateTubes(wkspName,calibrationfile):
   #
   # pixel by pixel efficiency correction for the linear detector
   #
   flood_wksp = "Vanadium_tube_calib_1to1_30062014"
   if  flood_wksp not in mtd:
		#Load("Vanadium_tube_calib_1to1_May14.nxs",OutputWorkspace="Vanadium_tube_calib_1to1_May14")
		#Load("Vanadium_tube_calib_1to1_25062014.nxs",OutputWorkspace="Vanadium_tube_calib_1to1_25062014")
		Load("Vanadium_tube_calib_1to1_30062014.nxs",OutputWorkspace="Vanadium_tube_calib_1to1_30062014")

   CopyInstrumentParameters("Vanadium_tube_calib_1to1_30062014",wkspName)

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
	e.g. createMapFile(1203,'w:\Users\Masks\Integrate_tubes_map_file.xml')
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
		print s
		f.write(s)
	s='</detector-grouping>\n'
	f.write(s)
	f.close()

def create1DPSD(rnum,lmin=0.8,lmax=13.0,Mapfile='w:\Users\Masks\Integrate_tubes_map_file.xml'):
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
