from mantid.simpleapi import *
from shutil import copyfile
import nxs as nxs
import os

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

def add_runs(runlist,pathout,instrument='LARMOR'):
	if(os.path.isdir(pathout)=='False'):
		return
	pfix=instrument
	runlist=parseRunList(runlist)
	runlist=runlit[0]
	added=Load(runlist[0])
	print "          uampHr = ", added.getRun().getProtonCharge()
	
	if(len(b)>1):
		for i in range(1,len(b)):
			wtemp=Load(runlist[i])
			print "          uampHr = ", wtemp.getRun().getProtonCharge()
			added=added+wtemp
		DeleteWorkspace("wtemp")
	
	nzeros=8-len(runlist[0])
	fpad=""
	for ii in range(nzeros):
		fpad+="0"
	print "writing file:   "+pathout+pfix+fpad+runlist[0]+"-add.nxs"
	SaveNexusProcessed("added",pathout+pfix+fpad+runlist[0]+"-add.nxs")  
	print "    total uampHr = ", added.getRun().getProtonCharge()
	DeleteWorkspace("added")

	