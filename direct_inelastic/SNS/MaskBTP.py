"""
Mask banks, tubes, and pixels on ARCS/SEQUOIA/CNCS/HYSPEC 
"""
from __future__ import print_function
from numpy import arange
from mantid import *
from mantid.simpleapi import *


def getEightPackHandle(inst,banknum):
	name=inst.getName()
	banknum=int(banknum)
	if name=="ARCS":
		if     (1<=banknum<= 38):
			return inst[3][banknum-1][0]
		elif  (39<=banknum<= 77):
			return inst[4][banknum-39][0]
		elif  (78<=banknum<=115):
			return inst[5][banknum-78][0]
		else: 
			raise ValueError("Out of range index for ARCS instrument bank numbers")
	elif name=="CNCS":
		if     (1<=banknum<= 50):
			return inst[3][banknum-1][0]
		else: 
			raise ValueError("Out of range index for CNCS instrument bank numbers")
	elif name=="HYSPEC":
		if     (1<=banknum<= 20):
			return inst[3][banknum-1][0]	
		else: 
			raise ValueError("Out of range index for HYSPEC instrument bank numbers")
	elif name=="SEQUOIA":
		if     (38<=banknum<= 74):
			return inst[3][banknum-38][0]
		elif  (75<=banknum<= 113):
			return inst[4][banknum-75][0]
		elif  (114<=banknum<=150):
			return inst[5][banknum-114][0]
		else: 
			raise ValueError("Out of range index for SEQUOIA instrument bank numbers")

		
def MaskBTP(**kwargs):

    instrument=kwargs.get('Instrument',None)
    banks=kwargs.get('Bank',None)
    tubes=kwargs.get('Tube',None)
    pixels=kwargs.get('Pixel',None)
    workspace=kwargs.get('Workspace',None)
    detlist=[]
    try:
#    if ((workspace!=None) and (mtd.workspaceExists(workspace) ==True)):
        w=mtd[workspace]
        inst=w.getInstrument()
        instrument =inst.getName()
    except:
        pass
    instrumentList=["ARCS","CNCS","HYSPEC","SEQUOIA"]
    try:
        instrumentList.index(instrument)
    except:
        print("Instrument not found")
        return detlist
    if (workspace==None):
        path=config["instrumentDefinition.directory"]
        LoadEmptyInstrument(Filename=path+instrument+"_Definition.xml",OutputWorkspace="temporaryWorkspaceForMasking")
        workspace="temporaryWorkspaceForMasking"
        w=mtd[workspace]
        inst=w.getInstrument()
    if (banks==None):
        if (instrument=="ARCS"):
            banks=arange(115)+1
        elif (instrument=="CNCS"):
            banks=arange(50)+1
        elif (instrument=="HYSPEC"):
            banks=arange(20)+1
        elif (instrument=="SEQUOIA"):
            banks=arange(113)+38
    else:
        # try to get the bank numbers in an array, even if the banks is string, array, or an integer
        banks=eval(str(banks))
        try:
            len(banks)
        except:
            banks=[banks]
    if(tubes==None):
        tubes=arange(8)+1
    else:
        tubes=eval(str(tubes))
        try:
            len(tubes)
        except:
            tubes=[tubes]
    if(pixels==None):
        pixels=arange(128)+1
    else:
        pixels=eval(str(pixels))
        try:
            len(pixels)
        except:
            pixels=[pixels]	
    for b in banks:
        ep=getEightPackHandle(inst,b)
        for t in tubes:
            if ((t<1) or (t>8)):
                raise ValueError("Out of range index for tube number")
            else:
                for p in pixels:
                    if ((p<1) or (p>128)):
                        raise ValueError("Out of range index for pixel number")
                    else:
                        pid=ep[int(t-1)][int(p-1)].getID()
                        detlist.append(pid)
    MaskDetectors(Workspace=workspace,DetectorList=detlist)
    return detlist

#MaskBTP(Instrument="SEQUOIA",Pixel="1,2,3,4,5,6,7,8,121,122,123,124,125,126,127,128")
#MaskBTP(Workspace="temporaryWorkspaceForMasking",Bank=110)
