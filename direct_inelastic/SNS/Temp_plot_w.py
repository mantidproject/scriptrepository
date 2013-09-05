"""
Temperature scan on a peak, using a mask
"""
from mantid.simpleapi import *
from numpy import *
import matplotlib as mplt

def ShiftTime(WName,lg_name):
	"""
	shift the time in a given log to match the time in the proton charge log"
	"""
	H_IN = mtd[WName]
	PC =  H_IN.getRun()['proton_charge'].firstTime()
	#print "P="+str(PC)+"\n"
	P =  H_IN.getRun()[lg_name].firstTime()
	#print "P="+str(P)+"\n"
	Tdiff = PC-P
	Tdiff_num = Tdiff.total_milliseconds()*1E-3
	#print "Tdiff="+str(Tdiff_num)+"\n"
	ChangeLogTime(InputWorkspace=WName, OutputWorkspace = WName, LogName = lg_name, TimeOffset = Tdiff_num)


IPTS_num='IPTS-6225'
run_num=30390

dat={}
Tstep=-1.0
Ts=arange(325,275,Tstep)
dcolors=256/len(Ts)
fname_str='/SNS/SEQ/%s/data/SEQ_%d_event.nxs'%(IPTS_num,run_num)
h_ws=Load(Filename=fname_str,OutputWorkspace='h_ws',LoadMonitors='1')
ShiftTime('h_ws','SampleTemp')
ConvertUnits(InputWorkspace='h_ws',OutputWorkspace='d_event',Target='dSpacing')
LoadMask(Instrument='SEQUOIA',InputFile=r'/SNS/SEQ/IPTS-6225/shared/19g/half_hlaf_0_group.xml',OutputWorkspace='roi')
#LoadMask(Instrument='SEQUOIA',InputFile=r'/SNS/SEQ/IPTS-6225/shared/19g/half_hlaf_0_bkg_group.xml',OutputWorkspace='roi2')
MaskDetectors(Workspace='d_event',MaskedWorkspace='roi')

#rgbvals=mplt.cm.jet(range(0,256,dcolors))
figure()
hold(True)
ydat_w=zeros(len(Ts))
for tidx,temp in enumerate(Ts):
     Tmin=temp-abs(Tstep/2.0)
     Tmax=temp+abs(Tstep/2.0)
     FilterByLogValue(InputWorkspace='d_event',OutputWorkspace='T_step',LogName='SampleTemp',MinimumValue=str(Tmin),MaximumValue=str(Tmax),LogBoundary='Left')
     Rebin(InputWorkspace='T_step',OutputWorkspace='T_step',Params='1,0.005,10',PreserveEvents='0')
     try:
       NormaliseByCurrent(InputWorkspace='T_step',OutputWorkspace='T_step')
       SumSpectra(InputWorkspace='T_step',OutputWorkspace='all_d',IncludeMonitors='0')
       Integration(InputWorkspace='all_d',OutputWorkspace='I',RangeLower='4.7999999999999998',RangeUpper='5.7999999999999998')
       h_I=mtd['I']
       ydat_w[tidx]=h_I.readY(0)
     except:
       pass
   
plot(Ts,ydatw,'bo')#color=tuple(rgbvals[tidx,:]))
show()

