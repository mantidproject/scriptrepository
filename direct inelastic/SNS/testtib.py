"""
Test if time independent background range is OK
Note - runs in MantidPlot
Note - uses pdftk to concatenate pdf files
"""
from numpy import *
from string import *
import os

def E2V(E):
# for energy in mev returns velocity in m/s
    return sqrt(E/5.227e-6)


def SpurionPromptPulse2(Ei, msd = 1800.0, tail_length_us = 3000.0, talk = False):
    #More sophisticated
    dist_mm = 39000.0 + msd + 4500.0
#    T0_moderator = 4.0 + 107.0 / (1.0 + (self._Ei / 31.0)*(self._Ei / 31.0)*(self._Ei / 31.0))
    T0_moderator = 0.0 
    t_focEle_us = 39000.0 / E2V(Ei) * 1000.0 + T0_moderator
    t_samp_us = (dist_mm - 4500.0) /E2V(Ei) * 1000.0 + T0_moderator
    t_det_us = dist_mm / E2V(Ei) * 1000 + T0_moderator
    frame_start_us = t_det_us - 16667/2
    frame_end_us = t_det_us + 16667/2
    index_under_frame = divide(int(t_det_us),16667)
    pre_lead_us = 16667 * index_under_frame
    pre_tail_us = pre_lead_us + tail_length_us
    post_lead_us = 16667 * (1+ index_under_frame)
    post_tail_us = post_lead_us + tail_length_us
    E_final_meV = -1
    E_transfer_meV = -1
    # finding an ok TIB range
    MinTIB_us = 2000.0
    slop_frac = 0.2
    #print t_focEle_us,pre_lead_us,frame_start_us,MinTIB_us,slop_frac
    if (t_focEle_us < pre_lead_us) and (t_focEle_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before focus element-1'
        TIB_high_us = t_focEle_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif (frame_start_us>pre_tail_us) and (t_focEle_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before focus element-2'
        TIB_high_us = t_focEle_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif t_focEle_us-pre_tail_us > MinTIB_us * (slop_frac + 1.0) and (t_focEle_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before focus element-3'
        TIB_high_us = t_focEle_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif t_samp_us-pre_tail_us > MinTIB_us * (slop_frac + 1.0) and (t_samp_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before sample-1'
        TIB_high_us = t_samp_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif t_samp_us-pre_tail_us > MinTIB_us / 1.5 * (slop_frac + 1.0) and (t_samp_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before sample-2'
        TIB_high_us = t_samp_us - MinTIB_us / 1.5 * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us / 1.5
    elif t_samp_us-pre_tail_us > MinTIB_us / 2.0 * (slop_frac + 1.0) and (t_samp_us-frame_start_us > MinTIB_us * (slop_frac + 1.0)):
        if talk:
            print 'choosing TIB just before sample-3'
        TIB_high_us = t_samp_us - MinTIB_us / 2.0 * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us / 2.0
    elif (pre_lead_us - frame_start_us > MinTIB_us * (slop_frac + 1.0)) and (t_focEle_us > pre_lead_us):
        if talk:
            print 'choosing TIB just before leading edge before elastic-1'
        TIB_high_us = pre_lead_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif (pre_lead_us - frame_start_us > MinTIB_us / 1.5 * (slop_frac + 1.0)) and (t_focEle_us > pre_lead_us):
        if talk:
            print 'choosing TIB just before leading edge before elastic-2'
        TIB_high_us = pre_lead_us - MinTIB_us / 1.5 * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us / 1.5
    elif (pre_lead_us - frame_start_us > MinTIB_us / 2.0 * (slop_frac + 1.0)) and (t_focEle_us > pre_lead_us):
        if talk:
            print 'choosing TIB just before leading edge before elastic-3'
        TIB_high_us = pre_lead_us - MinTIB_us / 2.0 * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us / 2.0
#    elif (pre_tail_us > frame_start_us) and (t_focEle_us - pre_tail_us > MinTIB_us * (slop_frac + 1.0)):
#        if talk:
#            print 'choosing TIB just before focus element'
#            print pre_tail_us, MinTIB_us, slop_frac
#        TIB_low_us = pre_tail_us + MinTIB_us * slop_frac / 2.0
#        TIB_high_us = TIB_low_us + MinTIB_us
    elif post_lead_us > frame_end_us:
        if talk:
            print 'choosing TIB at end of frame'
        TIB_high_us = frame_end_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    elif post_lead_us - t_det_us > MinTIB_us * (slop_frac + 1.0):
        if talk:
            print 'choosing TIB between elastic peak and later prompt pulse leading edge'
        TIB_high_us = post_lead_us - MinTIB_us * slop_frac / 2.0
        TIB_low_us = TIB_high_us - MinTIB_us
    else:
        if talk:
            print 'I cannot find a good TIB range'
        TIB_low_us = 0.0
        TIB_high_us = 0.0
    return [TIB_low_us, TIB_high_us]

runs=range(12312,12348)
runs.extend(range(12408,12434))

od='/SNS/users/3y9/Desktop/test/'
com='pdftk '
for run in runs:
	data=Load(Filename='/SNS/HYS/IPTS-8016/data/HYS_'+str(run)+'_event.nxs',OutputWorkspace='HYS_12312_event')
	sum=SumSpectra(data)
	Rebin(InputWorkspace='sum',OutputWorkspace='sum',Params='10')
	Ei=sum.run().getProperty('EnergyRequest').value[0]
	CreateWorkspace(OutputWorkspace='TIB',DataX=SpurionPromptPulse2(Ei),DataY='2,2',DataE='1,1')
	pl=plot(('sum','TIB'),0)
	l=pl.activeLayer()
	l.setScale(0,1,1000,1,5,5,1)
	l.setCurveLineWidth(0,4)
	l.setTitle("Ei = "+str(Ei)+" meV")
	filename=od+str(run)+'.pdf'
	l.exportVector(filename)
	pl.confirmClose(False)
	pl.close()
	com+=filename+' '
com+='cat output '+od+'tib.pdf'
os.system(com)
