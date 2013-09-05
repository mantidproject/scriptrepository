"""
Make some plots using matplotlib
"""
from mantid.simpleapi import*
outdir='/SNS/SEQ/IPTS-7672/shared/19g/tof_figs/'
indir='/SNS/SEQ/IPTS-7672/data/'
runs=range(33012,33193)

def plot_sum(x,y,err,titlestr='',xstr='',ystr=''):
	h_f=figure()
	errorbar(x,y,yerr=err,fmt='bo')
	title(titlestr)
	xlabel(xstr)
	ylabel(ystr)
	h_ax=gca()
	h_ax.set_yscale("log",nonposy='clip')

for runnum in runs:
	flname='SEQ_'+str(runnum)+'_event.nxs'
	Load(Filename=indir+flname,OutputWorkspace='OWS',LoadMonitors='1')
	Rebin(InputWorkspace='OWS',OutputWorkspace='OWS_h',Params='6000,10,16000',PreserveEvents='0')
	SumSpectra(InputWorkspace='OWS_h',OutputWorkspace='sum')
	h_sum=mtd['sum']
	tvec=h_sum.dataX(0)
	tvec=(tvec[:-1]+tvec[1:])/2.
	Ivec=h_sum.dataY(0)
	errvec=h_sum.dataE(0)
	plot_sum(tvec,Ivec,errvec,titlestr=str(runnum),xstr='TOF',ystr='I')
	savefig(outdir+'run'+str(runnum)+'.png')
