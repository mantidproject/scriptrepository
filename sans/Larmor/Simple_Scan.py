from __future__ import print_function
from matplotlib.pyplot import errorbar, plot, show
import numpy as np
import SANSroutines as lm
#import LSS.SANSroutines as lm
import matplotlib.pyplot as mpl
import genie_python.genie as gen
import math
from mantidplot import gui_cmd
from time import sleep

gen.set_instrument("IN:LARMOR")

#from genie_python.genie import *
def setupup():
    flipper1(1)
    flipper2(0)

def flipper1(state=0):
    if state == 0:
        #non flipping state
        set_pv("IN:LARMOR:KEPCO50_01:CURRENT:SP", 0.75)
    else:
        set_pv("IN:LARMOR:KEPCO50_01:CURRENT:SP", -0.75)
    sleep(3)        

def flipper2(state=0):
    if state == 0:
        #non flipping state
        set_pv("IN:LARMOR:KEPCO50_02:CURRENT:SP", 2.0)
    else:
        set_pv("IN:LARMOR:KEPCO50_02:CURRENT:SP", -2.0)
    sleep(3)

def cset_str(blockString, value):
        dict = {blockString: value}
        gen.cset(**dict)
    
def scan_axis(axis,startval,endval,npoints,frms,rtitle,usem4=0):
    lm.setuplarmor_nrscanning()

    gen.change(title=rtitle)
    gen.change(nperiods=npoints)
    
    gen.begin(paused=1)
    # setup the scan arrays and figure
    xval=np.zeros(npoints)
    yval=np.zeros(npoints)
    eval=np.zeros(npoints)

    stepsize=(endval-startval)/float(npoints-1)
    for i in range(npoints):
        xval[i]=(startval+i*stepsize)

    mpl.ion()
    fig1=mpl.figure(1)
    mpl.clf()
    ax = mpl.subplot(111)
    #ax.set_xlim((0,4))
    ax.set_xlabel(axis)
    ax.set_ylabel('Normalised Neutron counts')
    # reasonable x-Axis, necessary to get the full window from the first datapoint
    scanrange = np.absolute(endval - startval)
    mpl.xlim((startval-scanrange*0.05, endval+scanrange*0.05))
    mpl.draw()
    mpl.pause(0.001)
    
    for i in range(npoints):
        gen.change(period=i+1)
        cset_str(axis,xval[i])
        gen.waitfor(seconds=1)
        gen.waitfor_move()
        gfrm=gen.get_frames()
        gen.resume()
        gen.waitfor(frames=gfrm+frms)
        gen.pause()
        a1=gen.get_spectrum(1,i+1)
        msig=sum(a1['signal'])*100.0
        mesig=(math.sqrt(msig))
        print("msig="+str(msig)+" mesig="+str(mesig))
        # get the interesting monitor
        if usem4 < 1:
            a1=gen.get_spectrum(11,i+1)
            sig=sum(a1['signal'])*100.0
            a1=gen.get_spectrum(12,i+1)
            sig+=sum(a1['signal'])*100.0
            esig=math.sqrt(sig)
        else:
            a1=gen.get_spectrum(4,i+1)
            sig=sum(a1['signal'])*100.0
            esig=math.sqrt(sig)
        print("sig="+str(sig)+" esig="+str(esig))
        yval[i]=(sig/msig)
        eval[i]=(math.sqrt((sig/(msig*msig))+(sig*sig/(msig*msig*msig))))
        print("yval="+str(yval[i])+" esig="+str(eval[i]))
        ax.errorbar(xval[i], yval[i], eval[i], fmt = 'ko')
        fig1.canvas.draw()
        mpl.pause(0.001)
        f.open('u:/users/Larmor/lastscan.csv','w')
        s=str(xval[i])+','+str(yval[i])+','+str(eval[i])+'\n'
        f.write(s)
        f.close()
    gen.abort()
    #f.open('u:/users/Larmor/lastscan.csv','w')
    #for i in range(npoints):
    #    s=str(xval[i])+','+str(yval[i])+','+str(eval[i])+'\n'
    #    f.write(s)
    #f.close()
    '''
    for i in range(npoints):
        # get the upstream monitor
        a1=get_spectrum(i+1,1)
        msig=sum(a1['signal'])
        mesig=(sqrt(msig))
        print "msig="+str(msig)+" mesig="+str(mesig)
        # get the interesting monitor
        a1=get_spectrum(i+1,4)
        sig=sum(a1['signal'])
        esig=sqrt(sig)
        print "sig="+str(sig)+" esig="+str(esig)
        yval[i]=(sig/msig)
        eval[i]=(sqrt((sig/(msig*msig))+(sig*sig/(msig*msig*msig))))
        print "yval="+str(yval[i])+" esig="+str(eval[i])
    
    '''    
def scan_axis_mantid(axis,startval,endval,npoints,frms,rtitle,usem4=0):
    lm.setuplarmor_nrscanning()

    gen.change(title=rtitle)
    gen.change(nperiods=npoints)
    
    gen.begin(paused=1)
    # setup the scan arrays and figure
    xval=np.zeros(npoints)
    yval=np.zeros(npoints)
    eval=np.zeros(npoints)

    stepsize=(endval-startval)/float(npoints-1)
    for i in range(npoints):
        xval[i]=(startval+i*stepsize)

    gui_cmd(mpl.ion)
    fig1=gui_cmd(mpl.figure,1)
    gui_cmd(mpl.clf)
    ax=gui_cmd(mpl.subplot,111)
    #ax.set_xlim((0,4))
    gui_cmd(ax.set_xlabel,axis)
    gui_cmd(ax.set_ylabel,'Normalised Neutron counts')
    # reasonable x-Axis, necessary to get the full window from the first datapoint
    scanrange = np.absolute(endval - startval)
    gui_cmd(mpl.xlim,(startval-scanrange*0.05, endval+scanrange*0.05))
    gui_cmd(mpl.draw)
    gui_cmd(mpl.pause,0.001)
    
    for i in range(npoints):
        gen.change(period=i+1)
        cset_str(axis,xval[i])
        sleep(15)
        #gen.waitfor_move()
        
        gfrm=gen.get_frames()
        gen.resume()
        gen.waitfor(frames=gfrm+frms)
        gen.pause()
        a1=gen.get_spectrum(1,i+1)
        msig=sum(a1['signal'])*100.0
        mesig=(math.sqrt(msig))
        print("msig="+str(msig)+" mesig="+str(mesig))
        # get the interesting monitor
        if usem4 < 1:
            a1=gen.get_spectrum(11,i+1)
            sig=sum(a1['signal'])*100.0
            a1=gen.get_spectrum(12,i+1)
            sig+=sum(a1['signal'])*100.0
            esig=math.sqrt(sig)
        else:
            a1=gen.get_spectrum(4,i+1)
            sig=sum(a1['signal'])*100.0
            esig=math.sqrt(sig)
        print("sig="+str(sig)+" esig="+str(esig))
        yval[i]=(sig/msig)
        eval[i]=(math.sqrt((sig/(msig*msig))+(sig*sig/(msig*msig*msig))))
        print("yval="+str(yval[i])+" esig="+str(eval[i]))
        gui_cmd(ax.errorbar,xval[i], yval[i], eval[i], fmt = 'ko')
        gui_cmd(fig1.canvas.draw)
        gui_cmd(mpl.pause,0.001)
    
    gen.abort()

    
def scanloop():
    for i in range(1000):
        scan_axis(77,117,41,500,"500 frames 80tubes Scan Height 77-117 41 steps 3d printed aperture")
        scan_axis(117,77,41,500,"500 frames 80tubes Scan Height 117-77 41 steps 3d printed aperture ")
        
        
def polscan_axis(axis,startval,endval,npoints,frms,rtitle):
    lm.setuplarmor_nrscanning()

    gen.change(title=rtitle)
    gen.change(nperiods=npoints*2)
    
    gen.begin(paused=1)
    # setup the scan arrays and figure
    xval=np.zeros(npoints)
    yval=np.zeros(npoints)
    eval=np.zeros(npoints)

    stepsize=(endval-startval)/float(npoints-1)
    for i in range(npoints):
        xval[i]=(startval+i*stepsize)

    mpl.ion()
    fig1=mpl.figure(1)
    mpl.clf()
    ax = mpl.subplot(111)
    #ax.set_xlim((0,4))
    ax.set_xlabel(axis)
    ax.set_ylabel('Normalised Neutron counts')
    # reasonable x-Axis, necessary to get the full window from the first datapoint
    scanrange = np.absolute(endval - startval)
    mpl.xlim((startval-scanrange*0.05, endval+scanrange*0.05))
    mpl.draw()
    mpl.pause(0.001)
    flipper1(1)
    
    for i in range(npoints):
        gen.change(period=(i*2)+1)
        cset_str(axis,xval[i])
        flipper2(0)
        gen.waitfor_move()
        gfrm=gen.get_frames()
        resume()
        gen.waitfor(frames=gfrm+frms)
        pause()
        flipper2(1)
        gen.change(period=(i*2)+2)
        gfrm=gen.get_frames()
        resume()
        gen.waitfor(frames=gfrm+frms)
        pause()

        a1=gen.get_spectrum(1,(i*2)+1)
        msigup=sum(a1['signal'])*100.0
        mesigup=(sqrt(msigup))
        # get the interesting monitor
        a1=gen.get_spectrum(11,(i*2)+1)
        sigup=sum(a1['signal'])*100.0
        a1=gen.get_spectrum(12,(i*2)+1)
        sigup+=sum(a1['signal'])*100.0
        esigup=sqrt(sigup)

        a1=gen.get_spectrum(1,(i*2)+2)
        msigdo=sum(a1['signal'])*100.0
        mesigdo=(sqrt(msigdo))
        # get the interesting monitor
        a1=gen.get_spectrum(11,(i*2)+2)
        sigdo=sum(a1['signal'])*100.0
        a1=gen.get_spectrum(12,(i*2)+2)
        sigdo+=sum(a1['signal'])*100.0
        esigdo=sqrt(sigdo)
        
        yval[i]=(sigup-sigdo)/(sigup+sigdo)
        eval[i]=yval[i]*1e-3
        #eval[i]=(sqrt((sig/(msig*msig))+(sig*sig/(msig*msig*msig))))
        ax.errorbar(xval[i], yval[i], eval[i], fmt = 'ko')
        fig1.canvas.draw()
        mpl.pause(0.001)
    
    abort()
        