# -*- coding: utf-8 -*-
from __future__ import print_function

#from mantid.simpleapi import *
import numpy as np
import copy
from matplotlib import pyplot as plt
    
def findLine(chop_times,chopDist,moderator_limits):
    """
    calculates the lines on the limit of each chopper
    chop_times:a list of the opening and closing times of the chopper within the time frame
    chopDist: a list of the distance from moderator to chopper in meters
    moderator_limits: the earliest and latest times that neutrons can leave the
                    the moderator in microseconds
    """
    lines=[]
    for i in range(len(chop_times)):
        #final chopper openings
        leftM=(0-chopDist)/(moderator_limits[0]-chop_times[i][0])
        rightM=(0-chopDist)/(moderator_limits[1]-chop_times[i][1])
        leftC=0-leftM*moderator_limits[0]
        rightC=0-rightM*moderator_limits[1]
        lines.append([[leftM,leftC],[rightM,rightC]])
    return lines

def checkPath(chop_times, lines, chopDist,chop5Dist):
    """
    A recursive function to check for lines which can satisfy a window in the next chopper
    """
    if len(chop_times)>1:
        #recursive bit
        lines=checkPath(chop_times[1:],lines,chopDist[1:],chop5Dist)
    newLines=[]
    #print "Now looking at chopper ", (5-len(chop_times)), "distance ", chopDist[0]
    for line in lines:
        #for each line check to see if there is an opening in the right time window
        #fast first
        earlyT=(chopDist[0]-line[0][1])/line[0][0]
        #then slow
        lateT=(chopDist[0]-line[1][1])/line[1][0]
        
        #then compare this time window to when this chopper is open, keep the range if it is possible
        for i in range(len(chop_times[0])):
            if (chop_times[0][i][0]<earlyT) and ((chop_times[0][i][1]>lateT)):
                newLines.append(line)#the chopper window is larger than the maximum possible spread, change nothing
            elif (chop_times[0][i][0]>earlyT) and ((chop_times[0][i][1]<lateT)):
                #both are within the window, draw a new box
                chop5_open=(chop5Dist-line[0][1])/line[0][0]
                leftM=(chopDist[0]-chop5Dist)/(chop_times[0][i][0]-chop5_open)
                leftC=chop5Dist-leftM*chop5_open
                chop5_close=(chop5Dist-line[1][1])/line[1][0]
                rightM=(chopDist[0]-chop5Dist)/(chop_times[0][i][1]-chop5_close)
                rightC=chop5Dist-rightM*chop5_close
                newLines.append([[leftM,leftC],[rightM,rightC]])
            elif ((chop_times[0][i][1]< lateT)and(chop_times[0][i][1]>earlyT))and(chop_times[0][i][0]<earlyT):#(chop_times[0][i][1]<lateT) and ((chop_times[0][i][0]<lateT)and(chop_times[0][i][0]>earlyT)):
                #theleftmost range is fine but the rightmost is outside the window. Redefine it
                chop5_close=(chop5Dist-line[1][1])/line[1][0]
                rightM=(chopDist[0]-chop5Dist)/(chop_times[0][i][1]-chop5_close)
                rightC=chop5Dist-rightM*chop5_close
                newLines.append([line[0],[rightM,rightC]])
            elif (chop_times[0][i][1]>lateT)and((chop_times[0][i][0]>earlyT)and(chop_times[0][i][0]<lateT)):
                #the leftmost range is outside the chopper window
                chop5_open=(chop5Dist-line[0][1])/line[0][0]
                leftM=(chopDist[0]-chop5Dist)/(chop_times[0][i][0]-chop5_open)
                leftC=chop5Dist-leftM*chop5_open
                newLines.append([[leftM,leftC],line[1]])

    return newLines
    
def calcEnergy(lines,samDist):
    Ei = np.zeros(len(lines))
    massN=1.674929e-27
    for i in range(len(lines)):
        #look at the middle of the time window
        x1=((samDist-lines[i][0][1])/lines[i][0][0]+(samDist-lines[i][1][1])/lines[i][1][0])/2.
        m=samDist/(x1-lines[i][0][1])
        Ei[i]=(m*1e6)**2*massN/2.*6.242e21#mv^2/2
    return Ei

def calcRes(ei,chop_times,lastChopDist,samDist,detDist,ef=None,ef_type=None):
    #for each incident energy work out the moderator and chopper component of the resolution
    res=[]
    percent=[]
    #IMPORT POINT
    #The chopper opening times are the full opening, for the resolution we want FWHM
    #consequently divide each by a factor of 2 here
    #END IMPORTANT POINT
    chop_width=[(chop_times[0][1]-chop_times[0][0])/2.,(chop_times[1][1]-chop_times[1][0])/2.]
    
    for energy in ei:
        lamba=np.sqrt(81.81/energy)
        #this is the experimentally determined FWHM of moderator
        mod_FWHM=-3.143*lamba**2+49.28*lamba+0.535
        #the effective width at chopper 1
        mod_eff=0.6666*mod_FWHM
        #when running chopper 1 slowly the moderator is smaller than the chopper speed so use that
        if chop_width[0] >mod_eff:
            mod_width=mod_eff
            #print mod_width
        else:
            mod_width=chop_width[0]

        t_mod_chop=252.82*lastChopDist*lamba
        chopRes=(2*chop_width[1]/t_mod_chop)*((detDist+samDist+lastChopDist)/detDist)
        modRes =(2*mod_width/t_mod_chop)*(1+(samDist/detDist))
        res.append(np.sqrt(chopRes**2+modRes**2)*energy)
        percent.append(np.sqrt(chopRes**2+modRes**2))
    return res, percent
    
def plotFrame(lines, chop_times,dist,samDist,DetDist,fracEi,Eis):
    modSamDist=dist[-1]+samDist 
    totDist=modSamDist+DetDist
    for i in range(len(dist)):
        plt.plot([-20000,120000],[dist[i],dist[i]],c='k',linewidth=1.)
        for j in range(len(chop_times[i][:])):
            plt.plot(chop_times[i][j],[dist[i],dist[i]],c='white',linewidth=1.)
    
    plt.plot([-20000,120000],[totDist,totDist],c='k',linewidth=2.)    
    
    for i in range(len(lines)):
        x0=-lines[i][0][1]/lines[i][0][0]
        x1=(modSamDist-lines[i][0][1])/lines[i][0][0]
        plt.plot([x0,x1],[0,modSamDist],c='b')
        x2=(totDist-lines[i][0][1])/lines[i][0][0]
        plt.plot([x1,x2],[modSamDist,totDist],c='b')
        newline=[lines[i][0][0]*np.sqrt(1+fracEi),modSamDist-lines[i][0][0]*np.sqrt(1+fracEi)*x1]
        x3=(totDist-newline[1])/(newline[0])
        plt.plot([x1,x3],[modSamDist,totDist],c='r')
        
        newline=[lines[i][0][0]*np.sqrt(1-fracEi),modSamDist-lines[i][0][0]*np.sqrt(1-fracEi)*x1]
        x4=(totDist-newline[1])/(newline[0])
        plt.plot([x1,x4],[modSamDist,totDist],c='r')
        plt.text(x2,totDist+0.2,"{:3.1f}".format(Eis[i]))
        

    plt.xlabel('TOF ($\mu$sec)')  
    plt.ylabel('Distance (m)') 
    plt.xlim(0,100000)
    plt.show()

def calcFlux(Ei,freq1,percent):
    lamba=np.sqrt(81.81/Ei)
    
    #here are some constants (hahaha) relating to the instrument    
    intRef=0.885 #the flux at 5meV
    freqRef=150. # the frequency this corresponds to
    fluxProf=[0.0889,0.1003,0.1125,0.1213,0.1274,0.1358,0.1455,0.1562,0.1702,
              0.1902,0.2149,0.2496,0.2938,0.3537,0.4315,0.5244,0.6415,0.7856,
              0.9341,1.0551,1.1437,1.1955,1.2004,1.1903,1.1662,1.1428,1.1176,
              1.0875,1.0641,1.0562,1.0242,0.9876,0.9586,0.9415,0.924,0.8856,
              0.8865,0.8727,0.842,0.8125,0.7849,0.7596,0.7417,0.7143,0.6869,
              0.6608,0.6341,0.6073,0.581,0.5548,0.5304,0.507,0.4849,0.4639,
              0.4445,0.425,0.407,0.3902,0.3737,0.3579,0.3427,0.3274,0.3129,
              0.2989,0.2854,0.2724,0.2601,0.2483,0.2371,0.2267,0.2167,0.2072,
              0.1984,0.19,0.1821,0.1743,0.1669,0.1599,0.1532,0.1467,0.1404,
              0.1346,0.1291,0.1238,0.1189,0.1141,0.1097,0.1053,0.1014,0.0975,
              0.0938,0.0902,0.0866,0.0834,0.0801,0.077,0.0741,0.0712,0.0686,
              0.066,0.0637,0.0614,0.0593,0.0571,0.0551,0.0532,0.0512,0.0494,
              0.0477,0.0461,0.0445,0.043,0.0415,0.0401,0.0387]
    fluxLamba=np.linspace(0.5,11.9,num=len(fluxProf))
    flux=[]    
    
    for j in range(len(lamba)):
        i=(abs(fluxLamba-lamba[j])).argmin()
        intensity=fluxProf[i]
        if percent[j]<0.02:
            flux.append(5.6e4*intensity/intRef*(freqRef/freq1)**2)
        else:
            flux.append(5.6e4*intensity/intRef*(freqRef/freq1))
    return flux
    
def multi_rep(efocus, freq1, freqpr, chop2Phase=5, chop3=False):
    """
    A method to calculate the various possible incident energies with a given chopper setup on LET.
    The window of energy transfers plotted is 85% by default.
    efocus: The incident enrgy that all choppers are focussed on
    freq1: The frequency of the resolution choppers
    freqpr: frequency of the pulse removal chopper
    chop2Phase: the second choppers phase, adjustable to take the guessing out
    chop3: A boolean flag to describe if chopper 3 is running or not
    
    Original Matlab code R. Bewley STFC
    Rewritten in Python, D Voneshen STFC 2015
    """
    
    #conversion factors
    lam2TOF=252.82 #the conversion from wavelength to TOF at 1m, multiply by distance
    uSec=1e6 #seconds to microseconds
    #check that the chopperfrequencies make sense
    
    
    if chop3 == False:    
        print("no chop 3")
        dist=[7.83,8.4, 15.66, 23.5] #distance to each chopper in m choppers are in a logical order now
        freq = [freqpr/2., 10.0, freq1/2., freq1]# the frequencies of the choppers
        nslot=[6 ,1,6, 2]# number of slots in each chopper. assume that they are equally spaced
        slot_width = [40,890,52,10] # width of chopper slots in mm
        guide_width = [40,40,40,10] # width of the guide in mm
        radius = [290, 545, 290, 290] # radius in mm of each disk at centre of window
        numDisk = [2, 1, 1, 2]
    else:
        dist=[7.83,8.4,11.75, 15.66, 23.5] #distance to each chopper in m choppers are in a logical order now
        freq = [freq1, 10.0,freqpr, freq1/2., freq1]# the frequencies of the choppers
        #nslot=[6 ,1,2,6, 2]# number of slots in each chopper. assume that they are equally spaced
        nslot=[6 ,1,2,6, 1]
        slot_width = [40,890,56,52,10] # width of chopper slots in mm
        guide_width = [40,40,40,40,10] # width of the guide in mm
        radius = [290, 545, 290,290, 290] # radius in mm of each disk at centre of window
        numDisk = [2, 1,1, 1, 2]
    lam = np.sqrt(81.82/efocus) # convert from energy to wavelenth
    print(freq)
    samp_det = 3.5 #sample to detector distance
    chop_samp = 1.5 # final chopper to sample distance
    source_rep = 10 # rep rate of source
    tmod = 1000 # maximimum emmision window from moderator in us
    frac_ei = 0.90 #fraction of Ei to plot energy loss lines
    
    chop_times=[]#empty list to hold each chopper opening period
    
    #fig=plt.figure()
    
    #first we optimise on the main Ei
    for i in range(len(dist)):
        if i!=1: #chopper2 is a special case
            #loop over each chopper
            t_open=lam2TOF*lam*dist[i] #the opening time of the chopper so that it is open for the focus wavelength
            chopVel = 2*np.pi*radius[i]*numDisk[i]*freq[i] #effective chopper velocity (if 2 disks effective velocity is double)
            t_full_op=uSec*(slot_width[i]+guide_width[i])/chopVel #full opening time
            #t_tran100= uSec*(slot_width[i]-guide_width[i])/chopVel#opening time for 100% transmision
            
            #set the chopper phase to be as close to zero as possible
            next_win_t = uSec/(nslot[i]*freq[i])
            realTimeOp=np.array([(t_open-t_full_op/2.),(t_open+t_full_op/2.)])
            #realTransTime=np.array([(t_open-t_tran100/2.),(t_open+t_tran100/2.)])
        else:
            chopVel = 2*np.pi*radius[i]*numDisk[i]*freq[i] #effective chopper velocity (if 2 disks effective velocity is double)
            t_full_op=uSec*(slot_width[i]+guide_width[i])/chopVel #full opening time
            #t_tran100= uSec*(slot_width[i]-guide_width[i])/chopVel#opening time for 100% transmision
            next_win_t = uSec/(nslot[i]*freq[i])
            realTimeOp=np.array([chop2Phase,chop2Phase+t_full_op])
         #based on these optimsations, work out where all chopper openings are.
        realTimeOp-=next_win_t*np.ceil(realTimeOp[0]/next_win_t)
        #print realTimeOp
        chop_times.append([])
        while realTimeOp[0]<(uSec/source_rep+next_win_t):
            #there should be some plotting here too
            chop_times[i].append(copy.deepcopy(realTimeOp[:]))
            realTimeOp+=next_win_t
        #then we look for what else gets through
    #firstly calculate the bounding box for each window in final chopper
    lines = findLine(chop_times[-1],dist[-1],[0,tmod])
    lines = checkPath(chop_times[0:-1],lines,dist[:-1],dist[-1])
    #ok, now we know the possible neutron velocities. we now ned their energies
    Ei=calcEnergy(lines,(dist[-1]+chop_samp))
    #now calculate the resolution and flux.
    res,percent=calcRes(Ei,[chop_times[0][0],chop_times[-1][0]],dist[-1]-dist[0],chop_samp,samp_det)
    flux=calcFlux(Ei,freq1,percent)
    for idx in range(len(Ei)):
        print("For Ei {:3.2f} meV, the resolution is {:6.2f} ueV or {:3.1f}% and the flux {:1.0f} n/cm^2/s".format(Ei[idx],res[idx]*1000,percent[idx]*100,flux[idx]))
    plotFrame(lines, chop_times,dist,chop_samp,samp_det,frac_ei,Ei)    
    return lines, chop_times
    #return chop_times
    #return 'blah'
    
#multi_rep(6.27,160,80,chop2Phase=5,chop3=False)
