from __future__ import (absolute_import, division, print_function, unicode_literals)
from mantid.simpleapi import *
import numpy as np
import scipy
from scipy.optimize import fmin
from scipy.special import erfc
import matplotlib.pyplot as plt
import warnings
import time

start = time.time()

# script to calibrate the detector tubes on Merlin

def tube_calibrate_MER(run,tmin,tmax,*args):    
    #*args takes the form of 'bank' 'pack' 'tube'. If no arguments are given the whole detector array is fitted.
    #Load(Filename='MER'+str(run)+'.s02', OutputWorkspace='w1')
    Load(Filename='MER'+str(run)+'.raw', OutputWorkspace='w1')

    Rebin(InputWorkspace='w1', OutputWorkspace='w1', Params='1500,100,9000')
    mylen = len(args)

    
    all_ids = []
    w1 = mtd['w1']
    spec_num=np.zeros(w1.getNumberHistograms())
    for i in range(w1.getNumberHistograms()):
        myspectrum = w1.getSpectrum(i)
        spec_num[i]=myspectrum.getSpectrumNo()
        all_ids.extend(list(myspectrum.getDetectorIDs()))
    
    det_tubes = np.array(all_ids)/10000
    packstart = [1,1,1,1,1,1,1,1,1]
    packend = [4,4,5,4,4,4,4,3,3]
    tubestart = 1
    tubeend = 8
    
    print('File loaded')
    
    if mylen==0:
        doorstart = 1
        doorend = 9
        cal_info = 'doors-1_9'
    if mylen==1:
        bank = args[0]
        doorstart = bank
        doorend = bank
        cal_info = 'door-{0}'.format(bank)
    if mylen==2:
        bank = args[0]
        pack = args[1]
        doorstart = bank
        doorend = bank
        packstart[bank-1] = pack
        packend[bank-1] = pack
        cal_info = 'door-{0}_pack-{1}'.format(bank,pack)
    if mylen==3:
        bank = args[0]
        pack = args[1]
        tube = args[2]
        doorstart = bank
        doorend = bank
        packstart[bank-1] = pack
        packend[bank-1] = pack
        tubestart = tube
        tubeend = tube
        cal_info = 'door-{0}_pack-{1}_tube-{2}'.format(bank,pack,tube)
        print(doorstart)
        print(packstart)
        print(tubestart)

    run_dir =  config['defaultsave.directory']
    targ_file_name =  'calibration_Helen_192_{0}.csv'.format(cal_info)
    
    file = os.path.join(run_dir,targ_file_name )
    fid = open(file,'w')
    fid.write(' Tube_id, start, peak1, peak2, peak3, peak4, peak5, end\n')

    for bank in range(doorstart,doorend+1):
        for pack in range(packstart[np.array(bank)-1],packend[np.array(bank)-1]+1):
            for tube in range(tubestart,tubeend+1):
                tube_id = int(str(bank)+str(pack)+str(tube))
                myindex = np.nonzero(det_tubes.astype(int)==tube_id)
                spec_tube = np.array(spec_num)[[myindex[0]]]
                spec_min = int(min(spec_tube))
                spec_max = int(max(spec_tube))
                
                Rebin(InputWorkspace='w1', OutputWorkspace='w2', Params='1500,7500,9000')
                ExtractSpectra(InputWorkspace='w2',OutputWorkspace='w3',StartWorkspaceIndex=spec_min-1,EndWorkspaceIndex=spec_max-1)
                
                w3=mtd['w3']
                yval = w3.extractY()
                left_end = 0
                middle = 0
                right_end = 0
                if max(yval)!=0:
                    myout = myfit_data(bank,pack,tube,yval,mylen)
                if myout is not None:
                    print(tube_id,myout[0],myout[1],myout[2],myout[3],myout[4],myout[5],myout[6])
                    fid.write('{0:3.0f} {1:0.2f} {2:0.2f} {3:0.2f} {4:0.2f} {5:0.2f} {6:0.2f} {7:0.2f}\n'.format(tube_id,myout[0],myout[1],myout[2],myout[3],myout[4],myout[5],myout[6]))
                else:
                   plt.show(1)
                   print(' could not fit tube: ',tube_id)                    
                   fid.write('{0:3.0f} {1:0.2f} {2:0.2f} {3:0.2f} {4:0.2f} {5:0.2f} {6:0.2f} {7:0.2f}\n'.format(tube_id,0,0,0,0,0,0,0))                    
                    #raise RuntimeError('Could not fit')
    fid.close()
    print('***********************************************')
    print('*** Calibration info is written in file: {0}'.format(targ_file_name));
    print('*** Located in folder: {0}'.format(run_dir));    
    print('***********************************************')    
    
    
    

 
def myfit_data(bank,pack,tube,Intensity,mylen):
    error = np.sqrt(Intensity)
    position = range(1,513)
    left_end = 0
    middle = 0
    right_end = 0
    
    fac = 1

#fitting left hand end

    if np.logical_and(bank==3,pack==1)==False:
        if np.logical_and(bank==3,pack==2):
            fac=2
        maxint = max(np.array(Intensity[fac*20:fac*67]))
        if maxint==0:   #if there are no counts in the tube
            left_end = 0
            middle = 0
            right_end = 0
            return
        maxim = []
        maxim = np.where(np.array(Intensity)==maxint)[0]
        maxim = maxim+5
        if len(maxim)>1:
            maxim = maxim[np.where(np.logical_and(maxim>(fac*20)+1,maxim<(fac*67)+1))]
            maxim = maxim[0]
        else:
            maxim = maxim[0]
        if maxim>fac*66:
            left_end = 0
            middle = 0
            right_end = 0
            return
        minint = min(np.array(Intensity[maxim:fac*67]))
        minim = []
        minim = np.where(np.array(Intensity)==minint)[0]
        if len(minim)>1:
            minim = minim[np.where(np.logical_and(minim>maxim,minim<(fac*67)+1))]
            minim = minim[0]
        
        if minim<=maxim:    #spike at end of tube
            left_end = 0
            middle = 0
            right_end = 0
            return
            
        pos = position[maxim-10:minim]
        pos = np.reshape(np.array(pos),(-1,1))
        I = Intensity[maxim-10:minim]
        dI = error[maxim-10:minim]
        if mylen==3:
            fig=plt.figure()
            plt.get_current_fig_manager().window.setGeometry(5,5,1700,1000)
            s1=plt.subplot(2,4,1)
            s2=plt.subplot(2,4,2)
            s3=plt.subplot(2,4,3)
            s4=plt.subplot(2,4,4)
            s5=plt.subplot(2,4,5)
            s6=plt.subplot(2,4,6)
            s7=plt.subplot(2,4,7)
            plt.subplot(241)
            plt.errorbar(pos,I,yerr=dI,fmt='o')

    
        midint=(maxint-minint)/2.0
        v = abs((I/midint)-1.0)
        myindex = np.where(v==min(v))[0]
        pos_midint=pos[myindex[0]]
    
    #define starting parameters
        p0 = [-maxint/2, pos_midint, 6, 0]
    
    #error function to minimize
        eef = lambda p, pos, I: ((abs(endf(p,pos)-I))**2).sum()    
    
    # fitting the data with fmin
        p, fopt, iter, funcallas, warnflag = fmin(eef, p0, args=(pos,I),maxiter=2000,ftol=1e-06,maxfun=5000,full_output=True)
        if np.logical_or(warnflag==1,warnflag==2):
            left_end = 0
            middle = 0
            right_end = 0
            return
    
        left_end = p[1]
        print(left_end)
        if mylen==3:    
            myx = np.arange(min(pos),max(pos),0.5)
            myy = endf(p,myx)
            plt.plot(myx,myy)
            plt.draw()

#fitting right hand end
    if np.logical_and(bank==3,pack==2)==False:
        if np.logical_and(bank==3,pack==1):
            fac=2
        maxint = max(np.array(Intensity[512-fac*60:512-fac*20]))
        if maxint==0:   #if there are no counts in the tube
            left_end = 0
            middle = 0
            right_end = 0
            return
        maxim = []
        maxim = np.where(np.array(Intensity)==maxint)[0]
        maxim = maxim+5
        if len(maxim)>1:
            maxim = maxim[np.where(np.logical_and(maxim>(512-fac*60)+1,maxim<(512-fac*20)+1))]
            maxim = maxim[0]
        if maxim > 490:
            maxim = 490
        

        minint = min(np.array(Intensity[(512-fac*64):maxim]))
        minim = []
        minim = np.where(np.array(Intensity)==minint)[0]
        if len(minim)>1:
            minim = minim[np.where(np.logical_and(minim<maxim,minim>(512-fac*70)))]
            minim = minim[0]

        if maxim<=minim:    #spike at end of tube
            left_end = 0
            middle = 0
            right_end = 0
            return
    
        pos = position[minim:maxim]
        pos = np.reshape(np.array(pos),(-1,1))
        I = Intensity[minim:maxim]
        dI = error[minim:maxim]
        if mylen==3:
            plt.subplot(242)
            plt.errorbar(pos,I,yerr=dI,fmt='o')

    
        midint=(maxint-minint)/2.0
        v = abs((I/midint)-1.0)
        myindex = np.where(v==min(v))[0]
        pos_midint=pos[myindex[0]]
    
    #define starting parameters
        p0 = [maxint/2, pos_midint, 6, 0]
    
        eef = lambda p, pos, I: ((abs(endf(p,pos)-I))**2).sum()    
    
    # fitting the data with fmin
        p, fopt, iter, funcallas, warnflag = fmin(eef, p0, args=(pos,I),maxiter=2000,ftol=1e-06,maxfun=5000,full_output=True)
        if np.logical_or(warnflag==1,warnflag==2):
            left_end = 0
            middle = 0
            right_end = 0,
            return

        right_end = p[1]
        print(right_end)
        if mylen==3:
            myx = np.arange(min(pos),max(pos),0.5)
            myy = endf(p,myx)
            plt.plot(myx,myy)
            plt.draw()


#fitting stripes
    stripes = [] #create empty list
    if np.logical_and(bank==3,pack==1): #dealing with half length tubes
        vals=np.array([85,115,240,270])
        addval=3
    elif np.logical_and(bank==3,pack==2): 
        vals=np.array([310,340,390,420])
        addval=6
    else: #dealing with standard tubes
        vals=np.array([135,165,170,200,229,269,320,350,390,420])
        addval=3
    for i in range(0,int(len(vals)/2)):
        maxint = max(np.array(Intensity[vals[2*i]:vals[2*i+1]]))
        if maxint==0:   #if there are no counts in the tube
            left_end = 0
            middle = 0
            right_end = 0
            return
        maxim = []
        maxim = np.where(np.array(Intensity)==maxint)[0]
        if len(maxim)>1:
            maxim = maxim[np.where(np.logical_and(maxim>vals[2*i]-1,maxim<vals[2*i+1]+1))]
            maxim = maxim[0]
        
        lowx = maxim-11
        highx = maxim+11
        peak = position[maxim]
        pos = position[lowx:highx]
        pos = np.reshape(np.array(pos),(-1,1))
        I = Intensity[lowx:highx]
        dI = error[lowx:highx]
        if mylen==3:
            plt.subplot(2,4,i+addval)
            plt.errorbar(pos,I,yerr=dI,fmt='o')

    
        #define starting parameters
        p0 = [max(I), pos.mean(), 3.0, 1.0, (I[0]-I[len(pos)-1])/pos[0]-pos[len(pos)-1]]
        
        emf = lambda p, pos, I: ((abs(midf(p,pos)-I))**2).sum()    
        
        # fitting the data with fmin
        p, fopt, iter, funcallas, warnflag = fmin(emf, p0, args=(pos,I),maxiter=1000,ftol=1e-06,maxfun=5000,full_output=True)
        if np.logical_or(warnflag==1,warnflag==2):
            if np.logical_or(np.logical_and(bank!=3,pack!=1),np.logical_and(bank!=3,pack!=2)):
                left_end = 0
                middle = 0
                right_end = 0
                return
    
        if mylen==3:
            myx = np.arange(min(pos),max(pos),0.5)
            myy = midf(p,myx)
            plt.plot(myx,myy)
            plt.show()
            plt.draw()
        stripe1 = p[1]
        print(stripe1)
        stripes.append(stripe1)

    if np.logical_and(bank==3,pack==1): 
        return (0,0,0,0,stripes[0],stripes[1],right_end)
    elif np.logical_and(bank==3,pack==2): 
        return (left_end,stripes[0],stripes[1],0,0,0,0)
    elif np.logical_and(np.logical_and(bank==3,pack==3),1<=tube<=3): 
        return (left_end,stripes[0],stripes[1],256,stripes[3],stripes[4],right_end)
    else:
        return (left_end,stripes[0],stripes[1],stripes[2],stripes[3],stripes[4],right_end)
    
    
def endf(c,x):          #error function to fit end of tubes
    if c[3]<0:
        c[3]=0
    
    if c[0]>0:
        return c[0]*erfc((c[1]-x)/c[2]) + c[3]
    if c[0]<0:
        return -2*c[0] + c[0]*erfc((c[1]-x)/c[2]) + c[3]
        
def midf(c,x):          #gaussian function to fit stripes
    
    return c[3] + c[0]*np.exp(-(x-c[1])**2/(2*c[2]*c[2]))

if __name__ == "__main__":
        
    #####################################################################
    #This is the line to actually run the script
    tube_calibrate_MER(46470,1000,9000)   #In this example an optional argument is given to just look at door 3.
#####################################################################
