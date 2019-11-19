from __future__ import print_function
from genie_python.genie import *
import numpy as np
import matplotlib.pyplot as plt

import sys,os

def groupname(name=''):  # Currently Broken. Unsure Why.
    '''groupname() looks for a group name in the user directory 
    example : groupname('Fred')

    If the name is not found you are prompted to creat it. If it exist then
    it is made the current durectory 
    '''
    if name=='':
        # Ask for a name
        name=eval(input('What is your group name? '))
    if os.path.isdir('u:/%s'%name):
        os.chdir('u:/%s'%name)
    else:
        answer=eval(input(r'Group directory %s not found in u:\. Make a new one? '%name))
        if answer == 'y' or answer == 'Y':
            print('u:/%s'%name)
            #os.mkdir('u:/%s'%name)            

def waitformove():
    '''waitformove() waits till the table of motors reports that all motion is stopped
    '''
    waitfor_boolean(r'C:\LabVIEW Modules\Drivers\Galil DMC2280\Galil - System Functions.llb\Galil - Table of motor details.vi','In Motion',target=False)
#######################
## Polarization 
#######################
def flipper(state):
    '''flipper('up') or flipper('dn') 
    '''
    FlipperCurrent = 1.0
    
    PolarizationState={'up':0, 'dn':1, 'down':1}
        
    if(PolarizationState[state] == PolarizationState['up']): # UP    
        cset(flipper2i =  FlipperCurrent)
    elif(PolarizationState[state] == PolarizationState['dn']): # Down
        cset(flipper2i = -FlipperCurrent)
    
def pol_run(u=5000, d=5000, total=0, Title='', SampleName=''):
    '''pol_run is called in the following way
    pol_run(u=5000 d=5000 total=180000,Title='', SampleName='')
    '''
    polstate=[u,d]
    # Prep the run
    if Title=='' and SampleName=='':
        change(nperiods=2)
    else:
         change(nperiods=2,title=Title, sample_name=SampleName)
    framecount=0
    begin(paused=True)
    while framecount < total:
        # Loop over the polarization states
        for i in range(2):
            if i == 0 :
                change(period=1)
                flipper('up')
                print('counting on spin up for %d frames'%polstate[i], 'Framecount / Total ', framecount,'/ ', total)
            if i == 1 :
                change(period=2)
                flipper('dn')
                print('counting on spin dn for %d frames'%polstate[i], 'Framecount / Total ', framecount,'/ ', total)
        
            print("A small pause to stop cross contamination of polarisation states:")
            waitfor(seconds=3)  # CJK 11/10/2013 wait to stop cross contamination of pol states
            resume()
            waitfor(frames=(polstate[i] + framecount))
            framecount = framecount + polstate[i]  
            pause()
    end()
#######################
## scan
#######################
#def ascan(Block, Start, Stop, Npts, Frames, Io=3, I=4,Save=False):
def ascan(*args, **kwargs):

    BlockName = None
    Start=None
    Stop=None
    Npts=21
    Frames=10
    Io=3
    I=4
    Save=False
    
    # unpack ordered arguments
    nargs=len(args)
    if nargs==1:
        BlockName = args
    elif nargs==2:
        BlockName, Start = args 
    elif nargs==3:
        BlockName, Start, Stop = args
    elif nargs==4:
        BlockName, Start, Stop, Npts = args
    elif nargs==5:
        BlockName, Start, Stop, Npts, Frames = args
    elif nargs==6:
        BlockName, Start, Stop, Npts, Frames, Io = args
    elif nargs==7:
        BlockName, Start, Stop, Npts, Frames, Io, I = args
    elif nargs==8:
        BlockName, Start, Stop, Npts, Frames, Io, I, Save = args
    elif nargs >8:
        print('Too many arguments')
        return
               
    # Unpack any keyword arguments
    for k,v in kwargs.items():
        if k == 'Block':
            BlockName=v
        elif k == 'Start':
            Start=v
        elif k == 'Stop':
            Stop=v
        elif k == 'Npts':
            Npts=v
        elif k == 'Frames':
            Frames=v
        elif k == 'Io':
            Io=v
        elif k == 'I':
            I=v
        elif k == 'Save':
            Save=v
                
    print('#', BlockName, Start, Stop, Npts, Frames, Io, I, Save)
    
    # build array of axis values
    x_values = np.linspace(Start, Stop, Npts)
    y_values = x_values*0.0
    e_values = y_values*0.0

    #print np.shape(x_values), np.shape(y_values), np.shape(e_values)

    plt.clf()
    
    line, ( bottoms, tops), verts = plt.errorbar(x_values, y_values, yerr=e_values, fmt='o')

    CurrentState=get_runstate()
    if CurrentState[0] !=1 : # Here 1 = Setup
        print('Cant start a scan if instrument is "Running".')
        print(' Will abort current run in 30 sec. Use ctrl-c to stop action.')
        waitfor(seconds=30)
        abort()
    change(nperiods=Npts)
    begin(paused=True)
    print(('#%s\t\tI/Io\t\tErr\t\tI\t\tIo'%BlockName))
    for i,x in enumerate(x_values):
        #cset(stheta=x)
        exec('cset(%s=x)'%BlockName)
        waitformove()
        change(period=i+1)
        resume()
        waitfor(frames=Frames*(i+1))
        pause()
        # extract data from the dae
        spectraI=get_spectrum(I,i,dist=False)# dist=True => counts / time
        spectraIo=get_spectrum(Io,i,dist=False)# dist=True => counts / time
        CountsInI = np.array(spectraI['signal']).sum()
        CountsInIo = np.array(spectraIo['signal']).sum()
        if CountsInI > 0.0:
            y_values[i]=CountsInI/CountsInIo
            e_values[i]=abs(y_values[i])*np.sqrt( (np.sqrt(CountsInI)/CountsInI)**2 + (np.sqrt(CountsInIo)/CountsInIo)**2 )
        else:
            y_values[i]=0.0
            e_values[i]=0.0
        print(('%10.3f\t%10.3g\t%10.3g\t%10.3g\t%10.3g'%(x_values[i], y_values[i], e_values[i], CountsInI, CountsInIo)))
        # update the plot with the data
        line.set_ydata(y_values)
        bottoms.set_ydata(y_values - 0.5*e_values)
        tops.set_ydata(y_values + 0.5*e_values)
        #verts.set_ydata(e_values)
        plt.gca().relim()
        plt.gca().autoscale_view()
        plt.draw()
        #plt.flush_events()
    if Save==True:
        end()
    else:
        abort()
    return x_values, y_values, e_values
        
   
def scan(*args, **kwargs):
    '''scan(): scan an axis defined in a block to produce a Axis vs I/Io graph.
       This function is most often used to align the sample to the rest of the
       instrument.
       Examples:

       scan('stheta', start, stop, npts, frames, Io, I)
       
    '''

    xdata,ydata,edata = ascan(*args, **kwargs)  # Run the scan 

    ## Fit the data to a peak 
    from scipy.optimize import curve_fit
    # Define model function to be used to fit to the data above:
    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
    p0 = [ydata.max(),xdata.mean(), 1.0]# initial guess for the fitting coefficients (A, mu and sigma above)
    coeff, var_matrix = curve_fit(gauss, xdata, ydata, p0=p0)
    
    xscanlength=xdata.max()-xdata.min()
    xfit=np.linspace(xdata.min()-.1*xscanlength,xdata.max()+0.1*xscanlength)
    yfit = gauss(xfit, *coeff)

    plt.figure(1)
    plt.clf()
    plt.errorbar(xdata,ydata,yerr=edata,fmt='bo')
    plt.plot(xfit,yfit,'r-')
    plt.title('$A$=%g, $x_0$=%.5f, $\sigma$=%.5f'%(coeff[0],coeff[1],coeff[2]))
    plt.xlabel('%s'%args[0]); plt.ylabel('$I/I_0$')
    
    
         
#######################
## SE stuff
#######################
            
def pump(A=25, B=25, C=25, D=25,Flow=5, Wait=180):
    '''pump(A=25, B=25, C=25, D=25,Flow=5, Wait=180):
    '''
    if A+B+C+D != 100:
        print("Pump: Please make your channel values add up to 100")
        return
    if Wait == 180:
        print("Pump: Using default wait of 180 s")

    cset(concA=A, concB=B, concC=C, concD=D, hplcflow=Flow)
    waitfor(seconds=2)
    cset(hplcset=1)
    cset(pump_on_off=1)# turn the pump on
    print("Waiting for "+str(wait)+"s")
    waitfor(seconds=Wait)
    cset(pump_on_off=0)     # turn the pump off
