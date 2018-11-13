import os
import numpy as np
import pylab as py
from shutil import copyfile

# there is no need to create phx files as these are saved within the nxspe file

def create_MERInst_Files(run,one2onepar,one2oneold,one2onenew,ringmap):
    one_2one_par_file = create_one2onepar(run,one2onepar)
    create_one2onemap(one2oneold,one2onenew)
    create_ringmap(one_2one_par_file,ringmap)

# script to make one2one.par file - partially using Mantid script
def create_one2onepar(run,one2onepar):
    Load(Filename='MER'+str(run)+'.nxs', OutputWorkspace='w1')
    run_dir = os.path.dirname(os.path.realpath(__file__))
    tmp_file = os.path.join(run_dir,'test.par')
    SavePAR('w1',tmp_file)
    
    one2onepar= os.path.basename(one2onepar)    
    one2onepar = os.path.join(run_dir,one2onepar)    
    
    fid = open(one2onepar, 'w')
    
    data = np.genfromtxt(tmp_file,
                                        names="l2, 2theta, azi, pwid, phigh, crud",
                                        skip_header=1,
                                        dtype=(float, float, float, float, float, float))
                                        
    pixwid=0.0254
    pixhigh=0.0113
    
    fid.write('{0:8.0f}\n'.format(np.size(data)))
    for i in range(np.size(data)):
        fid.write('{0:8.3f} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f}\n'.format(data['l2'][i],data['2theta'][i],data['azi'][i],pixwid,pixhigh))

    fid.close()
    return one2onepar

#script to create rings mapping file. Inputs are the appropriate one2one.par file for the cycle (from which the 2theta values are taken
#and the other argument is the name of the rings.map file you want to create for this cycle
def create_ringmap(one2onepar,ringmap):
    run_dir = os.path.dirname(os.path.realpath(__file__))    
    ringmap= os.path.join(run_dir,os.path.basename(ringmap)) 
    
    fid = open(ringmap,'w')    
    det = np.genfromtxt(one2onepar,
                                    names="l2, 2theta, phi, pwid, phigh",
                                    skip_header=1,
                                    dtype =(float, float, float, float, float))
                                    
    ttheta=np.array(det['2theta'])
    group=0
    numspec_tot=0

    dtheta=0.63
    for angle in py.frange(2.83,136,dtheta):
        myindex=(ttheta>(angle-dtheta/2))*(ttheta<(angle+dtheta/2))
        spectra=np.asarray(np.where(myindex))
        spectra=spectra+1
        numspec=np.shape(spectra)[1]
        if np.shape(spectra)[1]>0:
            group=group+1
    
    fid.write('{0:4.0f}\n'.format(group))
    group=0
    for angle in py.frange(2.83,136,dtheta):
        myindex=(ttheta>(angle-dtheta/2))*(ttheta<(angle+dtheta/2))
        spectra=np.asarray(np.where(myindex))
        spectra=spectra+1
        numspec=np.shape(spectra)[1]
        if np.shape(spectra)[1]>0:
            group=group+1
        fid.write('{0:4.0f}\n'.format(group))
        fid.write('{0:5.0f}\n'.format(np.shape(spectra)[1]))
        for i in range(numspec):
            fid.write('{0:6.0f}\n'.format(spectra[0][i]))
    
    fid.close()
    
# since the one2onemap is just a list of spectra, with each group containing one file, it is unchanging from one cycle to the next. 
# Therefore just copy the previous cycle file across and rename for the new cycle. If the number of spectra changes, obviously this
# will not work.
def create_one2onemap(one2oneold,one2onenew):
    copyfile(one2oneold,one2onenew)

if __name__ == "__main__":	
	create_MERInst_Files(42385,'mynew_one2one.par','one2one_182.map','mynew_one2one.map','mynew_ring.map')    
	#create_one2onemap('one2one_174.map','mynew_one2one.map')
	#create_one2onepar(37394,'mynew_one2one.par')
	#create_ringmap('mynew_one2one.par','mynew_ring.map')    