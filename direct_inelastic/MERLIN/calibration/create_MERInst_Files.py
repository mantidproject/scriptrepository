import os
import numpy as np
import pylab as py
from shutil import copyfile
from sys import platform
from datetime import date
import subprocess

# there is no need to create phx files as these are saved within the nxspe file

def create_MERInst_Files(run,commit_changes,repository_path,cycle_index,det_dat_file,one2one_old_map):
    ''' create all Merlin calibration files and commit them into script repository if the repository is configured properly
    Inputs:
    run                      -- run number used for tube calibration
    commit_changes    -- if true, changes should be at least attempted to commit to the repository
    repository_path     -- the path where the repository with the instrument data are located.
    cycle_index           -- the string, which describes the current cycle and is used as the tag for marking
                                  the produced instrument files.
    one2one_old_map  -- old one-2-one map file, used as source for one-2-one file for current cycle.
    det_dat_file          -- the file the calibrated detectors, to commit to InstrumentRepository. empty if no commit is necessary.
    '''
    # working directory where the files shluld be saved
    working_dir = config['defaultsave.directory'];
    
    #
    print ('*************************************************************')
    print ('*** Creating MERLIN instrument files within working directory: {0}'.format(working_dir))    
    print ('*************************************************************')
    new_instrument_files_list = []
    #-----------------------------------------------------------------------------------------------------------------------
    # one2one. par file:
    par_file_template = 'one2one_';
    one2onepar = par_file_template+cycle_index+'.par'
    full_one2one_par_file = create_one2onepar(run,one2onepar) 
    print ('*** Have generated one2one par file: {0}'.format(full_one2one_par_file));
    new_instrument_files_list.append(full_one2one_par_file)
    #-----------------------------------------------------------------------------------------------------------------------    
    #rings map file:
    rings_file_template = 'ring_'
    ringmap = rings_file_template + cycle_index+'.map'

    ring_map_full_file = create_ringmap(full_one2one_par_file,ringmap) 
    print ('*** Have generated ring map file: {0}'.format(ring_map_full_file));    
    new_instrument_files_list.append(ring_map_full_file)    
    #-----------------------------------------------------------------------------------------------------------------------     
    # one2one map file:
    one2one_map_template = 'one2one_'
    one2one_map      = one2one_map_template+cycle_index+'.map'
    
    one2one_old_full_map = os.path.join(working_dir,one2one_old_map)
    if os.path.exists(one2one_old_full_map):           
       one2one_full_map = create_one2onemap(one2one_old_full_map,one2one_map)
       new_instrument_files_list.append(one2one_full_map)
    else:
        one2one_old_full_map = os.path.join(repository_path,one2one_old_map)        
        if  not os.path.exists(one2one_old_full_map) :
            print  '*** --> Can not generate one2one map file: {0} as the source file is not available'.format(one2one_map)
            print  '*** --> Copy source file: {0} into working directory: {1}'.format(one2one_old_map,working_dir)
            one2one_full_map = ''
        else:
            one2one_full_map = create_one2onemap(one2one_old_full_map,one2one_map)            
    if len(one2one_full_map)>0:
        print  '*** Have generated one2one map: {0}'.format(one2one_full_map)
        new_instrument_files_list.append(one2one_full_map)        
    #----------------------------------------------------------------------------------------------------------------------- 
    print ('*************************************************************')     
    print ('*** File generation completed  and new files are placed in folder: {0}'.format(working_dir))
    print ('*************************************************************')    
    if commit_changes:
        print '*** Adding new files to the repository in svn folder: {0}'.format(repository_path)
        # check if detector.dat file exist in working directory
        det_dat = os.path.join(working_dir,det_dat_file)
        if os.path.exists(det_dat):
            print '*** Preparing to commit detector calibration file: {0}'.format(det_dat_file)
            new_instrument_files_list.append(one2one_full_map)
        else:
            print '*** Detector calibration file: {0} does not exist in {1} so not committing it'.format(det_dat_file,working_dir) 
        commit_files_to_svn(repository_path,new_instrument_files_list,cycle_index)
    else:
        print '*** Changes will not be commited as folder: {0} does not exist'.format(repository_path)        
        
def commit_files_to_svn(svn_directory_path,fileslist,cycle_tag='current'):
     ''' routine copies the list of input files  into svn directory, 
      adds them to the svn repository and 
      commits them if proper commit righs are stored on the machine'''
     if not isinstance(fileslist,list):
          fileslist = [fileslist]
     n_added_files = 0
     
     subprocess.check_output(['svn','update'],cwd = svn_directory_path)
     
     for file in fileslist:
          fdir,fname = os.path.split(file)
          targ = os.path.join(svn_directory_path,fname)
          if os.path.exists(targ):
            copyfile(file,targ) 
          else:
             copyfile(file,targ) 
             n_added_files += 1            
             out = subprocess.check_output(['svn','add',targ],cwd = svn_directory_path,stderr=subprocess.STDOUT,shell=True)
             if len(out)>0:
                 print out
     
     if n_added_files > 0:
         print '*** Added {0} new files to svn script repository'.format(n_added_files)
     print '*** Updating svn repository with {0} existing and {1} new generated MERLIN insrument files'.format(len(fileslist),n_added_files)
         
     message = '-m\" Merlin autocommitted instrument files for cycle: '+cycle_tag+'\"'
     command = 'svn commit '+message+'; exit 0'
     out= subprocess.check_output(command ,cwd = svn_directory_path,stderr=subprocess.STDOUT,shell=True)
     print out
      

# script to make one2one.par file - partially using Mantid script
def create_one2onepar(run,one2onepar):
    Load(Filename='MER'+str(run)+'.nxs', OutputWorkspace='w1')
    #######run_dir = os.path.dirname(os.path.realpath(__file__))
      
    
    run_dir =  config['defaultsave.directory']
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
    run_dir = config['defaultsave.directory']
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
    return ringmap
    
def check_calibration_options(cycle_end_list):
    ''' Routine calculates the calibration options available as the function of existing machine setup.
      Inputs: 
        cycle_end_list -- the dictionary, defining the cycle index as the function of the cycle end date
       Ouput:
       repository_path -- the path to InstrumentFiles repository as the 
    '''
    # Check if we can immidiately commit our files to remote repository
    if platform == 'win32': # windows
        repository_path = r'c:\mprogs\ISIS\InstrumentFiles\merlin'
    else:                           # only linux is supported
       repository_path = '/usr/local/mprogs/InstrumentFiles/merlin'
    if os.path.exists(repository_path):
        commit_changes = True
    else:
        commit_changes = False
        repository_path = config['defaultsave.directory']
     # identify what cycle we are calibrating data for:
    cal_date = date.today()
    print ('*************************************************************')    
    print ('*** running MERLIN instrument files generation script')        
    print  '*** Current date: ',cal_date.isoformat()
    cycle_index = 'future'
    will_commit = 'will not be committed'
    # identify current cycle, when we are running the calibration
    for end_date,pref in sorted(cycle_end_list.items()):
         if cal_date < end_date : # we now in the cycle requested
            cycle_index = pref
            break
    if cycle_index == 'future':
         commit_changes = False

    print '*** for cycle: {0}'.format(cycle_index)  
    print('*** Instrument files will be build in the folder {0}'.format(config['defaultsave.directory']))                    
    if commit_changes:
         print('*** and will be committed to the repository at {0}'.format(repository_path))        
    else:
         print('*** but will not be committed to a repository')

    return (repository_path,commit_changes,cycle_index)

    
# since the one2onemap is just a list of spectra, with each group containing one file, it is unchanging from one cycle to the next. 
# Therefore just copy the previous cycle file across and rename for the new cycle. If the number of spectra changes, obviously this
# will not work.
def create_one2onemap(one2oneold,one2onenew):
    run_dir =  config['defaultsave.directory']
    new_one2one_full = os.path.join(run_dir,one2onenew)    
    copyfile(one2oneold,new_one2one_full)
    return new_one2one_full

if __name__ == "__main__":
     # list of the known cycles with approximate  end dates:
    cycle_end_list = {date(2019,03,29):'184',date(2019,07,19):'191',date(2019,10,25):'192',\
    date(2019,12,20):'193',date(2999,01,01):'future'}
    
    # old one2one map file used as source for one2one map file used in current cycle
    one2one_source_map = 'one2one_182.map' 
    # det.dat file with neutronic positions of the detectors, calculated using calibration scripts
    calibrated_det_dat     = 'det_corr_191_process_5.dat'
    
    # verify what calibration options are available as the function of the current machine settings
    repository_path,commit_changes,cycle_index = check_calibration_options(cycle_end_list)
    
      
    create_MERInst_Files(42385,commit_changes,repository_path,cycle_index,calibrated_det_dat,one2one_source_map)
    #create_one2onemap('one2one_174.map','mynew_one2one.map')
    #create_one2onepar(37394,'mynew_one2one.par')
    #create_ringmap('mynew_one2one.par','mynew_ring.map')
