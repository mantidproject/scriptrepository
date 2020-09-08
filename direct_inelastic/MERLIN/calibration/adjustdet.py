from __future__ import print_function
# import mantid algorithms
from mantid.simpleapi import *

import os
import numpy as np
import scipy
from scipy.optimize import fmin
from scipy.special import erfc
import matplotlib.pyplot as plt
import itertools as iter


#script to take fitted values and make adjustments to det_corr file
def adjust_detector_MER(fit_res,det_source_file):
    # Inputs:
    # fit_res --  the file produced by tubecalib.py file (in csv format)
    # det_source_file --  the ASCII file containing calibration information from the previous cycle.
    # Outputs: 
    # The file with name det_source_file +'_corrected.dat' containing the calibration information
    #
    # both source fit_res file and det_corr_file should be located in the script directory.
    #    
    fit_res = os.path.basename(fit_res)
    det_source_file  = os.path.basename(det_source_file)
    run_dir =  config['defaultsave.directory']
    fit_res = os.path.join(run_dir,fit_res)
    full_det_source_file = os.path.join(run_dir,det_source_file )

    finall_corr_file =   os.path.splitext(det_source_file)
    finall_corr_file  = finall_corr_file[0]+'_corrected.dat'
    finall_corr_file_fp = os.path.join(run_dir,finall_corr_file)
    
    tub_fit = np.genfromtxt(fit_res,
                                    names="Tube_id, tube_start, peak1, peak2, peak3, peak4, peak5, tube_end",
                                    dtype = (int, float, float, float, float, float, float, float),
                                    skip_header=True,delimiter=',')
    det = np.genfromtxt(full_det_source_file,
                                    names="detno, offset, l2, code, 2theta, phi, w_x, w_y, w_z, f_x, f_y, f_z, a_x, a_y, a_z, det1, det2, det3, det4",
                                    skip_header=12,
                                    dtype =(int, float, float, int, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float))
    print(' Fit tube shape: {0}'.format(tub_fit.shape))
    print(' Det shape: {0}'.format(det.shape))
    #take first 12 lines from original det_corr file and copy across into the new version
    with open(full_det_source_file) as myfile:
        head = list(iter.islice(myfile,12))

    fid = open(finall_corr_file_fp,'w')
    for item in head:
        fid.write(item)
    
    tube = np.floor(det['detno']/10000)
    num_tubes = tub_fit.size
    degrad_conv = 180/np.pi
    
    tube_len_short = 1.234              #active length of short tubes
    tube_len_long = 2.9                   #active length of long tubes 
    low_short_tube_pos = -1.5         #offset of the tubes below the beam stop
    high_short_tube_pos = 0.17      #offset of the tubes above the beam stop
    nDet = 512
    xbin_short = tube_len_short/nDet
    xbin_long = tube_len_long/nDet
    
    true_pos = [-0.48, -0.21, -0.15, 0, 0.15, 0.28, 0.48]       #fractional positions of stripes taken from calibration bin
    true_pos = np.array(true_pos)*tube_len_long
    
    for tub_num in range(0,num_tubes):
        tube_id = tub_fit['Tube_id'][tub_num]
        fit_par = [tub_fit['tube_start'][tub_num],tub_fit['peak1'][tub_num],tub_fit['peak2'][tub_num],tub_fit['peak3'][tub_num],tub_fit['peak4'][tub_num],
        tub_fit['peak5'][tub_num],tub_fit['tube_end'][tub_num]]

        if int(tube_id/10)==31:
            tube_type=1             #short upper tube
        elif int(tube_id/10)==32:
            tube_type=2             #short lower tube
        else:
            tube_type=0             #long tube
        
        index = np.nonzero(tube==tube_id)
        dist_old = det['l2'][index]
        th_old = det['2theta'][index]
        phi_old = det['phi'][index]
        [x_old,z_old,y_old] = mysph2cart(phi_old/degrad_conv,(th_old-90)/degrad_conv,dist_old)
        
        if tube_type==0:
            fit_pos = ((np.array(fit_par)-256)/512)*tube_len_long
            true_pos = [-0.48, -0.21, -0.15, 0, 0.15, 0.28, 0.48]       #fractional positions of stripes taken from calibration bin
            true_pos = np.array(true_pos)*tube_len_long
            ideal_det_pos = -0.5*tube_len_long + 0.5*xbin_long + np.array(list(range(0,nDet)))*xbin_long
        elif tube_type==1:
            fit_par = [tub_fit['peak4'][tub_num],tub_fit['peak5'][tub_num],tub_fit['tube_end'][tub_num]]
            true_pos = [0.15, 0.28, 0.48]       #fractional positions of stripes taken from calibration bin
            true_pos = np.array(true_pos)*tube_len_long 
            fit_pos = high_short_tube_pos + ((np.array(fit_par))/512)*tube_len_short
            ideal_det_pos = high_short_tube_pos + 0.5*xbin_short + np.array(list(range(0,nDet)))*xbin_short
        elif tube_type==2:
            fit_par = [tub_fit['tube_start'][tub_num],tub_fit['peak1'][tub_num],tub_fit['peak2'][tub_num]]
            true_pos = [-0.48, -0.21, -0.15]       #fractional positions of stripes taken from calibration bin
            true_pos = np.array(true_pos)*tube_len_long
            fit_pos = low_short_tube_pos + ((np.array(fit_par))/512)*tube_len_short
            ideal_det_pos = low_short_tube_pos + 0.5*xbin_short + np.array(list(range(0,nDet)))*xbin_short
        else:
            error('Adjust Detectors - Invalid argument')
        
        print('Tube N',tub_num,'With ID: ',tube_id,' has peaks fitted at: ',fit_par)
        print('True peaks positions: ',fit_pos)
        print('Stripes positions: '  ,true_pos)
        
        z_new = tube_correction(true_pos,fit_pos,ideal_det_pos)
        #
        # is this right that z is on place of Y coordinate?
        phi_new, theta_new, dist_new = mycart2sph(x_old,z_new,y_old)
        phi_new = phi_new*degrad_conv
        theta_new = theta_new*degrad_conv + 90
        
        
        det['l2'][index] = dist_new
        det['2theta'][index] = theta_new
        det['phi'][index] = phi_new
        
        for i in range(0,nDet):
            fid.write('{0:7.0f} {1:8.1f} {2:10.5f} {3:5.0f} {4:11.5f} {5:11.5f} {6:11.5f} {7:11.5f} {8:11.5f} {9:11.5f} {10:11.5f} {11:11.5f} {12:11.5f}' 
            '{13:11.5f} {14:11.5f} {15:3.0f} {16:3.0f} {17:3.0f} {18:3.0f}\n'
            .format(det['detno'][index][i], det['offset'][index][i], det['l2'][index][i], det['code'][index][i], det['2theta'][index][i], det['phi'][index][i], det['w_x'][index][i], det['w_y'][index][i], 
            det['w_z'][index][i],det['f_x'][index][i], det['f_y'][index][i], det['f_z'][index][i], det['a_x'][index][i], det['a_y'][index][i], det['a_z'][index][i], det['det1'][index][i], 
            det['det2'][index][i], det['det3'][index][i], det['det4'][index][i]))
    #------------------------------------------------------------------------------------
    print ('*************************************************************')
    print ('*** Detectors have been calibrated                        ***')    
    print('*** The name of the file with calibration information is: {0}'.format(finall_corr_file));
    print('*** The file is located in folder: {0}'.format(run_dir))
    print ('*************************************************************')
    
def mysph2cart(az,elev,r):
    z = r*np.sin(elev)
    x = r*np.cos(elev)*np.cos(az)
    y = r*np.cos(elev)*np.sin(az)
    return x, y, z
    
def mycart2sph(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)    
    elev = np.arctan2(z,np.sqrt(x**2+y**2))
    az = np.arctan2(y,x)
    return az, elev, r
    
def tube_correction(true_pos,meas_pos,ideal_det_pos):
    p = np.polyfit(meas_pos,true_pos,3)
    z_new = np.poly1d(p)(ideal_det_pos)
    return z_new

if __name__ == "__main__"  or __name__ == "__builtin__" or __name__ == "mantidqt.widgets.codeeditor.execution":
   # Inputs:
   # 1 the file produced by tubecalib.py file (in csv format)
   # 2 the file containing calibration information from the previous cycle.
   adjust_detector_MER('calibratioon_res_doors-1_9.csv','det_corr_194_process_5.dat')
   #Outputs:
   # The procedure writes calibrated file with the name {2}_corrected.dat