#--------------------------------------------------------------------
#                             anvred3.py
#
#
# Includes the option to perform a spherical absorption correction
# or calculate direction cosines for a polyhedral correction
# using the gaussian.f program.
#
# Includes a Tkinter gui interface.
#
#     A. J. Schultz, May 2014
#
#--------------------------------------------------------------------
#                             anvred2x.py
#
# Does not run from Isaw. User input read from anvred2x.inp file.
#   A. J. Schultz, February 2012
#--------------------------------------------------------------------
#
# Data reduction program:
#   Input is raw integrated intensities.
#   Output is relative Fsq's.
#
# Jython version:
#    A. J. Schultz, started December 2009
#
# anvred_py.py
#    Each spectrum is a separate file. The spectra files are created
#    by "TOPAZ_spectrum_multiple_banks.iss".
#
# anvred2.py:
#    This version reads one spectrum file containing spectra for
#    each detector. The spectra are created by "TOPAZ_spectrum.py".
#
# Modfications by Xiaoping Wang, April 2011
# Added Selection of neutron wavelengths limits wlMin, wlMax
# Omit zero intensity peaks in integrate file XP Wang 03/21/2011
# Changed to >=0 and used absolute value for minium I/sing(I) = 0  XP Wang 02/24/2011
# Added detector scale factors for vanadium/niobium spectrum XP Wang 09/24/2013
#
#
# 
# Comments from Fortran source:
# C**************************   ANVRED  ******************************
# C
# C ARGONNE NATIONAL LABORATORY VARIABLE WAVELENGTH DATA REDUCTION PROGRAM
# C
# C		Major contributions from:
# C			P. C. W. Leung
# C			A. J. Schultz
# C			R. G. Teller
# C			L. R. Falvello
# C
# C     The data output by this program are corrected  for  variations  in
# C  spectral distribution, variations in detector efficiency  across  the
# C  face  of  the  detector,  and the pertinent geometric factors such as
# C  (SIN(THETA))**2 and LAMBDA**4.
# C
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	Linux version:	A. Schultz   January, 2003                    !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# !	Version 4:		A. Schultz		June, 2003
# !		1. Process data from two detectors
# !		2. Does not use an x-file.
# !		3. Gets MONCNT from integrate file. If zero, sets CMONX = 1.0.
# !		4. Corrected ALPHAP for calculation of SPECT1.
#
# !	Version 5:		A. Schultz		July, 2003
# !		This version outputs a expnam.hkl file which can be input
# !		into SHELX with HKL 2.
# !	Version 5a:
# !		Cleaned-up and removed a lot of unused code.
# !		Added a test for dmin.
# !
# !	Version 6:		L. Falvello		January, 2004
# !		Polyhedral absorption correction with two detectors.
# !
# !	Version 7:		A. Schultz		2007
# !		Use spectrum obtained from each SCD detector
# !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	ANVRED_SNS:		A. Schultz		2008                                     !
# !		Process SNS data.                                                    !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	ANVRED_SNS_v2:		A. Schultz		June, 2008
# !		New spherical absorption correction. Removed all
# !		of the old correction code.
# !	ANVRED_SNS-v2.1: read detector parameters from integrate file.
# !	ANVRED_SNS-v2.2: get filename for spectrum file.
# !       ANVRED_SNS-v2.3: everything included in one file.  8/17/2008
# !       anvredSNS_2.4: made compatible with gfortran including removal
# !                      of FREIN3 and READ133.         10/8/2008
# !	anvredSNS_2.5: the datacom_SNS.inc file is no longer used. Cleaned
# !			up the code using ftnchek.    10/13/08
# !	anvredSNS_2.6: assign a common scale factor for each
# !			crystal setting, or for each detector. 1/29/09
# !
# !	4/13/09	Number of possible spectra increased from 2 to 100.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from __future__ import print_function

import os
import sys
import numpy
import math

if os.path.exists('/SNS/TOPAZ/shared/PythonPrograms/PythonLibrary'):
    sys.path.append('/SNS/TOPAZ/shared/PythonPrograms/PythonLibrary')
    sys.path.append('/SNS/software/ISAW/PythonSources/Lib')
else:
    sys.path.append('C:\ISAW_repo\PythonPrograms\PythonLibrary')
from readrefl_header import *
from readrefl_SNS import *
from readSpecCoef import *
from spectrumCalc import *
from absor_sphere import *

if sys.version_info > (3,):
    from tkinter import *
    import tkinter.simpledialog as tkSimpleDialog
else:
    from Tkinter import *
    import tkSimpleDialog
 
class MyDialog(tkSimpleDialog.Dialog):

    def body(self, master):
 
        parameters = []
        try:
            # Read user input parameters
            user_input = open('anvred3.inp', 'r')
            # while True:
            for i in range(17):
                lineString = user_input.readline()
                lineList = lineString.split()
                if len(lineList) == 0:  # check for blank record
                    parameters.append( '0' )
                    continue
                parameters.append(lineList[0])
            user_input.flush()
            user_input.close()
        except:
            for i in range(17):
                parameters.append('')
        
        Message(master, 
            text = (
                "If 'spherical' absorption correction is selected, output files are *.hkl and *.hkl_dn, "
                + "\nwhere the scale factors are based on crystal settings or detector number, respectively.  "
                + "\nIf 'polyhedral' absorption correction is selected, output files are *.hkl_no_abs and *.hkl_no_abs_dn, "
                + "\nwhere _nob_abs indicates the absorption correction will be performed in a later procedure."
                + "\n\nDefault user input values are obtained from anvred3.inp if it exists. A new anvred3.inp file will be written."), 
            anchor=E, aspect=1000, bg='yellow').grid(row=0, columnspan=2)
        
        self.label_text = []
        self.label_text.append( "Working directory" )
        self.label_text.append( "Experiment name" )
        self.label_text.append( "Total scattering linear absorption coefficent in units of cm^-1 (ignored if 'polyhedral')" )
        self.label_text.append( "True absorption linear absorption coefficent at 1.8 A in units of cm^-1 (ignored if 'polyhedral')" )
        self.label_text.append( "Radius of spherical crystal in units of cm (ignored if 'polyhedral')" )
        self.label_text.append( "Name of spectrum file" )
        self.label_text.append( "Normalize spectra to this wavelength" )
        self.label_text.append( "Minimum I/sig(I) ratio" )
        self.label_text.append( "Width of border in which peaks are rejected" )
        self.label_text.append( "Minimum integrated intensity" )
        self.label_text.append( "Minimum d-spacing" )
        self.label_text.append( "Multiply FSQ and sig(FSQ) by this factor" )
        self.label_text.append( "Minimum wavelength" )
        self.label_text.append( "Maximum wavelength" )
        self.label_text.append( "'spherical' or 'polyhedral' absorption correction" )
        self.label_text.append( "UB matrix file" )
        self.label_text.append( "'TOPAZ' or 'MaNDi' detector corrections" )
        
        for i in range( len( self.label_text) ):
            j = i+1
            Label(master, text = self.label_text[i]).grid(row = j, column = 1, sticky = W)
        
        # Populate entry boxes with default values
        self.entry = []
        number_of_params = len( parameters )
        for i in range(number_of_params):
            self.entry.append( Entry( master, width = 30) )
            self.entry[i].insert( 0, parameters[i] )
            j = i+1
            self.entry[i].grid( row = j, column = 0, sticky = W )
             
    def apply(self):
    
        self.result = []
        for i in range ( len(self.entry) ):
            if len( self.entry[i].get() ) == 0:  # check for blank input
                self.result.append( '0' )
                continue
            self.result.append( self.entry[i].get() )
                
def huq(h, k, l, UB):
    "Multiplies hkl times UB matrix to return q-vector"
    hh = [h, k, l]
    q = numpy.zeros((3))
    q = numpy.dot(hh, UB)
    return q

def rotate_matrix(UB, omega, chi, phi, SNS_or_IPNS):
    "Rotates UB matrix by setting angles"
    fmat = rotation_matrix(omega, chi, phi, SNS_or_IPNS)
    newmat = numpy.zeros((3,3))
    newmat = numpy.dot(UB, fmat)
    return newmat

def rotation_matrix(omega, chi, phi, SNS_or_IPNS):
    "Returns rotation matrix from setting angles"
    rad = 180. / math.pi

    ph = phi / rad
    cp = math.cos(ph)
    sp = math.sin(ph)
    R_phi = numpy.zeros((3,3))
    R_phi[0,0] = cp
    R_phi[0,1] = sp
    R_phi[1,0] = -sp
    R_phi[1,1] = cp
    R_phi[2,2] = 1.0

    ch = chi / rad        #changed -chi to chi, 8/23/07
    cc = math.cos(ch)
    sc = math.sin(ch)
    R_chi = numpy.zeros((3,3))
    R_chi[0,0] = 1.0
    R_chi[1,1] = cc
    R_chi[1,2] = sc
    R_chi[2,1] = -sc
    R_chi[2,2] = cc

    if SNS_or_IPNS == 'IPNS':
        om = -omega / rad   # for IPNS data set omega to -omega
    if SNS_or_IPNS == 'SNS':
        om = omega / rad      # for SNS data
    co = math.cos(om)
    so = math.sin(om)
    R_om = numpy.zeros((3,3))
    R_om[0,0] = co
    R_om[0,1] = so
    R_om[1,0] = -so
    R_om[1,1] = co
    R_om[2,2] = 1.0

    fmat = numpy.zeros((3,3))
    fmat = numpy.dot(R_phi, R_chi)
    fmat = numpy.dot(fmat, R_om)

    return fmat

def Rvec(twoth, az, L2):
    "Return R vectors for a peak in peaks file."
    
    R_IPNS = numpy.zeros((3))
    R_SNS = numpy.zeros((3))
    
    # IPNS axes with x as the beam direction and z is vertically upward
    R_IPNS[0] = math.cos(twoth) * L2
    R_IPNS[1] = math.cos(az) * math.sin(twoth) * L2
    R_IPNS[2] = math.sin(az) * math.sin(twoth) * L2
    
    # SNS axes with z as the beam direction and y is vertically upward
    R_SNS[0] = math.cos(az) * math.sin(twoth) * L2
    R_SNS[1] = math.sin(az) * math.sin(twoth) * L2
    R_SNS[2] = math.cos(twoth) * L2
    
    return R_IPNS, R_SNS
    
#--------------------------------------------------------
#               function spectrum2
#--------------------------------------------------------

#!  Obtain spectral correction from counts vs. time data
#!  in a Bankxx_spectrum.asc file.
#!  Fortran version: A. J. Schultz, July, 2009
#!  Jython version: A. J. Schultz, March, 2010

#  Also returns the relative sigma of the spectral correction.
#  A. J. Schultz, April, 2011

#  spectrum2 does not average over a +/- averageRange.
#  This is because TOPAZ_spectrum now includes
#  a Savitzky-Golay smoothing Filter.
#  A. J. Schultz, September, 2010

#  Parameters:
#  wavelength = wavelength in Angstroms
#  xtof = (L1 + detD)/hom; TOF = wl * xtof
#  spect1 = spectrum at normalization wavlength, usually 1 Angstrom
#  xtime = spectrum TOF array
#  xcounts = spectrum counts array

def spectrum2( wavelength, xtof, spect1, xtime, xcounts ):
    "Returns the relative spectrum and detector efficiency correction."
	
	
# TOF = WL * XTOF in units of microseconds
    TOF = wavelength * xtof    
	
    numTimeChannels = len( xtime )
    
    spect = 0.0
    spectx = 0.0
         
# begin determining the spectrum correction
    for j in range(numTimeChannels):
        
        if xtime[j] > TOF:
            deltaCounts = xcounts[j] - xcounts[j-1]
            deltaTime = xtime[j] - xtime[j-1]
            fraction = (TOF - xtime[j-1]) / deltaTime
            spectx = xcounts[j-1] + deltaCounts*fraction # interpolate
            break
    if spect1 == 0.0:
        spect = 0.0
        relSigSpect = 0.0
        return spect, relSigSpect
    elif spectx <= 0.0:
        spect = 0.0
        relSigSpect = 0.0
        return spect, relSigSpect    
    else:
        spect = spectx / spect1
    
    # relative sigma for spect
    # relSigSpect**2 = (sqrt(spectx)/spectx)**2 + (sqrt(spect1)/spect1)**2
    relSigSpect = sqrt((1.0/spectx) + (1.0/spect1))
    
    
    return spect, relSigSpect
    
    
#------------------ Begin ------------------
root = Tk()
root.withdraw()
root.title("ANVRED3 input")
d = MyDialog(root) 

directory_path = d.result[0]
expName = d.result[1]
smu = float( d.result[2] )
amu = float( d.result[3] )
radius = float( d.result[4] )
iSpec = 0
specCoeffFile = ' '
spectraFile = d.result[5]
normToWavelength = float( d.result[6] )
minIsigI = float( d.result[7] )
numBorderCh = int( d.result[8] )
intiMin = float( d.result[9] )
dMin = float( d.result[10] )
iIQ = 1
scaleFactor = float( d.result[11] )
wlMin = float( d.result[12] )
wlMax = float( d.result[13] )
abs_correc_type = d.result[14]
ub_matrix_file = d.result[15]
TOPAZ_or_MaNDi = d.result[16]

# Check that one of the two absorption correction types has been selected
if abs_correc_type != 'spherical' and abs_correc_type != 'polyhedral':
    print('')
    print('**********************************************************')
    print("Absorption correction type is not 'spherical' or 'polyhedral'.")
    print('Check your spelling.')
    print('**********************************************************')
    print('')
    exit()

# Write or over-write anvred3.inp file
user_input = open( 'anvred3.inp', 'w' )
for i in range( len( d.result ) ):
    user_input.write( d.result[i] + '     # ' + d.label_text[i] + '\n' )

# Read UB matrix if polyhedral absorption
# if abs_correc_type == 'polyhedral':
if True:
    # Open matrix file
    UB_input = open(ub_matrix_file,'r')

    # Initialize UB_IPNS matrix
    UB_IPNS = numpy.zeros((3,3))
    print('\n Input from matrix file ' + ub_matrix_file + ':\n')

    # Read matrix file into UB_IPNS matrix
    for i in range(3):
        linestring = UB_input.readline()
        print(linestring.strip('\n'))
        linelist = linestring.split()
        for j in range(3):
            UB_IPNS[i,j] = float(linelist[j])           
    # Read next 2 lines containing lattice constants and esd's
    for i in range(2):
        linestring = UB_input.readline()
        print(linestring.strip('\n'))
    print('\n')
    # End of reading and printing matrix file
                
# TOPAZ detector scale factors for vanadium/niobium spectrum
#Scolecite 2013B, XP Wang Sept 23, 2013 
detScale = {17:1.115862021, 18:0.87451341,\
      22:1.079102931, 26:1.087379072, 27:1.064563992, 28:0.878683269, \
      36:1.15493377, 37:1.010047685, 38:1.046416037, 39:0.83264528, \
      47:1.06806776, 48:0.872542083,\
      58:0.915242691}
#Scolecite 2014A, XP Wang March 11, 2014       
# detScale = {16:0.8857, 17:1.0649, 18:0.8836, 19:1.1039,\
    # 23:0.8976, 26:0.9996, 27:1.0727 ,28:0.9305, 29:1.0920,\
    # 33:1.0673, 36:1.0532, 37:0.9621, 38:1.0754, 39:0.9144,\
    # 46:0.9421, 47:1.0703, 48:0.9387, 49:1.1068,\
    # 58:0.9390}

# open the anvred.log file in the working directory
fileName = directory_path + 'anvred3.log'
logFile = open( fileName, 'w' )

# open the hkl file in the working directory
hklFileName = directory_path + expName + '.hkl'
if abs_correc_type == 'polyhedral':
    hklFileName = directory_path + expName + '.hkl_no_abs'
hklFile = open( hklFileName, 'w' )

# echo the input in the log file
logFile.write('\n********** anvred **********\n')
logFile.write('\nWorking directory: ' + directory_path)
logFile.write('\nExperiment name: ' + expName + '\n')

logFile.write('\nTotal scattering linear absorption coefficient: %6.3f cm^-1' % smu )
logFile.write('\nTrue absorption linear absorption coefficient: %6.3f cm^-1' % amu )
logFile.write('\nRadius of spherical crystal: %6.3f cm\n' % radius )

# logFile.write('\nIncident spectrum and detector efficiency correction.')
# logFile.write('\n    iSpec = 1. Spectrum fitted to 11 coefficient GSAS Type 2 function')
# logFile.write('\n    iSpec = 0. Spectrum data read from a spectrum file.')
# logFile.write('\niSpec: %i\n' % iSpec)

if iSpec == 1:   # spectrum is fitted to equation with 12 coefficients
    logFile.write('\nFile with spectrum coefficients: ' + specCoeffFile + '\n' )
    
if iSpec == 0:   # spectrum is read as TOF vs. counts
    logFile.write('\nFile with spectra: ' + spectraFile + '\n' )

logFile.write('\nNormalize spectra to a wavelength of %4.2f' % normToWavelength)
logFile.write('\nThe minimum I/sig(I) ratio: %i' % minIsigI )
logFile.write('\nWidth of border: %i channels' % numBorderCh )
logFile.write('\nMinimum integrated intensity: %i' % intiMin )
logFile.write('\nMinimum d-spacing : %4.2f Angstroms\n' % dMin )

# logFile.write('\nScale factor identifier:' )
# logFile.write('\n     IQ = 1. Scale factor per crystal setting.' )
# logFile.write('\n     IQ = 2. Scale factor for each detector in each setting.')
# logFile.write('\n     IQ = 3. Scale factor for each detector for all settings.')
# logFile.write('\nIQ: %i\n' % iIQ )

logFile.write('\nMultiply FSQ and sig(FSQ) by: %f\n' % scaleFactor )

logFile.write('\nMinimum wavelength: %f\n' % wlMin )
logFile.write('Maximum wavelength: %f\n' % wlMax )

logFile.write( '\n***** absorption correction type: ' + abs_correc_type + '\n' )
if abs_correc_type == 'polyhedral':
    logFile.write( '\nUB matrix file: ' + ub_matrix_file + '\n' )

# C
# C  CHECK ON THE EXISTANCE OF THE integrate FILE
# C
fileName = directory_path + expName + '.integrate'
integFile = open(fileName, 'r')

# !  Initial read of integrate file to get instrument and detectors calibration.
calibParam = readrefl_header( integFile )
L1 = float(calibParam[0])       # initial flight path length in cm
t0_shift = float(calibParam[1]) # t-zero offest in microseconds
nod = int(calibParam[2])    # number of detectors
print('********** nod = ', nod)

logFile.write('\nInitial flight path length: %10.4f cm' % L1 )
logFile.write('\nT-zero offset: %8.3f microseconds' % t0_shift )
logFile.write('\nNumber of detectors: %i' % nod )

# Initial values.
transmin = 1.0
transmax = 0.0
hom = 0.39559974    # Planck's constant divided by neutron mass

# Read spectrum coefficients if iSpec = 1
if iSpec == 1:
    specInput = open( specCoeffFile, 'r')
    # pj is a list of lists with dimensions (nod, 11)
    pj = readSpecCoef(specInput, logFile, nod)
    
# Read spectrum for each detector bank if iSpec = 0
if iSpec == 0:
    # spectra is an array of arrays containing the spectra in the
    # Spectrum_run1_run2.dat file.
    specInput = open( spectraFile, 'r' )
    
    for i in range(8):   # skip the first 8 lines
        lineString = specInput.readline()
        # print lineString
    
    # "spectra" is an array spectra[i][j] where i is the number
    # of the detector bank starting at zero, and j = 0 for
    # the array of times and j = 1 for the array of counts
    spectra = []
    
    lineString = specInput.readline()   # read "Bank 1" line
    
    for i in range( nod ):
        # set arrays to zero
        time = []
        counts = []
        
        print('Reading spectrum for ' + lineString[0:-1])
        while True:
            lineString = specInput.readline()
            lineList = lineString.split()
            if len(lineList) == 0: break
            if lineList[0] == 'Bank': break
            time.append( float( lineList[0] ) )
            counts.append( float( lineList[1] ) )
            
        spectra.append( [time, counts] )
    
    specInput.close()
    

# C-----------------------------------------------------------------------
# C  Calculate spectral correction at normToWavelength to normalize
# C  spectral correction factors later on.
spect1 = []     # spectrum value at normToWavelength for each detector
dist = []       # sample-to-detector distance
xtof = []       # = (L1+dist)/hom; TOF = wl * xtof

wavelength = normToWavelength
one = 1.0       # denominator in spectrum to calculate spect1

for id in range(nod):

    if iSpec == 1:  # The spectrum is calculated from coefficients
        
        spect = spectrumCalc(wavelength, calibParam, pj, id)
        spect1.append(spect)
        
    else:           # iSpec = 2           
        
        dist.append(calibParam[9][id])
        xtof.append((L1 + dist[id]) / hom)
        
        # spectra[id][0] are the times-of-flight
        # spectra[id][1] are the counts
        spectx = spectrum2( wavelength, xtof[id], \
            one, spectra[id][0], spectra[id][1] )
        spect = spectx[0]            # the spectral normalization parameter
        relSigSpect = spectx[1]      # the relative sigma of spect
        if spect == 0.0:
            print('*** Wavelength for normalizing to spectrum is out of range.')
        spect1.append(spect)
                          
# C-----------------------------------------------------------------------

# C
# C  SET THE CURRENT HISTOGRAM NUMBER TO 0 AND INITIALIZE THE MONITOR COUN
# C
curhst = 0
idet = 0
hstnum = 0
cmon = 100e+6
ncntr = 0      #!Number of processed reflections

nrun = 0
dn = 0
chi = 0.0
phi = 0.0
omega = 0.0
moncnt = 1000000.
eof = 999
hkllists =[]   # List of reflections, XP Wang, May 3013
fsqmax = 0.0
# C
# C   SET UP LOOP TO PROCESS THE REFLECTION DATA
# C
while True:

    peak = readrefl_SNS( integFile, eof, nrun, dn, chi, phi, omega,\
        moncnt)
    eof = peak[22]
    if eof == 0: break
    
    nrun = peak[0]
    dn = peak[1]
    chi = float( peak[2] )
    phi = float( peak[3] )
    omega = float( peak[4] )
    moncnt = peak[5]
    seqnum = peak[6]
    h = peak[7]
    k = peak[8]
    l = peak[9]
    col = peak[10]
    row = peak[11]
    chan = peak[12]
    L2 = peak[13]
    twoth = peak[14]  # radians
    az = peak[15]  # azimuthal angle in radians
    wl = peak[16]
    dsp = peak[17]
    ipkobs = peak[18]
    inti = peak[19]
    sigi = abs(peak[20])
    reflag = peak[21]
    
    if seqnum % 1000 == 0: print('seqnum =', seqnum)
    

    # set-up for new run or detector
    if nrun != curhst or dn != idet:
        if nrun != curhst:
            curhst = nrun
            if iIQ != 2: hstnum = hstnum + 1
            
            # Rotate UB matrix if polyhedral absorption correction
            # if abs_correc_type == 'polyhedral':
            if True:
                SNS_or_IPNS = 'SNS'
                # Using UB_IPNS with SNS rotation angles
                newmat = rotate_matrix(UB_IPNS, omega, chi, phi, SNS_or_IPNS)
                
                # Calculate direction cosines for reversed incident beam vector.
                # IPNS coordinates.
                R_reverse_incident = [ -L2, 0.0, 0.0 ]
                dir_cos_1 = [ 0, 0, 0 ]
                # Begin loop through a-star, b-star and c-star
                for i in range(3):
                    hkl = [ 0, 0, 0 ]
                    hkl[i] = 1
                    q_abc_star = huq( hkl[0], hkl[1], hkl[2], newmat )
                    length_q_abc_star = math.sqrt( numpy.dot( q_abc_star, q_abc_star) )                    
                    dir_cos_1[i] = ( numpy.dot( R_reverse_incident, q_abc_star ) 
                        / ( L2 * length_q_abc_star ) )
                
        idet = dn  #IDET and DN is the arbitrary detector number.
                   #ID is a sequential number in the order they are listed.
     
        for id in range(nod):
            detNum = calibParam[3][id]
            if detNum == dn: break
            
        if iIQ == 2: hstnum = hstnum + 1
        
        mnsum = moncnt
        
        if mnsum == 0:
            cmonx = 1.0
        else:
            cmonx = cmon / mnsum
            if cmonx == 0: cmonx = 1.0
        
        logFile.write('\n\nHISTOGRAM NUMBER %5d' % nrun)      
        logFile.write('\nDETECTOR BANK NUMBER %2d     DETECTOR SEQUENTIAL NUMBER %2d'\
            % (dn, id))
        logFile.write('\nANGLES ARE CHI =%7.2f   PHI =%7.2f   OMEGA=%7.2f\n'\
            % ( chi, phi, omega ))                        
        logFile.write('TOTAL MONITOR COUNTS ELAPSED%10d   CMONX =%10.4f\n'\
            % ( mnsum, cmonx ))
        logFile.write('* DATA SCALED TO 100 MILLION MONITOR COUNTS *\n')
        logFile.write('CORREC = SCALEFACTOR * CMONX * SINSQT /' + \
            '( SPECT * (DET EFF) * WL4 * ABTRANS )\n')
        if abs_correc_type == 'spherical':
            logFile.write('\n    H   K   L       FSQ     SIG     WL      INTI' + \
                '     SIG   SPECT  SINSQT   TRANS    TBAR\n')
        if abs_correc_type == 'polyhedral':
            logFile.write('\n    H   K   L       FSQ     SIG     WL      INTI' + \
                '     SIG   SPECT  SINSQT       l1       l2       m1       m2       n1       n2\n')
    # end of set-up for new run or detector
   
    # Omit zero intensity peaks from integrate file XP Wang 03/21/2011
    # Changed to >=0 and absolute value  XP Wang 02/24/2011

    # Omit peaks not indexed XP Wang, March, 2013
    if (h==0 and k==0 and l==0):
        logFile.write(' %4d *** Peak not indexed for run %4d det %4d   \n' \
            % (seqnum,nrun,dn))        
        continue  
    if inti == 0.0 :
        logFile.write(' %4d%4d%4d *** intI = 0.0 \n' \
            % (h, k, l))
        continue  

    if minIsigI >= 0 and inti < abs(minIsigI * sigi):
        logFile.write(' %4d%4d%4d *** inti < (minIsigI * sigi) \n' \
            % (h, k, l))
        continue
        
    if inti < intiMin:
        logFile.write(' %4d%4d%4d *** inti < intiMin \n' \
            % (h, k, l))
        continue

    # Set-up limits for neutron wavelentgh XP Wang 02/24/2011
    if wl < wlMin:
        logFile.write(' %4d%4d%4d *** wl < wlMin \n' \
            % (h, k, l))
        continue

    if wl > wlMax:
        logFile.write(' %4d%4d%4d *** wl > wlMax \n' \
            % (h, k, l))
        continue

    nRows = calibParam[4][id]
    nCols = calibParam[5][id]
    
    if col < numBorderCh:
        logFile.write(' %4d%4d%4d *** col < numBorderCh \n' \
            % (h, k, l))
        continue
        
    if col > (nCols - numBorderCh):
        logFile.write(' %4d%4d%4d *** col > (nCols - numBorderCh)\n' \
            % (h, k, l))
        continue
        
    if row < numBorderCh:
        logFile.write(' %4d%4d%4d *** row < numBorderCh \n' \
            % (h, k, l))
        continue
        
    if row > (nRows - numBorderCh):
        logFile.write(' %4d%4d%4d *** row > (nRows - numBorderCh)\n' \
            % (h, k, l))
        continue
                        
    if dsp < dMin:
        logFile.write(' %4d%4d%4d *** dsp < dMin \n' \
            % (h, k, l))
        continue
    
    ncntr = ncntr + 1
    
    if iSpec == 1:
        spect = spectrumCalc(wl, calibParam, pj, id)
        spect = spect / spect1[id]
    
    if iSpec == 0:
        spectx = spectrum2( wl, xtof[id], \
          spect1[id], spectra[id][0], spectra[id][1] )
        spect = spectx[0]
        relSigSpect = spectx[1]
    if spect == 0.0:
        logFile.write(' %4d%4d%4d *** spect == 0.0 \n' \
            % (h, k, l))
        continue
    
    # correct for the slant path throught the scintillator glass
    mu = (9.614 * wl) + 0.266    # mu for GS20 glass
    depth = calibParam[8][id]
    eff_center = 1.0 - exp(-mu * depth)  # efficiency at center of detector
    cosA = dist[id] / L2
    pathlength = depth / cosA
    eff_R = 1.0 - exp(-mu * pathlength)   # efficiency at point R
    sp_ratio = eff_center / eff_R  # slant path efficiency ratio
    
    sinsqt = ( wl / (2.0*dsp) )**2
    wl4 = wl**4
    
    if TOPAZ_or_MaNDi == 'TOPAZ':
        correc = scaleFactor * sinsqt * cmonx * sp_ratio / (wl4 * spect ) * detScale[detNum]
    if TOPAZ_or_MaNDi == 'MaNDi':
        correc = scaleFactor * sinsqt * cmonx * sp_ratio / (wl4 * spect )
        
    # absorption correction
    # trans[0] is the transmission
    # trans[1] is tbar
    if abs_correc_type == 'spherical':
        trans = absor_sphere(smu, amu, radius, twoth, wl)
        transmission = trans[0]
        if trans[0] < transmin: transmin = trans[0]
        if trans[0] > transmax: transmax = trans[0]
        
        correc = correc / trans[0]
    
    fsq = inti * correc

    sigfsq = sigi * correc
    
    #sigfsq = sqrt( sigfsq**2 + (relSigSpect*fsq)**2)  # not sure if last term is squared
    # Add instrument background constant to sigma, XP WAng June 2013
    sigfsq = sqrt( sigfsq**2 + (relSigSpect*fsq)**2 + 12.28/cmonx*scaleFactor)

    # Calculates direction cosines for scattered beam vector 
    R_IPNS, R_SNS = Rvec(twoth, az, L2)
    # Begin loop through a-star, b-star and c-star
    dir_cos_2 = [ 0, 0, 0 ]   # dir_cos_2 is the scattered beam with a*, b*, c*
    for i in range(3):
        abc = [ 0, 0, 0 ]
        abc[i] = 1
        q_abc_star = huq( abc[0], abc[1], abc[2], newmat )
        len_q_abc_star = math.sqrt( numpy.dot( q_abc_star, q_abc_star) )
        dir_cos_2[i] = numpy.dot( R_IPNS, q_abc_star ) / ( L2 * len_q_abc_star )
    
    # Write out to hkl file for spherical absorption correction.    
    if abs_correc_type == 'spherical':
        # tbar is the Coppen's tbar
        tbar = trans[1]
        
        # output reflection to log file and to hkl file
        logFile.write( ' '
            + 3*'%4d' % (h,k,l)
            + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
            + 4*'%8.4f' % (spect, sinsqt, trans[0], tbar)
            + '\n' )
        
        hklFile.write( 3*'%4d' % (h,k,l)
            + 2*'%8.2f' % (fsq, sigfsq) 
            + '%4d' % hstnum 
            + 2*'%8.5f' % (wl, tbar)
            + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
            + 2*'%6d' % (curhst, seqnum)
            + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )            

    # Write out to hkl file for polyhedral absorption correction.    
    if abs_correc_type == 'polyhedral':
    
        # output reflection to log file and to hkl file
        logFile.write(' '
            + 3*'%4d' % (h, k, l)
            + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
            + 2*'%8.4f' % (spect, sinsqt)
            + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
            + '\n' )
        
        tbar = 0.0
        transmission = 1.0
        hklFile.write( 3*'%4d' % (h, k, l)
            + 2*'%8.2f' % (fsq, sigfsq)
            + '%4d' % hstnum 
            + 2*'%8.5f' % (wl, tbar)
            + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
            + 2*'%6d' % (curhst, seqnum)
            + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )
            
    #Add to hkl list XP Wang May 2013
    hkllists.append([ h, k, l, fsq, sigfsq, hstnum, wl, tbar,
        dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], 
        dir_cos_1[2], dir_cos_2[2], curhst, seqnum, transmission,
        dn, twoth, dsp, col, row])
    if fsq > fsqmax: fsqmax = fsq        
            
if abs_correc_type == 'spherical':
    print('\nMinimum and maximum transmission = %6.4f, %6.4f' % (transmin, transmax))

logFile.write('\n\n***** Minimum and maximum transmission = %6.4f, %6.4f' \
    % (transmin, transmax))
    
if fsqmax > (10000.00 -1):
    print()
    print('################################################################')
    print('     Maximum FOSQ is', fsqmax)
    print('     Re-run anvred with a scale factor of', 0.1*scaleFactor)
    print('################################################################')

# last record all zeros for shelx
iz = 0
rz = 0.0
hklFile.write( 3*'%4d' % (iz, iz, iz)
    + 2*'%8.2f' % (rz, rz)
    + '%4d' % iz 
    + 2*'%8.5f' % (rz, rz)
    + 6*'%9.5f' % (rz, rz, rz, rz, rz, rz)
    + 2*'%6d' % (iz, iz)
    + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (rz, iz, rz, rz, rz, rz) )
    
# C-----------------------------------------------------------------------
logFile.close()
hklFile.close()


# Set scale ID equal to detector number.
# This code is from scale_by_detnum.py.
# Save reflections by DetNum  XP Wang May, 2013
from operator import itemgetter
from itertools import groupby

if iIQ == 1 or iIQ == 3:
    hklFileName1 = hklFileName + '_dn'
    hkl_output = open(hklFileName1, 'w')
    
    #Sort and save the result per module, XP Wang, March 2013
    hkllists.sort(key=itemgetter(17,0,1,2))  # sort by detector number
    nDet = 0
    for iDet, iGroup in groupby(hkllists, itemgetter(17)):
        nDet = nDet + 1
        for iHKL in iGroup:
            iHKL[5] = nDet                
            # output reflection sorted by detector number to hkl file                
            hkl_output.write( 3*'%4d' % (iHKL[0],iHKL[1],iHKL[2])
                + 2*'%8.2f' % (iHKL[3],iHKL[4])
                + '%4d' % iHKL[5] 
                + 2*'%8.5f' % (iHKL[6], iHKL[7])
                + 6*'%9.5f' % (iHKL[8],iHKL[9],iHKL[10],iHKL[11],iHKL[12],iHKL[13])
                + 2*'%6d' % (iHKL[14], iHKL[15])
                + '%7.4f%4d%9.5f%8.4f' % (iHKL[16],iHKL[17],iHKL[18],iHKL[19])
                + 2*'%7.2f' % (iHKL[20], iHKL[21])
                + '\n' )
                                
    # last record all zeros for shelx
    hkl_output.write( 3*'%4d' % (iz, iz, iz)
        + 2*'%8.2f' % (rz, rz)
        + '%4d' % iz 
        + 2*'%8.5f' % (rz, rz)
        + 6*'%9.5f' % (rz, rz, rz, rz, rz, rz)
        + 2*'%6d' % (iz, iz)
        + '%7.4f%4d%9.5f%8.4f' % ( rz, iz, rz, rz )
        + '%7.2f%7.2f' % ( rz, rz )
        + '\n' )

hkl_output.close()

print('All done!')

        


