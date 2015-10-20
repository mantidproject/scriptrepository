#!/usr/bin/env python
"""
Plot data in SCDcalib.log file from ISAW SCD Calibration.
A. J. Schultz, September 2015
"""

import pylab
import os
import math

# Make a ./plots subdirectory for the plot files.
if not os.path.exists('./plots'):
    os.mkdir('./plots')

input_fname = 'SCDcalib.log'
input = open(input_fname, 'r')

output_fname = 'SCDcalib_plot.log'
output = open(output_fname, 'w')
output.write('RMSD in mm units\n')
output.write(' ID  NumPeaks       Row    Column  Combined\n')

xcalc = [0,255]  
ycalc = [0,255]

# Begin reading and plotting.
# for line in input:
while True:

    line = input.readline()
    lineList = line.split()
    if len(lineList) == 0: continue
    if lineList[0] == '#': break
    
    # check for beginning of data
    if line[0:30] != 'Detector Row Number Comparison': continue
    
    title = line
    ID = lineList[-1]
    print 'ID ', ID

    x = []   # Theoretical
    y = []   # Measured
    chisq_row = 0.0
    
    # skip next 2 lines
    lineString = input.readline()
    lineString = input.readline()
    
    # plot row number comparison
    while True:
        lineString = input.readline()  
        lineList = lineString.split()
        if lineList[0] == 'Detector': break
        
        x.append( float( lineList[0] ) )
        y.append( float( lineList[1] ) )
        chisq_row = chisq_row + ( x[-1] - y[-1] )**2
    
    pylab.plot( x, y, 'r+' )
    pylab.plot( xcalc, ycalc )
    
    pylab.xlabel('Calculated Row Number')
    pylab.ylabel('Observed Row Number')
    pylab.grid(True)

    pylab.title(title)
    
    numPeaks = len(x)
    rmsd_row = math.sqrt( (1.0/numPeaks) * chisq_row )
    rmsd_row_mm = rmsd_row * 150 / 256
    reduced_chisq_row = chisq_row / ( numPeaks - 10 )
    textString = ('Number of peaks = %d \nreduced chisq = %.2f \nRMSD = %.2f ch (%.2f mm)' 
         % (numPeaks, reduced_chisq_row, rmsd_row, rmsd_row_mm))
    pylab.figtext( 0.5, 0.2, textString )
    
    filename = './plots/' + title[32:-1] + ' ' + title[9:13] + '.png'
    pylab.savefig(filename)
    pylab.clf()
    
    # ---------- plot column number comparison
    title = lineString
    
    # skip next 1 line
    lineString = input.readline()
    x = []   # Theoretical
    y = []   # Measured
    chisq_col = 0.0
    
    while True:
        lineString = input.readline()  
        lineList = lineString.split()
        if lineList[0] == 'Detector': break
        
        x.append( float( lineList[0] ) )
        y.append( float( lineList[1] ) )
        chisq_col = chisq_col + ( x[-1] - y[-1] )**2
            
    pylab.plot( x, y, 'r+' )
    pylab.plot( xcalc, ycalc )
    
    pylab.xlabel('Calculated Column Number')
    pylab.ylabel('Observed Column Number')
    pylab.grid(True)

    pylab.title(title)
    
    rmsd_col = math.sqrt( (1.0/numPeaks) * chisq_col )
    rmsd_col_mm = rmsd_col * 150 / 256
    reduced_chisq_col = chisq_col / ( numPeaks - 10 )
    textString = ('Number of peaks = %d \nreduced chisq = %.2f ch \nRMSD = %.2f ch (%.2f mm)' % 
        (numPeaks, reduced_chisq_col, rmsd_col, rmsd_col_mm))
    pylab.figtext( 0.5, 0.2, textString )
    
    filename = './plots/' + title[35:-1] + ' ' + title[9:15] + '.png'
    pylab.savefig(filename)
    pylab.clf()
    
    rmsd_combined = math.sqrt( (1.0/(2.0*numPeaks)) * ( chisq_col + chisq_row ) )
    rmsd_combined_mm = rmsd_combined * 150 / 256
    IDnum = int( ID )
    output.write(' %2d  %8d  %8.2f  %8.2f  %8.2f\n' % 
        (IDnum, numPeaks,rmsd_col_mm, rmsd_row_mm, rmsd_combined_mm) )
    
print '\nAll done!' 




    





