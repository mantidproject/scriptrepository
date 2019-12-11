# File: stats_vs_time.py
#
# This script will vary the elapsed time and output statistics.
#
# NOTE: All of the parameters that the user must specify are listed with 
# instructive comments in the sample configuration file: ReduceSCD.config.
#
from __future__ import print_function
import os
import sys
import pylab
import time
import math

if os.path.exists('/SNS/TOPAZ/shared/PythonPrograms/PythonLibrary'):
    sys.path.append('/SNS/TOPAZ/shared/PythonPrograms/PythonLibrary')
elif os.path.exists('/SNS/MANDI/shared/PythonPrograms/PythonLibrary'):
    sys.path.append('/SNS/MANDI/shared/PythonPrograms/PythonLibrary')
else:
    sys.path.append('C:\ISAW_repo\PythonPrograms\PythonLibrary')
    
import ReduceDictionary

if os.path.exists("/opt/Mantid/bin"):
    sys.path.append("/opt/Mantid/bin")
    # sys.path.append("/opt/mantidnightly/bin")
else:
    sys.path.append(os.environ['MANTIDPATH'])

from mantid.simpleapi import *

start_time = time.time()

#
# Load the parameter names and values from the specified configuration file 
# into a dictionary and set all the required parameters from the dictionary.
#
config_file_name = 'stats_vs_time.config'
params_dictionary = ReduceDictionary.LoadDictionary( config_file_name )

run                       = params_dictionary[ "run" ]
instrument_name           = params_dictionary[ "instrument_name" ]
calibration_file_1        = params_dictionary[ "calibration_file_1" ]
calibration_file_2        = params_dictionary[ "calibration_file_2" ]
data_directory            = params_dictionary[ "data_directory" ]
output_directory          = params_dictionary[ "output_directory" ]
min_tof                   = params_dictionary[ "min_tof" ] 
max_tof                   = params_dictionary[ "max_tof" ] 
min_pred_wl               = params_dictionary[ "min_pred_wl" ]
max_pred_wl               = params_dictionary[ "max_pred_wl" ]
min_pred_dspacing         = params_dictionary[ "min_pred_dspacing" ]
max_pred_dspacing         = params_dictionary[ "max_pred_dspacing" ]
peak_radius               = float( params_dictionary[ "peak_radius" ] )
bkg_inner_radius          = float( params_dictionary[ "bkg_inner_radius" ] )
bkg_outer_radius          = float( params_dictionary[ "bkg_outer_radius" ] )
integrate_if_edge_peak    = params_dictionary[ "integrate_if_edge_peak" ]
MultiplierBase            = float( params_dictionary[ "MultiplierBase" ] )

#
# Get the fully qualified input run file name, either from a specified data 
# directory or from findnexus
#
if data_directory is not None:
  full_name = data_directory + "/" + instrument_name + "_" + run + "_event.nxs"
else:
  temp_buffer = os.popen("findnexus --event -i "+instrument_name+" "+str(run) )
  full_name = temp_buffer.readline()
  full_name=full_name.strip()
  if not full_name.endswith('nxs'):
    print("Exiting since the data_directory was not specified and")
    print("findnexus failed for event NeXus file: " + instrument_name + " " + str(run))
    exit(0)

print("\nProcessing File: " + full_name + " ......\n")

# Name the files for this run
run_niggli_matrix_file = output_directory + "/" + run + "_Niggli.mat"
run_niggli_integrate_file = output_directory + run + "_Niggli.integrate"
temp_integrate_file = output_directory + "/temp.integrate"

# load integrate file and UB matrix into peaks_ws
peaks_ws = LoadIsawPeaks( Filename = run_niggli_integrate_file )
LoadIsawUB( InputWorkspace = peaks_ws, Filename = run_niggli_matrix_file )

# Load the monitor data and get the counting time in seconds
event_ws = LoadNexusMonitors( Filename=full_name )
total_time = event_ws.run()['duration'].value

number_of_steps = ( math.log(total_time) - math.log(60.0) ) / math.log(MultiplierBase)
number_of_steps = int(number_of_steps) + 1
print('\nnumber_of_steps = ', number_of_steps)
print('')
# ln_MultiplierBase = ( math.log(total_time) - math.log(60.0) ) / 10.0
# MultiplierBase = math.exp(ln_MultiplierBase)
# print 'For 10 steps, MultiplierBase =', MultiplierBase
                           
# Get complete list of peaks to be integrated and load the UB matrix into
# the predicted peaks workspace, so that information can be used by the
# PeakIntegration algorithm.
print("PREDICTING peaks to integrate....")
peaks_ws = PredictPeaks( InputWorkspace = peaks_ws,
        WavelengthMin = min_pred_wl, WavelengthMax = max_pred_wl,
        MinDSpacing = min_pred_dspacing, MaxDSpacing = max_pred_dspacing, 
        ReflectionCondition = 'Primitive' )
                
#
# Begin set-up for analysing statistics
#
statsFilename = output_directory + 'stats_vs_time.lst'
output = open(statsFilename, 'w')

output.write('\n        ***  SUMMARY OF INTENSITY STATISTICS  ***\n\n')
output.write('      Time     TOTAL      2SIG      3SIG      5SIG     10SIG\n')

# initialize arrays for plots
x = []
ytotal = []
y2 = []
y3 = []
y5 = []
y10 = []

#
# Begin loop to integrate peaks and analyze statistics
#
i = 0                
while True:
    time_stop = 60.0 * MultiplierBase**i
    if time_stop > total_time: 
        time_stop = total_time

    event_ws = LoadEventNexus( Filename = full_name, 
                    FilterByTofMin = min_tof, FilterByTofMax = max_tof,
                    FilterByTimeStop = time_stop )
                    
    #
    # Integrate predicted peaks in Q space using spheres, and save 
    # integrated intensities, with Niggli indexing.  First get an un-weighted 
    # workspace to do raw integration (we don't need high resolution or 
    # LorentzCorrection to do the raw sphere integration )
    #
    MDEW = ConvertToMD( InputWorkspace=event_ws, QDimensions="Q3D",
            dEAnalysisMode="Elastic", QConversionScales="Q in A^-1",
            LorentzCorrection='0', MinValues="-50,-50,-50", MaxValues="50,50,50",
            SplitInto='2', SplitThreshold='500',MaxRecursionDepth='10' )

    peaks_ws = IntegratePeaksMD( InputWorkspace=MDEW, PeakRadius=peak_radius,
            CoordinatesToUse="Q (sample frame)",
            BackgroundOuterRadius=bkg_outer_radius, 
            BackgroundInnerRadius=bkg_inner_radius,
            PeaksWorkspace=peaks_ws, 
            IntegrateIfOnEdge=integrate_if_edge_peak )

    # Save the final integrated peaks in a temporary file 
    SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, 
                   Filename = temp_integrate_file )

    # begin reading the integrate file
    sum = sig2 = sig3 = sig5 = sig10 = 0
    
    input = open(temp_integrate_file, 'r')
    lineString = input.readline()           # read first header line from integrate file
    while True:
        lineString = input.readline()
        lineList = lineString.split()
        if (len(lineList)) == 0: break
        
        formatFlag = int(lineList[0])       # test for line type
                        
        if formatFlag == 3:
            
            intI = float(lineList[14])  # intI is the integrated intensity
            sigI = float(lineList[15])  # sigI is the standard deviation
            sum = sum + 1
            if intI > (2.0 * sigI):
                sig2 = sig2 + 1
            if intI > (3.0 * sigI):
                sig3 = sig3 + 1
            if intI > (5.0 * sigI):
                sig5 = sig5 + 1
            if intI > (10.0 * sigI):
                sig10 = sig10 + 1
            
    # end reading the integrate file

    #
    # Print the results.
    #
    outString = '%10.1f%10d%10d%10d%10d%10d\n' % (time_stop, sum, sig2, sig3, sig5, sig10)
    print('')
    print(outString)
    output.write(outString)
    x.append( time_stop )
    ytotal.append( sum )
    y2.append( sig2 )
    y3.append( sig3 )
    y5.append( sig5 )
    y10.append( sig10 )

    input.close()
    
    if time_stop == total_time: break
    i = i + 1
        
os.remove( temp_integrate_file )                 
output.close()
end_time = time.time()
print('\nAnalyzed run ' + str(run) + ' in ' + str(end_time - start_time) + ' sec')

#
# Plot the results
#

pylab.plot( x, y2, marker = '^', linestyle = '-', label = '2sig' )
pylab.plot( x, y3, marker = '^', linestyle = '-', label = '3sig' )
pylab.plot( x, y5, marker = '^', linestyle = '-', label = '5sig' )
pylab.plot( x, y10, marker = '^', linestyle = '-', label = '10sig' )

pylab.legend( loc = 'lower right' )     # locate the legend in the lower right

pylab.xlabel('Elapsed time, seconds')
pylab.ylabel('Number of peaks')
plotTitle = 'stats_vs_epalsed_time'
pylab.title( plotTitle )
pylab.grid(True)
pylab.savefig( plotTitle )   # plot saved

# pylab.ion()                  # turn on interactive mode
pylab.show()

input('Type RETURN to continue.')
sys.exit(0)
