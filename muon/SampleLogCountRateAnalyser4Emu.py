"""*WIKI* 
Run this script on the Mantid Script window
Collates detector counts and various parameters as a function of run number
Outputs are in two tables, from which you can make graphs
Written by Koji Yokoyama on 23 Oct. 2018
Hint: just give run numbers first and try running the code
*WIKI*"""

import time
import numpy as np
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *

class GetPropertyFromRnum(object):
	def __init__(self,datapath,rnum,ring_dic):
		self.datapath = datapath
		self.rnum = rnum
		self.ring_dic = ring_dic

	def getproperty(self):
		# Get parameters log from specified run number
		thispath = self.datapath + rnum + ".nxs"
		wstuple = LoadMuonNexus(Filename=thispath)
		try:
			ws = wstuple[0][0] # first period if there are some
			nperiod = 2
		except:
			ws = wstuple[0] # for plain run
			nperiod = 1
		# Change/add this list depending on what log you want to check
		self.val1 = ws.run().getLogData("run_start").value
		self.val2 = ws.run().getLogData("Beamlog_Good_Frames_Total").times
		self.val3 = ws.run().getLogData("Beamlog_Good_Frames_Total").value
		self.val4 = ws.run().getLogData("Beamlog_Total_Counts").times
		self.val5 = ws.run().getLogData("Beamlog_Total_Counts").value
		self.val6 = ws.run().getLogData("Slits").times
		self.val7 = ws.run().getLogData("Slits").value
		self.val8 = ws.run().getLogData("Beamlog_TS1_Current").times
		self.val9 = ws.run().getLogData("Beamlog_TS1_Current").value
		self.val10 = ws.run().getLogData("field_danfysik").times
		self.val11 = ws.run().getLogData("field_danfysik").value
		self.val12 = ws.run().getLogData("field_t20").times
		self.val13 = ws.run().getLogData("field_t20").value
		# Get integrated count for each detector, and put them in a dictionary
		self.dicsum = {}
		for diclabel, data in ring_dic.items():
			integrated_counts = []
			for idx, ielem in enumerate(data):
				integrated_counts.append(sum(ws.readY(ielem-1)))
			self.dicsum[diclabel] = map(lambda x: int(x*nperiod), integrated_counts)

# About your runs
datapath = "\\\emu\data\EMU000"
rnum_start = 86090
rnum_end = 86317
# Table for detector array. This is for EMU
ring_dic = {
'ring_F1': [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
'ring_F2': [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
'ring_F3': [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48],
'ring_B1': [49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94],
'ring_B2': [50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 95],
'ring_B3': [51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96]
}

# Create output table for log params. Change/add this list depending on what log you want to check
log_table = CreateEmptyTableWorkspace()
log_table.setTitle("EventRates")
log_table.addColumn("int", "RunNumber")
log_table.addColumn("float", "MEvPerHour")
log_table.addColumn("float", "Slits")
log_table.addColumn("float", "TS1Current")
log_table.addColumn("float", "LF")
log_table.addColumn("float", "T20")

# Create output table for detectors.
det_table = CreateEmptyTableWorkspace()
det_table.setTitle("DetectorCountRates")
det_table.addColumn("int", "RunNumber")
for diclabel, data in ring_dic.items():
    for idx, ielem in enumerate(data):
		det_table.addColumn("float", "Det_"+str(ielem))
# ***Ignore from here***
# The first row is for detector numbers, starting with "0"
# DetNumbers = [0]
# for diclabel, data in ring_dic.items():
	# for idx, ielem in enumerate(data):
		# DetNumbers.append(float(ielem))
# det_table.addRow(DetNumbers)
# ***To here***

for n in range(rnum_start,rnum_end + 1):
	rnum = str(n)
	try:
		logdat = GetPropertyFromRnum(datapath,str(rnum),ring_dic)
		logdat.getproperty()
		time0 = time.mktime(time.strptime(logdat.val1,"%Y-%m-%dT%H:%M:%S")) # Calculates absolute time in s for time when run started
		
		# Number of frames
		date_time_array = map(str,np.array(logdat.val2))
		relative_times = [] # Makes a list of times (in absolute time) for the logged value.
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)	
		frames_array = map(float,np.array(logdat.val3))
		number_of_frames = frames_array[np.argmax(relative_times)] # Take the largest time for the final value
		# Total counts from Beamlog_Total_Counts
		date_time_array = map(str,np.array(logdat.val4))
		relative_times = []
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)		
		total_counts_array = map(float,np.array(logdat.val5))
		total_counts = total_counts_array[np.argmax(relative_times)]
		# Slits
		date_time_array = map(str,np.array(logdat.val6))
		relative_times = []
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)		
		slits_array = map(float,np.array(logdat.val7))
		relative_times_reduced = [] # Interested only in times (absolute time) after the start of run. Trims tables for time and log value
		slits_array_reduced = []
		for idx, ielem in enumerate(relative_times):
			if ielem >= 0:
				relative_times_reduced.append(ielem)
				slits_array_reduced.append(slits_array[idx])
		w = np.diff(relative_times_reduced, n=1)	# Parameters are recorded only when there is a change
		slits_avg = np.average(slits_array_reduced[:-1], weights = w) # Takes weighted average
		# TS1 current
		date_time_array = map(str,np.array(logdat.val8))
		relative_times = []
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)		
		ts1_current_array = map(float,np.array(logdat.val9))
		relative_times_reduced = []
		ts1_current_array_reduced = []
		for idx, ielem in enumerate(relative_times):
			if ielem >= 0:
				relative_times_reduced.append(ielem)
				ts1_current_array_reduced.append(ts1_current_array[idx])
		w = np.diff(relative_times_reduced, n=1)	
		ts1_current = np.average(ts1_current_array_reduced[:-1], weights = w)
		# LF field
		date_time_array = map(str,np.array(logdat.val10))
		relative_times = []
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)		
		field_array = map(float,np.array(logdat.val11))
		relative_times_reduced = []
		field_array_reduced = []
		for idx, ielem in enumerate(relative_times):
			if ielem >= 0:
				relative_times_reduced.append(ielem)
				field_array_reduced.append(field_array[idx])
		w = np.diff(relative_times_reduced, n=1)	
		field = np.average(field_array_reduced[:-1], weights = w)
		# TF field
		date_time_array = map(str,np.array(logdat.val12))
		relative_times = []
		for date_time in date_time_array:
			relative_times.append(time.mktime(time.strptime(date_time,"%Y-%m-%dT%H:%M:%S.000000000+0100")) - time0 -3600)		
		t20_array = map(float,np.array(logdat.val13))
		relative_times_reduced = []
		t20_array_reduced = []
		for idx, ielem in enumerate(relative_times):
			if ielem >= 0:
				relative_times_reduced.append(ielem)
				t20_array_reduced.append(t20_array[idx])
		w = np.diff(relative_times_reduced, n=1)	
		t20 = np.average(t20_array_reduced[:-1], weights = w)
		# Add obtained sample logs in the log table
		if number_of_frames == 0: number_of_frames=1
		log_table.addRow([n, (total_counts / number_of_frames)*40*3600, slits_avg, ts1_current, field, t20])
		# Integrated count for each detector. Normalised to frame. Add it in the table for count rate
		rnum_detcounts = [n]
		for diclabel in ring_dic:
			rnum_detcounts = rnum_detcounts + map(lambda x: x/number_of_frames, logdat.dicsum[diclabel])
		det_table.addRow(rnum_detcounts)
		
		print("Run ",n," processed")
	except:
		print("Run ",n," was not processed")
