"""
function to change initial time for logs to coincide with the first proton pulse
"""
def fixlogs(WS):
	run=mtd[WS].getRun()
	PC =  run['proton_charge'].firstTime()
	for log_name in list(run.keys()):
		if log_name not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
			try:
				P =  run[log_name].firstTime()
				Tdiff = (PC-P).total_milliseconds()*1E-3
				ChangeLogTime(InputWorkspace=WS, OutputWorkspace = WS, LogName = log_name, TimeOffset = Tdiff)
			except:
				pass
		    

