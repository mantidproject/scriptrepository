"""
Quicly look at powder data in NXSPE files
"""
# dat_runs= range(21796,21801)#Empty Can 34K
#dir="/SNS/SEQ/IPTS-6154/shared/empty_can/"

Temp = [34,150,250,320,430]
#Temp = [34,320]
f_pre="base_120meV"

# simple loop to load  data sets
for T in Temp:
	dir = "/SNS/SEQ/IPTS-6154/shared/{0}K/".format(T) 
	f_pre = "{0}K_120meV".format(T)
	fname = f_pre +'_120p0_Full'
	
	#load individual data sets
	LoadNXSPE(Filename= dir + fname +".nxspe",OutputWorkspace = "IWS")
	CloneWorkspace(InputWorkspace="IWS",OutputWorkspace=  "Data_Full_{0}K".format(T))
	#Bin Energy
	Rebin(InputWorkspace="IWS",OutputWorkspace="IWS",Params= "0,0.5,70")
	#Covert to SQW
	SofQW3(InputWorkspace = "IWS",OutputWorkspace = "IWS",QAxisBinning = "0.0,0.05,5",EMode ="Direct")
	#Tranpose (so Q is on x-axis in plots)
	Transpose(InputWorkspace = "IWS",OutputWorkspace = 'IWS')
	#SaveNXSPE(InputWorkspace = "IWS",Filename = "/SNS/SEQ/IPTS-6154/shared/SQW_{0}K.nxspe".format(T))
	#Bin Q
	CloneWorkspace(InputWorkspace="IWS",OutputWorkspace= "SQW_{0}K".format(T))
	#Apply smoothing (for nice looking data)
	SmoothData(InputWorkspace="IWS", OutputWorkspace="IWS", NPoints = 3)
	CloneWorkspace(InputWorkspace="IWS",OutputWorkspace= "SQW_{0}K_smooth".format(T))
	
	DeleteWorkspace(Workspace = "IWS")=
