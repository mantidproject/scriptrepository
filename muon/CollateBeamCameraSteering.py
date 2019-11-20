from __future__ import print_function
import time
import datetime
import math
import numpy

#RampLogs=CopyLogsToHist(FirstFile='//hifi/data/hifi00100175.nxs', LastFile='//hifi/data/hifi00100206.nxs', OutputWS='RampLogs', OmitUnmonitoredSpectra=True, ExtraLogs='Field_Main,Slits')
#FieldLog=ExtractSingleSpectrum(RampLogs,32)
#SlitLog=ExtractSingleSpectrum(RampLogs,33)

muonfile="//hifi/data/HIFI%08d.nxs"
camfile="//ndlt626/Users/Public/Documents/sxvh9usb/Tune Sept2016/IMG%d.FIT"
field=-4
if(field==0):
	# zero field
	piclist= list(range(2820,3013))
	runlist=list(range(108048,108068))
	resultname="steeringZF2"
	matrixname="steerMatZF2"
elif(field==2500):
	piclist= list(range(3013,3156))
	runlist=list(range(108068,108088))
	resultname="steering2500"
	matrixname="steerMat2500"
elif(field==5000):
	piclist= list(range(3170,3278))
	runlist=list(range(108088,108108))
	resultname="steering5000"
	matrixname="steerMat5000"
elif(field==7500):
	piclist= list(range(3278,3367))
	runlist=list(range(108108,108128))
	resultname="steering7500"
	matrixname="steerMat7500"
elif(field==10000):
	piclist= list(range(3380,3478))
	runlist=list(range(108128,108148))
	resultname="steering10000"
	matrixname="steerMat10000"
elif(field==12500):
	piclist= list(range(3478,3575))
	runlist=list(range(108148,108168))
	resultname="steering12500"
	matrixname="steerMat12500"
elif(field==15000):
	piclist= list(range(3585,3671))
	runlist=list(range(108168,108188))
	resultname="steering15000"
	matrixname="steerMat15000"
elif(field==-1):
	piclist= list(range(3722,3768))
	runlist=list(range(108188,108200))
	resultname="spotVsSlits"
	matrixname="spotVsSlitsMat"
elif(field==-2):
	piclist= list(range(3782,3844))
	runlist=list(range(108200,108204))
	resultname="spotVsTFX"
	matrixname="spotVsTFXMat"
elif(field==-3):
	piclist= list(range(3843,3883))
	runlist=list(range(108204,108208))
	resultname="spotVsTFY"
	matrixname="spotVsTFYMat"
elif(field==-4):
	piclist= list(range(0,0))
	runlist=list(range(108314,108323))
	resultname="FlypastSlits40SEP"
	matrixname="FlypastMat"
# full fitted results
# matrix of gradients
# rows = magnets
# columns =name, X0, Y0, dX/dMag,dY/dMag
results={}

for run in runlist:
	ww=LoadMuonNexus(Filename=muonfile % run)
	ws=ww[0]
	thisres={}
	if(run==runlist[0]):
		begin=datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:3])) # start of day
	events=0
	for i in range(ws.getNumberHistograms()):
		events+=numpy.sum(ws.dataY(i))
	thisres["rate"]=events/ws.getRun().getProperty("goodfrm").value
	events=0
	for i in range(ws.getNumberHistograms()/2):
		events+=numpy.sum(ws.dataY(i))
	thisres["rateF"]=events/ws.getRun().getProperty("goodfrm").value
	events=0
	for i in range(ws.getNumberHistograms()/2,ws.getNumberHistograms()):
		events+=numpy.sum(ws.dataY(i))
	thisres["rateB"]=events/ws.getRun().getProperty("goodfrm").value
	thisres["begin"]=(datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
	thisres["end"]=(datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_end").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
	com=ws.getComment()
	try:
		c1234=com.split("VSM1=")[1]
		c12,c34=c1234.split("HSM=")
		c1,c2=c12.split("VSM2=")
		c3,c4=c34.split("SEP=")
		thisres["VSM1"]=float(c1)
		thisres["VSM2"]=float(c2)
		thisres["HSM"]=float(c3)
		thisres["SEP"]=float(c4)
	except:
		thisres["VSM1"]=0.0
		thisres["VSM2"]=0.0
		thisres["HSM"]=0.0
		thisres["SEP"]=0.0		
	thisres["Run"]=run
	thisres["Field"]=ws.getRun().getProperty("sample_magn_field").value
	thisres["Slits"]=ws.getRun().getProperty("Slits").value[-1]
	results[run]=thisres
print(results)
#raise Exception("incomplete")

tt=WorkspaceFactory.createTable()
columns=["Run","Image","Field","VSM1","VSM2","HSM","SEP","Slits","rate","X0","eX0","Y0","eY0","XSig","eXSig","YSig","eYSig","Skew","eSkew","Background","eBackground","Intens","eIntens","Major","Minor"]
if(field==-4):
    columns=["Run","Field","VSM1","VSM2","HSM","SEP","Slits","rate","rateF","rateB"]
for co in columns:
	if(co=="Run" or co=="Image"):
		tt.addColumn("int",co,6)
	elif(co[0]=="e"):
		tt.addColumn("double",co,5)
	else:
		tt.addColumn("double",co,2)
		


for i in piclist: # range(5,405):
	img=LoadBeamCameraFile(SourceFile=(camfile % i), FilterNoiseLevel='50', BinSize='1', OutputWorkspace='IMG')
	starttime=(datetime.datetime(*(time.strptime(img.getRun().getProperty("start_time").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
	endtime=(datetime.datetime(*(time.strptime(img.getRun().getProperty("end_time").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
	
	for rn,r in list(results.items()):
		if(starttime>r["begin"] and endtime<r["end"]):
			print("image",i,"taken during run",rn)
		
			squash=MaskAndFlattenCameraImage(img,120,590,80,440)
			
			try:
				fplist=Fit(Function="name=TwoDimGaussianPlusBG",
					InputWorkspace='squash',
					Output='squash',
					OutputCompositeMembers=True,
					StartX=0,
					EndX=100000000,
					OutputNormalisedCovarianceMatrix='squash_NormalisedCovarianceMatrix',
					OutputParameters='squash_Parameters',
					OutputWorkspace='squash_Workspace',
					MaxIterations=20)
				
				#print fplist
				(YesNo, CostFn, CovarMat, Params,Workspace)=fplist
				(EllipseCurves,Major,Minor,Phi)=DrawGaussianEllipse(Params)
				(X0,Y0,XSig,YSig,Skew,Background,Intens,CostF2)=Params.column(1)
				r["X0"]=X0
				r["Y0"]=Y0
				r["XSig"]=XSig
				r["YSig"]=YSig
				r["Skew"]=Skew
				r["Background"]=Background
				r["Intens"]=Intens
				(eX0,eY0,eXSig,eYSig,eSkew,eBackground,eIntens,eCostF2)=Params.column(2)
				r["eX0"]=eX0
				r["eY0"]=eY0
				r["eXSig"]=eXSig
				r["eYSig"]=eYSig
				r["eSkew"]=eSkew
				r["eBackground"]=eBackground
				r["eIntens"]=eIntens
				r["Major"]=Major
				r["Minor"]=Minor
				r["Image"]=i
				SFac=math.sqrt(1-Skew**2/4)
				r["PeakInt"]=Intens/XSig/YSig*SFac
				r["XOverall"]=XSig/SFac
				r["YOverall"]=YSig/SFac
			except:
				print("fit failed for image ",i)
				# don't add to r, another image might be ok

#print results

for r in list(results.values()):
	if ("Image" in r) or (field==-4):
		row=[r[x] for x in columns]
		tt.addRow(row)

AnalysisDataService.addOrReplace(resultname,tt)

# to test: select from here down, Execute Selection
try:
	print(len(results))
except:
	mm=mtd["steeringZF2"]
	results={}
	for r in mm:
		print("Row->",r)
		results[r["Run"]]=r

print(results)

# find sets
if(field>=0):
	mags=["VSM1","VSM2","HSM","SEP"]
else:
	mags=["Field"]
tt2=WorkspaceFactory.createTable()
tt2.addColumn("str","Magnet")
tt2.addColumn("double","StartVal")
tt2.addColumn("double","X0")
tt2.addColumn("double","Y0")
tt2.addColumn("double","dXdMag")
tt2.addColumn("double","dYdMag")
for mag in mags:
	constmags=list(mags)
	del constmags[constmags.index(mag)]
	sorter={}
	common={}
	for r in list(results.values()):
		if "X0" in r:
			key=tuple([r[x] for x in constmags])
			if key not in sorter:
				sorter[key]={}
			sorter[key][r[mag]]=r
			if(r[mag] in common):
				common[r[mag]]+=1
			else:
				common[r[mag]]=1
		else:
			print("run",r["Run"],"doesn't have a successfully fitted image")
	magref=list(common.keys())[numpy.argmax(list(common.values()))]
	print(mag,"relative to ",magref)
	for (sv,s) in list(sorter.items()):
		if(len(s)>2):
			# got a set of points
			tmpM=numpy.array([x[mag]-magref for x in list(s.values())])
			print("found set of points for",constmags,"=",sv,"with",mag,"=",tmpM)
			tmpX=numpy.array([x["X0"] for x in list(s.values())])
			tmpXe=numpy.array([x["eX0"] for x in list(s.values())])
			tmpY=numpy.array([x["Y0"] for x in list(s.values())])
			tmpYe=numpy.array([x["eY0"] for x in list(s.values())])
			CreateWorkspace(DataX=tmpM,DataY=tmpX,DataE=tmpXe,NSpec=1,OutputWorkspace="tmpgrad")
			fr=Fit(Function="name=LinearBackground",InputWorkspace="tmpgrad",Output="Xvs"+mag)
			X0=fr[3].cell(0,1)
			gradX=fr[3].cell(1,1)
			DeleteWorkspace("tmpgrad")
			CreateWorkspace(DataX=tmpM,DataY=tmpY,DataE=tmpYe,NSpec=1,OutputWorkspace="tmpgrad")
			fr=Fit(Function="name=LinearBackground",InputWorkspace="tmpgrad",Output="Yvs"+mag)
			Y0=fr[3].cell(0,1)
			gradY=fr[3].cell(1,1)
			DeleteWorkspace("tmpgrad")
			print(mag,magref,X0,Y0,gradX,gradY))
			tt2.addRow((mag,magref,X0,Y0,gradX,gradY))
AnalysisDataService.addOrReplace(matrixname,tt2)

#del(results)
