# import mantid algorithms, numpy and matplotlib
import time
import datetime
import math
import numpy as np
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

class LostTimeWaiting(PythonAlgorithm):
    
    def PyInit(self):
        self.declareProperty(FileProperty(name="FirstFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
        self.declareProperty(FileProperty(name="LastFile",defaultValue="",action=FileAction.OptionalLoad,extensions = ["nxs"]))
        self.declareProperty("LogForRunStatus","run_status",doc="Block name")

    def category(self):
        return "Muon"
        
    def PyExec(self):
        
        #ws=self.getProperty("InputWS").value
        file1=self.getProperty("FirstFile").value
        file9=self.getProperty("LastFile").value
        if(file1 != file9):
            i1=file1.rindex('.')
            j1=i1-1
            while file1[j1-1].isdigit():
                j1=j1-1
            firstnum=int(file1[j1:i1])
            i9=file9.rindex('.')
            j9=i9-1
            while file9[j9-1].isdigit():
                j9=j9-1
            lastnum=int(file9[j9:i9])
            if(file1[:j9] != file9[:j9]):
                raise Exception("Files from different directories or instruments")
            if(file1[i1:] != file9[i9:]):
                raise Exception("Files of different types")
            if(i1-j1 != i9-j9):
                raise Exception("File numbering error")
            if(lastnum < firstnum):
                if(firstnum != 0):
                    raise Exception("Run numbers must increase")
                else:
                    # run zero: load temporary file instead
                    lastnum=10^(i1-j1)-1 # highest possible number as place holder
            #        runpaths.append(firstfn[:j]+str(n+firstnum).zfill(i-j)+firstfn[i:])
        else: # load any single file even if it's not numbered, eg Temporary File.
            firstnum=0
            lastnum=0
        # totals
        listWaiting=[]
        totalWaiting=0

        # loop and load files. Absolute numbers for now.
        gotToCurrentRun=False
        prog_reporter = Progress(self,start=0.0,end=1.0,nreports=lastnum+1-firstnum)
        for ff in range(firstnum,lastnum+1):
            if(file1 != file9):
                thispath=file1[:j1]+str(ff).zfill(i1-j1)+file1[i1:]
            else:
                thispath=file1
            # cope with Periods - have Workspace Group and extra return values from LoadMuonNexus
            try:
                wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
            except:
                thispath=file1[:j1]+"auto_A.tmp" # should really try all of them and pick newest!
                wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
                gotToCurrentRun=True
            try:
                ws=wstuple[0][0] # first period if there are some
            except:
                ws=wstuple[0] # plain run
            prog_reporter.report("Loading")
            if(gotToCurrentRun and int(ws.getRun().getProperty("run_number").value) != ff):
                print("temporary file is for the previous run, probably no new one started")
                break
            log=ws.getRun().getLogData(self.getProperty("LogForRunStatus").value)
            N=len(log.times)
            waits=0
            for i in range(1,N-1):
                if log.value[i]==4 and (log.times[i+1]-log.times[i-1])/1.E9<300: # differences in ns
                    waits=waits+(log.times[i+1]-log.times[i-1])/2.E9
            DeleteWorkspace(ws)
            if waits>0:
                totalWaiting += float(waits)
                listWaiting.append((ff,float(waits)))
        print ("Summary")
        for x in listWaiting:
            print ("Run",x[0],"waiting for",x[1],"seconds")
        print ("total waits = ",totalWaiting,"seconds, or",totalWaiting/3600.0,"hours")
        
AlgorithmFactory.subscribe(LostTimeWaiting)
