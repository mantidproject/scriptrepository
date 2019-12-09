from __future__ import print_function
import math
import numpy

class FitIndividualPhases(PythonAlgorithm):
    def PyInit(self):
        self.declareProperty(WorkspaceProperty("InputWorkspace","",direction=Direction.Input))
        self.declareProperty(FloatArrayProperty("FreqGuesses",values=[6.28*2*math.pi,6.76*2*math.pi]),doc="Angular frequency guess(es)")
        self.declareProperty("StartTime",0.1)
        self.declareProperty("EndTime",10.0)
        self.declareProperty(FileProperty("GroupingTable","M:/misc/beam/RIKEN/CHRONUS_Detector_Grouping_LFTF.xml",action=FileAction.Load,extensions=[".xml"]))
        self.declareProperty("ForwardGroupIndex",5)
        self.declareProperty("BackwardGroupIndex",6)
        self.declareProperty(FileProperty("InstrumentDefinition","M:/misc/beam/RIKEN/CHRONUS_Definition_JSL_2.xml",action=FileAction.Load,extensions=[".xml"]))
        self.declareProperty(WorkspaceProperty("OutputWorkspace","",direction=Direction.Output))
        self.declareProperty(FloatArrayProperty("FittedFrequencies",direction=Direction.Output))

    def category(self):
        return "Muon"

    def PyExec(self):
        inputWS=self.getProperty("InputWorkspace").value
        gw=LoadDetectorsGroupingFile(InputFile=self.getProperty("GroupingTable").value)
        fgi=self.getProperty("ForwardGroupIndex").value
        bgi=self.getProperty("BackwardGroupIndex").value
        fsl=[]
        bsl=[]
        if gw.getNumberHistograms() != inputWS.getNumberHistograms():
            print("oops, grouping for wrong instrument?")
        for i in range(gw.getNumberHistograms()):
            if gw.dataY(i)[0]==fgi:
                fsl.append(gw.getDetector(i).getID())
            elif gw.dataY(i)[0]==bgi:
                bsl.append(gw.getDetector(i).getID())
            #else:
            #    print "ignoring detector",i,"in group",gw.dataY(i)[0]
        FreqGuesses=self.getProperty("FreqGuesses").value # angular freqs please
        # initial fit to optimise frequencies
        NF=len(FreqGuesses)
        if NF==1:
            fitfn0="name=UserFunction,Formula=bg+a1*cos(w1*x+p1), a1=0.2,p1=0.1,bg=0.0,w1={0}".format(*FreqGuesses)
            fitfn="name=UserFunction,Formula=n0*(1+a1*cos({w1}*x+p1))*exp(-x/2.197), a1={a1},p1={p1},n0={n0}"
        elif NF==2:
            fitfn0="name=UserFunction,Formula=bg+a1*cos(w1*x+p1)+a2*cos(w2*x+p2), a1=0.1,a2=0.1,p1=0.1,p2=0.1,bg=0.0,w1={0},w2={1}".format(*FreqGuesses)
            fitfn="name=UserFunction,Formula=n0*(1+a1*cos({w1}*x+p1)+a2*cos({w2}*x+p2))*exp(-x/2.197), a1={a1},a2={a2},p1={p1},p2={p2},n0={n0}"
        elif NF==3:
            fitfn0="name=UserFunction,Formula=bg+a1*cos(w1*x+p1)+a2*cos(w2*x+p2)+a3*cos(w3*x+p3), a1=0.05,a2=0.05,a3=0.05,p1=0.1,p2=0.1,p3=0.1,bg=0.0,w1={w1},w2={w2},w3={w3}".format(*FreqGuesses)
            fitfn="name=UserFunction,Formula=n0*(1+a1*cos({w1}*x+p1)+a2*cos({w2}*x+p2)+a3*cos({w3}*x+p3))*exp(-x/2.197), a1={a1},a2={a2},a3={a3},p1={p1},p2={p2},p3={p3},n0={n0}"
        asy=AsymmetryCalc(InputWorkspace=inputWS,ForwardSpectra=fsl,BackwardSpectra=bsl,Alpha=1.0)

        t1=self.getProperty("StartTime").value
        t2=self.getProperty("EndTime").value

        fr0=Fit(Function=fitfn0,InputWorkspace=asy,StartX=t1,EndX=t2)
        #iv=dict(fr0.Function)
        iv={}
        fittedw=[]
        for i in range(1,NF+1):
            iv["a"+str(i)]=fr0.Function["a"+str(i)]
            iv["p"+str(i)]=fr0.Function["p"+str(i)]
            iv["w"+str(i)]=fr0.Function["w"+str(i)]
            fittedw.append(fr0.Function["w"+str(i)])
        for k in list(iv.keys()):
            if k[0]=="w":
                print("fitted freq",iv[k],"or",iv[k]/2/math.pi/0.01355,"G with ampl=",iv["a"+k[1:]])

        # bulk fit TF and re-generate workspace
        NS=inputWS.getNumberHistograms()
        
        #fitfn="name=UserFunction,Formula=n0*(1+a1*cos({w1}*x+p1)+a2*cos({w2}*x+p2))*exp(-x/2.197), a1={a1},a2={a2},p1={p1},p2={p2},n0={n0}"
        # w1g and w2g determined by prior fitting
        #iv={"w1":6.761*2*math.pi,"w2":6.286*2*math.pi,"a1":0.07,"a2":0.02,"p1":0.1,"p2":0.1}


        outWS=WorkspaceFactory.create("Workspace2D",NVectors=NS,XLength=NF*2+2,YLength=NF*2+1) # x bins n1,p1,n2,p2,nnr from 0 to 5
        empties=[]
        for i in range(NS):
            t1bin=int(numpy.searchsorted(inputWS.dataX(i),t1))
            hs=numpy.sum(inputWS.dataY(i)[t1bin:])
            print("histogram",i,"sum of good counts=",hs)
            if(hs>100):
                dt=inputWS.dataX(i)[1]-inputWS.dataX(i)[0]
                iv["n0"]=hs*dt/2.197*math.exp(t1/2.197)
                fn=fitfn.format(**iv)
                #print fn
                fr=Fit(Function=fn,InputWorkspace=inputWS,WorkspaceIndex=i,StartX=t1,EndX=t2,CreateOutput=True)
                #print fr
                ft=fr.OutputParameters
                ed={x["Name"]:x["Error"] for x in ft}
                for j in range(NF):
                    a=fr.Function["a"+str(j+1)]
                    p=fr.Function["p"+str(j+1)]
                    if(a<0):
                        a=-a
                        p=p+math.pi
                    p=math.fmod(p+8*math.pi,2*math.pi)
                    outWS.dataY(i)[j*2]=a
                    outWS.dataE(i)[j*2]=ed["a"+str(j+1)]
                    outWS.dataY(i)[j*2+1]=p
                    outWS.dataE(i)[j*2+1]=ed["p"+str(j+1)]
                outWS.dataX(i)[:]=list(range(2*NF+2))
                outWS.dataY(i)[2*NF]=fr.Function["n0"]
                outWS.dataE(i)[2*NF]=ed["n0"]
            else:
                outWS.dataX(i)[:]=list(range(2*NF+2))
                outWS.dataY(i)[:]=[0.0]*(2*NF+1)
                empties.append(i)
                

        AnalysisDataService.addOrReplace(self.getProperty("OutputWorkspace").valueAsStr,outWS)
        LoadInstrument(outWS,RewriteSpectraMap=True,Filename=self.getProperty("InstrumentDefinition").value)
        if(empties):
            MaskDetectors(outWS,WorkspaceIndexList=empties)
        # summary
        print("frequencies:")
        for i in range(NF):
            print(FreqGuesses[i]," optimised to ",iv["w"+str(i+1)]," with amplitude",iv["a"+str(i+1)])
        self.setProperty("OutputWorkspace",outWS)
        self.setProperty("FittedFrequencies",fittedw)

AlgorithmFactory.subscribe(FitIndividualPhases)
