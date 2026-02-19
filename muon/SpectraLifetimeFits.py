# import mantid algorithms, numpy and matplotlib
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

class SpectraLifetimeFits(PythonAlgorithm):

    def PyInit(self):
        # Declare properties
        self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input),doc="Raw loaded workspace e.g. from muon interface")
        self.declareProperty("FastLifetimeGuess",0.2,doc="Hint for fit")
        self.declareProperty("SlowLifetimeGuess",2.19,doc="Hint for fit")
        self.declareProperty("FirstGoodData",0.1,doc="Start time for fit, may be momentum dependent")
        self.declareProperty(WorkspaceProperty("AmplitudeRatios","",Direction.Output),doc="Ratio of slow to fast components")
        self.declareProperty(WorkspaceProperty("BackgroundRatios","",Direction.Output),doc="Ratio of background to fast component")

    def category(self):
        return 'Muon'

# lifetime analysis by spectrum

# expects two relaxing components (hign Z and low Z muonic atoms) and a "non relaxing" background

    def PyExec(self):
        lifetimeguesses=[self.getProperty("FastLifetimeGuess").value,self.getProperty("SlowLifetimeGuess").value]
        ampguesses=[30000,80]
        flatbg=True


        startX=self.getProperty("FirstGoodData").value
        endX=32.0

        N=len(lifetimeguesses)


        funcbits=["name=ExpDecay,Height={:f},Lifetime={:f}".format(ampguesses[i],lifetimeguesses[i]) for i in range(N)]
        if flatbg:
            funcbits.append("name=FlatBackground,A0={:f}".format(ampguesses[-1]))
        Func=";".join(funcbits)

        Ties=",".join(["f{:d}.Lifetime={:f}".format(i,lifetimeguesses[i]) for i in range(N)])

        #w=mtd["{inst:s}{run:d}_raw_data MA".format(inst=inst,run=run)]
        w=self.getProperty("InputWorkspace").value
        Nspec=w.getNumberHistograms()

        w2=WorkspaceFactory.create(w,XLength=2,YLength=1,NVectors=Nspec) # for amp1 / amp0
        w3=WorkspaceFactory.create(w,XLength=2,YLength=1,NVectors=Nspec) # for bg / amp0

        # stage 1, fit summed spectra
        # stage 1a, with fixed lifetimes to get amplitudes about right

        sum=SumSpectra(InputWorkspace=w,OutputWorkspace="sum")

        fr1=Fit(Function=Func,Ties=Ties,InputWorkspace="sum",WorkspaceIndex=0,Minimizer="Simplex",CostFunction="Poisson",StartX=startX,EndX=endX,CreateOutput=True,MaxIterations=5000,Output="sum1") # ,OutputParametersOnly=True

        #print (fr1.Function["f0.Height"])

        heightfits=[fr1.Function["f{:d}.Height".format(i)] for i in range(N)]
        self.log().notice("1st stage height fits were "+str(heightfits))
        if flatbg:
            self.log().notice("1st stage background was "+str(fr1.Function["f{:d}.A0".format(N)]))

        # stage 1b, full fit with optimised amplitudes

        func2bits=["name=ExpDecay,Height={:f},Lifetime={:f}".format(heightfits[i],lifetimeguesses[i]) for i in range(N)]
        if flatbg:
            bgfit=fr1.Function["f{:d}.A0".format(N)]
            func2bits.append("name=FlatBackground,A0={:f}".format(bgfit))
        Func2=";".join(func2bits)

        fr2=Fit(Function=Func2,InputWorkspace="sum",WorkspaceIndex=0,Minimizer="Simplex",CostFunction="Poisson",StartX=startX,EndX=endX,CreateOutput=True,MaxIterations=5000,Output="sum2") # ,OutputParametersOnly=True

        heightfits=[fr2.Function["f{:d}.Height".format(i)] for i in range(N)]
        lifetimefits=[fr2.Function["f{:d}.Lifetime".format(i)] for i in range(N)]
        self.log().notice("2nd stage height fits were "+str(heightfits))
        if flatbg:
            self.log().notice("2nd stage background was "+str(fr2.Function["f{:d}.A0".format(N)]))
        self.log().notice("2nd stage lifetime fits were "+str(lifetimefits))

        # stage 3, spectrum by spectrum

        func3bits=["name=ExpDecay,Height={:f},Lifetime={:f}".format(heightfits[i]/Nspec,lifetimefits[i]) for i in range(N)]
        if flatbg:
            bgfit=fr2.Function["f{:d}.A0".format(N)]
            func3bits.append("name=FlatBackground,A0={:f}".format(bgfit/Nspec))
        Func3=";".join(func3bits)
        Ties3=",".join(["f{:d}.Lifetime={:f}".format(i,lifetimefits[i]) for i in range(N)])

        for i in range(Nspec):
            if np.sum(w.dataY(i))>100:
                fr3=Fit(Function=Func3,InputWorkspace=w,Ties=Ties3,WorkspaceIndex=i,Minimizer="Simplex",CostFunction="Poisson",StartX=startX,EndX=endX,CreateOutput=True) # ,OutputParametersOnly=True
                heightfits=[fr3.Function["f{:d}.Height".format(i)] for i in range(N)]
                if bgfit:
                    heightfits.append(fr3.Function["f{:d}.A0".format(N)])
                self.log().debug(str(i)+" "+str(heightfits))
                w2.dataX(i)[0]=0.0
                w2.dataX(i)[1]=1.0
                w2.dataY(i)[0]=heightfits[1]/heightfits[0]
                w3.dataX(i)[0]=0.0
                w3.dataX(i)[1]=1.0
                w3.dataY(i)[0]=heightfits[-1]/heightfits[0]
            else:
                self.log().warning(str(i)+" is empty-ish")
                w2.dataX(i)[0]=0.0
                w2.dataX(i)[1]=1.0
                w2.dataY(i)[0]=0.0
                w3.dataX(i)[0]=0.0
                w3.dataX(i)[1]=1.0
                w3.dataY(i)[0]=0.0

        self.log().notice("lifetime guesses were "+str(lifetimeguesses))
        self.log().notice ("lifetime fits were "+str(lifetimefits))
        self.log().debug ("individual ties were "+str(Ties3))
        self.log().debug ("height fits were "+str(heightfits))
        self.log().debug ("last fit:")
        self.log().debug (str(fr3.Function))

        #AnalysisDataService.addOrReplace("AmpRatio{:d}".format(run),w2)
        #AnalysisDataService.addOrReplace("BGRatio{:d}".format(run),w3)
        self.setProperty("AmplitudeRatios",w2)
        self.setProperty("BackgroundRatios",w3)


AlgorithmFactory.subscribe(SpectraLifetimeFits)
