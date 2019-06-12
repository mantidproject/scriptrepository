# Calibrate muon dead times with "direct" fit
import numpy

class CalMuonDeadTime2(PythonAlgorithm):
    def PyInit(self):
        self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input),doc="Name of the input workspace")
        self.declareProperty(ITableWorkspaceProperty("DeadTimeTable","",Direction.Output),doc="The name of the TableWorkspace in which to store the dead times for each spectrum")
        self.declareProperty("FirstGoodData",0.5,doc="start time for fit in microseconds")
        self.declareProperty("LastGoodData",5.0,doc="end time for fit in microseconds")
        self.declareProperty("Minimiser","LeastSq",validator=StringListValidator(["LeastSq","Poisson"]))

    def category(self):
        return "Muon"
        
    def PyExec(self):
        w=self.getProperty("InputWorkspace").value
        t1=self.getProperty("FirstGoodData").value
        t2=self.getProperty("LastGoodData").value
        tmuon=2.19703 # could fetch from built in table of constants
        goodfrm=w.getRun().getProperty("goodfrm").value
        deltaT=w.dataX(0)[2]-w.dataX(0)[1] # assume uniform and the same for all spectra!
        mi=self.getProperty("Minimiser").value
        if mi=="LeastSq":
            costfn="Least squares"
            minim="Levenberg-MarquardtMD"
        elif mi=="Poisson":
            costfn="Poisson"
            minim="Levenberg-MarquardtMD"
        else:
            raise ValueError("Unknown minimiser choice")
        
        dtt=WorkspaceFactory.createTable()
        dtt.addColumn("int","Spectrum",1)
        dtt.addColumn("double","dead-time",2)
        
        for i in range(w.getNumberHistograms()):
            # dead time is defined by N=N0 exp(t/tau) = M/(1-M*(deadtime/tbin*frames)) where M is measured counts
            N0guess=numpy.amax(w.dataY(i))
            dguess=0.005/deltaT/goodfrm
            Func="name=UserFunction, Formula=1.0/(1.0/(N0*exp(-x/"+str(tmuon)+"))+d), N0="+str(N0guess)+", d="+str(dguess)
            fr=Fit(Function=Func,InputWorkspace=w,WorkspaceIndex=i,StartX=t1,EndX=t2,Minimizer=minim, CostFunction=costfn)
            N0=fr.Function["N0"]
            d=fr.Function["d"]
            deadtime=d*goodfrm*deltaT # independent of N0
            s=w.getSpectrum(i).getSpectrumNo()
            dtt.addRow([s,deadtime])
        self.setProperty("DeadTimeTable",dtt)

AlgorithmFactory.subscribe(CalMuonDeadTime2)
