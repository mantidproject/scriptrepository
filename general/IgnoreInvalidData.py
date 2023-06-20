# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import DeleteWorkspace, ConvertToPointData, ConvertToHistogram
from mantid.api import PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty
from mantid.kernel import Direction
import numpy as np


class IgnoreInvalidData(PythonAlgorithm):

    def name(self):
        return "IgnoreInvalidData"

    def summary(self):
        return "Ignore non-finite values in a workspace, by giving them a value of zero and an x value just outside the original data range."

    def PyInit(self):
        self.declareProperty(
            MatrixWorkspaceProperty("InputWorkspace", "", direction=Direction.Input),
            "input workspace to remove invalud values from",
        )
        self.declareProperty(
            MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
            "output workspace, will have no invalid values",
        )

    def PyExec(self):
        ws = self.getProperty("InputWorkspace").value  
        
        # need point data to do the check
        is_histogram = ws.isHistogramData()
        if is_histogram:
            ws = ConvertToPointData(ws)
            
        # loop over spectra
        for spec in range(ws.getNumberHistograms()):
            x = ws.dataX(spec)
            y = ws.dataY(spec)
            e = ws.dataE(spec)

            # move the NAN data to be at the end of the data range with a y value of 0
            large = np.max(x)
            dx = 1.e-6
            for j in range(len(x)-1,0, -1):
                if not np.isfinite(y[j]):
                    x[j] = large + dx
                    y[j] = 0.0
                    e[j] = 0.0
                    dx += 1.e-6
            
        # preserve the original distribution for the data
        if is_histogram:
            ws = ConvertToHistogram(ws)
        # output result and clean up
        self.setProperty("OutputWorkspace", ws)
        DeleteWorkspace(ws)

AlgorithmFactory.subscribe(IgnoreInvalidData)

