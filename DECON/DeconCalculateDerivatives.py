from mantid.simpleapi import (mtd, DeleteWorkspace, ExtractSingleSpectrum, Minus,
                              WienerSmooth, SplineSmoothing, AppendSpectra)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode, MatrixWorkspaceProperty,
                        TextAxis, Progress)
from mantid.kernel import (Direction, logger)

class DeconCalculateDerivatives(DataProcessorAlgorithm):
 
    _deriv_name	= None
    _input_ws  = None

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Smooth input in theta and calulate the derivatives."

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('DataWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Data workspace.")
        self.declareProperty(MatrixWorkspaceProperty('DerivativesWorkspace', 'Deriv',
                                                     direction = Direction.Output),
                             doc = "Name for the output Derivatives workspace.")
 
    def PyExec(self):
        self._setup()
        self._derivatives()
 
    def _setup(self):
        self._theta = self.getPropertyValue('DataWorkspace')
        self._deriv = self.getPropertyValue('DerivativesWorkspace')

    def _derivatives(self):

        #smooth data using FFT Wiener filter.
        smooth_prog = Progress(self, start=0.0, end=0.1, nreports=3)
        smooth_prog.report('Smoothing data ')

        self._weiner_smooth(self._theta, '__smooth')
        self._minus(self._theta, '__smooth', '__diff')
        self._append(self._theta, '__smooth', self._theta)
        self._append(self._theta, '__diff', self._theta)

        ax = TextAxis.create(3)
        for i, x in enumerate(['Input', 'Smooth', 'Diff']):
            ax.setLabel(i, x)
        mtd[self._theta].replaceAxis(1, ax)
        smooth_prog.report('Smoothing completed')

        #use smooth data to calculate derivatives using cubic spline.
        self._spline_smooth('__smooth', '__spline0', '__deriv1') 

        #use second derivative as input to calc higher derivatives
        self._extract('__deriv1_1', '__spline1', 1)
        self._spline_smooth('__spline1', '__spline2', '__deriv2') 

        # create the derivs WS
        create_prog = Progress(self, start=0.0, end=0.1, nreports=3)
        create_prog.report('Creating workspaces ')

        self._extract('__deriv1_1', self._deriv, 0)
        self._extract('__deriv1_1', '__deriv', 1)
        self._append(self._deriv, '__deriv', self._deriv)
        self._extract('__deriv2_1', '__deriv', 0)
        self._append(self._deriv, '__deriv', self._deriv)
        self._extract('__deriv2_1', '__deriv', 1)
        self._append(self._deriv, '__deriv', self._deriv)

        ax = TextAxis.create(4)
        for i, x in enumerate(['Deriv.1', 'Deriv.2', 'Deriv.3', 'Deriv.4']):
            ax.setLabel(i, x)
        mtd[self._deriv].replaceAxis(1, ax)

        workspaces = ['__smooth', '__diff', '__spline0', '__spline1', '__spline2', \
                      '__deriv1', '__deriv2', '__deriv']
        for ws in workspaces:
            self._delete_ws(ws)
        self.setProperty("DerivativesWorkspace", self._deriv)

    def _minus(self, lhs_ws, rhs_ws, output_ws):
        multiply_alg = self.createChildAlgorithm("Minus", enableLogging=False)
        multiply_alg.setProperty("LHSWorkspace", lhs_ws)
        multiply_alg.setProperty("RHSWorkspace", rhs_ws)
        multiply_alg.setProperty("OutputWorkspace", output_ws)
        multiply_alg.execute()
        mtd.addOrReplace(output_ws, multiply_alg.getProperty("OutputWorkspace").value)

    def _delete_ws(self, input_ws):
        delete_alg = self.createChildAlgorithm("DeleteWorkspace", enableLogging=False)
        delete_alg.setProperty("Workspace", input_ws)
        delete_alg.execute()

    def _weiner_smooth(self, input_ws, output_ws):
        wiener_smooth_alg = self.createChildAlgorithm("WienerSmooth", enableLogging = False)
        wiener_smooth_alg.setProperty("InputWorkspace", input_ws)
        wiener_smooth_alg.setProperty("WorkspaceIndexList", '')
        wiener_smooth_alg.setProperty("OutputWorkspace", output_ws)
        wiener_smooth_alg.execute()
        mtd.addOrReplace(output_ws, wiener_smooth_alg.getProperty("OutputWorkspace").value)
        mtd[output_ws].setDistribution(True)

    def _spline_smooth(self, input_ws, output_ws, deriv_ws):
        spline_smooth_alg = self.createChildAlgorithm("SplineSmoothing", enableLogging = False)
        spline_smooth_alg.setProperty("InputWorkspace", input_ws)
        spline_smooth_alg.setProperty("OutputWorkspace", output_ws)
        spline_smooth_alg.setProperty("OutputWorkspaceDeriv", deriv_ws)
        spline_smooth_alg.setProperty("DerivOrder", 2)
        spline_smooth_alg.setProperty("MaxNumberOfBreaks", 0)
        spline_smooth_alg.execute()
        mtd.addOrReplace(output_ws, spline_smooth_alg.getProperty("OutputWorkspace").value)
        mtd.addOrReplace(deriv_ws, spline_smooth_alg.getProperty("OutputWorkspaceDeriv").value)

    def _extract(self, input_ws, output_ws, index):
        extract_alg = self.createChildAlgorithm("ExtractSingleSpectrum", enableLogging = False)
        extract_alg.setProperty("InputWorkspace", input_ws)
        extract_alg.setProperty("OutputWorkspace", output_ws)
        extract_alg.setProperty("WorkspaceIndex", index)
        extract_alg.execute()
        mtd.addOrReplace(output_ws, extract_alg.getProperty("OutputWorkspace").value)

    def _append(self, input1_ws, input2_ws, output_ws):
        append_alg = self.createChildAlgorithm("AppendSpectra", enableLogging = False)
        append_alg.setProperty("InputWorkspace1", input1_ws)
        append_alg.setProperty("InputWorkspace2", input2_ws)
        append_alg.setProperty("OutputWorkspace", output_ws)
        append_alg.execute()
        mtd.addOrReplace(output_ws, append_alg.getProperty("OutputWorkspace").value)

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconCalculateDerivatives)
