from mantid.simpleapi import (mtd, CreateWorkspace, CloneWorkspace, DeleteWorkspace, ExtractSingleSpectrum,
                              Minus, Plus, Multiply, Divide, AddSampleLogMultiple, LoadNexus)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode, MatrixWorkspaceProperty,
                        WorkspaceGroupProperty, FileProperty, FileAction, Progress)
from mantid.kernel import (StringListValidator, StringMandatoryValidator, IntBoundedValidator,
                           FloatBoundedValidator, Direction, logger)
from mantid import config
import math, numpy as np


class DeconD4Result(DataProcessorAlgorithm):

    _input = None
    _stheta = None
    _theta_used = None
    _final_theta = None
    _name = None
    _lambda = 0.7
    _azero = 0.0
    _sofq = None
    _deriv = None
    _path = None
    _result = None
    _final_q = None
    _smooth = False
    _rebin_option = 'None'
    _rebin_qrange = None
    _rebin_qinc = None

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Perform deconvolution."

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('CorrectedWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input S(theta) workspace.")
        self.declareProperty(MatrixWorkspaceProperty('QWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Q workspace.")

        self.declareProperty(name = 'Wavelength', defaultValue = 0.7,
                             doc = 'Wavelength for output')
        self.declareProperty(name = 'ZeroAngle', defaultValue = 0.0,
                             doc = 'Zero angle for Q convertion')

        self.declareProperty(name = 'RebinOption', defaultValue = 'None',
                             validator = StringListValidator(['None','Interpolate','Spline']),
                             doc = 'Rebin options : None or Interpolate')
        self.declareProperty(name = 'RebinQrange', defaultValue = 'New',
                             validator = StringListValidator(['New','Snap']),
                             doc = 'Rebin Qrange : input New values or Snap to input values')
        self.declareProperty(name = 'RebinQinc', defaultValue = 0.1,
                             doc = 'Value of deltaQ for rebinning')

    def PyExec(self):
        self._setup()
        self._convert_result()
        self._rebin_result()

    def _setup(self):
        self._corr = self.getPropertyValue('CorrectedWorkspace')
        self._sofq = self.getPropertyValue('QWorkspace')
        self._final_q = self._sofq + '_corrected'

        self._lambda = self.getProperty('Wavelength').value
        self._azero = self.getProperty('ZeroAngle').value

        self._rebin_option = self.getProperty('RebinOption').value
        self._rebin_qrange = self.getProperty('RebinQrange').value
        self._rebin_qinc = self.getProperty('RebinQinc').value

    def _convert_result(self):                  #cutoff high angle data ouput _final_theta as *_theta_corrected
        convert_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        convert_prog.report('Converting result to Q')

        #convert *_theta_corrected to *_Q_corrected
        k0 = 4.0 * math.pi / self._lambda

        # Create a copy of the theta workspace and convert to Q
        self._clone_ws(self._corr, self._final_q)
        x_q = mtd[self._final_q].dataX(0)
        x_q = k0 * np.sin(0.5 * np.radians(x_q - self._azero))    #convert to Q after applying zero angle correction
        mtd[self._final_q].setX(0, x_q)
        unitx = mtd[self._final_q].getAxis(0).setUnit("MomentumTransfer")

        self._copy_log(self._corr, self._final_q, 'MergeReplaceExisting')
        convert_logs = [('lambda_out', self._lambda), ('zero_out', self._azero)]
        logger.information('Converting : %s ; from theta to Q as : %s' % (self._corr, self._final_q))
        logger.information('lambda = %f ; zero = %f' % (self._lambda, self._azero))
        log_names = [item[0] for item in convert_logs]
        log_values = [item[1] for item in convert_logs]
        self._add_sample_log_mult(self._final_q, log_names, log_values)

    def _rebin_result(self):                    #apply rebinning
        rebin_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        rebin_prog.report('Rebin result ')

        logger.information('Rebin option : ' + self._rebin_option)
        qrange = ''
        if mtd.doesExist(self._sofq):                  #check if S(Q) WS exists
		    logger.information('Sofq data from Workspace : %s' % self._sofq)
        else:                                    #read from nxs file
            sofq_path = FileFinder.getFullPath(self._sofq + '.nxs')
            LoadNexusProcessed(Filename=sofq_path,
                               OutputWorkspace=self._sofq,
                               EnableLogging=False)
            logger.information('Sq data from File : %s' % sofq_path)
        rebin_logs = [('rebin_option', self._rebin_option)]
        if self._rebin_option != 'None':          #rebin to be applied
            rebin_logs.append(('rebin_qrange', self._rebin_qrange))
            logger.information('Rebin qrange : %s' % self._rebin_qrange)
            if self._rebin_qrange == 'New':          #new Q range
                mtd[self._final_q].setDistribution(True)
                xs = mtd[self._final_q].readX(0)
                new_dq = float(self._rebin_qinc)       #increment in Q
                xmax = (int(xs[len(xs) -1 ] / new_dq) + 1) * new_dq    #find number of points & Q max
                qrange = '0.0, %f, %f' % (new_dq, xmax)   #create Q range
                self._rebin(self._final_q, self._final_q, qrange)
                x = mtd[self._final_q].readX(0)
                xshift = 0.5 * (x[0] - x[1])
                self._scale_x(self._final_q, self._final_q, xshift)
                logger.information('Output S(Q) rebinned for range : %s' % qrange)
            if self._rebin_qrange == 'Snap':         #use input Q range
                gR = mtd[self._sofq].getRun()      #input S(Q) WS
                stype = gR.getLogData('input_type').value
                logger.information('Rebin option : %s' % self._rebin_option)
                if stype != 'Q':             #check input was in Q
                    raise ValueError('Input type must be Q for Snap option')
                if self._rebin_option == 'Interpolate':
                    self._rebin_ws(self._final_q, self._sofq, self._final_q)
                    logger.information('Output S(Q) interpolated to input S(Q) : %s' % self._sofq)
                if self._rebin_option == 'Spline':
                    self._spline_interp(self._sofq, self._final_q, self._final_q, '', 2)
                    logger.information('Output S(Q) spline interpolated to input S(Q) :%s ' % self._sofq)
                    rebin_logs.append(('rebin_Q_file', self._sofq))
        log_names = [item[0] for item in rebin_logs]
        log_values = [item[1] for item in rebin_logs]
#        self._add_sample_log_mult(self._final_q, log_names, log_values)
        logger.information('Corrected WS created : %s' % self._final_q)

    def _create_ws(self, output_ws, x, y, e, nspec, names):
        create_alg = self.createChildAlgorithm("CreateWorkspace", enableLogging = False)
        create_alg.setProperty("OutputWorkspace", output_ws)
        create_alg.setProperty("DataX", x)
        create_alg.setProperty("DataY", y)
        create_alg.setProperty("DataE", e)
        create_alg.setProperty("Nspec", nspec)
        create_alg.setProperty("UnitX", 'MomentumTransfer')
        create_alg.setProperty("VerticalAxisUnit", 'Text')
        create_alg.setProperty("VerticalAxisValues", names)
        create_alg.setProperty("Distribution", True)
        create_alg.execute()
        mtd.addOrReplace(output_ws, create_alg.getProperty("OutputWorkspace").value)
        unitx = mtd[output_ws].getAxis(0).setUnit("Label")
        unitx.setLabel('2theta', 'deg')

    def _extract(self, input_ws, output_ws, index):
        extract_alg = self.createChildAlgorithm("ExtractSingleSpectrum", enableLogging = False)
        extract_alg.setProperty("InputWorkspace", input_ws)
        extract_alg.setProperty("OutputWorkspace", output_ws)
        extract_alg.setProperty("WorkspaceIndex", index)
        extract_alg.execute()
        mtd.addOrReplace(output_ws, extract_alg.getProperty("OutputWorkspace").value)

    def _clone_ws(self, input_ws, output_ws):
        clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
        clone_alg.setProperty("InputWorkspace", input_ws)
        clone_alg.setProperty("OutputWorkspace", output_ws)
        clone_alg.execute()
        mtd.addOrReplace(output_ws, clone_alg.getProperty("OutputWorkspace").value)

    def _rename_ws(self, input_ws, output_ws):
        rename_alg = self.createChildAlgorithm("RenameWorkspace", enableLogging=False)
        rename_alg.setProperty("InputWorkspace", input_ws)
        rename_alg.setProperty("OutputWorkspace", output_ws)
        rename_alg.execute()
        mtd.addOrReplace(output_ws, rename_alg.getProperty("OutputWorkspace").value)

    def _delete_ws(self, input_ws):
        delete_alg = self.createChildAlgorithm("DeleteWorkspace", enableLogging=False)
        delete_alg.setProperty("Workspace", input_ws)
        delete_alg.execute()

    def _scale(self, input_ws, output_ws, factor):
        scale_alg = self.createChildAlgorithm("Scale", enableLogging=False)
        scale_alg.setProperty("InputWorkspace", input_ws)
        scale_alg.setProperty("OutputWorkspace", output_ws)
        scale_alg.setProperty("Factor", factor)
        scale_alg.setProperty("Operation", 'Multiply')
        scale_alg.execute()
        mtd.addOrReplace(output_ws, scale_alg.getProperty("OutputWorkspace").value)

    def _scale_x(self, input_ws, output_ws, factor):
        scale_x_alg = self.createChildAlgorithm("ScaleX", enableLogging=False)
        scale_x_alg.setProperty("InputWorkspace", input_ws)
        scale_x_alg.setProperty("OutputWorkspace", output_ws)
        scale_x_alg.setProperty("Factor", factor)
        scale_x_alg.setProperty("Operation", 'Add')
        scale_x_alg.execute()
        mtd.addOrReplace(output_ws, scale_x_alg.getProperty("OutputWorkspace").value)

    def _minus(self, lhs_ws, rhs_ws, output_ws):
        minus_alg = self.createChildAlgorithm("Minus", enableLogging=False)
        minus_alg.setProperty("LHSWorkspace", lhs_ws)
        minus_alg.setProperty("RHSWorkspace", rhs_ws)
        minus_alg.setProperty("OutputWorkspace", output_ws)
        minus_alg.execute()
        mtd.addOrReplace(output_ws, minus_alg.getProperty("OutputWorkspace").value)

    def _plus(self, lhs_ws, rhs_ws, output_ws):
        plus_alg = self.createChildAlgorithm("Plus", enableLogging=False)
        plus_alg.setProperty("LHSWorkspace", lhs_ws)
        plus_alg.setProperty("RHSWorkspace", rhs_ws)
        plus_alg.setProperty("OutputWorkspace", output_ws)
        plus_alg.execute()
        mtd.addOrReplace(output_ws, plus_alg.getProperty("OutputWorkspace").value)

    def _multiply(self, lhs_ws, rhs_ws, output_ws):
        multiply_alg = self.createChildAlgorithm("Multiply", enableLogging=False)
        multiply_alg.setProperty("LHSWorkspace", lhs_ws)
        multiply_alg.setProperty("RHSWorkspace", rhs_ws)
        multiply_alg.setProperty("OutputWorkspace", output_ws)
        multiply_alg.execute()
        mtd.addOrReplace(output_ws, multiply_alg.getProperty("OutputWorkspace").value)

    def _divide(self, lhs_ws, rhs_ws, output_ws):
        divide_alg = self.createChildAlgorithm("Divide", enableLogging=False)
        divide_alg.setProperty("LHSWorkspace", lhs_ws)
        divide_alg.setProperty("RHSWorkspace", rhs_ws)
        divide_alg.setProperty("OutputWorkspace", output_ws)
        divide_alg.execute()
        mtd.addOrReplace(output_ws, divide_alg.getProperty("OutputWorkspace").value)

    def _copy_log(self, input_ws, output_ws, strategy):
        copy_log_alg = self.createChildAlgorithm("CopyLogs", enableLogging=False)
        copy_log_alg.setProperty("InputWorkspace", input_ws)
        copy_log_alg.setProperty("OutputWorkspace", output_ws)
        copy_log_alg.setProperty("MergeStrategy", strategy)
        copy_log_alg.execute()

    def _add_sample_log_mult(self, input_ws, log_names, log_values):
        sample_log_mult_alg = self.createChildAlgorithm("AddSampleLogMultiple", enableLogging=False)
        sample_log_mult_alg.setProperty("Workspace", input_ws)
        sample_log_mult_alg.setProperty("LogNames", log_names)
        sample_log_mult_alg.setProperty("LogValues", log_values)
        sample_log_mult_alg.execute()

    def _spline_interp(self, match_ws, interp_ws, output_ws, deriv_ws, order):
        spline_interp_alg = self.createChildAlgorithm("SplineInterpolation", enableLogging = False)
        spline_interp_alg.setProperty("WorkspaceToMatch", match_ws)
        spline_interp_alg.setProperty("WorkspaceToInterpolate", interp_ws)
        spline_interp_alg.setProperty("OutputWorkspace", output_ws)
        spline_interp_alg.setProperty("OutputWorkspaceDeriv", deriv_ws)
        spline_interp_alg.setProperty("DerivOrder", order)
        spline_interp_alg.execute()
        mtd.addOrReplace(output_ws, spline_interp_alg.getProperty("OutputWorkspace").value)

    def _rebin(self, input_ws, output_ws, params):
        rebin_alg = self.createChildAlgorithm("Rebin", enableLogging=False)
        rebin_alg.setProperty("InputWorkspace", input_ws)
        rebin_alg.setProperty("OutputWorkspace", output_ws)
        rebin_alg.setProperty("Params", params)
        rebin_alg.execute()
        mtd.addOrReplace(output_ws, rebin_alg.getProperty("OutputWorkspace").value)

    def _rebin_ws(self, rebin_ws, match_ws, output_ws):
        rebin_ws_alg = self.createChildAlgorithm("RebinToWorkspace", enableLogging=False)
        rebin_ws_alg.setProperty("WorkspaceToRebin", rebin_ws)
        rebin_ws_alg.setProperty("WorkspaceToMatch", match_ws)
        rebin_ws_alg.setProperty("OutputWorkspace", output_ws)
        rebin_ws_alg.execute()
        mtd.addOrReplace(output_ws, rebin_ws_alg.getProperty("OutputWorkspace").value)

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconD4Result)
#
