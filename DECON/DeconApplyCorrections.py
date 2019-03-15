from mantid.simpleapi import (mtd, CreateWorkspace, CloneWorkspace, DeleteWorkspace, ExtractSingleSpectrum,
                              Minus, Plus, Multiply, Divide, AddSampleLogMultiple, AppendSpectra)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode, MatrixWorkspaceProperty,
                        TextAxis, Progress)
from mantid.kernel import (StringListValidator, StringMandatoryValidator, IntBoundedValidator,
                           FloatBoundedValidator, Direction, logger)
from mantid import config
import math, numpy as np


class DeconApplyCorrections(DataProcessorAlgorithm):

    _data = None
    _data_used = None
    _final_theta = None
    _deriv = None
    _mome = None
    _coeff = None
    _result = None
    _nterms = 2
    _smooth = False
    _cutoff = False
    _cutoff_pt = 0.0

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Perform deconvolution."

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('DataWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Data workspace.")

        self.declareProperty(MatrixWorkspaceProperty('DerivativesWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Derivatives workspace.")

        self.declareProperty(MatrixWorkspaceProperty('MomentsWorkspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Moments workspace.")

        self.declareProperty(name = 'UseSmoothData', defaultValue = False,
                             doc = 'Switch smoothing Off/On')
        self.declareProperty(name = 'NumberTerms', defaultValue = 2,
                             validator=IntBoundedValidator(1),
                             doc = 'Number of corrections terms to use')
        self.declareProperty(name = 'Cutoff', defaultValue = False,
                             doc = 'Low cutoff option')
        self.declareProperty(name = 'CutoffPoint', defaultValue = 180.0,
                             doc = 'Low cutoff value')

    def PyExec(self):
        self._setup()

        # calc coefficients of M  into *_coeff
        self._calc_coeff()

        # calc corr terms  into *_corr
        self._corr_terms()

        # subtract corrections from S(Q) into *_result
        self._subtract_corr()

        # cutoff hi-theta & use chosen terms  into *_theta_corrected
        self._cut_result()

    def _setup(self):
        self._data = self.getPropertyValue('DataWorkspace')
        self._deriv = self.getPropertyValue('DerivativesWorkspace')
        self._mome = self.getPropertyValue('MomentsWorkspace')

        self._nterms = self.getProperty('NumberTerms').value
        self._smooth = self.getProperty('UseSmoothData').value
        self._cutoff = self.getProperty('Cutoff').value
        self._cutoff_pt = self.getProperty('CutoffPoint').value

    def _calc_coeff(self):         #creates output _coeff
        calc_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        calc_prog.report('Calculating coefficents ')

        smooth_ws = '__smooth'
        mome_ws = '__mome'
# start interpolating moments to S(theta) & normalising to m0
        logger.information('Calculating Coefficients')
        self._extract(self._data, smooth_ws, 1)
        self._spline_interp(smooth_ws, self._mome, mome_ws, '__deriv', 2)
        m0 = '__M0'
        self._extract(mome_ws, m0, 0)
        m1 = '__M1'
        self._extract(mome_ws, m1, 1)
        self._divide(m1, m0, m1)
        m2 = '__M2'
        self._extract(mome_ws, m2, 2)
        self._divide(m2, m0, m2)
        m3 = '__M3'
        self._extract(mome_ws, m3, 3)
        self._divide(m3, m0, m3)
        m4 = '__M4'
        self._extract(mome_ws, m4, 4)
        self._divide(m4, m0, m4)

### Henry 17 Sept 14:  From ADD2011 poster, here are the correct
###    signs for the coeffs, using the "positive" sign convention:
# c1 = -A1
# c2 = -A2 +A1*A1
# c3 = -A3 +2*A1*A2 -A1*A1*A1
# c4 = -A4 +2A1*A3 +A2*A2 -3A1*A1*A2  +A1*A1*A1*A1
### Note that the last two terms in c4 should have DIFFERENT signs,
### rather than the same signs as Spencer uses below:

# S(2q) = J +c1J1 +c2J2 +c3J3 +c4J4
# An = Mn/M0 normalised moments

### Henry 17 Sept 14:  Correct signs: c1 = -A1 (-m1)
        self._scale(m1, '__c1', -1.0)

### Henry 17 Sept 14:  Correct signs: c2 = A2 - A11 (-m2 +m1*m1)
        self._scale(m2, '__c2', -1.0)
        # c2 = -1.0 * A2
        self._multiply(m1, m1, '__a11')
        self._plus('__c2', '__a11', '__c2')
        # c2 += A1 * A1

### Henry 17 Sept 14:  Correct signs: c3 = A3 -2*A1*A2 +A11*A1
        self._scale(m3, '__c3', -1.0)
        # c3 = -1.0 * A3
        self._scale(m1, '__a12', 2.0)
        self._multiply('__a12', m2, '__a12')
        self._plus('__c3', '__a12', '__c3')
        # c3 += 2 * A1 * A2
        self._multiply('__a11', m1, '__a111')
        self._minus('__c3', '__a111', '__c3')
        # c3 += - A1 * A1 * A1

### Henry 17 Sept 14:  Correct signs: c4 = A4 -2*A1*A3 -A2*A2 +3*A11*A2 +A11*A11
        self._scale(m4, '__c4', -1.0)
        # c4 = -1.0 * A4
        self._scale(m1, '__a13', 2.0)
        self._multiply('__a13', m3, '__a13')
        self._plus('__c4', '__a13', '__c4')
        # c4 += 2 * A1 * A3
        self._multiply(m2, m2, '__a22')
        self._plus('__c4', '__a22', '__c4')
        # c4 += A2 * A2
        self._scale('__a11', '__a112', -3.0)
        self._multiply('__a112', m2, '__a112')
        self._plus('__c4', '__a112', '__c4')
        # c4 += - 3 * A1 * A1 * A2
        self._multiply('__a11', '__a11', '__a1111')
        self._plus('__c4', '__a112', '__c4')
        # c4 += A1 * A1 * A1 * A1

        self._coeff = self._data + '_coeff'             #coeffients WS
        self._clone_ws('__c1', self._coeff)
        self._append(self._coeff, '__c2', self._coeff)
        self._append(self._coeff, '__c3', self._coeff)
        self._append(self._coeff, '__c4', self._coeff)

        ax = TextAxis.create(4)
        for i, x in enumerate(['Coeff.1', 'Coeff.2', 'Coeff.3', 'Coeff.4']):
            ax.setLabel(i, x)
        mtd[self._coeff].replaceAxis(1, ax)

        self._copy_log(self._mome, self._coeff, 'MergeKeepExisting')
        workspaces = [m0, m1, m2, m3, m4, smooth_ws, mome_ws]
        for ws in workspaces:
            self._delete_ws(ws)
        workspaces = ['__a11', '__a111', '__a1111', '__a112', '__a12', '__a13', '__a22', \
                      '__c1', '__c2', '__c3', '__c4']
        for ws in workspaces:
            self._delete_ws(ws)
        logger.information('Coefficient WS created : %s' % self._coeff)
        calc_prog.report('Calculating coefficents completed')

    def _corr_terms(self):   #calculates the correction terms = coef*deriv as _corr
        calc_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        calc_prog.report('Correction terms ')

        logger.information('Calculating Correction terms')
        self._corr = self._data + '_corr'              #corrections WS
        self._extract(self._deriv, '__temp', 0)
        self._spline_interp('__temp', self._coeff, self._coeff, '', 2)
        self._multiply(self._coeff, self._deriv, self._corr)

        ax = TextAxis.create(4)
        for i, x in enumerate(['Corr.1', 'Corr.2', 'Corr.3', 'Corr.4']):
            ax.setLabel(i, x)
        mtd[self._corr].replaceAxis(1, ax)

        self._copy_log(self._mome, self._corr, 'MergeKeepExisting')
        self._delete_ws('__temp')
        logger.information('Correction terms WS created : %s' % self._corr)
        calc_prog.report('Correction terms completed')

    def _subtract_corr(self):    #subtract corrections from input to give _data_used & _result
        calc_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        calc_prog.report('Subtract corrections ')

        logger.information('Subtracting corrections')
        self._data_used = self._data + '_used'
        if self._smooth:                                #select which hist to use
            index = 1
        else:
            index= 0
        self._extract(self._data, self._data_used, index)

        self._extract(self._corr, '__wcr', 0)
        wsc1 = 'S-1C'    #1 term subtracted
        self._plus(self._data_used, '__wcr', wsc1)
        wsc2 = 'S-2C'    #2 terms subtracted
        self._extract(self._corr, '__wcr', 1)
        self._plus(wsc1, '__wcr', wsc2)
        wsc3 = 'S-3C'    #3 terms subtracted
        self._extract(self._corr, '__wcr', 2)
        self._plus(wsc2, '__wcr', wsc3)
        wsc4 = 'S-4C'    #4 terms subtracted
        self._extract(self._corr, '__wcr', 3)
        self._plus(wsc3, '__wcr', wsc4)

        self._result = self._data + '_result'             #results WS
        self._clone_ws(wsc1, self._result)
        self._append(self._result, wsc2, self._result)
        self._append(self._result, wsc3, self._result)
        self._append(self._result, wsc4, self._result)

        ax = TextAxis.create(4)
        for i, x in enumerate(['S-1C', 'S-2C', 'S-3C', 'S-4C']):
            ax.setLabel(i, x)
        mtd[self._result].replaceAxis(1, ax)

        subtract_logs = [('smooth', self._smooth)]
        log_names = [item[0] for item in subtract_logs]
        log_values = [item[1] for item in subtract_logs]
        self._add_sample_log_mult(self._result, log_names, log_values)
        workspaces = ['__wcr', wsc1, wsc2, wsc3, wsc4]
        for ws in workspaces:
            self._delete_ws(ws)
        logger.information('Results in WS %s' % self._result)

    def _cut_result(self):                  #cutoff high angle data ouput _final_theta as *_theta_corrected
        calc_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        calc_prog.report('Cut result ')

        logger.information('Cutting result')
        logger.information('Number of terms used : %i' % (self._nterms))
        temp_ws = '__final'
        self._extract(self._result, temp_ws, self._nterms -1)
        self._final_theta = self._data + '_corrected'   #theta corrected WS
        final_list = ['S-1C', 'S-2C', 'S-3C', 'S-4C']   #names for hist
        icut = 0
        cut_pt = 0
        cut_logs = [('correct_terms', self._nterms), ('cutoff', self._cutoff)]
        if self._cutoff:
            xs = mtd[self._data_used].readX(0)             #S(theta) data
            ys = mtd[self._data_used].readY(0)
            es = mtd[self._data_used].readE(0)
            xf = np.array(mtd[temp_ws].readX(0))
            xfa = np.array(xf)
            yf = np.array(mtd[temp_ws].readY(0))
            ef = np.array(mtd[temp_ws].readE(0))
            icut = np.where(self._cutoff_pt > xfa)[0][-1]
            cut_pt = xf[icut]                      #x-value for cutoff
            logger.information('Corrected data cutoff at : %f' % (cut_pt))
            xnew = np.array(xf[:icut])                #start off new array
            xnew = np.append(xnew,np.array(xs[icut:]))  #append input data
            ynew = np.array(yf[:icut])                #start off new array
            ynew = np.append(ynew,np.array(ys[icut:]))  #append input data
            enew = np.array(ef[:icut])
            enew = np.append(enew,np.array(es[icut:]))
            self._create_ws(self._final_theta, xnew, ynew, enew, 1, final_list[self._nterms - 1])
            cut_logs.append(('cutoff_point', cut_pt))
            self._delete_ws(temp_ws)
        else:
            self._rename_ws(temp_ws, self._final_theta)
            logger.information('Corrected data NOT cutoff')
        self._copy_log(self._result, self._final_theta, 'MergeReplaceExisting')
        log_names = [item[0] for item in cut_logs]
        log_values = [item[1] for item in cut_logs]
        self._add_sample_log_mult(self._final_theta, log_names, log_values)
        logger.information('Corrected WS created : ' + self._final_theta)

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

    def _append(self, input1_ws, input2_ws, output_ws):
        append_alg = self.createChildAlgorithm("AppendSpectra", enableLogging = False)
        append_alg.setProperty("InputWorkspace1", input1_ws)
        append_alg.setProperty("InputWorkspace2", input2_ws)
        append_alg.setProperty("OutputWorkspace", output_ws)
        append_alg.execute()
        mtd.addOrReplace(output_ws, append_alg.getProperty("OutputWorkspace").value)

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconApplyCorrections)
#
