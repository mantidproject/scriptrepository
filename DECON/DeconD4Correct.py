from mantid.simpleapi import (mtd, LoadNexus, CloneWorkspace, SaveNexusProcessed, CopyLogs,
                              AddSampleLogMultiple, DeconApplyCorrections, DeconD4Result)
from mantid.api import *
from mantid.kernel import *
from mantid import config
import math, os.path, numpy as np


class DeconD4Correct(DataProcessorAlgorithm):

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
    _mome = None
    _coeff = None
    _result = None
    _final_q = None
    _nterms = 2
    _smooth = False
    _cutoff = False
    _cutoff_pt = 0.0
    _rebin_option = 'None'
    _rebin_qrange = None
    _rebin_qinc = None
    _plot = False
    _saveNXS = False
    _saveAscii = False

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Perform deconvolution."

    def PyInit(self):
        self.declareProperty(name = 'S(theta) Input', defaultValue = 'Workspace',
                             validator = StringListValidator(['Workspace','File']),
                             doc = 'Sample input type')
        thetaWsCondition = VisibleWhenProperty('S(theta) Input', PropertyCriterion.IsEqualTo, 'Workspace')
        thetaFileCondition = VisibleWhenProperty('S(theta) Input', PropertyCriterion.IsEqualTo, 'File')

        self.declareProperty(MatrixWorkspaceProperty('S(theta) Workspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input S(theta) workspace.")
        self.setPropertySettings('S(theta) Workspace', thetaWsCondition)

        self.declareProperty(FileProperty('S(theta) File', '', action = FileAction.OptionalLoad,
	                                      extensions = ["_theta.nxs"]),
                             doc = 'File path for S(theta) file')
        self.setPropertySettings('S(theta) File', thetaFileCondition)

        self.declareProperty(name = 'Derivatives Input', defaultValue = 'Workspace',
                             validator = StringListValidator(['Workspace','File']),
                             doc = 'Derivatives input type')
        derivWsCondition = VisibleWhenProperty('Derivatives Input', PropertyCriterion.IsEqualTo, 'Workspace')
        derivFileCondition = VisibleWhenProperty('Derivatives Input', PropertyCriterion.IsEqualTo, 'File')

        self.declareProperty(MatrixWorkspaceProperty('Derivatives Workspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Derivatives workspace.")
        self.setPropertySettings('Derivatives Workspace', derivWsCondition)

        self.declareProperty(FileProperty('Derivatives File', '', action = FileAction.OptionalLoad,
                                          extensions = ["_deriv.nxs"]),
                             doc = 'File path for Derivatives file')
        self.setPropertySettings('Derivatives Workspace', derivWsCondition)

        self.declareProperty(name = 'Moments Input', defaultValue = 'Workspace',
                             validator = StringListValidator(['Workspace','File']),
                             doc = 'Sample input type')
        momeWsCondition = VisibleWhenProperty('Moments Input', PropertyCriterion.IsEqualTo, 'Workspace')
        momeFileCondition = VisibleWhenProperty('Moments Input', PropertyCriterion.IsEqualTo, 'File')

        self.declareProperty(MatrixWorkspaceProperty('Moments Workspace', '',
                                                     optional = PropertyMode.Optional, direction = Direction.Input),
                             doc = "Name for the input Moments workspace.")
        self.setPropertySettings('Moments Workspace', momeWsCondition)

        self.declareProperty(FileProperty('Moments File', '', action = FileAction.OptionalLoad,
                                          extensions = ["_mome.nxs"]),
                             doc = 'File path for Moments file')
        self.setPropertySettings('Moments File', momeFileCondition)

        self.declareProperty(MatrixWorkspaceProperty('QWorkspace', 'Q',
                                                     direction = Direction.Output),
                             doc = "Name for the output Derivatives workspace.")

        self.declareProperty(name = 'Use Smooth S(Q)', defaultValue = False,
                             doc = 'Switch smoothing Off/On')
        self.declareProperty(name = 'Number of terms to use', defaultValue = 2,
                             validator=IntBoundedValidator(1),
                             doc = 'Number of corrections terms to use')

        self.declareProperty(name = 'Cutoff', defaultValue = False,
                             doc = 'Low cutoff option')
        cutoffCondition = VisibleWhenProperty('Cutoff',PropertyCriterion.IsNotDefault)
        self.declareProperty(name = 'Cutoff point', defaultValue = 180.0,
                             doc = 'Low cutoff value')
        self.setPropertySettings('Cutoff point', cutoffCondition)

        self.declareProperty(name = 'Wavelength for output', defaultValue = 0.7,
                             doc = 'Wavelength for output')
        self.declareProperty(name = 'Zero angle for output', defaultValue = 0.0,
                             doc = 'Zero angle for Q convertion')

        self.declareProperty(name = 'Rebin option', defaultValue = 'None',
                             validator = StringListValidator(['None','Interpolate','Spline']),
                             doc = 'Rebin options : None or Interpolate')
        rebinOptionCondition = VisibleWhenProperty('Rebin option',PropertyCriterion.IsNotEqualTo, 'None')

        self.declareProperty(name = 'Rebin Qrange', defaultValue = 'New',
                             validator = StringListValidator(['New','Snap']),
                             doc = 'Rebin Qrange : input New values or Snap to input values')
        rebinQrangeCondition = VisibleWhenProperty('Rebin Qrange',PropertyCriterion.IsEqualTo, 'New')
        self.setPropertySettings('Rebin Qrange', rebinOptionCondition)

        self.declareProperty(name = 'Rebin Qinc', defaultValue = 0.1,
                             doc = 'Value of deltaQ for rebinning')
        self.setPropertySettings('Rebin Qinc', rebinQrangeCondition)

        self.declareProperty(name = 'Plot', defaultValue = False,
                             doc = 'Plot options')
        self.declareProperty(name = 'Save NXS format', defaultValue = False,
                             doc = 'Switch Save result to nxs file Off/On')
        self.declareProperty(name = 'Save ASCII format', defaultValue = False,
                             doc = 'Switch Save result to ascii file Off/On')

    def PyExec(self):
        import mantidplot as mp
        workdir = config['defaultsave.directory']
        self._setup()
        DeconApplyCorrections(DataWorkspace=self._stheta,
                              DerivativesWorkspace=self._deriv,
                              MomentsWorkspace=self._mome,
                              UseSmoothData=self._smooth,
                              NumberTerms=self._nterms,
                              Cutoff=self._cutoff,
                              CutoffPoint=self._cutoff_pt)
        DeconD4Result(CorrectedWorkspace=self._stheta + '_corrected',
                      QWorkspace=self._sofq,
                      RebinOption='None',
					  RebinQrange=self._rebin_qrange,
                      RebinQinc=self._rebin_qinc)

        if self._plot:
            mp.plotSpectrum(self._stheta + '_coeff', [0,1,2,3])
            mp.plotSpectrum(self._stheta + '_corr', [0,1,2,3])
            result_graph=mp.plotSpectrum(self._stheta + '_result', [0,1,2,3], False)
            mp.mergePlots(result_graph,mp.plotSpectrum(self._stheta + '_used', 0, False))
            self._plot_result([self._stheta + '_used', self._stheta + '_corrected'], 0)
            finalQ_graph=mp.plotSpectrum(self._sofq + '_corrected', 0, False)
            finalQ_layer = finalQ_graph.activeLayer()
            finalQ_layer.setAxisTitle(mp.Layer.Bottom, 'Q')
            mp.mergePlots(finalQ_graph,mp.plotSpectrum(self._sofq, 0, False))

        if self._saveNXS:
            save_nxs_prog = Progress(self, start=0.0, end=0.8, nreports=4)
            save_nxs_prog.report('Save NXS ')

            self._save_ws(self._stheta + '_coeff', 'Coefficients')
            self._save_ws(self._stheta + '_corr', 'Corrections')
            self._save_ws(self._stheta + '_result', 'Result')
            self._save_ws(self._stheta + '_corrected', 'Final theta corrected')
            self._save_ws(self._sofq + '_corrected', 'Final Q corrected')
            save_nxs_prog.report('Save NXS completed')

        if self._saveAscii:
            save_ascii_prog = Progress(self, start=0.0, end=0.8, nreports=3)
            save_ascii_prog.report('Save ascii ')

            self._save_Ascii(self._stheta + '_corrected', '.stc')
            self._save_Ascii(self._sofq + '_corrected', '.sqc')
            save_ascii_prog.report('Save ascii completed')

    def _setup(self):
        self._input = self.getPropertyValue('S(theta) Input')
        if self._input == 'Workspace':
            self._input_ws = self.getPropertyValue('S(theta) Workspace')
        else:
            self._input_ws = ''
        if self._input == 'File':
            self._path = self.getPropertyValue('S(theta) File')
        else:
            self._path = ''
        self._stheta = self._get_data()
        self._name = self._stheta[:-6]    #remove _theta

        #convert input S(theta) to Q
        input_prog = Progress(self, start=0.0, end=0.8, nreports=3)
        input_prog.report('Input data ')

        self._lambda = self.getProperty('Wavelength for output').value
        self._azero = self.getProperty('Zero angle for output').value
        self._sofq = self.getPropertyValue('QWorkspace')
        self._theta_to_q()

        input_prog.report('Input data completed')

        self._input = self.getPropertyValue('Derivatives Input')
        if self._input == 'Workspace':
            self._input_ws = self.getPropertyValue('Derivatives Workspace')
        else:
            self._input_ws = ''
        if self._input == 'File':
            self._path = self.getPropertyValue('Derivatives File')
        else:
            self._path = ''
        self._deriv = self._get_data()

        self._input = self.getPropertyValue('Moments Input')
        if self._input == 'Workspace':
            self._input_ws = self.getPropertyValue('Moments Workspace')
        else:
            self._input_ws = ''
        if self._input == 'File':
            self._path = self.getPropertyValue('Moments File')
        else:
            self._path = ''
        self._mome = self._get_data()

        self._nterms = int(self.getProperty('Number of terms to use').value)
        self._smooth = self.getProperty('Use Smooth S(Q)').value
        self._cutoff = self.getProperty('Cutoff').value
        self._cutoff_pt = float(self.getProperty('Cutoff point').value)
        self._rebin_option = self.getPropertyValue('Rebin option')
        self._rebin_qrange = self.getProperty('Rebin Qrange').value
        self._rebin_qinc = float(self.getProperty('Rebin Qinc').value)
        self._plot = self.getProperty('Plot').value
        self._saveNXS = self.getProperty('Save NXS format').value
        self._saveAscii = self.getProperty('Save ASCII format').value

    def _save_Ascii(self, inWS, ext):    #save result in ascii file
        workdir = config['defaultsave.directory']
        file = inWS + ext
        path = os.path.join(workdir, file)    #output file path
        logger.information('Creating Ascii file : %s' % path)
        handle = open(path, 'w')
        handle.write('# ' + file + " \n")
        handle.write("# \n")
        inGR = mtd[inWS].getRun()             #input WS
        len_h_in = inGR.getLogData('sofq_lines').value
        for n in range(0, int(len_h_in)):   #get header data
            h = inGR.getLogData('sofq_%i' % (n)).value
            handle.write('#' + h + " \n")
        handle.write("# \n")
        handle.write('# Block  2 - Decon' + " \n")
        if self._smooth:
            text = '# S(Q) smoothed'
        else:
            text = '# S(Q) unsmoothed'
        handle.write(text +" \n")
        handle.write('# moments file = ' + self._mome + '.mom' + " \n")
        if self._cutoff:
            text = '# Cutoff angle = %f' % (self._cutoff_pt)
        else:
            text = '# NO Cutoff'
        handle.write(text + " \n")
        handle.write(('# %i term corrected S(theta)' % (self._nterms)) +" \n")
        if ext == '.sqc':
            if self._rebin_option == 'None':
                handle.write('# NO rebinning' + " \n")
            else:
                handle.write('# rebin method : ' + self._rebin_option + " \n")
                if self._rebin_qrange == 'New':
                    handle.write(('# rebin Q-range is New with Qinc : %f' % self._rebin_qinc) + " \n")
                if self._rebin_qrange == 'Snap':
                    handle.write('# rebin Q-range is Snap to file : ' + self._sofq + " \n")
        handle.write('#   theta   S(theta)     sigma' + " \n")
        x = mtd[inWS].readX(0)            #data arrays
        y = mtd[inWS].readY(0)
        e = mtd[inWS].readE(0)
        for n in range(0,len(y)):
            handle.write(('%f   %f   %f' % (x[n], y[n], e[n])) + "\n")
        handle.close()

    def _plot_result(self, ws, list):                       #Plot
        import mantidplot as mp
        mp.plotSpectrum(ws, list)

    def _get_data(self):   #get data
        if self._input == 'Workspace':
            input_ws = self._input_ws
            name = self._input_ws
            logger.notice('Input from Workspace : %s' % input_ws)
        elif self._input == 'File':
            name = self._get_filename()
            input_ws = name
            LoadNexus(Filename = self._path,
                      OutputWorkspace = input_ws,
                      EnableLogging=False)
            unitx = mtd[input_ws].getAxis(0).setUnit("Label")
            unitx.setLabel('2theta', 'deg')
        else:
            raise ValueError('Input type not defined')
        self._add_sample_log_mult(input_ws, ['input_type'], ['theta'])
        return name

    def _get_filename(self):
        path = self._path
        if(os.path.isfile(path)):
            base = os.path.basename(path)
            name = os.path.splitext(base)[0]
            ext = os.path.splitext(base)[1]
            logger.information('Input file : %s' % path)
        else:
            raise ValueError('Could not find file: %s' % path)
        return name

    def _theta_to_q(self):
        k0 = 4.0 * math.pi / self._lambda
        self._clone_ws(self._stheta, self._sofq)
        x_q = mtd[self._sofq].dataX(0)
        x_q = k0 * np.sin(0.5 * np.radians(x_q - self._azero))    #convert to Q after applying zero angle correction
        mtd[self._sofq].setX(0, x_q)
        self._copy_log(self._stheta, self._sofq, 'MergeReplaceExisting')
        convert_logs = [('lambda_out', self._lambda), ('zero_out', self._azero), ('input_type', 'Q')]
        logger.information('Converting %s ; from theta to Q as  : %s' % (self._stheta, self._sofq))
        logger.information('lambda = %f ; zero = %f' % (self._lambda, self._azero))
        log_names = [item[0] for item in convert_logs]
        log_values = [item[1] for item in convert_logs]
        self._add_sample_log_mult(self._sofq, log_names, log_values)
        self.setProperty("QWorkspace", self._sofq)

    def _clone_ws(self, input_ws, output_ws):
        clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
        clone_alg.setProperty("InputWorkspace", input_ws)
        clone_alg.setProperty("OutputWorkspace", output_ws)
        clone_alg.execute()
        mtd.addOrReplace(output_ws, clone_alg.getProperty("OutputWorkspace").value)

    def _save_ws(self, input_ws, name):
        workdir = config['defaultsave.directory']
        path = os.path.join(workdir, input_ws + '.nxs')
        save_alg = self.createChildAlgorithm("SaveNexusProcessed", enableLogging = False)
        save_alg.setProperty("InputWorkspace", input_ws)
        save_alg.setProperty("Filename", path)
        save_alg.execute()
        logger.information('%s file : %s' % (name, path))

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

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconD4Correct)
#
