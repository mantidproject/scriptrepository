#pylint: disable=no-init,too-many-instance-attributes
from mantid.simpleapi import (mtd, ConvertToHistogram, CloneWorkspace, CropWorkspace, DeleteWorkspace,
                              RenameWorkspace, ExtractSingleSpectrum, GroupWorkspaces, CopyLogs,
                              AddSampleLogMultiple, ConvertSpectrumAxis, Fit, PlotPeakByLogValue,
                              ProcessIndirectFitParameters, ConvertTableToMatrixWorkspace)
from mantid.api import (PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, PropertyMode, Progress,
                        SpectraAxis, TextAxis)
from mantid.kernel import (Direction, StringListValidator, logger)
from mantid import config
import math, numpy as np
import os


class IndirectTwoPeakFit(PythonAlgorithm):

    _sample_ws = None
    _res_ws = None
    _e_min = None
    _e_max = None
    _hist_min = None
    _hist_max = None
    _bgd = None
    _elastic = None
    _parameter_table = None
    _output_workspace = None


    def category(self):
        return "Workflow\\Inelastic;PythonAlgorithms;Workflow\\MIDAS"


    def summary(self):
        return 'Convolution fit for 1 and 2 Lorentzians'


    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('SampleWorkspace', '',
                                                     optional=PropertyMode.Mandatory,
                                                     direction=Direction.Input),
                             doc="Name for the sample workspace.")

        self.declareProperty(name='EnergyMin', defaultValue=-0.5,
                             doc='Minimum energy for fit. Default=-0.5')

        self.declareProperty(name='EnergyMax', defaultValue=0.5,
                             doc='Maximum energy for fit. Default=0.5')

        self.declareProperty(name='Minimizer',defaultValue='Levenberg-Marquardt',
                             validator=StringListValidator(['Levenberg-Marquardt','FABADA']),
                             doc='Type of minimizer')

        self.declareProperty(name='MaxIterations', defaultValue=500,
                             doc='Max iterations. Default=500')

        self.declareProperty(name='OutputName', defaultValue='',
                             doc='Output workspace base name')


    def PyExec(self):

        setup_prog = Progress(self, start=0.05, end=0.95, nreports=3)

        self._tmp_fit_name = "__fit_ws"
        self._crop_ws(self._sample_ws, self._tmp_fit_name, self._e_min, self._e_max)

        convert_to_hist_alg = self.createChildAlgorithm("ConvertToHistogram", enableLogging=False)
        convert_to_hist_alg.setProperty("InputWorkspace", self._tmp_fit_name)
        convert_to_hist_alg.setProperty("OutputWorkspace", self._tmp_fit_name)
        convert_to_hist_alg.execute()
        mtd.addOrReplace(self._tmp_fit_name, convert_to_hist_alg.getProperty("OutputWorkspace").value)

        self._convert_to_elasticQ(self._tmp_fit_name)

        num_hist = self._sample_ws.getNumberHistograms()
        if self._hist_max is None:
            self._hist_max = num_hist - 1

        setup_prog.report('Fitting 1 peak')
        self._fit(1)
        setup_prog.report('Fitting 2 peaks')
        self._fit(2)
        self._delete_ws(self._tmp_fit_name)
		
        chi_group = self._output_name + '_ChiSq'
        chi_ws1 = self._output_name + '_1L_ChiSq'
        chi_ws2 = self._output_name + '_2L_ChiSq'
        self._clone_ws(chi_ws1, chi_group)
        self._append(chi_group, chi_ws2, chi_group)
        ws = mtd[chi_group]
        ax = TextAxis.create(2)
        for i, x in enumerate(['1 peak', '2 peaks']):
            ax.setLabel(i, x)
        ws.replaceAxis(1, ax)
        self._delete_ws(chi_ws1)
        self._delete_ws(chi_ws2)

        res_group = self._output_name + '_Result'
        res_ws1 = self._output_name + '_1L_Result'
        res_ws2 = self._output_name + '_2L_Result'
        self._extract(res_ws1, res_group, 1)
        self._extract(res_ws2, '__spectrum', 1)
        self._append(res_group, '__spectrum', res_group)
        self._extract(res_ws2, '__spectrum', 3)
        self._append(res_group, '__spectrum', res_group)
        ws = mtd[res_group]
        ax = TextAxis.create(3)
        for i, x in enumerate(['fwhm.1', 'fwhm.2.1', 'fwhm.2.2']):
            ax.setLabel(i, x)
        ws.replaceAxis(1, ax)
        self._delete_ws(res_ws1)
        self._delete_ws(res_ws2)
        self._delete_ws(self._output_name + '_1L_Parameters')
        self._delete_ws(self._output_name + '_2L_Parameters')

    def validateInputs(self):
        """
        Validate user input.
        """

        self._setup()
        issues = dict()
		
		# Energy max must be > min
        if self._e_max <= self._e_min:
            issues['EnergyMax'] = 'Energy maximum must be greater than minimum'

		# Spectrum max must be >= min
        if self._hist_max < self._hist_min:
            issues['HistogramMax'] = 'Histogram maximum cannot be less than minimum'

		# max iterations must be > 0
        if self._max_iterations == 0 or self._max_iterations == '':
            issues['MaxIterations'] = 'Maximum iterations must be greater than 0'
        return issues

    def _setup(self):
        """
        Gets algorithm properties.
        """

        self._sample_ws = self.getProperty('SampleWorkspace').value
        logger.information('Sample workspace = %s' % self._sample_ws)
 
        self._e_min = self.getProperty('EnergyMin').value
        self._e_max = self.getProperty('EnergyMax').value

        self._minimizer = self.getProperty('Minimizer').value
        self._max_iterations = self.getProperty('MaxIterations').value

        self._output_name = self.getProperty('OutputName').value
        if self._output_name == '':
            self._output_name = 'Two_peak'

    def _fit(self, numb):

        function = self._calc_function(numb)

    # Name stem for generated workspace
        self._output_workspace = self._output_name + '_%s' % (self._fit_type)
        logger.information('Fit workspace = '+ self._output_workspace)

    # Build input string for PlotPeakByLogValue
        num_hist = mtd[self._tmp_fit_name].getNumberHistograms()
    # _red file works with range(num_hist) BUT _sqw gives error when trying to do nhist+1 !!!
        input_str = [self._tmp_fit_name + ',i%s' % i for i in range(num_hist -1)]
        input_str = ';'.join(input_str)

        PlotPeakByLogValue(Input=input_str,
                           OutputWorkspace=self._output_workspace,
                           Function=function,
                           StartX=self._e_min,
                           EndX=self._e_max,
                           FitType='Sequential',
                           Minimizer=self._minimizer,
                           MaxIterations=self._max_iterations,
                           CreateOutput=True,
                           OutputCompositeMembers=True,
                           ConvolveMembers=True)

        # Remove unused workspaces
        self._delete_ws(self._output_workspace + '_NormalisedCovarianceMatrices')
        self._delete_ws(self._output_workspace + '_Parameters')

        # rename workspaces to match user input
        self._fit_group_name = self._output_workspace + '_Workspaces'
        if self._output_workspace + "_Workspaces" != self._fit_group_name:
            self._rename_ws(self._output_workspace + "_Workspaces", self._fit_group_name)
        self._parameter_name = self._output_workspace + '_Parameters'
        if self._output_workspace != self._parameter_name:
            self._rename_ws(self._output_workspace, self._parameter_name)

        # Create *_Result workspace
        if numb == 1:
            parameters = 'A0,Amplitude,FWHM'
            if self._elastic:
                parameters = 'A0,Amplitude,FWHM,Height'
        if numb == 2:
            parameters = 'A0,f0.Amplitude,f0.FWHM,f1.Amplitude,f1.FWHM'
            if self._elastic:
                parameters = 'A0,f0.Amplitude,f0.FWHM,f1.Amplitude,f1.FWHM,f2.Height'

        self._result_name = self._output_workspace + "_Result"
        pifp_alg = self.createChildAlgorithm("ProcessIndirectFitParameters", enableLogging=False)
        pifp_alg.setProperty("InputWorkspace", self._parameter_name)
        pifp_alg.setProperty("ColumnX", "axis-1")
        pifp_alg.setProperty("XAxisUnit", "MomentumTransfer")
        pifp_alg.setProperty("ParameterNames", parameters)
        pifp_alg.setProperty("OutputWorkspace", self._result_name)
        pifp_alg.execute()
        self._result_ws = pifp_alg.getProperty("OutputWorkspace").value
        mtd.addOrReplace(self._result_name, pifp_alg.getProperty("OutputWorkspace").value)
        self._transfer_sample_logs(self._result_name)

        chi_workspace = self._output_workspace + "_ChiSq"
        ctmw_alg = self.createChildAlgorithm("ConvertTableToMatrixWorkspace", enableLogging=False)
        ctmw_alg.setProperty("InputWorkspace", self._parameter_name)
        ctmw_alg.setProperty("OutputWorkspace", chi_workspace)
        ctmw_alg.setProperty("ColumnX", "axis-1")
        ctmw_alg.setProperty("ColumnY", 'Chi_squared')
        ctmw_alg.execute()
        mtd.addOrReplace(chi_workspace, ctmw_alg.getProperty("OutputWorkspace").value)
        self._transfer_sample_logs(chi_workspace)

        # Process generated workspaces
        wsnames = mtd[self._fit_group_name].getNames()
        for i, workspace in enumerate(wsnames):
            output_ws = self._output_workspace + '_%d_Workspace' % i
            self._rename_ws(workspace, output_ws)
            self._transfer_sample_logs(output_ws)

    def _transfer_sample_logs(self, ws):
        """
        Copy the sample logs from the input workspace and add them to the output workspaces
        """

        sample_logs  = {'e_min': self._e_min, 'e_max': self._e_max, 'fit_type': self._fit_type}

        self._copy_log(self._sample_ws, ws)

        log_names = [item for item in sample_logs]
        log_values = [sample_logs[item] for item in sample_logs]

        self._add_sample_log_mult(ws, log_names, log_values)

    def _calc_function(self, numb):
        import math, numpy as np
        
        if numb == 1:
            bgd = self._sample_ws.readY(0)[10]
            self._fit_type = '1L'
            self._extract(self._sample_ws, '__spectrum', 0)
            amp, FWHM = self._fit_lor('__spectrum')
            lor_fun = 'name=Lorentzian,Amplitude=%f,PeakCentre=0.0,FWHM=%f' % (amp, FWHM)
            lor_fun += ',constraint=(Amplitude>0.0,FWHM>0.0)'
            function = lor_fun
            logger.information('Function is %s' % function)
            self._fit_1 = self._fit_type
            self._delete_ws('__spectrum')
            return function

        if numb == 2:
            a0 = mtd[self._result_name].readY(0)[0]
            self._fit_type = '2L'
            amp0 = mtd[self._result_name].readY(0)[0]
            FWHM0 = mtd[self._result_name].readY(1)[0]
            self._extract(self._fit_group_name[:-10] + '0_Workspace', '__diff_ws', 2)
            y_array = mtd['__diff_ws'].readY(0)
            y_array = -y_array
            for n in range(len(y_array)):
                if y_array[n] < 0.0:
                    y_array[n] = -y_array[n]
            mtd['__diff_ws'].setY(0, np.array(y_array))
            amp1, FWHM1 = self._fit_lor('__diff_ws')
            lor_fun = 'name=Lorentzian,Amplitude=%f,PeakCentre=0.0,FWHM=%f' % (amp0, FWHM0)
            lor_fun += ',constraint=(Amplitude>0.0,FWHM>0.0)'
            lor_fun += ';name=Lorentzian,Amplitude=%f,PeakCentre=0.0,FWHM=%f' % (amp1, FWHM1)
            lor_fun += ',constraint=(Amplitude>0.0,FWHM>0.0)'
            tie_fun = ';ties=(f1.PeakCentre=f0.PeakCentre)'
            function = lor_fun + tie_fun
            logger.information('Function is %s' % function)
            self._delete_ws('__diff_ws')
            return function

    def _fit_lor(self, ws):
        fun = 'name=Lorentzian, Amplitude=1.0, PeakCentre=0.0, FWHM=%f' % self._resolution
        Fit(InputWorkspace=ws,
                Function=fun,
                Output='__peak',
                OutputParametersOnly=True,
                EnableLogging=False)
        params_table = '__peak_Parameters'
        para_y = np.asarray(mtd[params_table].column('Value'))
        height = para_y[0]
        FWHM = para_y[2]
        logger.information('Width is %f' % FWHM)
        self._delete_ws('__peak_Parameters')
        self._delete_ws('__peak_NormalisedCovarianceMatrix')
        return height, FWHM

    def _get_Efixed(self, workspace):
        inst = mtd[workspace].getInstrument()

        if inst.hasParameter('Efixed'):
            self._e_fixed = inst.getNumberParameter('EFixed')[0]

        if inst.hasParameter('analyser'):
            analyser_name = inst.getStringParameter('analyser')[0]
            analyser_comp = inst.getComponentByName(analyser_name)

            if analyser_comp.hasParameter('Efixed'):
                self._e_fixed = analyser_comp.getNumberParameter('EFixed')[0]
                logger.information('Efixed = %f' % self._e_fixed)

                if analyser_comp.hasParameter('resolution'):
                    self._resolution = analyser_comp.getNumberParameter('resolution')[0]
                else:
                    logger.information('Resolution not found')
                    self._resolution = 0.01
                logger.information('Resolution = %f' % self._resolution)
            else:
                raise ValueError('Efixed parameter not found')
        else:
            raise ValueError('Analyser parameter not found')

    def _convert_to_elasticQ(self, input_ws, output_ws=None):
        """
        Helper function to convert the spectrum axis of a sample to ElasticQ.

        @param input_ws - the name of the workspace to convert from
        @param output_ws - the name to call the converted workspace
        """

        if output_ws is None:
            output_ws = input_ws

        self._get_Efixed(input_ws)
        axis = mtd[input_ws].getAxis(1)
        if axis.isSpectra():
            convert_axis_alg = self.createChildAlgorithm("ConvertSpectrumAxis", enableLogging=False)
            convert_axis_alg.setProperty("InputWorkspace", input_ws)
            convert_axis_alg.setProperty("Target", 'ElasticQ')
            convert_axis_alg.setProperty("EMode", 'Indirect')
            convert_axis_alg.setProperty("EFixed", self._e_fixed)
            convert_axis_alg.setProperty("OutputWorkspace", output_ws)
            convert_axis_alg.execute()
            mtd.addOrReplace(output_ws, convert_axis_alg.getProperty("OutputWorkspace").value)

        elif axis.isNumeric():
        # Check that units are Momentum Transfer
            if axis.getUnit().unitID() != 'MomentumTransfer':
                raise RuntimeError('Input must have axis values of Q')

            self._clone_ws(input_ws, output_ws)

        else:
            raise RuntimeError('Input workspace must have either spectra or numeric axis.')

    def _clone_ws(self, input_ws, output_ws):
        clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
        clone_alg.setProperty("InputWorkspace", input_ws)
        clone_alg.setProperty("OutputWorkspace", output_ws)
        clone_alg.execute()
        mtd.addOrReplace(output_ws, clone_alg.getProperty("OutputWorkspace").value)

    def _crop_ws(self, input_ws, output_ws, xmin, xmax):
        crop_alg = self.createChildAlgorithm("CropWorkspace", enableLogging=False)
        crop_alg.setProperty("InputWorkspace", input_ws)
        crop_alg.setProperty("OutputWorkspace", output_ws)
        crop_alg.setProperty("XMin", xmin)
        crop_alg.setProperty("XMax", xmax)
        crop_alg.execute()
        mtd.addOrReplace(output_ws, crop_alg.getProperty("OutputWorkspace").value)

    def _delete_ws(self, input_ws):
        delete_alg = self.createChildAlgorithm("DeleteWorkspace", enableLogging=False)
        delete_alg.setProperty("Workspace", input_ws)
        delete_alg.execute()

    def _rename_ws(self, input_ws, output_ws):
        rename_alg = self.createChildAlgorithm("RenameWorkspace", enableLogging=False)
        rename_alg.setProperty("InputWorkspace", input_ws)
        rename_alg.setProperty("OutputWorkspace", output_ws)
        rename_alg.execute()
        mtd.addOrReplace(output_ws, rename_alg.getProperty("OutputWorkspace").value)

    def _extract(self, input_ws, output_ws, index):
        extract_alg = self.createChildAlgorithm("ExtractSingleSpectrum", enableLogging = False)
        extract_alg.setProperty("InputWorkspace", input_ws)
        extract_alg.setProperty("WorkspaceIndex", index)
        extract_alg.setProperty("OutputWorkspace", output_ws)
        extract_alg.execute()
        mtd.addOrReplace(output_ws, extract_alg.getProperty("OutputWorkspace").value)

    def _append(self, input1_ws, input2_ws, output_ws):
        append_alg = self.createChildAlgorithm("AppendSpectra", enableLogging = False)
        append_alg.setProperty("InputWorkspace1", input1_ws)
        append_alg.setProperty("InputWorkspace2", input2_ws)
        append_alg.setProperty("OutputWorkspace", output_ws)
        append_alg.execute()
        mtd.addOrReplace(output_ws, append_alg.getProperty("OutputWorkspace").value)

    def _group_ws(self, input_ws, output_ws):
        group_alg = self.createChildAlgorithm("GroupWorkspaces", enableLogging=False)
        group_alg.setProperty("InputWorkspaces", input_ws)
        group_alg.setProperty("OutputWorkspace", output_ws)
        group_alg.execute()
        mtd.addOrReplace(output_ws, group_alg.getProperty("OutputWorkspace").value)

    def _copy_log(self, input_ws, output_ws):
        copy_log_alg = self.createChildAlgorithm("CopyLogs", enableLogging=False)
        copy_log_alg.setProperty("InputWorkspace", input_ws)
        copy_log_alg.setProperty("OutputWorkspace", output_ws)
        copy_log_alg.execute()
        mtd.addOrReplace(output_ws, copy_log_alg.getProperty("OutputWorkspace").value)

    def _add_sample_log_mult(self, input_ws, log_names, log_values):
        sample_log_mult_alg = self.createChildAlgorithm("AddSampleLogMultiple", enableLogging=False)
        sample_log_mult_alg.setProperty("Workspace", input_ws)
        sample_log_mult_alg.setProperty("LogNames", log_names)
        sample_log_mult_alg.setProperty("LogValues", log_values)
        sample_log_mult_alg.execute()


# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectTwoPeakFit)
