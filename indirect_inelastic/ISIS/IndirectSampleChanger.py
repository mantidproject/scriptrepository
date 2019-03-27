# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import os.path

from mantid import config
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
import mantidplot as mp


class IndirectSampleChanger(DataProcessorAlgorithm):
    _plot = False
    _save = False
    _instrument = 'IRIS'
    _number_runs = 1
    _runs = None
    _temperature = None
    _spectra_range = None
    _elastic_range = None
    _inelastic_range = None
    _total_range = None
    _sample_log_name = None
    _sample_log_value = None
    _msdfit = False
    _widthfit = False
    _q1_workspaces = None

    def category(self):
        return "Workflow\\MIDAS"

    def summary(self):
        return "Create elastic window scans for sample changer"

    def PyInit(self):
        self.declareProperty(name='Instrument', defaultValue='IRIS',
                             validator=StringListValidator(['IRIS', 'OSIRIS']),
                             doc='Instrument name')
        self.declareProperty(name='Analyser', defaultValue='',
                             validator=StringListValidator(['graphite', 'mica', 'fmica']),
                             doc='Analyser bank used during run.')
        self.declareProperty(name='Reflection', defaultValue='',
                             validator=StringListValidator(['002', '004', '006']),
                             doc='Reflection number for instrument setup during run.')

        self.declareProperty(name="FirstRun", defaultValue=-1,
                             validator=IntBoundedValidator(lower=0),
                             doc="First Sample run-number.")
        self.declareProperty(name='LastRun', defaultValue=-1,
                             validator=IntBoundedValidator(lower=0),
                             doc="Last Sample run-number.")
        self.declareProperty(name='NumberSamples', defaultValue=-1,
                             validator=IntBoundedValidator(lower=0),
                             doc="Increment for run-number.")

        self.declareProperty(IntArrayProperty(name='SpectraRange', values=[0, 1],
                                              validator=IntArrayLengthValidator(2)),
                             doc='Comma separated range of spectra number to use.')
        self.declareProperty(FloatArrayProperty(name='ElasticRange',
                                                validator=FloatArrayLengthValidator(2)),
                             doc='Range of background to subtract from raw data in time of flight.')
        self.declareProperty(FloatArrayProperty(name='InelasticRange',
                                                validator=FloatArrayLengthValidator(2)),
                             doc='Range of background to subtract from raw data in time of flight.')
        self.declareProperty(FloatArrayProperty(name='TotalRange'),
                             doc='Energy range for the total energy component.')

        self.declareProperty(name='MsdFit', defaultValue=False,
                             doc='Run msd fit')

        self.declareProperty(name='WidthFit', defaultValue=False,
                             doc='Perform a 2 peak width Fit, do not use with GroupingMethod as "All"')

        self.declareProperty(name='Plot', defaultValue=False,
                             doc='Switch Plot Off/On')
        self.declareProperty(name='Save', defaultValue=False,
                             doc='Switch Save result to nxs file Off/On')

    def PyExec(self):
        from IndirectImport import import_mantidplot
        mp = import_mantidplot()
        workdir = config['defaultsave.directory']
        # self._setup()

        q2_workspaces = []
        scan_alg = self.createChildAlgorithm("EnergyWindowScan", 0.05, 0.95)
        for numb in range(self._number_samples):
            run_numbers = []
            run_names = []
            first_run = self._run_first + numb
            for idx in range(int(self._number_runs)):
                run = str(first_run + idx * self._number_samples)
                run_numbers.append(run)
                run_names.append(self._instrument + run)
            q0 = self._instrument.lower() + run_numbers[0] + '_to_' + run_numbers[-1] + '_s' + str(numb)
            output_ws = q0 + '_red'
            scan_ws = q0 + '_scan'
            scan_alg.setProperty('InputFiles', run_names)
            scan_alg.setProperty('LoadLogFiles', True)
            scan_alg.setProperty('CalibrationWorkspace', '')
            scan_alg.setProperty('Instrument', self._instrument_name)
            scan_alg.setProperty('Analyser', self._analyser)
            scan_alg.setProperty('Reflection', self._reflection)
            scan_alg.setProperty('SpectraRange', self._spectra_range)
            scan_alg.setProperty('ElasticRange', self._elastic_range)
            scan_alg.setProperty('InelasticRange', self._inelastic_range)
            scan_alg.setProperty('TotalRange', self._total_range)
            scan_alg.setProperty('DetailedBalance', Property.EMPTY_DBL)
            scan_alg.setProperty('GroupingMethod', 'Individual')
            scan_alg.setProperty('SampleEnvironmentLogName', self._sample_log_name)
            scan_alg.setProperty('SampleEnvironmentLogValue', self._sample_log_value)
            scan_alg.setProperty('msdFit', self._msdfit)
            scan_alg.setProperty('ReducedWorkspace', output_ws)
            scan_alg.setProperty('ScanWorkspace', scan_ws)
            scan_alg.execute()

            logger.information('OutputWorkspace : %s' % output_ws)
            logger.information('ScanWorkspace : %s' % scan_ws)

            q1_ws = scan_ws + '_el_eq1'
            q2_ws = scan_ws + '_el_eq2'
            q2_workspaces.append(q2_ws)
            eisf_ws = scan_ws + '_eisf'
            el_elt_ws = scan_ws + '_el_elt'
            inel_elt_ws = scan_ws + '_inel_elt'
            tot_elt_ws = scan_ws + '_total_elt'
#            output_workspaces = [q1_ws, q2_ws, eisf_ws, el_elt_ws, inel_elt_ws, tot_elt_ws]
            output_workspaces = [q1_ws, eisf_ws, el_elt_ws, inel_elt_ws, tot_elt_ws]

            if self._plot:
                for ws in output_workspaces:
                    mp.plotSpectrum(ws, 0, error_bars=True)
                if self._msdfit:
                    mp.plotSpectrum(scan_ws + '_msd', 0, error_bars=True)

        if self._widthfit:
            result_workspaces = list()
            chi_workspaces = list()
            temperatures = list()
        # Get input workspaces
            fit_progress = Progress(self, 0.0, 0.05, 3)
            input_workspace_names = mtd[output_ws].getNames()
            x = mtd[input_workspace_names[0]].readX(0)
            xmin = x[0]
            xmax = x[len(x) - 1]
            for input_ws in input_workspace_names:

                red_ws = input_ws[:-3] + 'red'
                # Get the sample temperature
                temp = self._get_temperature(red_ws)
                if temp is not None:
                    temperatures.append(temp)
                else:
                # Get the run number
                    run_no = self._get_InstrRun(input_ws)[1]
                    run_numbers.append(run_no)

                num_hist = mtd[input_ws].getNumberHistograms()
                logger.information('Reduced histograms : %i' % num_hist)
                result = input_ws[:-3] + 'fit'
                func = 'name=Lorentzian,Amplitude=1.0,PeakCentre=0.0,FWHM=0.01'
                func += ',constraint=(Amplitude>0.0,FWHM>0.0)'
                for idx in range(num_hist):
                    fit_progress.report('Fitting workspace: %s ; spectrum %i' % (input_ws, idx))
                    IndirectTwoPeakFit(SampleWorkspace=input_ws,
                                       EnergyMin=xmin,
                                       EnergyMax=xmax,
                                       Minimizer='Levenberg-Marquardt',
                                       MaxIterations=500,
                                       OutputName=result)
                result_workspaces.append(result + '_Result')
                chi_workspaces.append(result + '_ChiSq')


    def _setup(self):
        self._run_first = self.getProperty('FirstRun').value
        self._run_last = self.getProperty('LastRun').value
        self._number_samples = self.getProperty('NumberSamples').value
        self._number_runs = (self._run_last - self._run_first + 1) / self._number_samples
        logger.information('Number of runs : %i' % self._number_runs)
        logger.information('Number of scans : %i' % self._number_samples)

        self._instrument_name = self.getPropertyValue('Instrument')
        self._analyser = self.getPropertyValue('Analyser')
        self._reflection = self.getPropertyValue('Reflection')

        self._spectra_range = self.getProperty('SpectraRange').value
        self._elastic_range = self.getProperty('ElasticRange').value
        self._inelastic_range = self.getProperty('InelasticRange').value
        self._total_range = self.getProperty('TotalRange').value
 
        self._sample_log_name = 'Position'
        self._sample_log_value = 'last_value'

        self._msdfit = self.getProperty('MsdFit').value

        self._widthfit = self.getProperty('WidthFit').value

        self._plot = self.getProperty('Plot').value
        self._save = self.getProperty('Save').value

    def validateInputs(self):
        self._setup()
        issues = dict()

        if self._run_first > self._run_last:
            issues["FirstRun"] = 'First run must be before last run'

        if self._number_runs < self._number_samples:
            issues["NumberSamples"] = 'There must be at least 1 run per sample'

        if self._spectra_range[0] > self._spectra_range[1]:
            issues['SpectraRange'] = 'Range must be in format: lower,upper'

        return issues

    def _get_temperature(self, ws_name):
        """
        Gets the sample temperature for a given workspace.

        @param ws_name Name of workspace
        @returns Temperature in Kelvin or None if not found
        """
        instr, run_number = self._get_InstrRun(ws_name)

        facility = config.getFacility()
        pad_num = facility.instrument(instr).zeroPadding(int(run_number))
        zero_padding = '0' * (pad_num - len(run_number))

        run_name = instr + zero_padding + run_number
        log_filename = run_name.upper() + '.log'

        run = mtd[ws_name].getRun()

        if self._sample_log_name in run:
            # Look for temperature in logs in workspace
            tmp = run[self._sample_log_name].value
            value_action = {'last_value': lambda x: x[len(x) - 1],
                            'average': lambda x: x.mean()
                            }
            temp = value_action[self._sample_log_value](tmp)
            logger.debug('Temperature %s K found for run: %s' % (temp, run_name))
            return temp

        else:
            # Logs not in workspace, try loading from file
            logger.information('Log parameter not found in workspace. Searching for log file.')
            log_path = FileFinder.getFullPath(log_filename)

            if log_path != '':
                # Get temperature from log file
                LoadLog(Workspace=ws_name, Filename=log_path)
                run_logs = mtd[ws_name].getRun()
                if self._sample_log_name in run_logs:
                    tmp = run_logs[self._sample_log_name].value
                    temp = tmp[len(tmp) - 1]
                    logger.debug('Temperature %d K found for run: %s' % (temp, run_name))
                    return temp
                else:
                    logger.warning('Log entry %s for run %s not found' % (self._sample_log_name, run_name))
            else:
                logger.warning('Log file for run %s not found' % run_name)

        # Can't find log file
        logger.warning('No temperature found for run: %s' % run_name)
        return None

    def _get_InstrRun(self, ws_name):
        """
        Get the instrument name and run number from a workspace.

        @param ws_name - name of the workspace
        @return tuple of form (instrument, run number)
        """

        run_number = str(mtd[ws_name].getRunNumber())
        if run_number == '0':
            # Attempt to parse run number off of name
            match = re.match(r'([a-zA-Z]+)([0-9]+)', ws_name)
            if match:
                run_number = match.group(2)
            else:
                raise RuntimeError("Could not find run number associated with workspace.")

        instrument = mtd[ws_name].getInstrument().getName()
        if instrument != '':
            for facility in config.getFacilities():
                try:
                    instrument = facility.instrument(instrument).filePrefix(int(run_number))
                    instrument = instrument.lower()
                    break
                except RuntimeError:
                    continue

        return instrument, run_number

    def _save_output(self, input_ws):
        from mantid.simpleapi import SaveNexusProcessed
        workdir = config['defaultsave.directory']
        el_eq1_path = os.path.join(workdir, input_ws + '_el_eq1.nxs')
        logger.information('Creating file : %s' % el_eq1_path)
        self._save_ws(input_ws + '_el_eq1', el_eq1_path)
        el_eq2_path = os.path.join(workdir, input_ws + '_el_eq2.nxs')
        logger.information('Creating file : %s' % el_eq2_path)
        self._save_ws(input_ws + '_el_eq2', el_eq2_path)

        inel_eq1_path = os.path.join(workdir, input_ws + '_inel_eq1.nxs')
        logger.information('Creating file : %s' % inel_eq1_path)
        self._save_ws(input_ws + '_inel_eq1', inel_eq1_path)
        inel_eq2_path = os.path.join(workdir, input_ws + '_inel_eq2.nxs')
        logger.information('Creating file : %s' % inel_eq2_path)
        self._save_ws(input_ws + '_inel_eq2', inel_eq2_path)

        eisf_path = os.path.join(workdir, input_ws + '_eisf.nxs')
        logger.information('Creating file : %s' % eisf_path)
        self._save_ws(input_ws + '_eisf', eisf_path)

        if self._msdfit:
            msd_path = os.path.join(workdir, input_ws + '_msd.nxs')
            logger.information('Creating file : %s' % msd_path)
            self._save_ws(input_ws + '_msd', msd_path)
            msd_fit_path = os.path.join(workdir, input_ws + '_msd_fit.nxs')
            logger.information('Creating file : %s' % msd_fit_path)
            self._save_ws(input_ws + '_msd_fit', msd_fit_path)

#        if self._widthfit:
#            mp.plotSpectrum(self._output_ws + '_Diffusion', 0, error_bars=True)


AlgorithmFactory.subscribe(IndirectSampleChanger)  # Register algorithm with Mantid
