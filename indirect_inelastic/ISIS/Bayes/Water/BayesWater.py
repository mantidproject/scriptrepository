#pylint: disable=invalid-name,too-many-instance-attributes,too-many-branches,no-init
from IndirectImport import *

from mantid.api import (PythonAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, PropertyMode,
                        WorkspaceGroupProperty, Progress)
from mantid.kernel import StringListValidator, Direction
from mantid.simpleapi import *
from mantid import config, logger
from IndirectCommon import *
import os
import numpy as np

if is_supported_f2py_platform():
    QLwat     = import_f2py("QLwat")

class BayesWater(PythonAlgorithm):

    _program = 'QLwat'
    _sam_ws = None
    _res_ws = None
    _resnorm_ws = None
    _e_min = None
    _e_max = None
    _sam_bins = None
    _res_bins = None
    _elastic = None
    _background = None
    _res_norm = None
    _width = False
    _wfile = None
    _loop = None
 
    def category(self):
        return "PythonAlgorithms;Workflow\\MIDAS"

    def summary(self):
        return "This algorithm runs the Fortran QWater program which fits Water function. The"+\
               " whole function is then convoled with the resolution function."

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('SampleWorkspace', '', direction=Direction.Input),
                             doc='Name of the Sample input Workspace')

        self.declareProperty(MatrixWorkspaceProperty('ResolutionWorkspace', '', direction=Direction.Input),
                             doc='Name of the resolution input Workspace')

        self.declareProperty(name='MinRange', defaultValue=-0.2,
                             doc='The start of the fit range. Default=-0.2')

        self.declareProperty(name='MaxRange', defaultValue=0.2,
                             doc='The end of the fit range. Default=0.2')

        self.declareProperty(name='SampleBins', defaultValue=1,
                             doc='The number of sample bins')

        self.declareProperty(name='ResolutionBins', defaultValue=1,
                             doc='The number of resolution bins')

        self.declareProperty(name='Elastic', defaultValue=True,
                             doc='Fit option for using the elastic peak')

        self.declareProperty(name='Background', defaultValue='Flat',
                             validator=StringListValidator(['Sloping','Flat','Zero']),
                             doc='Fit option for the type of background')

        self.declareProperty(name='Loop', defaultValue=True,
                             doc='Switch Sequential fit On/Off')

        self.declareProperty(WorkspaceGroupProperty('OutputWorkspaceFit', '', direction=Direction.Output),
                             doc='The name of the fit output workspaces')

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspaceResult', '', direction=Direction.Output),
                             doc='The name of the result output workspaces')


    def validateInputs(self):
        self._get_properties()
        issues = dict()

        # Validate fitting range in energy
        if self._e_min > self._e_max:
            issues['MaxRange'] = 'Must be less than EnergyMin'

        return issues

    def _get_properties(self):
        self._program = 'QLwat'
        self._sam_ws = self.getPropertyValue('SampleWorkspace')
        self._res_ws = self.getPropertyValue('ResolutionWorkspace')
        self._resnorm_ws = ''
        self._res_norm = False
        self._width = False
        self._wfile = ''

        self._e_min = self.getProperty('MinRange').value
        self._e_max = self.getProperty('MaxRange').value
        self._sam_bins = self.getPropertyValue('SampleBins')
        self._res_bins = self.getPropertyValue('ResolutionBins')
        self._elastic = self.getProperty('Elastic').value
        self._background = self.getPropertyValue('Background')

        self._loop = self.getProperty('Loop').value

    def PyExec(self):
        run_f2py_compatibility_test()

        from IndirectBayes import (CalcErange, GetXYE, C2Fwater)
        from IndirectCommon import (CheckXrange, CheckAnalysers, getEfixed, GetThetaQ,
                                    CheckHistZero, CheckHistSame)

        erange = [self._e_min, self._e_max]
        nbins = [self._sam_bins, self._res_bins]
        #convert true/false to 1/0 for fortran
        o_el = 1 if self._elastic else 0
        o_w1 = 1 if self._width else 0
        o_res = 1 if self._res_norm else 0

        #fortran code uses background choices defined using the following numbers
        if self._background == 'Sloping':
            o_bgd = 2
        elif self._background == 'Flat':
            o_bgd = 1
        elif self._background == 'Zero':
            o_bgd = 0

        fitOp = [o_el, o_bgd, o_w1, o_res]

        workdir = config['defaultsave.directory']
        if not os.path.isdir(workdir):
            workdir = os.getcwd()
            logger.information('Default Save directory is not set. Defaulting to current working Directory: ' + workdir)

        array_len = 4096                           # length of array in Fortran
        CheckXrange(erange, 'Energy')

        nbin,nrbin = nbins[0], nbins[1]

        logger.information('Sample is ' + self._sam_ws)
        logger.information('Resolution is ' + self._res_ws)

        CheckAnalysers(self._sam_ws, self._res_ws)
        efix = getEfixed(self._sam_ws)
        theta, Q = GetThetaQ(self._sam_ws)

        nsam,ntc = CheckHistZero(self._sam_ws)
        totalNoSam = nsam

        #check if we're performing a sequential fit
        if self._loop != True:
            nsam = 1

        prog = 'QLw'                        # res file
        logger.information('Program is ' + prog)
        logger.information(' Number of spectra = %i' % nsam)
        logger.information(' Erange : %f to %f' % (erange[0], erange[1]))

        dtn, xsc = self._read_width_file(self._res_norm, self._resnorm_ws, totalNoSam)
        Wy, We = self._read_norm_file(self._width,self._wfile,totalNoSam)

        fname = self._sam_ws[:-4] + '_'+ prog
        fit_ws = fname + '_Fit'
        wrks=os.path.join(workdir, self._sam_ws[:-4])
        logger.information(' lptfile : ' + wrks + '_' + prog + '.lpt')
        lwrk=len(wrks)
        wrks.ljust(140, ' ')
        wrkr=self._res_ws
        wrkr.ljust(140, ' ')

        group = ''
        workflow_prog = Progress(self, start=0.3, end=0.7, nreports=nsam*3)
        for m in range(nsam):
            logger.information('Group %i at angle %f' % (m, theta[m]))
            nsp = m + 1
            nout, bnorm, Xdat, Xv, Yv, Ev = CalcErange(self._sam_ws, m, erange, nbin)
            Ndat = nout[0]
            Imin = nout[1]
            Imax = nout[2]
            if prog == 'QLd':
                mm = m
            else:
                mm = 0
            Nb, Xb, Yb, Eb = GetXYE(self._res_ws, mm, array_len)     # get resolution data
            numb = [nsam, nsp, ntc, Ndat, nbin, Imin, Imax, Nb, nrbin]
            rscl = 1.0
            reals = [efix, theta[m], rscl, bnorm]

            workflow_prog.report('Processing spectrum number %i as Water' % m)
            nd,xout,yout,eout,yfit,yprob=QLwat.qlwat(numb, Xv, Yv, Ev, reals, fitOp,
													 Xdat, Xb, Yb, Wy, We, dtn, xsc,
													 wrks, wrkr, lwrk)

            dataX = xout[:nd]
            dataX = np.append(dataX, 2*xout[nd - 1] - xout[nd - 2])
            yfit_list = np.split(yfit[:4*nd], 4)
            dataF1 = yfit_list[1]
            workflow_prog.report('Processing data')
            dataG = np.zeros(nd)
            datX = dataX
            datY = yout[:nd]
            datE = eout[:nd]
            datX = np.append(datX, dataX)
            datY = np.append(datY, dataF1[:nd])
            datE = np.append(datE, dataG)
            res1 = dataF1[:nd] - yout[:nd]
            datX = np.append(datX, dataX)
            datY = np.append(datY, res1)
            datE = np.append(datE, dataG)
            nsp = 3
            names = 'data, fit.1, diff.1'
            res_plot = [0, 1, 2]

            # create result workspace
            fit_ws = fname + '_Workspaces'
            fout = fname + '_Workspace_' + str(m)

            workflow_prog.report('Creating OutputWorkspace')
            CreateWorkspace(OutputWorkspace=fout,
                            DataX=datX,
                            DataY=datY,
                            DataE=datE,
                            Nspec=nsp,
                            UnitX='DeltaE',
                            VerticalAxisUnit='Text',
                            VerticalAxisValues=names)

            # append workspace to list of results
            group += fout + ','

        comp_prog = Progress(self, start=0.7, end=0.8, nreports=2)
        comp_prog.report('Creating Group Workspace')
        GroupWorkspaces(InputWorkspaces=group,
                        OutputWorkspace=fit_ws)

        out_ws = C2Fwater(fname)

        #Add some sample logs to the output workspaces
        CopyLogs(InputWorkspace=self._sam_ws,
                 OutputWorkspace=out_ws)
        self._add_sample_logs(out_ws, prog, erange, nbins)
        CopyLogs(InputWorkspace=self._sam_ws,
                 OutputWorkspace=fit_ws)
        self._add_sample_logs(fit_ws, prog, erange, nbins)

        self.setProperty('OutputWorkspaceFit', fit_ws)
        self.setProperty('OutputWorkspaceResult', out_ws)

    def _add_sample_logs(self, workspace, fit_program, e_range, binning):

        sample_binning, res_binning = binning
        energy_min, energy_max = e_range

        sample_logs = [('res_workspace', self._res_ws),
                       ('fit_program', fit_program),
                       ('background', self._background),
                       ('elastic_peak', self._elastic),
                       ('energy_min', energy_min),
                       ('energy_max', energy_max),
                       ('sample_binning', sample_binning),
                       ('resolution_binning', res_binning)]

        resnorm_used = (self._resnorm_ws != '')
        sample_logs.append(('resnorm', str(resnorm_used)))
        if resnorm_used:
            sample_logs.append(('resnorm_file', str(self._resnorm_ws)))

        width_file_used = (self._wfile != '')
        sample_logs.append(('width', str(width_file_used)))
        if width_file_used:
            sample_logs.append(('width_file', str(self._wfile)))

        log_alg = self.createChildAlgorithm('AddSampleLogMultiple', 0.9, 1.0, False)
        log_alg.setProperty('Workspace', workspace)
        log_alg.setProperty('LogNames', [log[0] for log in sample_logs])
        log_alg.setProperty('LogValues', [log[1] for log in sample_logs])
        log_alg.execute()

    def _read_norm_file(self, readRes, resnormWS, nsam):  # get norm & scale values
        resnorm_root = resnormWS
        # Obtain root of resnorm group name
        if '_Intensity' in resnormWS:
            resnorm_root = resnormWS[:-10]
        if '_Stretch' in resnormWS:
            resnorm_root = resnormWS[:-8]

        if readRes:  # use ResNorm file option=o_res
            Xin = s_api.mtd[resnorm_root + '_Intensity'].readX(0)
            nrm = len(Xin)  # no. points from length of x array
            if nrm == 0:
                raise ValueError('ResNorm file has no Intensity points')
            Xin = s_api.mtd[resnorm_root + '_Stretch'].readX(0)  # no. points from length of x array
            if len(Xin) == 0:
                raise ValueError('ResNorm file has no xscale points')
            if nrm != nsam:  # check that no. groups are the same
                raise ValueError('ResNorm groups (' + str(nrm) + ') not = Sample (' + str(nsam) + ')')
            else:
                dtn, xsc = self._get_res_norm(resnorm_root, 0)
        else:
            # do not use ResNorm file
            dtn, xsc = self._get_res_norm(resnorm_root, nsam)
        return dtn, xsc

    def _get_res_norm(self, resnormWS, ngrp):
        if ngrp == 0:  # read values from WS
            dtnorm = s_api.mtd[resnormWS + '_Intensity'].readY(0)
            xscale = s_api.mtd[resnormWS + '_Stretch'].readY(0)
        else:  # constant values
            dtnorm = []
            xscale = []
            for _ in range(0, ngrp):
                dtnorm.append(1.0)
                xscale.append(1.0)
        dtn = PadArray(dtnorm, 51)  # pad for Fortran call
        xsc = PadArray(xscale, 51)
        return dtn, xsc

    # Reads in a width ASCII file
    def _read_width_file(self, readWidth, widthFile, numSampleGroups):
        widthY, widthE = [], []
        if readWidth:
            logger.information('Width file is ' + widthFile)
            # read ascii based width file
            try:
                wfPath = s_api.FileFinder.getFullPath(widthFile)
                handle = open(wfPath, 'r')
                asc = []
                for line in handle:
                    line = line.rstrip()
                    asc.append(line)
                handle.close()
            except Exception:
                raise ValueError('Failed to read width file')
            numLines = len(asc)
            if numLines == 0:
                raise ValueError('No groups in width file')
            if numLines != numSampleGroups:  # check that no. groups are the same
                raise ValueError('Width groups (' + str(numLines) + ') not = Sample (' + str(numSampleGroups) + ')')
        else:
            # no file: just use constant values
            widthY = np.zeros(numSampleGroups)
            widthE = np.zeros(numSampleGroups)
        # pad for Fortran call
        widthY = PadArray(widthY, 51)
        widthE = PadArray(widthE, 51)

        return widthY, widthE

if is_supported_f2py_platform():
    # Register algorithm with Mantid
    AlgorithmFactory.subscribe(BayesWater)         # Register algorithm with Mantid
