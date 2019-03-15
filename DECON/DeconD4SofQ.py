from mantid.simpleapi import (mtd, SaveNexusProcessed, DeconD4LoadSofQ, DeconCalculateDerivatives)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction)
from mantid.kernel import (StringListValidator, StringMandatoryValidator, IntBoundedValidator,
                           FloatBoundedValidator, Direction, logger)
from mantid import config
import math, os.path, numpy as np

class DeconD4SofQ(DataProcessorAlgorithm):
 
    _input_path = None
    _type = None
    _input_ws = None
    _theta_name	= None
    _deriv_name	= None
    _lambda = None
    _azero = None
    _plot = False
    _save = False

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Read S(Q) file in ascii format."

    def PyInit(self):
        self.declareProperty(FileProperty('File', '', action = FileAction.OptionalLoad,
                                          extensions = [".sot",".soq"]),
                             doc = 'File path for structure')
        self.declareProperty(name = 'Wavelength for input', defaultValue=0.7,
                             doc = 'Wavelength for output')
        self.declareProperty(name='Zero angle for input', defaultValue=0.0,
                             doc = 'Zero angle for Q convertion')

        self.declareProperty(name='Plot', defaultValue=False,
                             doc = 'Plot options')
        self.declareProperty(name='Save', defaultValue=False,
                             doc = 'Switch Save result to nxs file Off/On')
 
    def PyExec(self):
        self._setup()
        DeconD4LoadSofQ(FilePath=self._input_path,
                        Wavelength=self._lambda,
                        ZeroAngle=self._azero)

        self._deriv_name = self._input_name + '_deriv'
        logger.information('Derivatives workspace : %s' % self._deriv_name)
        DeconCalculateDerivatives(DataWorkspace=self._theta_name,
                                  DerivativesWorkspace=self._deriv_name)

        if self._save:                        #Save option
            self._save_result()
        if self._plot:                       #Plot option
            self._plot_result()

    def _setup(self):
        self._input_path = self.getPropertyValue('File')
        self._get_filename()
        self._lambda = self.getProperty('Wavelength for input').value
        self._azero = self.getProperty('Zero angle for input').value

        self._plot = self.getProperty('Plot').value
        self._save = self.getProperty('Save').value

    def _save_result(self):                        #Save 
        self._save_ws(self._theta_name, 'S(theta)')
        self._save_ws(self._deriv_name, 'Deriv')

        if self._type == 'Q':
            self._save_ws(self._q_name, 'S(Q)')
			
    def _plot_result(self):                       #Plot 
        import mantidplot as mp
        mp.plotSpectrum(self._theta_name, [0,1,2])
        mp.plotSpectrum(self._deriv_name, [0,1,2,3])

    def _get_filename(self):
        if(os.path.isfile(self._input_path)): 
            base = os.path.basename(self._input_path)
            self._input_name = os.path.splitext(base)[0]
            ext = os.path.splitext(base)[1]
        else:
            raise ValueError('Could not find file: %s' % self._input_path)

        logger.information('SofQ extension : %s' % (ext))
        self._theta_name = self._input_name + '_theta'              #S(theta) WS
        if ext == '.sot':
            self._type = 'angle'
        elif ext == '.soq':
            self._type = 'Q'
            self._q_name = self._input_name + '_Q'                 #S(Q) WS
        else:
            raise ValueError('File type not angle or Q')
        logger.information('SofQ input data is in %s ; file : %s' % (self._type, self._input_path))

    def _save_ws(self, input_ws, text):
        workdir = config['defaultsave.directory']
        path = os.path.join(workdir, input_ws + '.nxs')
        save_alg = self.createChildAlgorithm("SaveNexusProcessed", enableLogging = True)
        save_alg.setProperty("InputWorkspace", input_ws)
        save_alg.setProperty("Filename", path)
        save_alg.execute()
        logger.information('%s file saved as %s' % (text, path))

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconD4SofQ)
#
