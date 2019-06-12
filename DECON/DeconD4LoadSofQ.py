from mantid.simpleapi import (mtd, CreateWorkspace, CloneWorkspace, DeleteWorkspace, AddSampleLogMultiple)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction, Progress)
from mantid.kernel import (StringListValidator, StringMandatoryValidator, IntBoundedValidator,
                           FloatBoundedValidator, Direction, logger)
from mantid import config
import math, os.path, numpy as np

class DeconD4LoadSofQ(DataProcessorAlgorithm):
 
    _input_path = None
    _input_name = None
    _temp = None
    _type  = None
    _theta_ws  = None
    _q_ws  = None
    _lambda = None
    _azero = None

    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Read S(Q) file in ascii format."

    def PyInit(self):
        self.declareProperty(name='FilePath', defaultValue='',
                             doc = 'File path for structure')
        self.declareProperty(name = 'Wavelength', defaultValue=0.7,
                             doc = 'Wavelength for output')
        self.declareProperty(name='ZeroAngle', defaultValue=0.0,
                             doc = 'Zero angle for Q convertion')
 
    def PyExec(self):
        self._setup()
        self._read_sofq()
        self._delete_ws(self._temp)

    def _setup(self):
        self._input_path = self.getProperty('FilePath').value
        self._get_filename()
        self._lambda = self.getProperty('Wavelength').value
        self._azero = self.getProperty('ZeroAngle').value

    def _get_filename(self):
        if(os.path.isfile(self._input_path)): 
            base = os.path.basename(self._input_path)
            self._input_name = os.path.splitext(base)[0]
            ext = os.path.splitext(base)[1]
        else:
            raise ValueError('Could not find file: %s' % self._input_path)

        logger.information('SofQ extension : %s' % (ext))
        self._theta_ws = self._input_name + '_theta'              #S(theta) WS
        if ext == '.sot':
            self._type = 'angle'
        elif ext == '.soq':
            self._type = 'Q'
            self._q_ws = self._input_name + '_Q'                 #S(Q) WS
        else:
            raise ValueError('File type not angle or Q')
        logger.information('SofQ input data is in %s ; file : %s' % (self._type, self._input_path))

    def _read_sofq(self):  #reads the structure data
        read_prog = Progress(self, start=0.0, end=0.1, nreports=3)
        read_prog.report('Reading data ')

        handle = open(self._input_path, 'r')
        asc = []
        for line in handle:                    #read lines into list 'asc'
            line = line.rstrip()
            asc.append(line)
        handle.close()
        len_asc = len(asc)
        len_head, head = self._read_header(asc)        #find header block
        self._sofq_header  = [('input_type', self._type), ('sofq_lines', len_head)] 
        for m in range(0, len_head):          #number of header lines
            self._sofq_header.append(('sofq_%i' % (m), head[m]))

        #get data from list
        x, y, e = self._read_sofq_data(asc, len_head, len_asc)
        dataX = np.array(x)
        dataY = np.array(y)
        dataE = np.array(e)
        self._temp = '__temp'
        self._create_ws(self._temp, dataX, dataY, dataE)
        log_names = [item[0] for item in self._sofq_header]
        log_values = [item[1] for item in self._sofq_header]
        self._add_sample_log_mult(self._temp, log_names, log_values)
        if self._type == 'angle':       #_input_ws = _theta_temp
            self._clone_ws(self._temp, self._theta_ws)
            unitx = mtd[self._theta_ws].getAxis(0).setUnit("Label")
            unitx.setLabel('2theta', 'deg')

        if self._type == 'Q':
            self._clone_ws(self._temp, self._q_ws)
            self._q_to_theta()       #converts _input_ws from Q to theta in _theta_tmp

        read_prog.report('Reading data completed')
#        logger.information('Theta workspace : %s' % self._theta_ws)

    def _q_to_theta(self):   #convert from Q to 2theta
        k0 = 4.0 * math.pi / float(self._lambda)
        self._clone_ws(self._temp, self._theta_ws)
        x_th = mtd[self._theta_ws].dataX(0)
        x_th = 2.0 * np.degrees(np.arcsin(x_th / k0))    #convert to angle
        x_th = x_th - self._azero                   #apply zero angle correction
        mtd[self._theta_ws].setX(0, x_th)
        unitx = mtd[self._theta_ws].getAxis(0).setUnit("Label")
        unitx.setLabel('2theta', 'deg')
        self._sofq_logs = [('lambda_in', self._lambda), ('zero_in', self._azero)]
        log_names = [item[0] for item in self._sofq_logs]
        log_values = [item[1] for item in self._sofq_logs]
        self._add_sample_log_mult(self._theta_ws, log_names, log_values)
        logger.information('Convert Q to 2theta : lambda = %f ; zero= %f' %
                            (self._lambda, self._azero))
        logger.information('2Theta : %f to %f' % (x_th[0], x_th[len(x_th) - 1]))

    def _read_sofq_data(self, asc, l_head, l_asc):    #reads the S(theta) data
        x = []
        y = []
        e = []
        n = 0
        for m in range(l_head + 1, l_asc):
            val = self._read_line(asc[m])          #splits line into 3 entries
            x.append(val[0])
            y.append(val[1])
            e.append(val[2])
            n += 1
        l_data = n-1
        logger.information('Data :')
        logger.information('Point   1 : %f %f %f' % (x[0], y[0], e[0]))
        logger.information('Point %i : %f %f %f' %
                           (l_data, x[l_data], y[l_data], e[l_data]))
        return x, y, e
	
    def _read_line(self, asc):
        extracted = []
        elements = asc.split()							#split line on spaces
        for n in elements:
            extracted.append(float(n))
        return extracted

    def _read_header(self, asc):
        head = []
        lines = 0
        for m in range(30):
            char = asc[m]
            if char.startswith('#'):        #check if line begins with a #
                head.append(asc[m])         #list of lines
                lines = m                   #number of lines
        logger.information('Data Header : ')
        for m in range(0, lines - 2):
            logger.information(head[m])
        return lines, head

    def _create_ws(self, output_ws, x, y, e):
        create_alg = self.createChildAlgorithm("CreateWorkspace", enableLogging = False)
        create_alg.setProperty("OutputWorkspace", output_ws)
        create_alg.setProperty("DataX", x)
        create_alg.setProperty("DataY", y)
        create_alg.setProperty("DataE", e)
        create_alg.setProperty("Nspec", 1)
        create_alg.setProperty("UnitX", 'MomentumTransfer')
        create_alg.setProperty("Distribution", True)
        create_alg.execute()
        mtd.addOrReplace(output_ws, create_alg.getProperty("OutputWorkspace").value)
        if self._type == 'angle':       #_input_ws = _theta_temp
            unitx = mtd[output_ws].getAxis(0).setUnit("Label")
            unitx.setLabel('2theta', 'deg')

    def _add_sample_log_mult(self, input_ws, log_names, log_values):
        sample_log_mult_alg = self.createChildAlgorithm("AddSampleLogMultiple", enableLogging=False)
        sample_log_mult_alg.setProperty("Workspace", input_ws)
        sample_log_mult_alg.setProperty("LogNames", log_names)
        sample_log_mult_alg.setProperty("LogValues", log_values)
        sample_log_mult_alg.execute()

    def _clone_ws(self, input_ws, output_ws):
        clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
        clone_alg.setProperty("InputWorkspace", input_ws)
        clone_alg.setProperty("OutputWorkspace", output_ws)
        clone_alg.execute()
        mtd.addOrReplace(output_ws, clone_alg.getProperty("OutputWorkspace").value)

    def _delete_ws(self, input_ws):
        delete_alg = self.createChildAlgorithm("DeleteWorkspace", enableLogging=False)
        delete_alg.setProperty("Workspace", input_ws)
        delete_alg.execute()

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconD4LoadSofQ)
#
