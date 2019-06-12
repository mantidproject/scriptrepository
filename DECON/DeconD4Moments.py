from mantid.simpleapi import (mtd, CreateWorkspace, DeleteWorkspace,
                              AddSampleLogMultiple, SaveNexusProcessed)
from mantid.api import (DataProcessorAlgorithm, AlgorithmFactory, PropertyMode, MatrixWorkspaceProperty,
                        WorkspaceGroupProperty, FileProperty, FileAction, Progress)
from mantid.kernel import *
from mantid import config
import os.path, numpy as np

class DeconD4Moments(DataProcessorAlgorithm):

    _mome_path = None
    _mome_name = None
    _freeze = False
    _freeze_pt = 0.0
    _plot = False
    _save = False
 
    def category(self):
        return "Workflow\\DECON;PythonAlgorithms"

    def summary(self):
        return "Read moments file in ascii format."

    def PyInit(self):
        self.declareProperty(FileProperty('Moments File', '', action = FileAction.Load,
                                          extensions = [".mom"]),
                             doc = 'File path for structure in angle')

        self.declareProperty(name = 'Freeze', defaultValue=False,
                             doc = 'Low cutoff option')
        freezeCondition = VisibleWhenProperty('Freeze',PropertyCriterion.IsNotDefault)
        self.declareProperty(name = 'Freeze point', defaultValue=0.0,
                             doc = 'Low cutoff value')
        self.setPropertySettings('Freeze point', freezeCondition)

        self.declareProperty(name = 'Plot', defaultValue=False,
                             doc = 'Plot options')
        self.declareProperty(name = 'Save', defaultValue=False,
                             doc = 'Switch Save result to nxs file Off/On')
 
    def PyExec(self):
        import mantidplot as mp
        self._setup()
        self._read_mome()

        if self._save:                        #Save option
            self._save_ws()
        if self._plot:                       #Plot option
            mp.plotSpectrum(self._mome_name, [0])
            mp.plotSpectrum(self._mome_name, [1,2,3,4])
 
    def _setup(self):
        self._mome_path = self.getPropertyValue('Moments File')
        self._get_filename()
        self._freeze = self.getProperty('Freeze').value
        self._freeze_pt = self.getProperty('Freeze point').value

        self._plot = self.getProperty('Plot').value
        self._save = self.getProperty('Save').value

    def _read_mome(self):
        read_prog = Progress(self, start=0.0, end=0.1, nreports=3)
        read_prog.report('Reading data ')

        handle = open(self._mome_path, 'r')
        asc = []
        for line in handle:                    #read lines into list 'asc'
            line = line.rstrip()
            asc.append(line)
        handle.close()
        len_asc = len(asc)
        len_head, head = self._read_header(asc)        #find header block
        mome_header  = [('mome_lines', len_head)] 
        for m in range(0, len_head):          #number of header lines
            mome_header.append(('mome_%i' % (m), head[m]))

        x, m0, m1, m2, m3, m4 = self._read_mome_data(asc,len_head,len_asc)   #get data from list
        xQ = np.array(x)
        read_prog.report('Reading data completed')

        create_prog = Progress(self, start=0.0, end=0.1, nreports=3)
        create_prog.report('Creating workspaces ')

        ifr = 0
        logger.information('Freeze : %s' % self._freeze)
        mome_logs = [('freeze', self._freeze)]
        if self._freeze:                             #Freeze option
            ifr = np.where(self._freeze_pt > xQ)[0][-1]
        logger.information('Freeze point : %f' % (x[ifr]))
        mome_logs.append(('freeze_point', x[ifr]))
        e0 = np.zeros(len(x))                  # set errors to zero
        dataX = xQ
        dataY = self._mome_cut(m0, ifr)      #apply freeze
        dataE = e0
        dataX = np.append(dataX,xQ)
        Ytmp = self._mome_cut(m1, ifr)
        dataY = np.append(dataY,Ytmp)
        dataE = np.append(dataE,e0)
        dataX = np.append(dataX,xQ)
        Ytmp = self._mome_cut(m2, ifr)
        dataY = np.append(dataY,Ytmp)
        dataE = np.append(dataE,e0)
        dataX = np.append(dataX,xQ)
        Ytmp = self._mome_cut(m3, ifr)
        dataY = np.append(dataY,Ytmp)
        dataE = np.append(dataE,e0)
        dataX = np.append(dataX,xQ)
        Ytmp = self._mome_cut(m4, ifr)
        dataY = np.append(dataY,Ytmp)
        dataE = np.append(dataE,e0)
        names = 'M0, M1, M2, M3, M4'              #names for hists
        self._create_ws(self._mome_name, dataX, dataY, dataE, names)
        log_names = [item[0] for item in mome_header]
        log_values = [item[1] for item in mome_header]
        self._add_sample_log_mult(self._mome_name, log_names, log_values)
        log_names = [item[0] for item in mome_logs]
        log_values = [item[1] for item in mome_logs]
        self._add_sample_log_mult(self._mome_name, log_names, log_values)
        create_prog.report('Creating workspaces completed')

    def _mome_cut(self, mome, ifr):   #applies the freeze option to moments
        Ytmp = np.array(mome)
        dataY = []
        if self._freeze:
            fr_array = np.empty(ifr)
            fr_array.fill(Ytmp[ifr])   #fill list with value at freeze point
            dataY = fr_array           #start off Yarray with freeze value
            dataY = np.append(dataY, Ytmp[ifr:])   #append values from data
        else:
            dataY = Ytmp               #copy all original
        return dataY

    def _save_ws(self):                        #Save 
        workdir = config['defaultsave.directory']
        mome_path = os.path.join(workdir, self._mome_name + '.nxs')
        save_alg = self.createChildAlgorithm("SaveNexusProcessed", enableLogging = False)
        save_alg.setProperty("InputWorkspace", self._mome_name)
        save_alg.setProperty("Filename", mome_path)
        save_alg.execute()
        logger.information('Moments file is %s' % mome_path)

			
    def _get_filename(self):
        path = self._mome_path
        if(os.path.isfile(path)): 
            base = os.path.basename(path)
            self._mome_name = os.path.splitext(base)[0] + '_mome'
            ext = os.path.splitext(base)[1]
            logger.information('Mome input file : %s' % path)
        else:
            raise ValueError('Could not find file: %s' % path)

    def _read_mome_data(self, asc, l_head, l_asc):    #reads the S(theta) data
        x = []
        m0 = []
        m1= []
        m2 = []
        m3 = []
        m4 = []
        n = 0
        for m in range(l_head + 1, l_asc):
            val = self._read_line(asc[m])          #splits the line into 6 entries
            x.append(val[0])
            m0.append(val[1])
            m1.append(val[2])
            m2.append(val[3])
            m3.append(val[4])
            m4.append(val[5])
            n += 1
        l_data = n - 1
        logger.information('Data :')
        logger.information('Point   1 : %f %f %f' % (x[0], m1[0], m2[0]))
        logger.information('Point %i : %f %f %f' %
                           (l_data, x[l_data], m1[l_data], m2[l_data]))
        return x, m0, m1, m2, m3, m4
	
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

    def _create_ws(self, output_ws, x, y, e, names):
        create_alg = self.createChildAlgorithm("CreateWorkspace", enableLogging = False)
        create_alg.setProperty("OutputWorkspace", output_ws)
        create_alg.setProperty("DataX", x)
        create_alg.setProperty("DataY", y)
        create_alg.setProperty("DataE", e)
        create_alg.setProperty("Nspec", 5)
        create_alg.setProperty("UnitX", 'MomentumTransfer')
        create_alg.setProperty("VerticalAxisUnit", 'Text')
        create_alg.setProperty("VerticalAxisValues", names)
        create_alg.setProperty("Distribution", True)
        create_alg.execute()
        mtd.addOrReplace(output_ws, create_alg.getProperty("OutputWorkspace").value)
        unitx = mtd[output_ws].getAxis(0).setUnit("Label")
        unitx.setLabel('2theta', 'deg')

    def _add_sample_log_mult(self, input_ws, log_names, log_values):
        sample_log_mult_alg = self.createChildAlgorithm("AddSampleLogMultiple", enableLogging=False)
        sample_log_mult_alg.setProperty("Workspace", input_ws)
        sample_log_mult_alg.setProperty("LogNames", log_names)
        sample_log_mult_alg.setProperty("LogValues", log_values)
        sample_log_mult_alg.execute()

# Register algorithm with Mantid
AlgorithmFactory.subscribe(DeconD4Moments)
#
