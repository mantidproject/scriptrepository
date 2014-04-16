# Algorithm to start IN16reduce
from mantid.simpleapi import *
from mantid.kernel import StringListValidator, StringMandatoryValidator
from mantid.api import PythonAlgorithm, AlgorithmFactory
from mantid import config, logger, mtd
import sys, os.path, math

class IndirectILLreduction(PythonAlgorithm):
 
	def category(self):
		return "Workflow\\MIDAS;PythonAlgorithms"

	def PyInit(self):
		self.declareProperty(name='Analyser', defaultValue='silicon', validator=StringListValidator(['silicon']), doc = 'Analyser crystal')
		self.declareProperty(name='Reflection', defaultValue='111', validator=StringListValidator(['111']), doc = 'Analyser reflection')
		self.declareProperty(name='RunName', defaultValue='', validator=StringMandatoryValidator(), doc = 'Run name (after <Instr>_)')
		self.declareProperty(name='Map', defaultValue='default', validator=StringListValidator(['default','user']), doc = 'Use detector map')
		self.declareProperty(name='Mirror', defaultValue=False, doc = 'Mirror mode')
		self.declareProperty(name='Verbose', defaultValue=True, doc = 'Switch Verbose Off/On')
		self.declareProperty(name='Save', defaultValue=False, doc = 'Switch Save result to nxs file Off/On')
		self.declareProperty(name='Plot', defaultValue=False, doc = 'Plot options')
 
	def PyExec(self):

		self.log().information('IndirectILLreduction input')
		workdir = config['defaultsave.directory']
		instr = 'IN16B'
		ana = self.getPropertyValue('Analyser')
		refl = self.getPropertyValue('Reflection')
		run = self.getPropertyValue('RunName')
		map = self.getProperty('Map').value
		mirror = self.getProperty('Mirror').value
		Verbose = self.getProperty('Verbose').value
		Save = self.getProperty('Save').value
		Plot = self.getProperty('Plot').value

		run_name = instr+'_'+run
		raw_name = run_name+'_'+ana+refl+'_raw'
		run_path = os.path.join(workdir, run+'.nxs')		# path name for sample nxs file
		LoadILLIndirect(FileName=run_path, OutputWorkspace=raw_name)
		AddSampleLog(Workspace=raw_name, LogName="mirror_sense", LogType="String", LogText=str(mirror))
		if Verbose:
			logger.notice('Nxs file : '+run_path)
		from IndirectNeutron import ILLindirect
		ILLindirect(raw_name,instr,ana,refl,map,Verbose,Plot,Save)
 
AlgorithmFactory.subscribe(IndirectILLreduction)      # Register algorithm with Mantid
