from mantid.simpleapi import *
from mantid.api import PythonAlgorithm, AlgorithmFactory
from mantid.kernel import StringListValidator, StringMandatoryValidator
from mantid import config
import os.path, math

class IndirectTransmission(PythonAlgorithm):
 
	def category(self):
		return "Workflow\\MIDAS;PythonAlgorithms"

	def PyInit(self):
		self.declareProperty(name='Instrument',defaultValue='IRIS',validator=StringListValidator(['IRIS','OSIRIS']), doc='Instrument')
		self.declareProperty(name='Analyser',defaultValue='graphite',validator=StringListValidator(['graphite','fmica']), doc='Analyser')
		self.declareProperty(name='Reflection',defaultValue='002',validator=StringListValidator(['002','004']), doc='Reflection')
		self.declareProperty(name='Chemical Formula',defaultValue='',validator=StringMandatoryValidator(), doc='Sample chemical formula')
		self.declareProperty(name='Number Density', defaultValue=0.1, doc='Number denisty. Default=0.1')
		self.declareProperty(name='Thickness', defaultValue=0.1, doc='Sample thickness. Default=0.1')
 
	def PyExec(self):

		self.log().information('Indirect transmission')
		instr = self.getPropertyValue('Instrument')
		ana = self.getPropertyValue('Analyser')
		refl = self.getPropertyValue('Reflection')
		formula = self.getPropertyValue('Chemical Formula')
		dens = self.getPropertyValue('Number Density')
		thick = self.getPropertyValue('Thickness')
		
		idf_dir = config['instrumentDefinition.directory']
		idf = idf_dir + instr + '_Definition.xml'
		mtWS = '__empty_'+instr
		LoadEmptyInstrument(OutputWorkspace=mtWS, Filename=idf)
		ipf = idf_dir + instr + '_' + ana + '_' + refl + '_Parameters.xml'
		LoadParameterFile(Workspace=mtWS, Filename=ipf)
		ins = mtd[mtWS]
		efixed = ins.getInstrument().getNumberParameter('efixed-val')[0]
		wave=1.8*math.sqrt(25.2429/efixed)
		logger.notice('Analyser : ' +ana+refl +' with energy = ' + str(efixed))
		result = SetSampleMaterial(InputWorkspace=mtWS,ChemicalFormula=formula)
		mat = ins.sample().getMaterial()
		absXS = result[5]*wave/1.7982
		logger.notice('Absorption Xsection at wavelength ' + str(wave) +
			' A = '+str(absXS))
		cohXS = result[4]
		logger.notice('Coherent Xsection = ' + str(cohXS))
		incXS = result[3]
		logger.notice('Incoherent Xsection = ' + str(incXS))
		scatXS = result[3]
		logger.notice('Total scattering Xsection = ' + str(scatXS))
		totalXS = absXS + scatXS
		logger.notice('Number density = ' + dens)
		dens = float(dens)
		logger.notice('Thickness = ' + thick)
		thick = float(thick)
		trans = math.exp(-dens*totalXS*thick)
		scatter = 1.0 - math.exp(-dens*scatXS*thick)
		logger.notice('Transmission (abs+scatt) = ' +str(trans))
		logger.notice('Total scattering = ' +str(scatter))

AlgorithmFactory.subscribe(IndirectTransmission)         # Register algorithm with Mantid
