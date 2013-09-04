from mantid.api import PythonAlgorithm, AlgorithmFactory
from mantid.kernel import StringListValidator, StringMandatoryValidator
from qtiGenie import *
class c2et(PythonAlgorithm):

	def PyInit(self):
		self.declareProperty("instrument","Full",StringListValidator(["mar","maps"]))
		self.declareProperty("Runfile","", validator = StringMandatoryValidator(), doc = "runfile")
		self.declareProperty("dVan", "", validator = StringMandatoryValidator(), doc = "Detector Vanadium")
		#self.declareWorkspaceProperty("Run file","", Direction = Direction.Input, doc = "A value for the end of the range(inclusive)")
		#self.declareWorkspaceProperty("Whitebeam", "", Direction = Direction.Input, doc = "Optional preamble")
		self.declareProperty("Ei", "",StringMandatoryValidator(),doc='Ei guess')
		self.declareProperty("Rebin","", StringMandatoryValidator(),doc='Rebin string')
		self.declareProperty("keywords","",doc='Key words')
		self.declareProperty("mapfile", "", StringMandatoryValidator(), doc = 'mapfile')
		
	def PyExec(self):
		run = self.getProperty("Runfile")
		dVan = self.getProperty("dVan")
		Ei = float(self.getProperty("Ei"))
		rebin = str(self.getProperty("Rebin"))
		print rebin
		mapfile=self.getProperty("mapfile")
		print mapfile
		inst=self.getProperty("instrument")
		kw=self.getProperty("keywords")
		
		kw=kw.split(',')
		iliad_setup(inst)
		dd={}
		for item in kw:
			key,value = item.split('=')
			dd[key] = value
		
		iliad(dVan,run,Ei,rebin,mapfile,det_cal_file=wb,**dd)
		
		
#############################################################################################

# Register algorithm with Mantid
AlgorithmFactory.subscribe(c2et)
