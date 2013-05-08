from MantidFramework import *
from qtiGenie import *
class c2et(PythonAlgorithm):

	def PyInit(self):
		self.declareProperty("instrument","Full",ListValidator(["mar","maps"]))
		self.declareProperty("Runfile","", Validator = MandatoryValidator(), Description = "runfile")
		self.declareProperty("dVan", "", Validator = MandatoryValidator(), Description = "Detector Vanadium")
		#self.declareWorkspaceProperty("Run file","", Direction = Direction.Input, Description = "A value for the end of the range(inclusive)")
		#self.declareWorkspaceProperty("Whitebeam", "", Direction = Direction.Input, Description = "Optional preamble")
		self.declareProperty("Ei", "",MandatoryValidator(),Description='Ei guess')
		self.declareProperty("Rebin","", MandatoryValidator(),Description='Rebin string')
		self.declareProperty("keywords","",Description='Key words')
		self.declareProperty("mapfile", "", MandatoryValidator(), Description = 'mapfile')
		
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
mantid.registerPyAlgorithm(c2et())
