## Wizard scripts
## version 1.02 December 2015

from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *

class QuantumWizard1(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'
		
	def PyInit(self):
		PeriodicTable=[
		'','Mu','e',
		'H','D','3He',
		'6Li','7Li','Li','Be','10B','11B','13C','14N','N','15N','17O','F','21Ne',
		'Na','25Mg','Al','29Si','P','33S','35Cl','37Cl',
		'39K','K','40K','41K','43Ca','Sc','47Ti','49Ti','50V','51V','V','53Cr','Mn','57Fe','Co','61Ni','63Cu','65Cu','67Zn',
		'69Ga','71Ga','73Ge','As','77Se','79Br','81Br','83Kr',
		'85Rb','87Rb','87Sr','Y','91Zr','Nb','95Mo','97Mo','99Ru','101Ru','Rh','105Pd','107Ag','109Ag','111Cd','113Cd',
		'113In','115In','In','115Sn','117Sn','119Sn','121Sb','123Sb','123Te','125Te','I','129Xe','131Xe',
		'Cs','135Ba','137Ba','138La','139La','La',
		'Pr','143Nd','145Nd','147Sm','149Sm','151Eu','153Eu','155Gd','157Gd','Tb','161Dy','163Dy','Ho','167Er','Tm','171Yb','173Yb','175Lu','Lu','176Lu',
		'177Hf','179Hf','180Ta','181Ta','Ta','183W','185Re','187Re','187Os','189Os','191Ir','193Ir','195Pt','Au','199Hg','201Hg','203Tl','205Tl','207Pb','Bi','235U' ]
		self.declareProperty('Spin0','Mu',StringListValidator(PeriodicTable))
		self.declareProperty('Spin1','e',StringListValidator(PeriodicTable))
		self.declareProperty('Spin2','',StringListValidator(PeriodicTable))
		self.declareProperty('Spin3','',StringListValidator(PeriodicTable))
		self.declareProperty('Spin4','',StringListValidator(PeriodicTable))
		self.declareProperty('Spin5','',StringListValidator(PeriodicTable))
		self.declareProperty('HyperfineOptions','To One Electron Only',StringListValidator(['None (nuclear)','To One Electron Only','All pairs']))

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard1)

class QuantumWizard2(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'
		
	def PyInit(self):
		self.declareProperty('CalculationType','Plain',StringListValidator(['Plain','RF resonance','Relaxation and Conversion']))
		self.declareProperty('States_for_conversion',1,IntBoundedValidator(1,6))

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard2)

class QuantumWizard3(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('Vary_position_between_sites',True)
		self.declareProperty('Vary_hyperfine_between_sites',True)
		self.declareProperty('Vary_quadrupole_between_sites',True)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard3)

class QuantumWizard4(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('Use_dipolar',False)
		self.declareProperty('x',0.0)
		self.declareProperty('y',0.0)
		self.declareProperty('z',0.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard4)

class QuantumWizard5(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('A',0.0,doc='Isotropic part')
		self.declareProperty('D',0.0,doc='Axial anisotropy')
		self.declareProperty('AxisX',0.0,doc='Symmetry axis of axial HFC')
		self.declareProperty('AxisY',0.0)
		self.declareProperty('AxisZ',1.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard5)

class QuantumWizard6(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('nuQ',0.0,doc='Axial hyperfine splitting')
		self.declareProperty('AxisX',0.0,doc='Symmetry axis of axial EFG')
		self.declareProperty('AxisY',0.0)
		self.declareProperty('AxisZ',1.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard6)

class QuantumWizard7(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('RelaxSpin0',0.0)
		self.declareProperty('RelaxSpin1',0.0)
		self.declareProperty('RelaxSpin2',0.0)
		self.declareProperty('RelaxSpin3',0.0)
		self.declareProperty('RelaxSpin4',0.0)
		self.declareProperty('RelaxSpin5',0.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard7)

class QuantumWizard8(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('ConversionOption','AllEqual',StringListValidator(['AllEqual','TwoLevels']))
		self.declareProperty('BetweenAllStates',1.0)
		self.declareProperty('From_0_To_1',1.0)
		self.declareProperty('From_1_To_0',1.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard8)

class QuantumWizard9(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('PopulationState0',0.1)
		self.declareProperty('PopulationState1',0.1)
		self.declareProperty('PopulationState2',0.1)
		self.declareProperty('PopulationState3',0.1)
		self.declareProperty('PopulationState4',0.1)
		self.declareProperty('PopulationState5',0.1)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard9)

class QuantumWizard10(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('RF_Field_strength_G',1.0)
		self.declareProperty('RF_Frequency_MHz',20.0)
		self.declareProperty('RF_Phase_Deg',0.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard10)

class QuantumWizard11(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('Averaging','LFuniform',StringListValidator(['LFuniform','LFrandom','TFuniform','TFrandom','Fixed axis LF']))
		self.declareProperty('NumberOfAverages',100)
		self.declareProperty('AxisX',0.0)
		self.declareProperty('AxisY',0.0)
		self.declareProperty('AxisZ',1.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard11)

class QuantumWizard12(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('MagneticFieldTesla',0.0)

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard12)

class QuantumWizard13(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('MeasurementType','Time Spectra',StringListValidator(['Time Spectra','Integral Asymmetry']))

	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard13)

class QuantumWizard14(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum\Wizard'

	def PyInit(self):
		self.declareProperty('LoopVariable','Magnetic Field',StringListValidator(['<none>','Magnetic Field','Relaxation rate of Electron','Conversion rates','RF Frequency','Muon isotropic hyperfine constant']))
		self.declareProperty('StartValue',0.0)
		self.declareProperty('EndValue',1.0)
		self.declareProperty('NumberOfPoints',100,IntBoundedValidator(lower=2))
		self.declareProperty('LogScale',False)
		
	def PyExec(self):
		pass

AlgorithmFactory.subscribe(QuantumWizard14)

		