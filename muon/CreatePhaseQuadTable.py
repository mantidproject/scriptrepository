## CreatePhaseQuadTable - make a blank Mantid table for Phase Quad algorithm
## Author: James Lord
## Version 2.0, September 2015
## now sets Efficiency column (to expected asymmetry) and sets dead time=0
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
import math

class CreatePhaseQuadTable(PythonAlgorithm):
	def category(self):
		return 'Muon'

	def PyInit(self):
		self.declareProperty("Instrument","MUSR")
		self.declareProperty("FieldAxis","z",StringListValidator(["x","y","z"]))
		self.declareProperty(ITableWorkspaceProperty('TableName','Table1',Direction.Output))
		
	def PyExec(self):

		inst=self.getProperty("Instrument").value

		cdws=self.createChildAlgorithm("CreateSimulationWorkspace")
		cdws.setProperty("Instrument",inst)
		cdws.setProperty("BinParams","0,1,32")
		cdws.execute()
		dws=cdws.getProperty("OutputWorkspace").value
		#dws=CreateSimulationWorkspace(inst,"0,1,32")
		
		nr=dws.getNumberHistograms()

		cetw=self.createChildAlgorithm("CreateEmptyTableWorkspace")
		cetw.execute()
		Tab3=cetw.getProperty("OutputWorkspace").value
		#Tab2=CreateEmptyTableWorkspace()
		Tab3.addColumn("bool","ok")
		Tab3.addColumn("double","asym")
		Tab3.addColumn("double","phase")
		Tab3.addColumn("double","dead")

		ax=self.getProperty("FieldAxis").value

		for i in range(nr):
			det=dws.getDetector(i).getPos()-dws.getInstrument().getSample().getPos()
			r=math.sqrt(det.X()**2+det.Y()**2+det.Z()**2)
			if(ax=="x"):
				phi=math.atan2(det.Z(),det.Y())
				ampl=math.sqrt(det.Z()**2+det.Y()**2)/r
			elif(ax=="y"):
				phi=math.atan2(det.X(),det.Z())
				ampl=math.sqrt(det.X()**2+det.Z()**2)/r
			else: # z
				phi=math.atan2(det.Y(),det.X())
				ampl=math.sqrt(det.Y()**2+det.X()**2)/r
			Tab3.addRow([True,ampl,phi,0.0])

		self.setProperty("TableName",Tab3)
		del Tab3
		del dws

AlgorithmFactory.subscribe(CreatePhaseQuadTable)
