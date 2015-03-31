## CreatePhaseQuadTable - make a blank Mantid table for Phase Quad algorithm
## Author: James Lord
## Version 1.0, January 2015
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *

class CreatePhaseQuadTable(PythonAlgorithm):
	def category(self):
		return 'Muon'

	def PyInit(self):
		self.declareProperty("Instrument","MUSR")
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

		for i in range(nr):
			phi=dws.getDetector(i).getPhi()
			Tab3.addRow([True,0.2,phi,0.01])

		self.setProperty("TableName",Tab3)
		del Tab3
		del dws

AlgorithmFactory.subscribe(CreatePhaseQuadTable)
