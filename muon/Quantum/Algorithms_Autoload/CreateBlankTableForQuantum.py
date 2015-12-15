## CreateBlankTableForQuantum - make a blank Mantid table with 2 text columns and about 30 rows
## Author: James Lord
## Version 1.02, December 2015
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *

class CreateBlankTableForQuantum(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum'

	def PyInit(self):
		self.declareProperty("BlankRows",30)
		self.declareProperty(ITableWorkspaceProperty('TableName','Table1',Direction.Output))
		
	def PyExec(self):

		nr=self.getProperty("BlankRows").value

		cetw=self.createChildAlgorithm("CreateEmptyTableWorkspace")
		cetw.execute()
		Tab3=cetw.getProperty("OutputWorkspace").value
		#Tab2=CreateEmptyTableWorkspace()
		Tab3.addColumn("str","Code")
		Tab3.addColumn("str","Value")

		for i in range(nr):
			Tab3.addRow(["",""])

		self.setProperty("TableName",Tab3)
		del Tab3

AlgorithmFactory.subscribe(CreateBlankTableForQuantum)
