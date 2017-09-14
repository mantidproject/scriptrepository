# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from runtypes import HeData

class HeliumLog(PythonAlgorithm):
    def category(self):
        return "SESANS"

    def PyInit(self):
        self.declareProperty(FileProperty(name="LogFile", defaultValue="",
                                          action=FileAction.Load,
                                          extensions=["tsv"]))
        self.declareProperty(
            ITableWorkspaceProperty(name="OutputWorkspace",
                                    defaultValue="",
                                    direction=Direction.Output))

    def _load_helium_file(self, helium_file):
        with open(helium_file, "r") as infile:
            infile.readline()  # read header
            return [self._convert_he(line.split("\t"))
                    for line in infile]

    @staticmethod
    def _convert_he(i):
        """
        Convert a line from the ³He log into an HeData object

        Parameters
        ----------
        i
        A list of strings, representing each cell in the row of the
        ³He log file

        Returns
        -------
        A HeData object containing the information from this row of the log
        """
        import datetime
        run = int(i[0])
        cell = i[1]
        pl = float(i[2])
        phe = float(i[3]) / 100.0
        dt = datetime.datetime.strptime(i[4] + " " + i[5], "%m/%d/%Y %H:%M")
        fid = float(i[6])
        t1 = float(i[7])

        return HeData(run, cell, 0.0733 * pl * phe, dt, fid, t1)

    def PyExec(self):
        path = self.getProperty("LogFile").value
        hetemp = self._load_helium_file(path)
        my_table = WorkspaceFactory.createTable()
        my_table.addColumn("int", "Number")
        my_table.addColumn("str", "Cell")
        my_table.addColumn("float", "scale")
        my_table.addColumn("str", "Start time")
        my_table.addColumn("float", "fid")
        my_table.addColumn("float", "Time Constant")

        for run in hetemp:
            my_table.addRow([run.run, run.cell,
                             run.scale,
                             run.dt.isoformat(),
                             run.fid, run.t1])
        self.setProperty("OutputWorkspace", my_table)

AlgorithmFactory.subscribe(HeliumLog)
