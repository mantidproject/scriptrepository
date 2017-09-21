# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from runtypes import HeData
import numpy as np
import datetime


class HeliumPolarisation(PythonAlgorithm):
    def category(self):
        return "SESANS"

    def PyInit(self):
        self.declareProperty(
            WorkspaceProperty(name="HeliumLogTable", defaultValue="",
                              direction=Direction.Input))
        self.declareProperty(
            WorkspaceProperty(name="LogTable", defaultValue="",
                              direction=Direction.Input))
        self.declareProperty(
            WorkspaceProperty(name="BaseLine", defaultValue="",
                              direction=Direction.Input))
        self.declareProperty(name="RunNumber", defaultValue=0)
        self.declareProperty(name="SimSteps", defaultValue=10)
        self.declareProperty(
            WorkspaceProperty(name="OutputWorkspace", defaultValue="",
                              direction=Direction.Output))



    @staticmethod
    def _he3pol(scale, time):
        """
        Create a ³He polarisation compensator

        The polarisation is given by

        P = tanh(scale × exp(-time) × λ)

        Parameters
        ----------
        scale : float
          The initial constant for the he3 polarisation at time=0
        time : float
          The amount of time that has passed since polarisation
          in units of the cell's time constants

        Returns
        -------
        A function which takes a wavelength and gives a
        ³He polarising efficiency.
        """
        def pol(wavelength):
            """
            Calculate the polarisation efficiency as a given wavelength

            Parameters
            ----------
            wavelength
              The wavelength of the neutron in Å

            Return
            ------
            The polarising efficiency of the ³He analyser
            """
            return np.tanh(scale * np.exp(-time) * wavelength)
        return pol


    @staticmethod
    def _he3_stats(timecode, helog, table):
        """Return the information about the ³He analyser's state during this
        run.  This function requires that get_he3_log has already been
        run, creating the "helium_log" table.

        Parameters
        ----------
        run
          A RunData object containing the run whose analyser
          statistics are needed

        Return
        ------
        A HeData object describing the state of the analyser during the run.

        """
        stats = [x for x in helog
                 if timecode > x["Start time"]][-1]
        return HeData(stats["Number"], stats["Cell"], stats["scale"],
                      datetime.datetime.strptime(stats["Start time"],
                                                 "%Y-%m-%dT%H:%M:%S"),
                      stats["fid"], stats["Time Constant"])

    def PyExec(self):
        helog = self.getProperty("HeliumLogTable").value
        runs = self.getProperty("LogTable").value
        run = self.getProperty("RunNumber").value
        base = self.getProperty("BaseLine").value
        nbins = base.blocksize()
        dx = base.extractX()[0]
        dx = (dx[1:]+dx[:-1])/2

        runline = [x for x in runs if x["Run Number"] == run][0]
        start = datetime.datetime.strptime(runline["Start time"],
            "%Y-%m-%dT%H:%M:%S")
        end = datetime.datetime.strptime(runline["End time"],
            "%Y-%m-%dT%H:%M:%S")
        
        he_stat = self._he3_stats(runline["Start time"], helog, runs)
        start = (start-he_stat.dt).seconds/3600/he_stat.t1
        end = (end-he_stat.dt).seconds/3600/he_stat.t1
        pols = []
        for time in np.linspace(start, end, 10):
            temp = self._he3pol(he_stat.scale, time)(dx)
            self.log().notice(repr(temp))
            pols.append(temp)
        wpol = WorkspaceFactory.create(base, NVectors=1,
                                       XLength=nbins + 1,
                                       YLength=nbins)
        wpol.setY(0, np.mean(pols, axis=0))
        wpol.setX(0, base.dataX(0))

        self.setProperty("OutputWorkspace", wpol)

AlgorithmFactory.subscribe(HeliumPolarisation)
