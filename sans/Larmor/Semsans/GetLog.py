# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
import datetime
import numpy as np
import re
import xml.etree.cElementTree
from .runtypes import QuickData, RunData

class GetLog(PythonAlgorithm):

    JPATH = r'\\isis\inst$\NDXLARMOR\Instrument\logs\journal'

    def category(self):
        return "SESANS"

    def PyInit(self):
        self.declareProperty("start", defaultValue=0,
                             doc="Starting run number")
        self.declareProperty("stop", defaultValue=0,
                             doc="Ending run number")
        self.declareProperty("run", defaultValue=r'(.+)_run:_(\d+)_SANS',
                             doc="Ending run number")
        self.declareProperty("trans", defaultValue=r".+_TRANS",
                             doc="Ending run number")
        self.declareProperty("can_sans", defaultValue=r"D2._2mm_SANS",
                             doc="Ending run number")
        self.declareProperty("can_trans", defaultValue=r"D2._2mm_TRANS",
                             doc="Ending run number")
        self.declareProperty("direct_sans", defaultValue="MT_Beam_SANS",
                             doc="Ending run number")
        self.declareProperty("direct_trans", defaultValue="MT_Beam_TRANS",
                             doc="Ending run number")
        self.declareProperty(WorkspaceProperty(name="LogTable",
                                               direction=Direction.Output,
                                               defaultValue="LogTable"))
        self.declareProperty(WorkspaceProperty(name="SamplesTable",
                                               direction=Direction.Output,
                                               defaultValue="SampleTable"))


    def PyExec(self):
        start = self.getProperty("start").value
        stop = self.getProperty("stop").value
        self.run = self.getProperty("run").value
        self.trans = self.getProperty("trans").value
        self.can_sans = self.getProperty("can_sans").value
        self.can_trans = self.getProperty("can_trans").value
        self.direct_sans = self.getProperty("direct_sans").value
        self.direct_trans = self.getProperty("direct_trans").value
        runs = np.arange(start, stop+1)

        log_file = self.JPATH + "\\" + self._get_relevant_log(min(runs))
        results = []
        with open(log_file, "r") as infile:
            journal = xml.etree.cElementTree.iterparse(infile)
            for _, child in journal:
                if "NXentry" in child.tag:
                    num = self._get_xml_run_number(child)
                    if num in runs:
                        for param in child:
                            if "title" in param.tag:
                                sample = param.text.replace(" ", "_")
                            elif "start_time" in param.tag:
                                start = datetime.datetime.strptime(
                                    param.text,
                                    "%Y-%m-%dT%H:%M:%S")
                            elif "end_time" in param.tag:
                                stop = datetime.datetime.strptime(
                                    param.text,
                                    "%Y-%m-%dT%H:%M:%S")
                            elif "duration" in param.tag:
                                duration = datetime.timedelta(
                                    seconds=int(param.text))
                            elif "proton_charge" in param.tag:
                                proton_charge = float(param.text)
                        results.append(
                            QuickData(num, sample, start, stop, duration,
                                      proton_charge))
                    child.clear()
                    if num > max(runs):
                        break
        trans = [run for run in results
                 if re.match(self.trans, run[1])]
        csans = [run for run in results
                 if re.match(self.can_sans, run[1])]
        ctrans = [run for run in results
                  if re.match(self.can_trans, run[1])]
        dtrans = [run for run in results
                  if re.match(self.direct_trans, run[1])]
        temp = [self.convert_run(run, trans, csans, ctrans, dtrans)
                for run in results
                if (re.match(self.run, run.sample) or
                    re.match(self.can_sans, run.sample) or
                    re.match(self.direct_sans, run.sample))
                and run.charge/run.duration.seconds > 0.005]

        d = {}
        for run in temp:
            if run.sample in list(d.keys()):
                d[run.sample].append(run)
            else:
                d[run.sample] = [run]

        my_table = WorkspaceFactory.createTable()
        my_table.addColumn("int", "Run Number")
        my_table.addColumn("str", "Sample")
        my_table.addColumn("str", "Start time")
        my_table.addColumn("str", "End time")
        my_table.addColumn("int", "Trans run")
        my_table.addColumn("int", "Can Sans run")
        my_table.addColumn("int", "Can Trans run")
        my_table.addColumn("int", "Direct Trans run")

        for k, v in list(d.items()):
            for run in v:
                my_table.addRow(
                    [run.number, run.sample,
                     run.start.isoformat(),
                     run.end.isoformat(),
                     run.trans, run.csans, run.ctrans,
                     run.direct])
        self.setProperty("LogTable", my_table)

        my_table = WorkspaceFactory.createTable()
        my_table.addColumn("str", "Name")
        for k, _ in list(d.items()):
            my_table.addRow([k])
        self.setProperty("SamplesTable", my_table)

    def _get_relevant_log(self, run):
        """

        Find the correct journal log for the run in question.  This allows
        our code to work over multiple run cycles.

        Parameters
        ----------
        run
          The run number being analysed

        Returns
        -------
        A string containing the path to the journal file containing that run.

        """
        with open(self.JPATH+r"\journal_main.xml", "r") as infile:
            journals = xml.etree.cElementTree.parse(infile)
        for child in journals.getroot():
            if int(child.attrib["first_run"]) <= run and \
               run <= int(child.attrib["last_run"]):
                return child.attrib["name"]

    @staticmethod
    def _get_xml_run_number(node):
        """

        Find the run number of an xml node

        Parameters
        ----------
        node
        The xml node for the run

        Returns
        -------
        The run number
        """
        for child in node:
            if "run_number" in child.tag:
                return int(child.text)

    def convert_run(self, run, trans, csans, ctrans, dtrans):
        """
        Fully populate the data for a run.

        Parameters
        ----------
        run
          A QuickData object about the run to be converted
        trans
          A list QuickData objects for the sample transmission runs
        csans
          A list of QuickData objects for the can sans runs
        ctrans
          A list of QuickData objects for the can transmission runs
        dtrans
          A list of QuickData objects for the direct beam measurements

        Returns
        -------
        A fully populated RunData object for the requested run.
        """
        if re.match(self.run, run.sample):
            sample = re.match(self.run, run.sample).group(1)
        else:
            sample = "Full_Blank"

        tr = -1
        for tran in trans:
            if sample in tran.sample:
                tr = tran.number
        csans.sort(key=lambda x: x.start-run.start)
        ctrans.sort(key=lambda x: x.start-run.start)
        dtrans.sort(key=lambda x: x.start-run.start)
        return RunData(run.number, sample, run.start, run.end, tr,
                       int(csans[0].number), int(ctrans[0].number),
                       int(dtrans[0].number))

AlgorithmFactory.subscribe(GetLog)
