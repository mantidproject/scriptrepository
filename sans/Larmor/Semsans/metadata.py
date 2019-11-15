# -*- coding: utf-8 -*-
import re
import datetime
import xml.etree
from mantid.simpleapi import CreateEmptyTableWorkspace, RenameWorkspace
from .runtypes import HeData, RunData, QuickData

RUN_IDENTIFIERS = {
    "run": r'(.+) run: (\d+)_SANS',
    "trans": r".+_TRANS",
    "can_sans": r"D2. 2mm_SANS",
    "can_trans": r"D2. 2mm_TRANS",
    "direct_sans": "MT Beam_SANS",
    "direct_trans": "MT Beam_TRANS"}


def convert_he(i):
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
    run = int(i[0])
    cell = i[1]
    pl = float(i[2])
    phe = float(i[3]) / 100.0
    dt = datetime.datetime.strptime(i[4] + " " + i[5], "%m/%d/%Y %H:%M")
    fid = float(i[6])
    t1 = float(i[7])

    return HeData(run, cell, 0.0733 * pl * phe, dt, fid, t1)


def load_helium_file(helium_file):
    """

    Turn the ³He log file into a list of HeData

    Parameters
    ----------
    helium_file
      The path to the ³He log file

    Returns
    -------
    A list of HeData
    """
    with open(helium_file, "r") as infile:
        infile.readline()  # read header
        return [convert_he(line.split("\t"))
                for line in infile]


def convert_run(run, trans, csans, ctrans, dtrans):
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
    if re.match(RUN_IDENTIFIERS["run"], run.sample):
        sample = re.match(RUN_IDENTIFIERS["run"], run.sample).group(1)
    else:
        sample = "Full Blank"

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


JPATH = r'\\isis\inst$\NDXLARMOR\Instrument\logs\journal'


def get_relevant_log(run):
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
    with open(JPATH+r"\journal_main.xml", "r") as infile:
        journals = xml.etree.ElementTree.parse(infile)
    for child in journals.getroot():
        if int(child.attrib["first_run"]) <= run and \
           run <= int(child.attrib["last_run"]):
            return child.attrib["name"]


def get_xml_run_number(node):
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


def get_he3_log(path):
    """
    Load the ³He log data into Mantid Table named "helium_log"

    Parameters
    ----------
    path
      A string with the path to the ³He as a tsv file
    """
    hetemp = load_helium_file(path)
    my_table = CreateEmptyTableWorkspace()
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
    RenameWorkspace(my_table, "helium_log")


def get_log(runs):
    """
    Uses the run journal to identify which run numbers are associated with
    which samples and create a table for each sample, containing all the
    information needed to analyse each run.

    Parameters
    ----------
    runs
      A list of integer run numbers
    """
    log_file = JPATH + "\\" + get_relevant_log(min(runs))
    results = []
    with open(log_file, "r") as infile:
        journal = xml.etree.cElementTree.iterparse(infile)
        for _, child in journal:
            if "NXentry" in child.tag:
                num = get_xml_run_number(child)
                if num in runs:
                    for param in child:
                        if "title" in param.tag:
                            sample = param.text
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
             if re.match(RUN_IDENTIFIERS["trans"], run[1])]
    csans = [run for run in results
             if re.match(RUN_IDENTIFIERS["can_sans"], run[1])]
    ctrans = [run for run in results
              if re.match(RUN_IDENTIFIERS["can_trans"], run[1])]
    dtrans = [run for run in results
              if re.match(RUN_IDENTIFIERS["direct_trans"], run[1])]
    temp = [convert_run(run, trans, csans, ctrans, dtrans)
            for run in results
            if (re.match(RUN_IDENTIFIERS["run"], run.sample) or
                re.match(RUN_IDENTIFIERS["can_sans"], run.sample) or
                re.match(RUN_IDENTIFIERS["direct_sans"], run.sample))
            and run.charge/run.duration.seconds > 0.005]

    d = {}
    for run in temp:
        if run.sample in list(d.keys()):
            d[run.sample].append(run)
        else:
            d[run.sample] = [run]

    for k, v in list(d.items()):
        my_table = CreateEmptyTableWorkspace()
        my_table.addColumn("int", "Run Number")
        my_table.addColumn("str", "Sample")
        my_table.addColumn("str", "Start time")
        my_table.addColumn("str", "End time")
        my_table.addColumn("int", "Trans run")
        my_table.addColumn("int", "Can Sans run")
        my_table.addColumn("int", "Can Trans run")
        my_table.addColumn("int", "Direct Trans run")

        for run in v:
            my_table.addRow(
                [run.number, run.sample,
                 run.start.isoformat(),
                 run.end.isoformat(),
                 run.trans, run.csans, run.ctrans,
                 run.direct])
        RenameWorkspace(my_table, k+"_runs")
