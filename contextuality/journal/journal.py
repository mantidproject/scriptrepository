from lxml import etree
import collections
import logging
import numpy as np

from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *

SCHEMA = '{http://definition.nexusformat.org/schema/3.0}'
BASE_PATH = r"\\isis\inst$\NDX{}\Instrument\logs\journal\journal_{}.xml"
LOG_PATH = r"\\isis\inst$\NDX{inst:}\Instrument\data\cycle_{cyc:}\{inst:}{num:08d}.log"


def load_journal(instrument, cycle):
    """Load the instrument journal"""
    with open(BASE_PATH.format(instrument, cycle), "r") as infile:
        return etree.parse(infile)


def safe_float(x):
    """Turn a string into a float if possible"""
    try:
        return float(x.strip())
    except ValueError:
        return x


def load_log(instrument, cycle, number):
    """Get the motor log for a run"""
    with open(LOG_PATH.format(inst=instrument, cyc=cycle, num=number),
              "r") as infile:
        data = [x.split("\t") for x in infile.readlines()]
    d = collections.defaultdict(list)
    for x in data:
        d[x[1]].append({"time": x[0], "value": safe_float(x[2])})
    return d


def meets(instrument, cycle, start, stop, motors, values, prog=None):
    results = []
    for i in range(start, stop):
        prog.report("Filtering Runs")
        d = load_log(instrument, cycle, i)
        for m, v in zip(motors, values):
            if not (np.abs(d[m] and d[m][-1]["value"] - float(v)) < 1e-3):
                break
        else:
            results.append(i)
    return results
