"""

This module contains the datatypes used to analysing semsans runs

"""
import datetime
from collections import namedtuple

HeData = namedtuple("HeData", "run cell scale dt fid t1")
RunData = namedtuple("RunData", "number sample start end "
                     "trans csans ctrans direct")
QuickData = namedtuple("QuickData", "number sample start end duration charge")


def table_to_run(table):
    """
    Turns a mantid table of sample runs into a list of RunData.  The table
    should have been generated get get_log to ensure that it is in the correct
    format.

    Parameters
    ----------
    table:
      A mantid table object containing the runs of interest

    Returns
    -------
    A list of RunData objects
    """
    return [
        RunData(x["Run Number"], x["Sample"],
                datetime.datetime.strptime(x["Start time"],
                                            "%Y-%m-%dT%H:%M:%S"),
                datetime.datetime.strptime(x["End time"],
                                           "%Y-%m-%dT%H:%M:%S"),
                x["Trans run"], x["Can Sans run"],
                x["Can Trans run"], x["Direct Trans run"])
        for x in table]
