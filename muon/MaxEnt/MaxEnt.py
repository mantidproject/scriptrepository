""" A thin Mantid algorithm wrapper around an F2Py-compiled version of the Fortran
"MaxEnt" program, which is currently being used from within OpenGENIE.

The existing functionality of the OpenGENIE version was to specify a run number
and several options, and have the program do its calculations and plot the result.
There were also three files given to the program containing groupings as well as
phase and deadtime data.  These were updated after the program had run its course.

We try to do the same here, except that we can work out a few of the options
ourselves without having to ask the user for them.  We also supply the plotted
data as an OutputWorkspace for users to play with.

Another difference is that in most cases we are removing the responsiblity of
the Fortran to write to and read from files; doing this via Python is easier
and a lot more robust (doesn't crash if it can't get access rights).

We're supplying this via the Script Repository to keep F2Py out of the main one.
The main disadvantage to this is that we will have to run the unit tests manually
after every change to make sure things haven't regressed.  Another potential
downside is that this script requires users to add the containing directory to
their "Python extensions" in Mantid preferences.  This will mean that the .pyd
will be used (locked) by the Mantid process and so updating it via the script
repository will not be possible without first removing the directory again from
the preferences.

Example usage (SUBJECT TO CHANGE WITH TESTING AND FEEDBACK FROM USERS):

# Execute with default input files, for the run "48033" on "MUSR".
result = MaxEnt(RunNumber="48033", Instrument="MUSR")
"""

from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import (LoadMuonNexus, DeleteWorkspace,
                              LoadDetectorsGroupingFile, CreateWorkspace,
                              GroupWorkspaces, CreateEmptyTableWorkspace)

import numpy as np
import os
import string
import random
import sys
import itertools

# Nasty variables needed to match the magic numbers found in the Fortran or in
# the original OpenGENIE version of MaxEnt.
MAX_HISTOS = 96 # (Originally 64 but changed to 96 to support EMU.)
MAX_INPUT_DATA_SIZE = 262144
DEFAULT_LEVEL_VALUE = 0.1
DEFAULT_SIGMA_LOOSENESS_VALUE = 1.02
FRAMES = 12894
RES = 16000
POINTS_TO_FIT = 4096

# Add this script's directory to the path so that we may find the necessary
# files/module.
path = os.path.dirname(os.path.realpath(__file__))
if not path in sys.path:
    sys.path.insert(0,path)

# If we're on a supported platform, import the Fortran module needed for MaxEnt.
import platform
os_env = platform.system() + platform.architecture()[0]
if os_env == "Windows64bit":
    import mantid_maxent_win64 as maxent
else:
    raise RuntimeError("Currently, MaxEnt will only run on 64-bit Windows " \
                       "installations of Mantid.  Your current installation" \
                       "is %s." % os_env)

# Check whether we can plot by seeing if we are inside Mantid.
try:
    import mantidplot
    INSIDE_MANTIDPLOT = True
except ImportError:
    # Probably running the unit tests.
    INSIDE_MANTIDPLOT = False

RUN_NUM_PROP          = "RunNumber"
INSTRUMENT_PROP       = "Instrument"
FIT_DEADTIME_PROP     = "FitDeadTimes"
FIX_PHASES_PROP       = "FixPhases"
DEFAULT_LEVEL         = "DefaultLevel"
SIGMA_LOOSENESS_PROP  = "SigmaLooseness"
GROUPINGS_PROP        = "GroupingsFile"
PHASES_PROP           = "PhasesInFile"
DEADTIMES_PROP        = "DeadTimesInFile"
PHASES_RESULT_PROP    = "PhasesOutFile"
DEADTIMES_RESULT_PROP = "DeadTimesOutFile"

OUT_WS_PROP = "OutputWorkspace"

SUPPORT_INSTRUMENTS = ["EMU", "HIFI", "MUSR"]

# Have the default input files be the ones we provide in the same directory
# as this script.
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_GROUPINGS_FILENAME = os.path.join(SCRIPT_DIR, "%s_Detector_Grouping_MaxEnt.xml")
DEFAULT_PHASES_FILENAME    = os.path.join(SCRIPT_DIR, "phase.dat")
DEFAULT_DEADTIMES_FILENAME = os.path.join(SCRIPT_DIR, "taud.dat")

def run_from_script_dir(func):
    '''Decorator that changes the working directory to the directory of this
    script for the duration of the decorated function.
    '''
    def change_dir_wrapper(*args, **kwargs):
        cwd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = func(*args, **kwargs)
        os.chdir(cwd)
        return result

    return change_dir_wrapper

class MaxEnt(PythonAlgorithm):
 
    def category(self):
        return 'Muon'
    
    def summary(self):
        return "Runs an underlying \"MaxEnt\" Fortran routine.  Please enter a run "\
               "number.  All other fields are optional.\n\n"\
               "Currently only single-period data is officially supported.  When "\
               "given more than one period, this routine will simply take the first one."

    def PyInit(self):
        """ Alg initialisation. """

        self.declareProperty(name=RUN_NUM_PROP,
                             defaultValue="",
                             direction=Direction.Input,
                             validator=StringMandatoryValidator(),
                             doc="The run number of the data.")

        self.declareProperty(name=INSTRUMENT_PROP,
                             defaultValue="MUSR",
                             direction=Direction.Input,
                             validator=StringListValidator(SUPPORT_INSTRUMENTS),
                             doc="The Instrument from which the data was taken.")

        self.declareProperty(name=FIT_DEADTIME_PROP,
                             defaultValue=True,
                             direction=Direction.Input,
                             doc="Fit dead times.")

        self.declareProperty(name=FIX_PHASES_PROP,
                             defaultValue=False,
                             direction=Direction.Input,
                             doc="Fix phases.")

        self.declareProperty(name=DEFAULT_LEVEL,
                             defaultValue=DEFAULT_LEVEL_VALUE,
                             direction=Direction.Input,
                             doc="Default level.")

        self.declareProperty(name=SIGMA_LOOSENESS_PROP,
                             defaultValue=DEFAULT_SIGMA_LOOSENESS_VALUE,
                             direction=Direction.Input,
                             doc="Sigma losseness.")

        self.declareProperty(FileProperty(name=GROUPINGS_PROP,
                                          defaultValue="",
                                          action=FileAction.OptionalLoad,
                                          extensions=["xml", "map"]),
                             doc="Detector groupings file.  Note that you do " \
                                 "not have to provide this -- INST_Detector_G" \
                                 "rouping_MaxEnt.xml will be used as a default")

        self.declareProperty(FileProperty(name=PHASES_PROP,
                                          defaultValue=DEFAULT_PHASES_FILENAME,
                                          action=FileAction.Load,
                                          extensions = ["dat"]),
                             doc="Input phases file.")

        self.declareProperty(FileProperty(name=DEADTIMES_PROP,
                                          defaultValue=DEFAULT_DEADTIMES_FILENAME,
                                          action=FileAction.Load,
                                          extensions = ["dat"]),
                             doc="Input dead times file.")

        self.declareProperty(FileProperty(name=PHASES_RESULT_PROP,
                                          defaultValue=DEFAULT_PHASES_FILENAME,
                                          action=FileAction.Save,
                                          extensions = ["dat"]),
                             doc="Output phases file.")

        self.declareProperty(FileProperty(name=DEADTIMES_RESULT_PROP,
                                          defaultValue=DEFAULT_DEADTIMES_FILENAME,
                                          action=FileAction.Save,
                                          extensions = ["dat"]),
                             doc="Output dead times file.")
        
        self.declareProperty(WorkspaceProperty(name=OUT_WS_PROP,
                                               defaultValue="",
                                               direction=Direction.Output,
                                               optional=PropertyMode.Optional),
                             doc="Output data.")

    @run_from_script_dir
    def PyExec(self):
        """ Alg execution. """
        instrument         = self.getProperty(INSTRUMENT_PROP).value
        run_number         = self.getProperty(RUN_NUM_PROP).value
        fit_deadtime       = self.getProperty(FIT_DEADTIME_PROP).value
        fix_phases         = self.getProperty(FIX_PHASES_PROP).value
        default_level      = self.getProperty(DEFAULT_LEVEL).value
        sigma_looseness    = self.getProperty(SIGMA_LOOSENESS_PROP).value
        groupings_file     = self.getProperty(GROUPINGS_PROP).value
        in_phases_file     = self.getProperty(PHASES_PROP).value
        in_deadtimes_file  = self.getProperty(DEADTIMES_PROP).value
        out_phases_file    = self.getProperty(PHASES_RESULT_PROP).value
        out_deadtimes_file = self.getProperty(DEADTIMES_RESULT_PROP).value

        isis = config.getFacility('ISIS')
        padding = isis.instrument(instrument).zeroPadding(0)
        run_name = instrument + str(run_number).zfill(padding)

        try:
            run_number = int(run_number)
        except:
            raise RuntimeError("'%s' is not an integer run number." % run_number)
        try:
            run_file_path = FileFinder.findRuns(run_name)[0]
        except:
            raise RuntimeError("Unable to find file for run %i" % run_number)

        if groupings_file == "":
            groupings_file = DEFAULT_GROUPINGS_FILENAME % instrument

        # Load data and other info from input files.

        def temp_hidden_ws_name():
            """Generate a unique name for a temporary, hidden workspace."""
            selection = string.ascii_lowercase + string.ascii_uppercase + string.digits
            return '__temp_MaxEnt_' + ''.join(random.choice(selection) for _ in range(20))

        input_data_ws_name = temp_hidden_ws_name()
        LoadMuonNexus(Filename=run_file_path, OutputWorkspace=input_data_ws_name)
        input_data_ws = mtd[input_data_ws_name]
        
        if isinstance(input_data_ws, WorkspaceGroup):
            Logger.get("MaxEnt").warning("Multi-period data is not currently supported.  Just using first period.")
            input_data_ws = input_data_ws[0]

        groupings_ws_name = temp_hidden_ws_name()
        LoadDetectorsGroupingFile(InputFile=groupings_file, OutputWorkspace=groupings_ws_name)
        groupings_ws = mtd[groupings_ws_name]

        def yield_floats_from_file(path):
            """Given a path to a file with a float on each line, will return
            the floats one at a time.  Throws otherwise.  Strips whitespace
            and ignores empty lines."""
            with open(path, 'r') as f:
                for i, line in enumerate(line.strip() for line in f):
                    if line == "":
                        continue
                    try:
                        yield float(line)
                    except:
                        raise RuntimeError("Parsing error in '%s': Line %d: '%s'." % 
                                           (path, i, line))

        input_phases         = np.array(list(yield_floats_from_file(in_phases_file)))
        input_phases_size    = len(input_phases)
        input_deadtimes      = np.array(list(yield_floats_from_file(in_deadtimes_file)))
        input_deadtimes_size = len(input_deadtimes)

        n_bins      = input_data_ws.blocksize()
        n_detectors = input_data_ws.getNumberHistograms()

        def time_value_to_time_channel_index(value):
            """Given a time value, will return the index of the time channel in
            which the value falls."""
            bin_width = input_data_ws.readX(0)[1] - input_data_ws.readX(0)[0]
            diff = value - input_data_ws.readX(0)[0]
            return int(diff / bin_width)

        # Mantid corrects for time zero on loading, so we want to find the actual channels
        # where 0.0 occurs, and where we have values of 0.1 onwards.
        time_zero_channel  = time_value_to_time_channel_index(0.0)
        first_good_channel = time_value_to_time_channel_index(0.1)

        input_data = np.concatenate([input_data_ws.readY(i) for i in range(n_detectors)])

        groupings = [groupings_ws.readY(row)[0] for row in range(groupings_ws.getNumberHistograms())]
        groupings = map(int, groupings)
        n_groups = len(set(groupings))

        # Cleanup.

        input_data_ws.delete()
        groupings_ws.delete()

        # We're faced with the problem of providing more than a dozen parameters to
        # the Fortran, which can be a bit messy (especially on the Fortran side of
        # things where we need to make "Cf2py" declarations).  A cleaner way of
        # doing this is to simply pass in a few callbacks -- one for each input
        # type -- and have the Fortran provide the name of the variable it wants
        # to the callback.  The callback will then look up the corresponding value
        # and feed it back to the Fortran.
        #
        # We also have a callback for printing to the results log.

        self.int_vars = {
            "RunNo"       : run_number,
            "frames"      : FRAMES,
            "res"         : RES,
            "Tzeroch"     : time_zero_channel,
            "firstgoodch" : first_good_channel,
            "ptstofit"    : POINTS_TO_FIT,
            "histolen"    : n_bins,
            "nhisto"      : n_detectors,
            "n_groups"    : n_groups,
        }

        self.float_vars = {
            "deflevel" : default_level,
            "sigloose" : sigma_looseness,
        }

        self.bool_vars = {
            "fixphase" : fix_phases,
            "fitdt"    : fit_deadtime,
        }

        self._assert_map_values_are_of_expected_type()

        def lookup(par_name, par_map, default):
            """The basis of the callbacks passed to the Fortran.  Given a parameter
            name it will consult the appropriate variable map, and return the
            corresponding value of the parameter.  Else return a default and log a
            warning if a parameter with the name does not exist."""
            par_name = par_name.strip()
            if par_name in par_map:
                return par_map[par_name]
            msg = """WARNING: tried to find a value for parameter with name %s but
            could not find one.  Default of \"%s\" provided.""" % (par_name, default)
            Logger.get("MaxEnt").warning(msg)
            return default

        def log(priority, message):
            """Log the given message with given priority."""
            try:
                logger = getattr(Logger.get("MaxEnt"), priority.lower())
            except AttributeError:
                # If we don't recognise the priority, use warning() as a default.
                logger = getattr(Logger.get("MaxEnt"), "warning")
            logger(message)
            return True

        # The Fortran expects arrays to be of a certain size, so any arrays that
        # aren't big enough need to be padded.
        input_phases    = self._pad_to_length_with_zeros(input_phases, MAX_HISTOS)
        input_deadtimes = self._pad_to_length_with_zeros(input_deadtimes, MAX_HISTOS)
        input_data      = self._pad_to_length_with_zeros(input_data, MAX_INPUT_DATA_SIZE)
        groupings       = self._pad_to_length_with_zeros(groupings, MAX_HISTOS)

        # TODO: Return the contents of "NNNNN.max", instead of writing to file.
        f_out, fchan_out, output_deadtimes, output_phases, chi_sq = maxent.mantid_maxent(
            # Input data and other info:
            input_data,
            groupings,
            input_deadtimes,
            input_phases,
            # Variable-lookup callbacks:
            lambda par_name: lookup(par_name, self.int_vars,   0),
            lambda par_name: lookup(par_name, self.float_vars, 0.0),
            lambda par_name: lookup(par_name, self.bool_vars,  False),
            # Callback for logging:
            log
        )

        def write_items_to_file(path, items):
            """Given a path to a file and a list of items, will write the items
            to the file, one on each line."""
            with open(path, 'w') as f:
                for item in items:
                    f.write(str(item) + "\n")

        # Chop the padded outputs back down to the correct size.
        output_phases    = output_phases[:input_phases_size]
        output_deadtimes = output_deadtimes[:input_deadtimes_size]
        input_phases     = input_phases[:input_phases_size]
        input_deadtimes  = input_deadtimes[:input_deadtimes_size]
        fchan_out        = fchan_out[:n_bins]
        f_out            = f_out[:n_bins]

        write_items_to_file(out_phases_file,    output_phases)
        write_items_to_file(out_deadtimes_file, output_deadtimes)
                 
        log_output = "\nDead times in:\n" +  str(input_deadtimes) + "\n" +\
                     "\nDead times out:\n" + str(output_deadtimes) + "\n" +\
                     "\nPhases in:\n" +      str(input_phases) + "\n" +\
                     "\nPhases out:\n" +     str(output_phases) + "\n" + \
                     "\nGroupings:\n" +      str(groupings) + "\n" +\
                     "\nChi Squared:\n" +    str(chi_sq) + "\n" +\
                     "\nInput variables:\n"

        for type_map in self.int_vars, self.float_vars, self.bool_vars:
            for name, value in type_map.items():
                log_output += str(name) + " = " + str(value) + "\n"

        Logger.get("MaxEnt").notice(log_output)

        # Generate our own output ws name if the user has not provided one.
        out_ws_name = self.getPropertyValue(OUT_WS_PROP)
        if out_ws_name == "":
            out_ws_name = run_name + "; MaxEnt"
            self.setPropertyValue(OUT_WS_PROP, out_ws_name)

        out_ws = CreateWorkspace(OutputWorkspace=out_ws_name,
                                 DataX=fchan_out[:n_bins],
                                 DataY=f_out[:n_bins])
        self.setProperty(OUT_WS_PROP, out_ws)

        # MaxEnt inputs table.
        input_table_name = run_name + "; MaxEnt Input"
        input_table = CreateEmptyTableWorkspace(OutputWorkspace = input_table_name)
        input_table.addColumn("str", "Name")
        input_table.addColumn("str", "Value")
        inputs = itertools.chain(self.int_vars.items(), 
                                 self.float_vars.items(),
                                 self.bool_vars.items())
        for name, value in inputs:
            input_table.addRow([str(name), str(value)])

        # Deadtimes and phases input/output table.
        dead_phases_table_name = run_name + "; MaxEnt Deadtimes & Phases"
        dead_phases_table = CreateEmptyTableWorkspace(OutputWorkspace = dead_phases_table_name)
        for column_name in "Deadtimes In", "Deadtimes Out", "Phases In", "Phases Out":
          dead_phases_table.addColumn("double", column_name)
        for row in zip(input_deadtimes, output_deadtimes, input_phases, output_phases):
            dead_phases_table.addRow(list(map(float, row)))

        # Chi-squared output table.
        chisq_table_name = run_name + "; MaxEnt Chi^2"
        chisq_table = CreateEmptyTableWorkspace(OutputWorkspace = chisq_table_name)
        chisq_table.addColumn("int", "Cycle")
        for iteration in range(10):
          chisq_table.addColumn("double", "Iter " + str(iteration + 1))
        for cycle, data in enumerate(chi_sq):
            chisq_table.addRow([cycle + 1] + list(map(float,data)))

        all_output_ws = [input_table_name,
                         dead_phases_table_name,
                         chisq_table_name,
                         out_ws_name]

        # The output workspaces of this algorithm belong in the same groups
        # that are created by the muon interface.  If the appropriate group
        # doesn't exist already then it needs to be created.
        if not run_name in mtd:
            GroupWorkspaces(InputWorkspaces = all_output_ws,
                            OutputWorkspace = run_name)
        else:
            group = mtd[run_name]
            for output_ws in all_output_ws:
              if not group.contains(output_ws):
                group.add(output_ws)

        out_ws.getAxis(0).getUnit().setLabel("Field", "G")
        out_ws.setYUnitLabel("P(B)")

        if INSIDE_MANTIDPLOT:
            mantidplot.plotSpectrum(out_ws, 0)

    def _pad_to_length_with_zeros(self, array, length):
        """Pad out the given array to the given length, with zeros.  Note that
        we can't use numpy.pad here as numpy version 1.7.0 is required and
        we're only using 1.6.2 in Mantid at present.
        """
        pad_width = length - len(array)
        if pad_width <= 0:
            return array
        return np.concatenate((array, np.zeros(pad_width)))

    def _assert_map_values_are_of_expected_type(self):
        """Just a quick programming check to make sure we're supplying values to the
        Fortran of the types we expect."""
        type_map = {
            int : self.int_vars,
            float : self.float_vars,
            bool : self.bool_vars
        }

        for var_type, var_map in type_map.items():
            for name, value in var_map.items():
                assert isinstance(value, var_type), (
                    "The variable '%s' is supposed to be of type %s but is of type %s."
                    % (name, str(var_type), str(type(value))))

AlgorithmFactory.subscribe(MaxEnt)
