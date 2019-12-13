import numpy as np
from collections import Sequence
from mantid.api import WorkspaceGroup, Workspace  # type: ignore
from mantid.simpleapi import Load, Divide, Rebin, SaveSESANS, RenameWorkspace, ConvertUnits, SumSpectra, ExtractSingleSpectrum, MaskDetectors, DeleteWorkspace, Logarithm, CloneWorkspace, Scale  # type: ignore

POL_SCALE=10

def patterson(wstemp, const):
    # type: (Workspace, float) -> Workspace
    """Convert workspace into spin echo form"""
    temp = CloneWorkspace(wstemp)
    x = temp.extractX()
    x = (x[:, 1:]+x[:, :-1])/2
    temp = Logarithm(temp)
    for i in range(x.shape[0]):
        temp.setY(i, temp.readY(i)/x[i]**2/temp.sample().getThickness())
        temp.setE(i, temp.readE(i)/x[i]**2/temp.sample().getThickness())
        temp = ConvertUnits(temp, "SpinEchoLength", EFixed=const)
    return temp


def echo_cal(angle):
    # type: (float) -> float
    # return 1/np.tan(np.array(np.abs(angle))*np.pi/180)*0.05647
    # return np.polyval([-4.00540061e-07, 8.88410387e-05, -
    #                    7.57804680e-03, 2.54373296e-01], np.abs(angle))
    return np.polyval([5.10018568e-09, -1.51816309e-6, 1.7463217e-04, -
                       1.02919899e-02, 2.84070156e-01], np.abs(angle))


def get_angle(run):
    # type: (MetaData) -> float
    """Find the angle that a run was measured at"""
    # takes a single run
    return run["SimTheta"][-1]
    # return np.median(run.getRun().getProperty("Mag1_Theta").value)


class RunSet(object):
    """A group of individual runs that form a complete dataset"""

    def __init__(self, sample, p0=None, trans=None, p0trans=None, title=None,
                 binning="1.8,0.05,12.0", mask="MaskWorkspace", progress=None):
        # type: (List[int], Optional[List[int]], Optional[List[int]], Optional[List[int]], Optional[str], str, str) -> None
        if p0 is None:
            p0 = []
        if trans is None:
            trans = []
        if p0trans is None:
            p0trans = []
        if progress:
            progress.setNumSteps(POL_SCALE*len(sample)+POL_SCALE*len(p0)+len(trans)+len(p0trans))
        self._sample = map(lambda x: MetaData(x, binning, mask, progress), sample)
        self._p0 = map(lambda x: MetaData(x, binning, mask, progress), p0)
        self._p0trans = get_counts(p0trans, binning, progress)
        self._trans = get_counts(trans, binning, progress)
        self._title = title

    @property
    def blank(self):
        """The background polarisation run."""
        return self._p0

    @property
    def blank_trans(self):
        """The background transmission run."""
        return self._p0trans

    @property
    def sample(self):
        """The sample polarisation run."""
        return self._sample

    @property
    def trans(self):
        """The sample transmission run."""
        return self._trans

    @property
    def title(self):
        """The title of the measurement group"""
        if self._title:
            return self._title
        return self._sample[0].title

    def reduction(self):
        """Perform a reduction on the RunSet"""
        if not self.blank:
            return calculate_polarisation(self.sample)
        return self.transmission_normalised_sesans()

    def transmission_normalised_sesans(self):
        """Calculate the Sesans measurement of a run, correcting for
        instrument polarisation and sample transmission."""

        angle = get_angle(self.sample[0])
        wstemp = Divide(calculate_polarisation(self.sample),
                        calculate_polarisation(self.blank),
                        OutputWorkspace=self.title)
        if self.trans and self.blank_trans:
            samp = calculate_counts(self.sample)
            blnk = calculate_counts(self.blank)
            # Losses on an unscattered beam
            lost = (samp/self.trans)/(blnk/self.blank_trans)

            wstemp *= lost
            DeleteWorkspace(lost)
        SaveSESANS(wstemp, (self.title + ".ses").replace(":", "colon"),
                   ThetaZMax=0.1, ThetaYMax=0.1,
                   EchoConstant=echo_cal(angle)*10000, Sample=self.sample[0].title)
        wstemp = patterson(wstemp, echo_cal(angle))
        # DeleteWorkspace(wstemp)
        RenameWorkspace(wstemp, "{}_sel".format(self.title))
        return wstemp


class MetaData(object):
    """The additional information associated with a run"""

    def __init__(self, x, binning, mask, progress):
        # type: (int, str, str) -> None
        self._run_number = x
        temp = _cheap_load(x, binning, mask)
        self._data = temp
        self._title = temp[0].getTitle()
        self._env = {}  # type: Dict[str, Union[str, float]]
        # Must make a copy of the run data before the
        # temp variable gets garbage collected
        for k in temp[0].getRun().keys():
            result = temp[0].getRun()[k].value
            if isinstance(result, Sequence) and not isinstance(result, str):
                self._env[k] = result[0]
            else:
                self._env[k] = result
        if progress:
            progress.reportIncrement(POL_SCALE,
                                     "Loading run {}".format(x))

    def __del__(self):
        DeleteWorkspace(self._data)

    @property
    def is_transmission(self):
        # type: () -> bool
        """Does this run hold transmission data"""
        return self["nperiods"] == 1

    @property
    def is_blank(self):
        # type: () -> bool
        """Was this run performed on a blank?"""
        return "Blank" in self._title

    @property
    def is_echo(self):
        # type: () -> bool
        """Is this an echo scan?"""
        return self["nperiods"] > 2

    @property
    def is_sample(self):
        # type: () -> bool
        """Was this run performed on a sample in Sesans mode?"""
        return not self.is_blank and not self.is_transmission

    @property
    def run_number(self):
        # type: () -> int
        """What is the index of the run in the run journal?"""
        return self._run_number

    @property
    def title(self):
        # type: () -> str
        """What was the title of the run?"""
        return self._title

    @property
    def data(self):
        # type: () -> Workspace
        """The actual measurement data"""
        return self._data

    def is_appropriate(self, run):
        # type: (MetaData) -> bool
        """Were these two runs performed under the same conditions?"""
        trial = [round(run[x][-1], 2) == round(self[x][-1], 2)
                 for x in ["Mag2_Theta", "DCMagField1", "Echo_Coil"]]
        return self != run and all(trial)

    def __getitem__(self, key):
        return self._env[key]

    def find_nearest_appropriate(self, values, match_sample=False,
                                 match_current=True):
        """Given a list of metadata, find the last run before this one with
        was performed under the same instrument conditions."""
        options = []
        for meta in values:
            if (not match_current or self.is_appropriate(meta)) \
               and (not match_sample or self.title in meta.title):
                options.append(meta)
        if options:
            return min(options,
                       key=lambda x: abs(self.run_number - x.run_number))
        return None


def _cheap_load(run, binning, mask):
    # type: (int, str, str) -> Workspace
    """Load a run with proper memoisation"""
    temp = Load(str(run), LoadMonitors=True)
    temp = ConvertUnits("temp", "Wavelength", AlignBins=1)
    temp = Rebin(temp, binning, PreserveEvents=False)
    temp_monitors = ExtractSingleSpectrum("temp_monitors", 0)
    temp_monitors = ConvertUnits(temp_monitors, "Wavelength", AlignBins=1)
    temp_monitors = Rebin(temp_monitors, binning, PreserveEvents=False)
    temp /= temp_monitors
    DeleteWorkspace("temp_monitors")
    MaskDetectors(temp, MaskedWorkspace=mask)
    temp = SumSpectra(temp)
    RenameWorkspace(temp, "{}_raw".format(run))
    return temp


def calculate_polarisation(runs):
    # type: (List[MetaData]) -> Workspace
    """Calculate the combined polarisation for a list of run numbers"""
    title = "{}_pol_combined".format(
        "_".join(map(lambda x: str(x.run_number), runs)))

    w1 = CloneWorkspace(runs[0].data)
    for run in runs[1:]:
        w1 += run.data
    # RenameWorkspace(w1, title)

    up = w1[0]
    dn = w1[1]
    pol = (up - dn) / (up + dn)
    pol = RenameWorkspace(pol, OutputWorkspace=title)
    return pol

def calculate_counts(runs):
    # type: (List[MetaData]) -> Workspace
    """Calculate the combined intensity for a list of run numbers"""
    title = "{}_int_combined".format(
        "_".join(map(lambda x: str(x.run_number), runs)))

    w1 = CloneWorkspace(runs[0].data)
    for run in runs[1:]:
        w1 += run.data
    # RenameWorkspace(w1, title)

    up = w1[0]
    dn = w1[1]
    intensity = Scale(up + dn, 1.0/len(runs))
    intensity = RenameWorkspace(intensity, OutputWorkspace=title)
    return intensity

def get_counts(runs, binning, progress):
    # type: (List[int], str) -> Optional[Workspace]
    """Calculate the combined intensity for a group of transmission runs"""
    if not len(runs):
        return None
    title = "{}_int_combined".format(
        "_".join(map(lambda x: str(x), runs)))

    w1 = Load(str(runs[0]))
    for run in runs[1:]:
        progress.report("Loading Transmission {}".format(run))
        wtemp = Load(str(run))
        w1 += wtemp
    if len(runs) > 1:
        DeleteWorkspace(wtemp)
    # RenameWorkspace(w1, title)

    w1 = ConvertUnits("w1", "Wavelength", AlignBins=1)
    w1 = Rebin(w1, binning, PreserveEvents=False)
    w1_base = ExtractSingleSpectrum(w1, 0)
    w1 = ExtractSingleSpectrum(w1, 3)
    w1 /= w1_base
    w1 = RenameWorkspace(w1, OutputWorkspace=title)
    return w1
