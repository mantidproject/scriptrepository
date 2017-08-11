# -*- coding: utf-8 -*-
# %matplotlib qt
# import matplotlib.pyplot as plt
import datetime
import os.path
import numpy as np
from scipy.optimize import curve_fit
from mantid.simpleapi import (mtd, ConjoinWorkspaces, Load, ConvertUnits,
                              ExtractSingleSpectrum, SaveNexusProcessed,
                              DeleteWorkspace, RebinToWorkspace,
                              RenameWorkspace, Rebin, LoadMask, CloneWorkspace,
                              SumSpectra, MaskDetectors, GroupWorkspaces,
                              CreateWorkspace, DeleteWorkspaces, WeightedMean,
                              Fit)
from .runtypes import HeData


BASE = r"LARMOR{:08d}.nxs"
# BASE = "{}"

RANGE = "1.8,0.1,8"

def sumToShim(rnum, output_dir=None):
    """
    Combine both spin states into a single workspace

    Parameters
    ----------
    rnum : int
      The run number to be shimmed
    output_dir : string
      If given, the folder where the workspace should be saved

    """
    try:
        wtemp = Load(BASE.format(rnum), LoadMonitors=True)
        RebinToWorkspace('wtemp_1', 'wtemp_monitors_1', PreserveEvents=False,
                         OutputWorkspace='wtemp_1')
        RebinToWorkspace('wtemp_2', 'wtemp_monitors_1', PreserveEvents=False,
                         OutputWorkspace='wtemp_2')
        wtemp_1 = ConjoinWorkspaces('wtemp_monitors_1', 'wtemp_1')
        wtemp_2 = ConjoinWorkspaces('wtemp_monitors_2', 'wtemp_2')
    except:
        wtemp_monitors = Load(BASE.format(rnum))
    wtempShim = mtd['wtemp_monitors_1'] + mtd['wtemp_monitors_2']
    RenameWorkspace(wtempShim, 'LARMOR{:08d}'.format(rnum))
    if output_dir:
        SaveNexusProcessed('LARMOR{:08d}'.format(rnum),
                           os.path.join(output_dir,
                                        "LARMOR{:08d}-add.nxs".format(rnum)))
    RenameWorkspace('LARMOR{:08d}'.format(rnum),
                    'LARMOR{:08d}-add'.format(rnum))


def he3pol(scale, time):
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


def he3_stats(run):
    """
    Return the information about the ³He analyser's state during this
    run.  This function requires that get_he3_log has already been
    run, creating the "helium_log" table.

    Parameters
    ----------
    run
      A RunData object containing the run whose analyser statistics are needed

    Return
    ------
    A HeData object describing the state of the analyser during the run.

    """
    timecode = run.start.isoformat()
    stats = [x for x in mtd["helium_log"]
             if timecode > x["Start time"]][-1]
    return HeData(stats["Number"], stats["Cell"], stats["scale"],
                  datetime.datetime.strptime(stats["Start time"],
                                             "%Y-%m-%dT%H:%M:%S"),
                  stats["fid"], stats["Time Constant"])


def int3samples(runs, name, masks, binning='0.5, 0.05, 8.0'):
    """
    Finds the polarisation versus wavelength for a set of detector tubes.

    Parameters
    ----------
    runs: list of RunData objects
      The runs whose polarisation we are interested in.

    name: string
      The name of this set of runs

    masks: list of string
      The file names of the masks for the sequential tubes that are being used
      for the SEMSANS measurements.

    binning: string
      The binning values to use for the wavelength bins.  The default value is
      '0.5, 0.025, 10.0'
    """
    for tube, _ in enumerate(masks):
        for i in [1, 2]:
            final_state = "{}_{}_{}".format(name, tube, i)
            if final_state in mtd.getObjectNames():
                DeleteWorkspace(final_state)

    for rnum in runs:
        w1 = Load(BASE.format(rnum.number), LoadMonitors=True)
        w1mon = ExtractSingleSpectrum('w1_monitors', 0)
        w1 = ConvertUnits('w1', 'Wavelength', AlignBins=1)
        w1mon = ConvertUnits(w1mon, 'Wavelength')
        w1 = Rebin(w1, binning, PreserveEvents=False)
        w1mon = Rebin(w1mon, binning)
        w1 = w1 / w1mon
        for tube, mask in enumerate(masks):
            Mask_Tube = LoadMask('LARMOR', mask)
            w1temp = CloneWorkspace(w1)
            MaskDetectors(w1temp, MaskedWorkspace="Mask_Tube")
            Tube_Sum = SumSpectra(w1temp)
            for i in [1, 2]:
                final_state = "{}_{}_{}".format(name, tube, i)
                if final_state in mtd.getObjectNames():
                    mtd[final_state] += mtd["Tube_Sum_{}".format(i)]
                else:
                    mtd[final_state] = mtd["Tube_Sum_{}".format(i)]

    x = mtd["{}_0_1".format(name)].extractX()[0]
    dx = (x[1:] + x[:-1]) / 2
    pols = []

    for run in runs:
        he_stat = he3_stats(run)
        start = (run.start-he_stat.dt).seconds/3600/he_stat.t1
        end = (run.end-he_stat.dt).seconds/3600/he_stat.t1
        for time in np.linspace(start, end, 10):
            temp = he3pol(he_stat.scale, time)(dx)
            pols.append(temp)
    wpol = CreateWorkspace(x, np.mean(pols, axis=0),
                           # and the blank
                           UnitX="Wavelength",
                           YUnitLabel="Counts")

    for tube, _ in enumerate(masks):
        up = mtd["{}_{}_2".format(name, tube)]
        dn = mtd["{}_{}_1".format(name, tube)]
        pol = (up - dn) / (up + dn)
        pol /= wpol
        DeleteWorkspaces(["{}_{}_{}".format(name, tube, i)
                          for i in range(1, 3)])
        RenameWorkspace("pol",
                        OutputWorkspace="{}_{}".format(name, tube))
    DeleteWorkspaces(["Tube_Sum_1", "Tube_Sum_2"])

    GroupWorkspaces(["{}_{}".format(name, tube)
                     for tube, _ in enumerate(masks)
                     for i in range(1, 3)],
                    OutputWorkspace=str(name))


def norm(sample, blank, masks):
    """
    Normalise the neutron polarisation on a tube by tube basis

    Parameters
    ----------
    sample: string
      The name of the sample to be normalised.  The individual tubes are
      assumed to be in workspaces with names like sample_2
    blank: string
      The name of the blank to be normalised agains.  The individual tubes are
      assumed to be in workspaces with names like blank_2
    masks: list of string
      The file names for the masks used for the individual tubes.
    """
    for t, _ in enumerate(masks):
        wtemp = mtd[sample + "_{}".format(t)] / \
            mtd[blank + "_{}".format(t)]

        y = mtd[blank + "_{}".format(t)].extractY()
        e = wtemp.extractE()

        e[np.abs(y) < 0.2] *= 1e9
        wtemp.setE(0, e[0])

        RenameWorkspace("wtemp", sample + "_{}_norm".format(t))
    wtemp = WeightedMean(sample + "_0_norm",
                         sample + "_1_norm")
    for tube in range(2, len(masks)):
        wtemp = WeightedMean(wtemp, sample + "_{}_norm".format(tube))
    RenameWorkspace(wtemp, OutputWorkspace=sample + "_Norm")
    DeleteWorkspaces(["{}_{}_norm".format(sample, tube)
                      for tube, _ in enumerate(masks)])


def sel_const(runs, dist=4.0, thickness=5e-3,
              show_fits=False, show_quality=False):
    """Calculate the spin echo length of the instrument

    Parameters
    ----------
    runs: list of Workspaces
      A list of the workspecaes containing the polarisation versus wavelength
      for consecutive detector tubes
    dist: float
      The distances from the sample to the detector in meters.  The default
      is 4.0
    thickness: float
      The distance between detector tubes in meters.  Defaults to 5mm.
    show_fits: bool
      If true, plots the sinusoid fits used to calculate the frequency
    show_quality: bool
      If true, plots the frequency versus tube position to confirm that the
      frequency grows linearly with position.

    Returns
    -------
    A float containing the spin echo length, in nanometers, of a one angstrom
    neutron.
    """
    freqs = []
    for run in runs:
        x = run.extractX()[0]
        x = (x[1:] + x[:-1]) / 2
        p = run.extractY()[0]
        p[np.isnan(p)] = 0
        fp = np.fft.fft(p)

        conv = len(x) / (np.max(x) - np.min(x))
        max_arg = (np.nanargmax(np.abs(fp[1:int(len(p) / 2)])) + 1)
        amp = np.abs(fp[max_arg]) / len(x)
        max_arg = max_arg * conv / len(x)

        model = "name=UserFunction,Formula=e*cos(x*f),f={},e={}"
        Fit(Function=model.format(max_arg*2*np.pi, amp),
            InputWorkspace=run, StartX=3, EndX=7, CreateOutput=True)

        result = mtd[run.getName()+"_Parameters"]

        freqs.append(np.abs(result.column(1)[1]/2/np.pi))
        if not show_fits:
            DeleteWorkspace(run.getName()+"_Parameters")
            DeleteWorkspace(run.getName()+"_NormalisedCovarianceMatrix")
            DeleteWorkspace(run.getName()+"_Workspace")

    def model(x, m, b):
        """A simple linear model"""
        return m * x + b
    xs = np.arange(len(runs))*thickness
    fit, _ = curve_fit(model, xs, freqs)

    sel_fits = CreateWorkspace(xs, freqs)
    Fit(Function="name=LinearBackground",
        InputWorkspace=sel_fits, CreateOutput=True)
    result = 0.1 * mtd["sel_fits_Parameters"].column(1)[1] * dist
    if not show_quality:
        DeleteWorkspace("sel_fits_Parameters")
        DeleteWorkspace("sel_fits_NormalisedCovarianceMatrix")
        DeleteWorkspace("sel_fits_Workspace")
        DeleteWorkspace("sel_fits")

    return result


def sel(sample, const):
    """
    Convert workspace into normalised spin echo units

    Parameters
    ----------
    sample : str
      The name of the workspace to be converted
    const : float
      The spin echo length of a one angstrom neutron in nanometers

    Returns
    -------
    The name of a workspace with the x coordinate in proper spin echo
    units and the y coordinate in log(P/P₀)/λ²
    """
    wtemp = mtd[sample]
    x = wtemp.extractX()
    y = wtemp.extractY()
    e = wtemp.extractE()

    dx = (x[:, 1:] + x[:, :-1]) / 2

    wtemp = ConvertUnits(sample, "SpinEchoLength", EFixed=const)

    wenv = CreateWorkspace(wtemp.extractX(), np.log(y) / dx**2, e / y / dx**2,
                           UnitX="SpinEchoLength",
                           YUnitLabel="Counts")
    RenameWorkspace(wenv, OutputWorkspace=sample + "_sel")
