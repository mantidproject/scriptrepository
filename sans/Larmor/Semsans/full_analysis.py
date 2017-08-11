"""

This module provides helper functions for automatinc the semsans
analysis.  The most important function from the user perspective is
analyse.

"""

from os.path import join, exists, basedir
import sans.command_interface.ISISCommandInterface as ici
from mantid.api import mtd
from mantid.simpleapi import DeleteWorkspaces
from .Semsans import sumToShim, sel_const, int3samples, norm, sel
from .runtypes import table_to_run


def get_shimed(run, path, refresh=False):
    """

    Gives the path for the detector image for a run with both spin
    states merged.  If the spin state combination has already been
    performed, then the file path is immediately returned.  If the
    file does not already exist, then it can be created.

    Parameters
    ----------
    run
      The run number whose detector image is being requested
    path
      The folder where the detector images should be saved
    reload
      Whether to always recreate the shimmed file.  The default is False

    Return
    ------
    A string containing the location of the shimmed file

    """
    f = join(path, "LARMOR{:08d}-add.nxs".format(run))
    if refresh or not exists(f):
        sumToShim(run, path)
    return f


def analyse(data_table, masks, output_file,
            show_fits=False, show_quality=False):
    """

    Analyse a set of combined Sans and Semsans measurements.
    The result will be a set of workspaces with the Semsans data
    and a csv file that can be loaded into the batch sans interface
    to load the sans data.

    Parameters
    ----------
    data_table
      The name of WorkspaceTable with the runs that need to be analysed.  This
      table should be generated by the get_log function.
    masks
      A list of strings containing the filenames of the masks for the
      individual tubes used in the semsans polarisation calculation
    output_file
      Where to save the CSV file that will generate the SANS runs
    show_fits
      Whether to show the individual tube fits used to calculate the spin echo
      length.  Defaults to False
    show_quality
      Whether to show the quality of the linear fit of the precession
      frequencies used to calculate the spin echo length.  Defaults to False

    """
    if "Full Blank" not in mtd.getObjectNames():
        int3samples(table_to_run(mtd["Full Blank_runs"]), "Full Blank", masks)
    const = sel_const([mtd["Full Blank_{}".format(tube)]
                       for tube, _ in enumerate(masks)],
                      show_fits=show_fits, show_quality=show_quality)

    k = data_table[:-5]
    runs = table_to_run(mtd[data_table])
    for idx, run in enumerate(runs):
        get_shimed(run.number, basedir(output_file))
        int3samples([run], "{}_run{:02d}".format(k, idx), masks)
        # DeleteWorkspace("{}_sans_nxs".format(run.number))
        semsans_ws = "{}_run{:02d}".format(k, idx)
        norm(semsans_ws, "Full Blank", masks)
        sel(semsans_ws+"_norm", const)
        DeleteWorkspaces([semsans_ws,
                          semsans_ws+"_Norm"])
    with open(output_file, "w") as outfile:
        framework = "sample_sans,{}-add," \
                    "sample_trans,{}-add," \
                    "sample_direct_beam,{}-add,"\
                    "can_sans,{}-add,"\
                    "can_trans,{}-add,"\
                    "can_direct_beam,{}-add,"\
                    "output_as,hours_{:0.2f}\n"
        for idx, run in enumerate(runs):
            outfile.write(
                framework.format(run.number, run.trans, run.direct,
                                 run.csans, run.ctrans, run.direct,
                                 (run.start-runs[0].start).seconds/3600))
