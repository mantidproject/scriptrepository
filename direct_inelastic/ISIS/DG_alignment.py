from mantid.simpleapi import *
from mantidqt.utils.qt import import_qt
import numpy as np
from mantidqt.widgets.instrumentview.api import get_instrumentview
from mantid.geometry import SpaceGroupFactory, CrystalStructure, ReflectionGenerator, ReflectionConditionFilter, SpaceGroup, UnitCell
from scipy.spatial.transform import Rotation as R
from enum import Enum
from scipy.optimize import minimize
from datetime import datetime
import sys, time
from pathlib import Path

#====================================================#
#Script for determining UB Matrices for samples on machines at ISIS
#Users also have options for peak prediction and reorientation
#====================================================#
#     User input:                                    #
#====================================================#
#     Lattice Parameters                             #
lattice_parameters = {
    "a": 5.4, "b": 5.4, "c": 5.4,
    "alpha": 90, "beta": 90, "gamma": 90
}
#==================================================#
#           Runs + Peaks                              #
run_numbers = [52146, 52146]
corresponding_hkls = [[2,2,0], [1,1,1]] # please write in same order the runs are written in
# e.g[[h1,k1,l1],[h2,k2,l2]]

peak_picking_method = 'ROI'  # Either 'MANUAL' or 'ROI'

UB_method = 'TWO_PEAKS' # Available methods: ONE_PEAK, TWO_PEAKS, THREE_PEAKS, FIVEPLUS_PEAKS
#ONE_PEAK = Method for known scattering plane + 1 peak
#TWO_PEAK = Method for if you are able to index 2 non-collinear peaks
#THREE_PEAK = 3 Peaks, lattice parameters are not required
#FIVEPLUS_PEAKS = First automatically index peaks + find UB algorithm
#for method = ONE PEAK, vectors defining the scattering plane are required
vector_1 = []
vector_2 = []

#    Peak to Orient to u vector
Reorient = True # Enter True or False
hkl_of_desired_alignment= [2,2,0] # please input as a single peak: [h,k,l]

Predict_peaks = True # Enter true or false
min_dspacing = 1.5 # in A^-1
max_dspacing = 9.0 # in A^-1

unit_cell_space_group = "F d -3 m" # this either needs to be a string of the symbol e.g "F d -3 m" or and integer with no "",
# if you don't know your space group number enter "N/A"
unit_cell_symmetry = 'PRIMITIVE'
#Use from available list and enter in the format reflection_condition.{CELL_SYMMETRY} list and enter in the format reflection_condition.{CELL_SYMMETRY},
# PRIMITIVE by default change by replacing PRIMITIVE with desired symmetry condition
# which are any of: PRIMITIVE, BODY_CENTERED, C_FACE_CENTERED, B_FACE_CENTERED, A_FACE_CENTERED, ALL_FACE_CENTERED, RHOMBO_CENTERED_REV, HEX_CENTERED_REV

#===============================#
#   Staff Input
instrument_name = 'INSTRUMENT_NAME'
cycle = 'CYCLE_ID'
rbnum = 'USER_RB_FOLDER'

cycle_shortform = cycle[2:] if cycle.startswith('20') else cycle
savedir = f'/data/analysis/{instrument_name}/RBNumber/{rbnum}/'
datadir = f'/data/instrument/{instrument_name}/CYCLE20{cycle_shortform.replace("_","")}/{rbnum}/'

d_tolerance = 0.1 # Tolerenace on peak position, in d-spacing in A^-1
hkl_tolerance = 0.2 # Used in UB quality validation, for how close indexes calculated have to be to whole
gonio_axes = {"Axis0":"Rot,0,1,0,1"} # See https://docs.mantidproject.org/nightly/algorithms/SetGoniometer-v1.html

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #Output Details#
config['default.instrument'] = instrument_name

##Checking quality of inputs, this checks if there are enough listed hkls for peaks given
if len(run_numbers) != len(corresponding_hkls):
    raise ValueError("Number of hkls = number of peaks")

#this block gets a list of all the space group symbols in mantid, and their corresponding numbers
spgr_symbols = SpaceGroupFactory.getAllSpaceGroupSymbols()
spgr_numbers = [SpaceGroupFactory.createSpaceGroup(symbol).getNumber() for symbol in spgr_symbols]
space_group = {}

#dictionary with space group symbol, & corresponding number for dictionary lookup
for x in range(len(spgr_symbols)):
    space_group[spgr_symbols[x]] = spgr_numbers[x]

#working out mantid allowed space group from input
keys = lattice_parameters.keys()
lattice_string = ""
for x in keys:
    lattice_string += str(lattice_parameters[x]) + " "
match unit_cell_space_group:
    case str():
        if unit_cell_space_group != "N/A":
            xtal = CrystalStructure(lattice_string,unit_cell_space_group, "")
    case int():
        keys = [key for key, val in space_group.items() if val == unit_cell_space_group]
        general_space_group = keys[0]
        xtal = CrystalStructure(lattice_string, general_space_group, "")
    case _:
        logger.warning("Invalid Space group input")
        unit_cell_space_group = "N/A"
        if Predict_peaks:
            print(f"The reflection condition {unit_cell_symmetry} will be used instead")


##calculating d_spacing for each hkl, this is to set limits on IV
sample_cell = UnitCell(*lattice_parameters.values())
corresponding_dspacing = []
for x in range(len(run_numbers)):
    h, k, l = float(corresponding_hkls[x][0]), float(corresponding_hkls[x][1]), float(corresponding_hkls[x][2])
    d_spacing = sample_cell.d(h, k, l)
    corresponding_dspacing.append(d_spacing)

#Manual peak picking method
match peak_picking_method:
    case "MANUAL":
        for x in range(len(run_numbers)):
            ws = Load(f'{datadir}/{instrument_name[:3]}{run_numbers[x]}', EnableLogging=False)
            if x == 0:
                CreatePeaksWorkspace(ws, 0, OutputWorkspace="SingleCrystalPeakTable",  EnableLogging=False)
            ws = ConvertUnits(ws, Target='dSpacing', EMode='Elastic', EFixed=None, EnableLogging=False) #converting  TOF to dspacing

            #if Ei exists delete it as it causes peak pick to error
            if ws.getRun().hasProperty("Ei"):
                DeleteLog('ws', 'Ei', EnableLogging=False)

            SetGoniometer(ws, **gonio_axes, EnableLogging=False)

            #displays IV for each run
            iv = get_instrumentview(ws)
            iv.select_pick_tab()
            #sets range on IV
            bin_range_lower= corresponding_dspacing[x] - d_tolerance -0.1
            bin_range_upper= corresponding_dspacing[x] + d_tolerance +0.1
            iv.set_bin_range(bin_range_lower, bin_range_upper)
            iv.show_view()
            print(f"The Peak you are searching for is {corresponding_hkls[x]}")
            while True:
                time.sleep(1)
                try:
                    iv.get_pick_tab()
                except:
                    break
            #sets hkl on the peaksworkspace
            SingleCrystalPeakTable = mtd["SingleCrystalPeakTable"]
            SingleCrystalPeakTable.getPeak(x).setHKL(corresponding_hkls[x][0], corresponding_hkls[x][1],
                                                 corresponding_hkls[x][2])
        SingleCrystalPeakTable = mtd["SingleCrystalPeakTable"]



##ROI + pick best pixel
    case "ROI":
        for x in range(len(run_numbers)):
            ws = Load(f'{datadir}/{instrument_name[:3]}{run_numbers[x]}', EnableLogging=False)
            #creates peak table on first iteration
            if x == 0:
                CreatePeaksWorkspace(ws, 0, OutputWorkspace="SingleCrystalPeakTable",  EnableLogging=False)
            ws = ConvertUnits(ws, Target='dSpacing', EMode='Elastic', EFixed=None,  EnableLogging=False)#converting TOF to dspacing
            SetGoniometer(ws, **gonio_axes, EnableLogging=False)

            if ws.getRun().hasProperty("Ei"):
                DeleteLog('ws', 'Ei', EnableLogging=False)

            #opens instrument view for each run so user can draw mask
            iv = get_instrumentview(ws)
            iv.select_tab(2)
            bin_range_lower= corresponding_dspacing[x] - d_tolerance -0.1
            bin_range_upper= corresponding_dspacing[x] + d_tolerance +0.1
            iv.set_bin_range(bin_range_lower, bin_range_upper)
            iv.show_view()
            print(f"The Peak you are searching for is {corresponding_hkls[x]}")
            while True:
                time.sleep(1)
                try:
                    iv.get_pick_tab()
                except:
                    break

            #get mask workspace which user has gotten
            mask = mtd["MaskWorkspace"]
            ispec = np.flatnonzero(~mask.extractY().astype(bool))
            #getting corresponding spectra
            ws_out = ExtractSpectra(ws, OutputWorkspace="ws_out", WorkspaceIndexList=ispec,  EnableLogging=False)
            lower_limit = corresponding_dspacing[x] - d_tolerance
            upper_limit = corresponding_dspacing[x] + d_tolerance
            #integrate
            ws_int = Integration(InputWorkspace=ws_out, Outputworkspace="ws_int", RangeLower=lower_limit, RangeUpper=upper_limit,  EnableLogging=False)
            ispec_max = ispec[np.argmax(ws_int.extractY())]
            #find the corresponding spec number
            #get det id
            detids = ws.getSpectrum(int(ispec_max)).getDetectorIDs()
            detid = detids[int(round(len(detids)/2,0))]

            SingleCrystalPeakTable = mtd["SingleCrystalPeakTable"]
            AddPeak(SingleCrystalPeakTable, ws, corresponding_dspacing[x], detid)
            SingleCrystalPeakTable.getPeak(x).setHKL(corresponding_hkls[x][0], corresponding_hkls[x][1],
                                                 corresponding_hkls[x][2])
            DeleteWorkspace(ws_out,  EnableLogging=False)
            DeleteWorkspace(ws_int, EnableLogging=False)
            DeleteWorkspace(ws,  EnableLogging=False)
            DeleteWorkspace(mask,  EnableLogging=False)
            # deleting for good file management
        SingleCrystalPeakTable = mtd["SingleCrystalPeakTable"]

    case _:
        raise RuntimeError(f'Unknown peak picking method {peak_picking_method}')

Number_of_peaks = SingleCrystalPeakTable.getNumberPeaks()
used_hkls = {}
for x in range(len(corresponding_hkls)):
    hkl_provided = np.array(SingleCrystalPeakTable.getPeak(x).getHKL())
    peak_predicted_d=sample_cell.d(hkl_provided[0],hkl_provided[1],hkl_provided[2])
    if not np.isclose(corresponding_dspacing[x], peak_predicted_d, atol=d_tolerance):
        logger.warning(f"The {hkl_provided} peak may not be indexed correctly, the dspacing calculated from hkl on workspace for this peak was {peak_predicted_d}, d spacing expected for this peak is {corresponding_d[x]}")

match UB_method:
    case "ONE_PEAK":
        if Number_of_peaks > 1:
            logger.information("This method gives a rough estimate, CalculateUMatrix or FindUBUsingLatticeParameters might be more appropriate")
        FindUBFromScatteringPlane(Vector1=vector_1, Vector2=vector_2, **lattice_parameters, PeaksWorkspace="SingleCrystalPeakTable",  EnableLogging=False)
        used_hkls[run_numbers[0]]= np.array(corresponding_hkls[0])
    case "TWO_PEAKS":
        if Number_of_peaks < 2:
            raise ValueError("Not Enough Peaks to use Calculate U matrix method")
        else:
            CalculateUMatrix(PeaksWorkspace="SingleCrystalPeakTable", **lattice_parameters,  EnableLogging=False)
        ##checking output UB is reasonable
    case "THREE_PEAKS":
        if Number_of_peaks < 3:
            raise ValueError("Not Enough Peaks to use FindUBUsingIndexedPeaks")
        FindUBUsingIndexedPeaks("SingleCrystalPeakTable", tolerance=tolerance,  EnableLogging=False)

    case "FIVEPLUS_PEAKS":
        if Number_of_peaks < 5:
            raise ValueError("Not enough peaks to use FindUBUsingLatticeParameters")
        else:
            FindUBUsingLatticeParameters(PeaksWorkspace="SingleCrystalPeakTable", **lattice_parameters, NumInitial=Number_initial_peaks,  EnableLogging=False)
    case _:
        raise ValueError(f"UB Method {UB_method} doesn't exist")

if UB_method != 'ONE_PEAK':
    Number_of_peaks = SingleCrystalPeakTable.getNumberPeaks()
    temporary_workspace = CloneWorkspace("SingleCrystalPeakTable", "SingleCrystalPeakTable_temp", EnableLogging=False)
    Indexed=IndexPeaks("temporary_workspace", tolerance=hkl_tolerance, RoundHKLs=1,  EnableLogging=False)
    percentage_indexed=int(Indexed[0])/Number_of_peaks*100
    logger.information("Percentage of peaks provided indexed within set tolerance: {}%".format(percentage_indexed))
    if percentage_indexed < 80:
        logger.warning("Relatively low percentage were indexed within tolerance, you may wish to check your inputted hkl or adjust the set tolerance")
    logger.information(f"Average error in hkl: {Indexed[1]}, given hkl tolerance: {hkl_tolerance}")
    DeleteWorkspace(temporary_workspace,  EnableLogging=False)

    for x in range(len(corresponding_hkls)):
        used_hkls[run_numbers[x]] = np.array(corresponding_hkls[x])
#doing some auto index to check validity of UB

##Creates predict peaks table,
if Predict_peaks == True:
    if unit_cell_space_group != "N/A":
        PredictPeaks(SingleCrystalPeakTable, MinDSpacing=min_dspacing, MaxDSpacing=max_dspacing, OutputWorkspace = "PredictedPeaks")
    else:
        PredictPeaks(SingleCrystalPeakTable,ReflectionCondition=unit_cell_symmetry.value,  MinDSpacing=min_dspacing, MaxDSpacing=max_dspacing, OutputWorkspace = "PredictedPeaks")

#defining rotation functions to optimise
def rotation(theta, vertical_dir):
    axis = np.array(vertical_dir) / np.linalg.norm(vertical_dir)
    rotation_object = R.from_rotvec(theta * axis)
    return rotation_object


def objective(angle, vertical_dir, target_vector, initial_vector):
    rotation_object = rotation(angle, vertical_dir)
    rotated_vector = rotation_object.apply(initial_vector)
    cost = (np.linalg.norm(np.array(target_vector) - rotated_vector)) ** 2
    return cost

#getting the UB that's been calculated
OrientedLattice = SingleCrystalPeakTable.sample().getOrientedLattice()

#getting up and horizontal for rotations
Inst = SingleCrystalPeakTable.getInstrument()
lab_orientation = Inst.getReferenceFrame()
up = lab_orientation.vecPointingUp()
beam_direction = np.array(OrientedLattice.getuVector())
beam_dir_q = np.array(OrientedLattice.qFromHKL(beam_direction))
horizontal = np.cross(beam_dir_q, np.array(up))
UB_unaligned = OrientedLattice.getUB()
#getting initial uv before re alignment
u_unaligned = OrientedLattice.getuVector()
print(u_unaligned)
u_max = max(abs(np.array(u_unaligned)))
v_unaligned = OrientedLattice.getvVector()
v_max = max(abs(np.array(v_unaligned)))

#getting peaks used for UB calculations
used_peaks = ""
for x in range(len(used_hkls)):
    used_peaks += f"{list(used_hkls.values())[x]} for Run Number: {list(used_hkls.keys())[x]} \n"

#date to be used in filename
ct = datetime.now()
ct = ct.strftime('%Y-%m-%d %H:%M:%S')
print(ct)
#performing the reorientation via optimisation
if Reorient == True:
    hkl_q=np.array(OrientedLattice.qFromHKL([hkl_of_desired_alignment[0], hkl_of_desired_alignment[1], hkl_of_desired_alignment[2]]))

    result = minimize(objective, 0, args=(np.array(up), beam_dir_q, hkl_q), bounds=[(-np.pi, np.pi)],
                      method="Nelder-Mead")
    theta_new = result.x[0]
    theta_user = round(np.degrees(theta_new),2)
    rotation_object = R.from_rotvec(theta_new * np.array(up) / np.linalg.norm(np.array(up)))
    rotation_matrix = rotation_object.as_matrix()
    UB_final = rotation_matrix @ UB_unaligned
    temporary_workspace = CloneWorkspace("SingleCrystalPeakTable", "SingleCrystalPeakTable_temp", EnableLogging=False)
    SetUB(Workspace =temporary_workspace,UB=UB_final, EnableLogging=False)
    temporary_workspace = mtd["temporary_workspace"]
    OrientedLattice_aligned= temporary_workspace.sample().getOrientedLattice()
    u_aligned = OrientedLattice_aligned.getuVector()
    v_aligned = OrientedLattice_aligned.getvVector()
    #setting new UB on a different workspace to get u,v

    u_unaligned = u_unaligned / u_max
    v_unaligned = v_unaligned / v_max

    u_al_max = max(abs(np.array(u_aligned)))
    v_al_max = max(abs(np.array(v_aligned)))

    u_aligned = u_aligned/u_al_max
    v_aligned = v_aligned/v_al_max

    #getting phi rotation
    if len(gonio_axes) > 1:
        match_1 = re.findall(r"\d", Axes["Axis1"])
        match_1.pop()
        Axis1 = np.array([float(x) for x in match_1])
    else:
        Axis1 = np.array(horizontal)
    beam_dir_new_q = np.array(OrientedLattice_aligned.qFromHKL(u_aligned))
    print(beam_dir_new_q)
    hkl_new_q = np.array(OrientedLattice_aligned.qFromHKL([hkl_of_desired_alignment[0], hkl_of_desired_alignment[1], hkl_of_desired_alignment[2]]))
    print(hkl_new_q)
    result_phi = minimize(objective, 0, args=(Axis1, beam_dir_new_q, hkl_new_q), bounds=[(-np.pi/2, np.pi/2)],
                      method="Nelder-Mead")
    Phi = round(np.degrees(result_phi.x[0]),4)
    #writing the file
    results = f"{ct}\n UB Calculation Method: {UB_method}\n\n" \
              f"Used Peaks: \n{used_peaks}\n" \
              f"UB From Peak Data:\n {np.round(UB_unaligned,4)},\n\n" \
              f" u vector: {np.round(u_unaligned,4)},\n" \
              f" v vector:{np.round(v_unaligned,4)},\n\n" \
              f"To align q = {hkl_of_desired_alignment} along beam direction:\n\n" \
              f"Rotate to Psi0 {theta_user},\n\n" \
              f" Using {Axis1} as second Axis \n\n" \
              f"The offset from the plane of the provided q is: {Phi} \n\n" \
              f"UB After Psi Rotation:\n {np.round(UB_final,4)},\n\n" \
              f" u vector: {np.round(u_aligned,4)},\n" \
              f" v_vector: {np.round(v_aligned,4)} "
    print(results)
    file_loc = savedir + "DG_Alignment_output" + str(ct)
    try:
        data_file = open(file_loc, "x")
    except FileNotFoundError:
        print('Could not write output data')
    else:
        data_file.write(results)
        data_file.close()

    names = ["SingleCrystalPeakTable", "temporary_workspace"]
    WorkspaceNames = ["Processed Peak Table Unaligned", "Processed Peak Table Aligned"]

elif Reorient == False:
    u_unaligned = u_unaligned / u_max
    v_unaligned = v_unaligned / v_max
    results = f"{ct}\n\n UB Calculation Method: {UB_method}\n\nUsed Peaks: \n{used_peaks}\nUB From Peak Data:\n {UB_unaligned},\n\n u vector: {u_unaligned},\n v vector:{v_unaligned}"
    print(results)
    file_loc = datadir + "DG_Alignment_output" + str(ct)
    try:
        data_file = open(file_loc, "x")
    except FileNotFoundError:
        print('Could not write output data')
    else:
        data_file.write(results)
        data_file.close()

    names = ["SingleCrystalPeakTable"]
    WorkspaceNames = ["Processed Peak Table"]

RenameWorkspaces(names, WorkspaceNames)

if Predict_peaks == True:
    InstrumentWidget = import_qt("._instrumentview", "mantidqt.widgets.instrumentview", "InstrumentWidget")
    peaks = mtd["PredictedPeaks"]
    ws = Load(f'{datadir}/{instrument_name[:3]}{run_numbers[x]}', EnableLogging=False)
    SetGoniometer(ws, **gonio_axes, EnableLogging=False)
    ws = ConvertUnits(ws, Target='dSpacing', EMode='Elastic', EFixed=None,  EnableLogging=False)
    iv = get_instrumentview(ws)
    iv.container.widget.overlay(peaks.name())
    bin_range_lower = min_dspacing
    bin_range_upper = max_dspacing
    iv.set_bin_range(bin_range_lower, bin_range_upper)

    iv.show_view()
