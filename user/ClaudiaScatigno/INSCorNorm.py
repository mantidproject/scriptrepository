##### INSCorNorm (PythonAlgorithm) Author: Claudia Scatigno <claudia.scatigno@cref.it>
#####  How to cite: Scatigno, Claudia, et al. "A python algorithm to analyze inelastic neutron scattering spectra based on the y-scale formalism." Journal of Chemical Theory and Computation 16.12 (2020): 7671-7680

import mantid
import numpy as np
import scipy.special as sps
np.set_printoptions(precision=4, linewidth=150)

# The value of the hydrogen bound scattering cross section to be used for the normalization and for the evaluation of the energy-dependent Cross Section for Free Hydrogen Gas or Capelli Romanelli models
hydrogen_XS = 82.03  # barn
# instrument specific
# The final energy of the INS spectrometer, TOSCA in this case
final_energy = 3.51  # meV
# The Avogadro constant
avog_const = 6.022E+23 # mol^-1
# The forward and backward scattering angles of the INS spectrometer in degrees.
scattering_angles = [42.55, 137.73]
# Conversion to the cosine values of the angles, used later
scattering_cost = np.abs(np.cos(np.array(scattering_angles) / 180. * np.pi))
# The conversion factor from cm^-1 to meV
inv_cm_to_meV = 0.12398
# The two values of the secant for F(orward) and B(ackward) angles
secF=1./(np.cos(np.radians(scattering_angles[0])))
secB=1./(np.cos(np.radians(scattering_angles[1])))


# Function to import the data and to convert them in meV if they are in cm-1
def check_and_import_data(data_name):
    # Load the original data
    data = mtd[data_name]
    # Check if the data are in cm-1
    if str(data.getAxis(0).getUnit().symbol()) == "cm^-1":
        input_in_inv_cm = 'True'
        # conversion from cm-1 to meV
        ConvertUnits(InputWorkspace=data_name, OutputWorkspace=data_name + '_scaled', Target='DeltaE', EMode='Indirect', EFixed=final_energy, AlignBins=True)
        # conversion from energy transfer to incident energy
        ScaleX(InputWorkspace=data_name + "_scaled", OutputWorkspace=data_name + '_scaled', Operation="Add", Factor=final_energy)
    # If the data are in meV, add the final energy
    elif str(data.getAxis(0).getUnit().symbol()) == "meV":
        input_in_inv_cm = 'False'
        # conversion from energy transfer to incident energy
        ScaleX(InputWorkspace=data_name, OutputWorkspace=data_name + '_scaled', Operation="Add", Factor=final_energy)
    else:
        print("The units of the workspace are neither in cm^-1 nor in meV!")
        return "Data should be provided as energy transfer either in meV or in cm^-1"
    return data, input_in_inv_cm
 # The algorithm INSCorNorm
class INSCorNorm(PythonAlgorithm):
    
    # Brief explanation of the algorithm
    def summary(self):
        return ("""This algorithm corrects INS spectra by empty container and sample self-shielding.
Total cross section can be provided as experimental Transmission data or it can be modeled with two different methods.
In addition, once the hydrogen mean kinetic energy is provided, the algorithm normalizes the data using the Y-scaling formalism.
The description of the corrections and normalization procedures can be found in Colognesi J. Neutron Research 19 (2017).
Details about this algorithm are in Scatigno, Claudia, et al. "A python algorithm to analyze inelastic neutron scattering spectra based on the y-scale formalism." Journal of Chemical Theory and Computation 16.12 (2020): 7671-7680.""")
 
    def category(self):
        return 'MyTools'
    
# Workspaces and parameters that must be provided by the user for the execution of the algorithm
    def PyInit(self):
        
    # INS SPECTRA INPUT 
    # The Workspace containing the INS data of the sample
       self.declareProperty(WorkspaceProperty('INSSampleWorkspace', '', direction=Direction.Input), 
       doc='The sample INS workspace to correct and/or normalize.')
       
    # CAN CORRECTION INPUTS
    # If the sample was in a can, select this option
       self.declareProperty("Sample in Can?", True, doc="If there is a can to subtract or the sample was in a can.")
    # The transmission of the empty can, a value ranging from 0 to 1
       self.declareProperty('EmptyCanTransmission', "0.99", direction=Direction.Input, doc="The transmission value of the empty can. 1-mm-thick Al containers correspond to 0.99")
    # The Workspace containing the INS data of the can
       self.declareProperty(WorkspaceProperty('INS Can Workspace', '', direction=Direction.Input), 
       doc='The can INS workspace to subtract.')
 
    # SELF-SHIELDING CORRECTION INPUTS
       self.declareProperty("Cross Section", "From Transmission Data", StringListValidator(["From Transmission Data", "Model: Free Hydrogen Gas", "Model: Capelli Romanelli"]),
       doc="""Method to handle the cross section:\n"From Transmission Data" uses an 
       experimental sample transmission spectrum;\n"Model: Free Hydrogen Gas" uses the hydrogen cross 
       section from a free gas model;\n"Capelli Romanelli" uses the average hydrogen cross section for CH, NH, 
       and OH from J. Appl. Cryst. (2019). 52, 1233?1237""")

    # If the option 'From Transmission Data' is selected, the user must select the Workspace containing the data
      self.declareProperty(WorkspaceProperty('Transmission Workspace', '', direction=Direction.Input), 
      doc="The Transmission workspace to use" if cross_section == "From Transmission Data")
    # Numbers of Hydrogen atoms per unit cell
      self.declareProperty('Number of Hydrogen atoms', "2", direction=Direction),Input, doc="The number of Hydrogen Atoms from Stoichiometry.")
    # The thickness of the sample, expressed in cm
      self.declareProperty('Sample thickness', "0.1", direction=Direction), Input, doc="The sample thickness in cm, used if Cross Section must be calculated from the models.")
    # Density of the sample experimentally measured
      self.declareProperty('Sample density', "3.35", direction=Direction), Input, doc="The experimental density of the sample in g/cm^3, given by the mass of the sample and the volume of the can.")
    # Facot rof the inner volume of the can filled with the sample. A value ranging from 0 to 1
      self.declareProperty('Can filling factor', "1.00", direction=Direction),Input, doc="The filling factor of the can with the sample, given by the percentage of can volume filled with the sample (always <= 1)")
     # The molar mass of the sample
    self.declareProperty('Sample molar mass', "93.24", direction=Direction.Input, doc="The molar mass of the sample in g/mol")
     # The proton Mean Kinetic Energy in meV
    self.declareProperty('Proton Kinetic Energy', "110", direction=Direction.Input, doc="The Hydrogen Mean Kinetic Energy in meV,either from DINS measurements or PhonDOS.")
    # If selected, this option allow the normalization of the data
    self.declareProperty("Normalize", True, doc="Wether to normalize the corrected data using the Yscaling.")
    # Numbers of Hydrogen atoms per unit cell
    self.declareProperty('Hydrogens Per Unit Cell', "2", direction=Direction.Input, doc="The number of Hydrogen Atoms per unit cell.")
    # The unit of measurement of the energy region used for the normalization
    self.declareProperty("Normalization Region Units", "cm-1", StringListValidator(["cm-1", "meV"]))
    # The lower bound of the energies for the normalization
    self.declareProperty('Normalization Region Min', "6000", direction=Direction.Input, doc="Lower bound of the energy region where the spectrum is normalized.")
    # The upper bound of the energies for the normalization
    self.declareProperty('Normalization Region Max', "8000", direction=Direction.Input, doc="Upper bound of the energy region where the spectrum is normalized.")
 
    # If selected, this option makes the script producing an output in cm^-1 instead of meV
    self.declareProperty("Output in cm-1", True, doc="Output in cm^-1 or meV")
    def PyExec(self):
     
     
 # INITIALIZATION
    # All corrections and evaluations are performed in meV. Workspaces provided in cm^-1 will be 
    converted and processed in meV.
    # The string names of the INS workspaces for:
    # - The sample
    global suffix, final_energy_bin, hydrogens_per_unit_cell, Emin, Emax ins_name = str(self.getProperty("INS Sample Workspace").value)
    # - The can
    ins_can_name = str(self.getProperty("INS Can Workspace").value)
    # The sample INS workspaces is first converted to meV (if it was in cm^-1), 
    # then transformed from energy transfer to incident neutron energy 
    check_and_import_data(ins_name)
    ws_ins = mtd[ins_name + '_scaled']
    # The same procedure is made for the can INS workspace
    check_and_import_data(ins_can_name)
    ws_can = mtd[ins_can_name + '_scaled']
    # If the sample is inside in a can the value (non energy dependent) of the transmission of the empty 
    can is collected
    incan = self.getProperty("Sample in Can?").value
    if incan:
    empty_can_transmission = float(self.getProperty("EmptyCanTransmission").value)
    # Otherwise, the transmission is 1 
    else:
    empty_can_transmission = 1.00
    # Depending on the user's choice, the cross section is obtained from measurements or from one of the 
    two models proposed 
    cross_section = self.getProperty("Cross Section").value
 
    # Cross section from transmission data
    if cross_section == "From Transmission Data":
    trans_name = str(self.getProperty("Transmission Workspace").value)
    # Rebin the transmission + cross section according to the INS workspace binning
    RebinToWorkspace(WorkspaceToRebin=trans_name,
    WorkspaceToMatch=ins_name + '_scaled', OutputWorkspace=trans_name + '_reb')
    # The transmission contributions from experimental data
    ws_trans = mtd[trans_name + '_reb']
    # The suffix to indicate transmission from measurements
    suffix = '_meas_trans'
    # The workspace with the transmission data is created
    CreateWorkspace(OutputWorkspace='Trans_measured_' + ins_name, DataX=ws_trans.readX(0), 
    DataY=ws_trans.readY(0), UnitX='meV', VerticalAxisUnit='Text', VerticalAxisValues='1', 
    YUnitLabel='Transmission', WorkspaceTitle='Transmission', Distribution=True)
 
    # Cross section from models
    else:
    # Take the values of the parameters necessary to the theoretical models to produce the cross 
    sections value
    number_hydrogens = float(self.getProperty("Number of Hydrogen atoms").value)
    sample_thickness = float(self.getProperty("Sample thickness").value)
    sample_density = float(self.getProperty("Sample density").value)
    filling_factor = float(self.getProperty("Can filling factor").value)
    sample_molar_mass = float(self.getProperty("Sample molar mass").value)
    effective_thickness = sample_density * avog_const * number_hydrogens * sample_thickness / sample_molar_mass / filling_factor
    # The scaled data are taken from the workspace
    data_x = mtd[ins_name + "_scaled"].dataX(0)
    data_y = mtd[ins_name + "_scaled"].dataY(0)
    # Cross section from the Free Hydrogen Gas model
    if cross_section == "Model: Free Hydrogen Gas":
    # From Eq. 41 in Colognesi and from Turchin, Valentin Fedorovich. Slow neutrons. Vol. 2112. Israel 
    Program for Scientific Translations, 1965.
    hydrogen_kinetic_energy = float(self.getProperty("Proton Kinetic Energy").value) xn = np.sqrt(3. * data_x / 2. / hydrogen_kinetic_energy) model_XS = hydrogen_XS / 4. / xn * ((xn + 1. / 2. / xn) * sps.erf(xn) + np.pi ** (-0.5) * np.exp(-xn ** 2))
    # The information about the model used are saved for naming the ouput
    trans_model_name = 'Trans_model_FHG_' + ins_name
    suffix = '_FHG_trans'
    Cross section from Capelli Romanelli model
    elif cross_section == "Model: Capelli Romanelli":
    # Average hydrogen XS from Capelli, Romanelli, J. Appl. Cryst. (2019). 52, 1233-1237
    a, b, x0, y0 = 73.4, -1.99, 1.86, 19.79
    # conversion from energy in meV to wavelength in Angstrom to obtain the cross sections according to the Capelli Romanelli model
    x = np.sqrt(81.8047 / data_x)
    model_XS = y0 + a / (1. + (x / x0) ** b)
    # The information about the model used are saved for naming the ouput
    trans_model_name = 'Trans_model_CR_' + ins_name
    suffix = '_CR_trans'
 
    # The transmission data from model are populated
    for bin_num in range(len(data_y)):
    # Parameter in barn^-1 (= cm^-24) when [sample_thickness] = cm; [unit_cell_volume] = A^3 (= cm^24)
    meanXS = (model_XS[bin_num] + model_XS[bin_num + 1]) / 2.
    data_y[bin_num] = np.exp(- effective_thickness * meanXS * 1E-24)
    # The data of the transmission from models are saved in a workspace 
    ws_trans = CreateWorkspace(OutputWorkspace=trans_model_name, DataX=data_x, DataY=data_y, 
    UnitX='meV', VerticalAxisUnit='Text', VerticalAxisValues='1', YUnitLabel='Transmission', 
    WorkspaceTitle='Transmission', Distribution=True)
    # A new worskapce with the data, corrected with the transmission data (from models or measurement) is created
    CloneWorkspace(InputWorkspace=ins_name + '_scaled', OutputWorkspace=ins_name + suffix + '_corrected')
    ws = mtd[ins_name + suffix + '_corrected']
 
    # The bin corresponging to the final energy is found
    for k in range(ws_trans.blocksize() - 1): 
    if (ws_trans.readX(0)[k] <= final_energy and ws_trans.readX(0)[k + 1] > final_energy):
    final_energy_bin = k
 
    # If the user requires the normalization this part of the script is executed
    normalize = self.getProperty("Normalize").value 
    # The transmission of the first wall of the can is obtained
    incident_wall_transmission = np.sqrt(empty_can_transmission)
    # This is the transmission values obtained from the measurement or one of the two models
    incident_sample_transmission = ws_trans.readY(0)
    for direction in range(2): # runs over forward and backward directions
    # Empty can scattered wall transmission
    scattered_wall_transmission = np.exp(np.log(incident_wall_transmission) / scattering_cost[direction])
    # Corrected contribution from empty can
    empty_can_corrected = ws_can.readY(direction) / (incident_wall_transmission + scattered_wall_transmission)
    # Sample corrections
    scattered_sample_transmission = np.exp(np.log(ws_trans.readY(0)[final_energy_bin]) / scattering_cost[direction])
    # The ample self shielding is initialized
    sample_self_shielding = 0. * incident_sample_transmission + 1.
    for bin_num in range(ws.blocksize()):
    if incident_sample_transmission[bin_num] <= scattered_sample_transmission: sample_self_shielding[bin_num] = 1.
    else:
    # For forward detectors
    if direction==0:
    FTheta= (np.exp(-hydrogen_XS*sample_thickness)-np.exp(-hydrogen_XS*sample_thickness*secF))/(hydrogen_XS*secF-hydrogen_XS)
    # For backward detectors
    else:
    FTheta= (np.exp(-hydrogen_XS*sample_thickness)-np.exp(-hydrogen_XS*sample_thickness*secB))/(hydrogen_XS*secB-hydrogen_XS)
    FTheta*= -np.exp(hydrogen_XS*sample_thickness*secB)
 sample_self_shielding[bin_num] = incident_sample_transmission[bin_num] -FTheta*scattered_sample_transmission
 sample_self_shielding[bin_num] /= np.log(incident_sample_transmission[bin_num]) -np.log(scattered_sample_transmission)
    if incan:
    ws.dataY(direction)[bin_num] = ws_ins.readY(direction)[bin_num] -empty_can_corrected[bin_num] * (incident_sample_transmission[bin_num] * incident_wall_transmission + scattered_sample_transmission * scattered_wall_transmission)
    ws.dataY(direction)[bin_num] /= incident_wall_transmission * sample_self_shielding[bin_num] * scattered_wall_transmission
 
    for bin_num in range(ws.blocksize()):
    ws.dataY(2)[bin_num] = (ws.readY(0)[bin_num] + ws.dataY(1)[bin_num]) / 2.
    ws.dataE(2)[bin_num] = np.sqrt(ws.readE(0)[bin_num] ** 2 + ws.dataE(1)[bin_num] ** 2) / 2.
    # Normalization for the self scatterting contribution following the procedure and the aproximation 
    made in [1]
    if normalize:
    # A new workspace is created
    CloneWorkspace(InputWorkspace=ins_name + suffix + '_corrected', OutputWorkspace=ins_name + suffix + '_corrected_normalized')
    ws_nor = mtd[ins_name + suffix + '_corrected_normalized']
 
    # The parameters are provided by the user
    proton_kinetic_energy = float(self.getProperty("Proton Kinetic Energy").value)
    hydrogens_per_unit_cell = float(self.getProperty("Hydrogens Per Unit Cell").value)
    normalization_region_units = self.getProperty("Normalization Region Units").value)
    Emin = float(self.getProperty("Normalization Region Min").value)
    Emax = float(self.getProperty("Normalization Region Max").value)
    E = ws_nor.dataX(0)
 
    # Constants for the calculations
    hbar = 2.0446
    Mn, Mh = 1.0086, 1.0079 # Mn neutron mass, Mh hydrogen mass
    # Momentum transfer
    Q = np.sqrt(2. * Mn * (E + final_energy - 2. * np.sqrt(E * final_energy) * scattering_cost[direction])) / hbar # Eq. 2
    # y-scaling
    y = Mh / Q / hbar ** 2 * (E - final_energy - hbar ** 2 * Q ** 2 / 2. / Mh) # Eq. 21
    # Sigma momentum
    sigP = np.sqrt(2. / 3. * proton_kinetic_energy * Mh) / hbar # Eq. 26
    # First factor the formula for the self scattering 
    ss = Mh / hbar ** 2 / Q / (np.sqrt(2. * np.pi) * sigP) * np.exp(-y ** 2 / 2. / sigP ** 2) # Eq. 29
    # In the low-temperature approximation of the final state effects, the isotropically-approximated 
    coefficients of Eq.30 becomes (from Eq.35):
    b3 = np.sqrt(2.) / 12. * sigP
    b4 = sigP ** 2 / 24.
    b6 = sigP ** 2 / 144.
    x = y / np.sqrt(2.) / sigP
    # Hermite polynomials for the contributions in Eq. 29
    H3 = 8. * x ** 3 - 12. * x
    H4 = 16. * x ** 4 - 48. * x ** 2 + 12.
    H6 = 64. * x ** 6 - 480. * x ** 4 + 720. * x ** 2 - 120.
    # Final calculation of the self scattering Eq. 22
    ss *= (1. + b3 / Q * H3 + b4 / Q ** 2 * H4 + b6 / Q ** 2 * H6)
    # Contribution to the corrected measured spectrum due to the Hydrogen self scattering Eq. 4
    ss *= hydrogen_XS / 4. / np.pi * hydrogens_per_unit_cell spectrum_integral, normalization_integral = 0., 0.
 
    # The enerrgy range for normalization is taken from the user's choice
    if normalization_region_units != "meV": E /= inv_cm_to_meV
    for bin_num in range(ws_nor.blocksize()):
    if E[bin_num] > Emin and E[bin_num] < Emax: spectrum_integral += ws_nor.dataY(direction)[bin_num] normalization_integral += ss[bin_num]
 
    # The data are normalized for the normalization integral
    for bin_num in range(ws_nor.blocksize()):
    ws_nor.dataY(direction)[bin_num] = ws_nor.readY(direction)[bin_num] * normalization_integral / spectrum_integral
    ws_nor.dataE(direction)[bin_num] = ws_nor.readE(direction)[bin_num] * normalization_integral / spectrum_integral
 
    # The average spectrum of forward and backscattering data is calculated
    for bin_num in range(ws_nor.blocksize()):
    ws_nor.dataY(2)[bin_num] = (ws_nor.readY(0)[bin_num] + ws_nor.dataY(1)[bin_num]) / 2.
    ws_nor.dataE(2)[bin_num] = np.sqrt(ws_nor.readE(0)[bin_num] ** 2 + 
    ws_nor.dataE(1)[bin_num] ** 2) / 2.
    ScaleX(InputWorkspace=ins_name + suffix + '_corrected_normalized', 
    OutputWorkspace=ins_name + suffix + '_corrected_normalized', Operation="Add", Factor=-final_energy)
    output = self.getProperty("Output in cm-1").value
    if output:
    ConvertUnits(ins_name+suffix+'_corrected', OutputWorkspace=ins_name+suffix+'_corrected', 
    Target='DeltaE_inWavenumber', EMode='Indirect', EFixed=final_energy, AlignBins=True)
    if normalize:
    ConvertUnits(ins_name+suffix + '_corrected_normalized', OutputWorkspace=ins_name+suffix + '_corrected_normalized', Target='DeltaE_inWavenumber', EMode='Indirect', EFixed=final_energy, AlignBins=True)
    # Some useless workspace are deleted
    DeleteWorkspace(ins_name + '_scaled')
    DeleteWorkspace(ins_can_name + '_scaled')
    if cross_section == "From Transmission Data":
    DeleteWorkspace(trans_name + '_reb')
 
    AlgorithmFactory.subscribe(INSCorNorm)