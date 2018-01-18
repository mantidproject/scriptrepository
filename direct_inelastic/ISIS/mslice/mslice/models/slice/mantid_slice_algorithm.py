from __future__ import (absolute_import, division, print_function)
import numpy as np

from mantid.simpleapi import BinMD, LoadCIF
from mantid.api import IMDEventWorkspace
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
from scipy import constants

from .slice_algorithm import SliceAlgorithm
from mslice.models.alg_workspace_ops import AlgWorkspaceOps
from mslice.models.workspacemanager.mantid_workspace_provider import MantidWorkspaceProvider

KB_MEV = constants.value('Boltzmann constant in eV/K') * 1000
HBAR_MEV = constants.value('Planck constant over 2 pi in eV s') * 1000
E_TO_K = np.sqrt(2*constants.neutron_mass)/HBAR_MEV
E2L = 1.e23 * constants.h**2 / (2 * constants.m_n * constants.e)  # energy to wavelength conversion E = h^2/(2*m_n*l^2)
crystal_structure = {'Copper': ['3.6149 3.6149 3.6149', 'F m -3 m', 'Cu 0 0 0 1.0 0.05'],
                     'Aluminium': ['4.0495 4.0495 4.0495', 'F m -3 m', 'Al 0 0 0 1.0 0.05'],
                     'Niobium': ['3.3004 3.3004 3.3004', 'I m -3 m', 'Nb 0 0 0 1.0 0.05'],
                     'Tantalum': ['3.3013 3.3013 3.3013', 'I m -3 m', 'Ta 0 0 0 1.0 0.05']}


class MantidSliceAlgorithm(AlgWorkspaceOps, SliceAlgorithm):
    def __init__(self):
        self._workspace_provider = MantidWorkspaceProvider()

    def compute_slice(self, selected_workspace, x_axis, y_axis, smoothing, norm_to_one):
        workspace = self._workspace_provider.get_workspace_handle(selected_workspace)
        assert isinstance(workspace, IMDEventWorkspace)
        self._fill_in_missing_input(x_axis, workspace)
        self._fill_in_missing_input(y_axis, workspace)
        n_x_bins = self._get_number_of_steps(x_axis)
        n_y_bins = self._get_number_of_steps(y_axis)
        x_dim_id = workspace.getDimensionIndexByName(x_axis.units)
        y_dim_id = workspace.getDimensionIndexByName(y_axis.units)
        x_dim = workspace.getDimension(x_dim_id)
        y_dim = workspace.getDimension(y_dim_id)
        xbinning = x_dim.getName() + "," + str(x_axis.start) + "," + str(x_axis.end) + "," + str(n_x_bins)
        ybinning = y_dim.getName() + "," + str(y_axis.start) + "," + str(y_axis.end) + "," + str(n_y_bins)
        thisslice = BinMD(InputWorkspace=workspace, AxisAligned="1", AlignedDim0=xbinning, AlignedDim1=ybinning)
        # perform number of events normalization
        with np.errstate(invalid='ignore'):
            plot_data = thisslice.getSignalArray() / thisslice.getNumEventsArray()
        # rot90 switches the x and y axis to to plot what user expected.
        plot_data = np.rot90(plot_data)
        self._workspace_provider.delete_workspace(thisslice)
        boundaries = [x_axis.start, x_axis.end, y_axis.start, y_axis.end]
        if norm_to_one:
            plot_data = self._norm_to_one(plot_data)
        plot = [plot_data, None, None, None, None, None]
        return plot, boundaries

    def compute_boltzmann_dist(self, sample_temp, e_axis):
        '''calculates exp(-E/kBT), a common factor in intensity corrections'''
        if sample_temp is None:
            return None
        kBT = sample_temp * KB_MEV
        energy_transfer = np.linspace(e_axis.end, e_axis.start, self._get_number_of_steps(e_axis))
        return np.exp(-energy_transfer / kBT)

    def compute_chi(self, scattering_data, boltzmann_dist, e_axis):
        energy_transfer = np.linspace(e_axis.end, e_axis.start, self._get_number_of_steps(e_axis))
        signs = np.sign(energy_transfer)
        signs[signs == 0] = 1
        chi = (signs + (boltzmann_dist * -signs))[:, None]
        chi = np.pi * chi * scattering_data
        return chi

    def compute_chi_magnetic(self, chi):
        if chi is None:
            return None
        # 291 milibarns is the total neutron cross-section for a moment of one bohr magneton
        chi_magnetic = chi / 291
        return chi_magnetic

    def compute_d2sigma(self, scattering_data, workspace, e_axis):
        Ei = self._workspace_provider.get_EFixed(self._workspace_provider.get_parent_by_name(workspace))
        if Ei is None:
            return None
        ki = np.sqrt(Ei) * E_TO_K
        energy_transfer = np.linspace(e_axis.end, e_axis.start, self._get_number_of_steps(e_axis))
        kf = (np.sqrt(Ei - energy_transfer)*E_TO_K)[:, None]
        d2sigma = scattering_data * kf / ki
        return d2sigma

    def compute_symmetrised(self, scattering_data, boltzmann_dist, e_axis):
        energy_transfer = np.arange(e_axis.end, 0, -e_axis.step)
        negatives = scattering_data[len(energy_transfer):] * boltzmann_dist[len(energy_transfer):,None]
        return np.concatenate((scattering_data[:len(energy_transfer)], negatives))

    def compute_gdos(self, scattering_data, boltzmann_dist, q_axis, e_axis):
        energy_transfer = np.linspace(e_axis.end, e_axis.start, self._get_number_of_steps(e_axis))
        momentum_transfer = np.linspace(q_axis.start, q_axis.end, self._get_number_of_steps(q_axis))
        momentum_transfer = np.square(momentum_transfer, out=momentum_transfer)
        gdos = scattering_data / momentum_transfer
        gdos *= energy_transfer[:,None]
        gdos *= (1 - boltzmann_dist)[:,None]
        return gdos

    def sample_temperature(self, ws_name, sample_temp_fields):
        ws = self._workspace_provider.get_parent_by_name(ws_name)
        # mantid drops log data during projection, need unprojected workspace.
        sample_temp = None
        for field_name in sample_temp_fields:
            try:
                sample_temp = ws.run().getLogData(field_name).value
            except RuntimeError:
                pass
        try:
            float(sample_temp)
        except (ValueError, TypeError):
            pass
        else:
            return sample_temp
        if isinstance(sample_temp, str):
            sample_temp = self.get_sample_temperature_from_string(sample_temp)
        if isinstance(sample_temp, np.ndarray) or isinstance(sample_temp, list):
            sample_temp = np.mean(sample_temp)
        return sample_temp

    def get_sample_temperature_from_string(self, string):
        pos_k = string.find('K')
        if pos_k == -1:
            return None
        k_string = string[pos_k - 3:pos_k]
        sample_temp = float(''.join(c for c in k_string if c.isdigit()))
        return sample_temp

    def compute_recoil_line(self, axis, relative_mass=1):
        momentum_transfer = np.arange(axis.start, axis.end, axis.step)
        line = np.square(momentum_transfer * 1.e10 * constants.hbar) / (2 * relative_mass * constants.neutron_mass) /\
            (constants.elementary_charge / 1000)
        return momentum_transfer, line

    def compute_powder_line(self, ws_name, axis, element, cif_file=False):
        efixed = self._workspace_provider.get_EFixed(self._workspace_provider.get_parent_by_name(ws_name))
        if axis.units == 'MomentumTransfer':
            x0 = self._compute_powder_line_momentum(ws_name, axis, element, cif_file)
        elif axis.units == 'Degrees':
            x0 = self._compute_powder_line_degrees(ws_name, axis, element, efixed, cif_file)
        else:
            raise RuntimeError("units of axis not recognised")
        x = sum([[xv, xv, np.nan] for xv in x0], [])
        y = sum([[efixed / 20,  -efixed / 20, np.nan] for xv in x0], [])
        return x, y

    def _compute_powder_line_momentum(self, ws_name, q_axis, element, cif_file):
        d_min = (2 * np.pi) / q_axis.end
        d_max = (2 * np.pi) / q_axis.start
        structure = self._crystal_structure(ws_name, element, cif_file)
        dvalues = self.compute_dvalues(d_min, d_max, structure)
        x = (2 * np.pi) / dvalues
        return x

    def _crystal_structure(self, ws_name, element, cif_file):
        if cif_file:
            ws = self._workspace_provider.get_parent_by_name(ws_name)
            LoadCIF(ws, cif_file)
            return ws.sample().getCrystalStructure()
        else:
            return CrystalStructure(crystal_structure[element][0], crystal_structure[element][1],
                                    crystal_structure[element][2])

    def _compute_powder_line_degrees(self, ws_name, theta_axis, element, efixed, cif_file):
        wavelength = np.sqrt(E2L / efixed)
        d_min = wavelength / (2 * np.sin(np.deg2rad(theta_axis.end / 2)))
        d_max = wavelength / (2 * np.sin(np.deg2rad(theta_axis.start / 2)))
        structure = self._crystal_structure(ws_name, element, cif_file)
        dvalues = self.compute_dvalues(d_min, d_max, structure)
        x = 2 * np.arcsin(wavelength / 2 / dvalues) * 180 / np.pi
        return x

    def compute_dvalues(self, d_min, d_max, structure):
        generator = ReflectionGenerator(structure)
        hkls = generator.getUniqueHKLsUsingFilter(d_min, d_max, ReflectionConditionFilter.StructureFactor)
        dvalues = np.sort(np.array(generator.getDValues(hkls)))[::-1]
        return dvalues

    def _norm_to_one(self, data):
        data_range = data.max() - data.min()
        return (data - data.min())/data_range
