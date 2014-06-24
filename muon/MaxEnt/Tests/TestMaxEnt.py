import unittest
import os

from mantid.simpleapi import *

SCRIPT_DIR = str(os.path.dirname(os.path.realpath(__file__)))

class TestMaxEntMUSR(unittest.TestCase):

    def setUp(self):
        self.max_ent_alg = AlgorithmManager.createUnmanaged('MaxEnt')
        self.max_ent_alg.initialize()
        self.max_ent_alg.setChild(True)
        self.max_ent_alg.setProperty("RunNumber",        "48033")
        self.max_ent_alg.setProperty("Instrument",       "MUSR")
        self.max_ent_alg.setProperty("PhasesInFile",     os.path.join(SCRIPT_DIR, "input_phase.dat"))
        self.max_ent_alg.setProperty("DeadTimesInFile",  os.path.join(SCRIPT_DIR, "input_taud.dat"))
        self.max_ent_alg.setProperty("PhasesOutFile",    os.path.join(SCRIPT_DIR, "output_phase.dat"))
        self.max_ent_alg.setProperty("DeadTimesOutFile", os.path.join(SCRIPT_DIR, "output_taud.dat"))

    def tearDown(self):
        os.remove(os.path.join(SCRIPT_DIR, "output_phase.dat"))
        os.remove(os.path.join(SCRIPT_DIR, "output_taud.dat"))

    def check_actual_ws_and_expected_file_match(self, actual_ws_name, expected_file_name):
        actual_ws = mtd[actual_ws_name]
        expected_ws = LoadNexus(os.path.join(SCRIPT_DIR, expected_file_name))
        result = CheckWorkspacesMatch(Workspace1=actual_ws, Workspace2=expected_ws)
        self.assertEquals(result, "Success!", "\"%s\"" % result)

    def test_fit_and_fix_checked(self):
        self.max_ent_alg.setProperty("FitDeadTimes", True)
        self.max_ent_alg.setProperty("FixPhases", True)
        self.max_ent_alg.execute()

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt",
            "musr_fit_and_fix_checked_result.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Input",
            "musr_fit_and_fix_checked_input.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Deadtimes & Phases",
            "musr_fit_and_fix_checked_deadtimes_and_phases.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Chi^2",
            "musr_fit_and_fix_checked_chi_squared.nxs")

    def test_fit_checked(self):
        self.max_ent_alg.setProperty("FitDeadTimes", True)
        self.max_ent_alg.setProperty("FixPhases", False)
        self.max_ent_alg.execute()

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt",
            "musr_fit_checked_result.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Input",
            "musr_fit_checked_input.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Deadtimes & Phases",
            "musr_fit_checked_deadtimes_and_phases.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Chi^2",
            "musr_fit_checked_chi_squared.nxs")

    def test_fix_checked(self):
        self.max_ent_alg.setProperty("FitDeadTimes", False)
        self.max_ent_alg.setProperty("FixPhases", True)
        self.max_ent_alg.execute()

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt",
            "musr_fix_checked_result.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Input",
            "musr_fix_checked_input.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Deadtimes & Phases",
            "musr_fix_checked_deadtimes_and_phases.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Chi^2",
            "musr_fix_checked_chi_squared.nxs")

    def test_none_checked(self):
        self.max_ent_alg.setProperty("FitDeadTimes", False)
        self.max_ent_alg.setProperty("FixPhases", False)
        self.max_ent_alg.execute()

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt",
            "musr_none_checked_result.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Input",
            "musr_none_checked_input.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Deadtimes & Phases",
            "musr_none_checked_deadtimes_and_phases.nxs")

        self.check_actual_ws_and_expected_file_match(
            "MUSR00048033; MaxEnt Chi^2",
            "musr_none_checked_chi_squared.nxs")

if __name__ == '__main__':
    unittest.main()