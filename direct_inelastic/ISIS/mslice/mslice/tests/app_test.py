import mock
import unittest
from mock import patch
import PyQt4

from mslice.app.mainwindow import MainWindow
from mslice.app.mainwindow_ui import Ui_MainWindow
from mslice.presenters.powder_projection_presenter import PowderProjectionPresenter
from mslice.presenters.workspace_manager_presenter import WorkspaceManagerPresenter
from mslice.presenters.slice_plotter_presenter import SlicePlotterPresenter
from mslice.presenters.cut_presenter import CutPresenter
from mslice.views.mainview import MainView
from mslice.views.powder_projection_view import PowderView
from mslice.views.slice_plotter_view import SlicePlotterView
from mslice.views.cut_view import CutView
from mslice.views.workspace_view import WorkspaceView
from mslice.models.projection.powder.projection_calculator import ProjectionCalculator
from mslice.models.slice.slice_plotter import SlicePlotter
from mslice.models.workspacemanager.workspace_provider import WorkspaceProvider
from mslice.models.cut.cut_algorithm import CutAlgorithm
from mslice.models.cut.cut_plotter import CutPlotter

mainview = mock.create_autospec(MainView)
slice_view = mock.create_autospec(SlicePlotterView)
powder_view = mock.create_autospec(PowderView)
cut_view = mock.create_autospec(CutView)
workspace_view = mock.create_autospec(WorkspaceView)

class AppTests(unittest.TestCase):
    def setUp(self):
        self.workspace_provider = mock.create_autospec(WorkspaceProvider)
        self.slice_plotter = mock.create_autospec(SlicePlotter)
        self.projection_calculator = mock.create_autospec(ProjectionCalculator)
        self.cut_algorithm = mock.create_autospec(CutAlgorithm)
        self.cut_plotter = mock.create_autospec(CutPlotter)
        # using globals is a bit ugly but I can't think of another way to access the required variables
        # because they are only defined in the scope of MainWindow / setupUI
        global workspace_view
        global slice_view
        global powder_view
        global cut_view
        self.slice_presenter = SlicePlotterPresenter(slice_view, self.slice_plotter)
        self.workspace_presenter = WorkspaceManagerPresenter(workspace_view, self.workspace_provider)
        self.powder_presenter = PowderProjectionPresenter(powder_view, self.projection_calculator)
        self.cut_presenter = CutPresenter(cut_view, self.cut_algorithm, self.cut_plotter)
        workspace_view.get_presenter = mock.Mock(return_value=self.workspace_presenter)
        slice_view.get_presenter = mock.Mock(return_value=self.slice_presenter)
        powder_view.get_presenter = mock.Mock(return_value=self.powder_presenter)
        cut_view.get_presenter = mock.Mock(return_value=self.cut_presenter)

    def mock_setupUi(self, mock_setup):
        self.wgtWorkspacemanager = workspace_view
        self.wgtSlice = slice_view
        self.wgtPowder = powder_view
        self.wgtCut = cut_view
        global mainview
        mainview = mock_setup

    @patch.object(Ui_MainWindow, 'setupUi', mock_setupUi)
    @patch.object(PyQt4.QtGui.QMainWindow, '__init__', lambda x: None)
    def test_mainwindow(self):
        """Test the MainWindow initialises correctly"""
        MainWindow()
        self.projection_calculator.set_workspace_provider.assert_called_with(self.workspace_provider)
        self.slice_plotter.set_workspace_provider.assert_called_with(self.workspace_provider)
        self.cut_algorithm.set_workspace_provider.assert_called_with(self.workspace_provider)
