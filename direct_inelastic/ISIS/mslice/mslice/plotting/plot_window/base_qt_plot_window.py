from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt4 import QtGui

from .base_plot_window import BasePlotWindow

class MatplotlibCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class BaseQtPlotWindow(BasePlotWindow, QtGui.QMainWindow):
    """Inherit from this and a Ui_MainWindow from QT Designer to get a working PlotWindow

    The central widget will be replaced by the canvas"""
    def __init__(self, number, manager):
        super(BaseQtPlotWindow,self).__init__(number,manager)
        QtGui.QMainWindow.__init__(self)
        self.setupUi(self)
        self.canvas = MatplotlibCanvas(self)
        self.setCentralWidget(self.canvas)
        self.setWindowTitle('Figure %i'%number)


    def closeEvent(self, event):
        self._manager.figure_closed(self.number)
        event.accept()
