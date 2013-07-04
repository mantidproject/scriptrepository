"""
    Script used to start the MAGGOT GUI from MantidPlot
"""

import FunctionalMAGGOT as FunctionalMAGGOT
MainWindow = QtGui.QMainWindow()
ui = FunctionalMAGGOT.Ui_MAGGOTWindow()
ui.setupUi(MainWindow)
MainWindow.show()
