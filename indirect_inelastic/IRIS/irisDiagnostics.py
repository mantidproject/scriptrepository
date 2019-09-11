""" Plots heatmap of raw & normalised counts in each detector for a 
    specified run range to quickly identify dead detectors """


# import mantid algorithms
from mantid.simpleapi import *

from mantid.kernel import *
from mantid.api import *

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import numpy as np

class irisDiagnostics(PythonAlgorithm):

    def category(self):
        return 'Diagnostics\\Inelastic'

    def summary(self):
        return ("Plots a heatmap of the raw number of counts in each " 
                "spectroscopy detector within specified time of flight "  
                "boundaries for the contiguous run range given.") 


    def PyInit(self):
        # Declare properties
        self.declareProperty('Start run', 90000, 
        IntBoundedValidator(lower=41545),
        doc = 'First run to integrate')
        
        self.declareProperty('End run', 90050,
        IntBoundedValidator(lower=41545),
        doc = 'Last run to integrate')
        
        self.declareProperty('Start ToF integration', 56000,
        IntBoundedValidator(lower=18000),
        doc = "Time of flight for fastest neutron to count (microseconds)")

        self.declareProperty('End ToF integration', 76000,
        IntBoundedValidator(lower=18000), 
        doc = "Time of flight for slowest neutron to count (microseconds)")
        
    def PyExec(self):
        # Run the algorithm

        start_wsi_mon = 0 # workspace index: start counting at 0!
        end_wsi_mon = 1
        dets_mon = np.array(range(start_wsi_mon+1, end_wsi_mon+2))
        
        start_wsi_pg = 2 # workspace index: start counting at 0!
        end_wsi_pg = 52
        dets_pg = np.array(range(start_wsi_pg+1, end_wsi_pg+2))

        start_wsi_mi = 53 # workspace index: start counting at 0!
        end_wsi_mi = 103
        dets_mi = np.array(range(start_wsi_mi+1, end_wsi_mi+2))

        start_wsi_diff = 104 # workspace index: start counting at 0!
        end_wsi_diff = 111
        dets_diff = np.array(range(start_wsi_diff+1, end_wsi_diff+2))

        start_run = self.getProperty('Start run').value
        end_run = self.getProperty('End run').value
        raw_run = [i for i in range(start_run, end_run + 1)]
        
        raw_data_mon = []
        raw_data_pg = []
        raw_data_mi = []
        raw_data_diff = []
        binned_mon =[]
        binned_pg =[]
        binned_mi =[]
        binned_diff =[]
        bin_mon_summed = []
        bin_pg_summed = []
        bin_mi_summed = []
        bin_diff_summed = []

        for i in raw_run:
            # open file
            ws = Load('IRIS000%s.raw' % i)
            # integrate time bins
            ws_mon = Integration(ws, StartWorkspaceIndex=start_wsi_mon, EndWorkspaceIndex=end_wsi_mon)
            ws_pg = Integration(ws, StartWorkspaceIndex=start_wsi_pg, EndWorkspaceIndex=end_wsi_pg)
            ws_mi = Integration(ws, StartWorkspaceIndex=start_wsi_mi, EndWorkspaceIndex=end_wsi_mi)
            ws_diff = Integration(ws, StartWorkspaceIndex=start_wsi_diff, EndWorkspaceIndex=end_wsi_diff)
            # extract data
            binned_mon.append(ws_mon.extractY())
            binned_pg.append(ws_pg.extractY())
            binned_mi.append(ws_mi.extractY())
            binned_diff.append(ws_diff.extractY())
            # sum all detectors
            bin_mon_summed.append(np.sum(ws_mon.extractY()))
            bin_pg_summed.append(np.sum(ws_pg.extractY()))
            bin_mi_summed.append(np.sum(ws_mi.extractY()))
            bin_diff_summed.append(np.sum(ws_diff.extractY()))
        t_resolved_mon = np.squeeze(np.stack(binned_mon))
        t_resolved_pg = np.squeeze(np.stack(binned_pg))
        t_resolved_mi = np.squeeze(np.stack(binned_mi))
        t_resolved_diff = np.squeeze(np.stack(binned_diff))

        # plot
        fig = plt.figure()
#        ax0 = fig.add_subplot(241, title = 'Mon')
#        plt.imshow(t_resolved_mon, cmap = 'magma', extent=[min(dets_mon)-0.5, max(dets_mon)+0.5,
#        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto' )
#        plt.xlabel('Detector')
#        plt.ylabel('Run')

#        ax1 = fig.add_subplot(349)
#        plt.scatter(raw_run, bin_mon_summed)

#        ax2 = fig.add_subplot(245, title = 'Mon norm')
#        norm_mon = np.expand_dims(bin_pg_summed,1)
#        plt.imshow(t_resolved_mon / norm_mon, cmap = 'magma', extent=[min(dets_mon)-0.5, max(dets_mon)+0.5,
#        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto' )
        
#        plt.xlabel('Detector')
#        plt.ylabel('Run')

        ax3 = fig.add_subplot(231, title = 'PG')
        plt.imshow(t_resolved_pg, cmap = 'magma', extent=[min(dets_pg)-0.5, max(dets_pg)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto'  )
        plt.setp(ax3.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')

#        ax4 = fig.add_subplot(3,4,10)
#        plt.scatter(raw_run, bin_pg_summed)
                
        ax5 = fig.add_subplot(234, title = 'PG norm')
        norm_pg = np.expand_dims(bin_pg_summed,1)
        plt.imshow(t_resolved_pg / norm_pg, cmap = 'magma', extent=[min(dets_pg)-0.5, max(dets_pg)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto' )
        plt.setp(ax5.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')
#        plt.ylabel('Run')

        ax6 = fig.add_subplot(232, sharey=ax3, title = 'Mica')
        plt.imshow(t_resolved_mi, cmap = 'magma', extent=[min(dets_mi)-0.5, max(dets_mi)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto'  )
        plt.setp(ax6.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')

#        ax7 = fig.add_subplot(3,4,11)
#        plt.scatter(raw_run, bin_mi_summed)
                
        ax8 = fig.add_subplot(235, sharey=ax5, title = 'Mica norm')
        norm_mi = np.expand_dims(bin_mi_summed,1)
        plt.imshow(t_resolved_mi / norm_mi, cmap = 'magma', extent=[min(dets_mi)-0.5, max(dets_mi)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto' )
        plt.setp(ax8.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')
#        plt.ylabel('Run')
                
        ax9 = fig.add_subplot(233, sharey=ax3, title = 'Diff')
        plt.imshow(t_resolved_diff, cmap = 'magma', extent=[min(dets_diff)-0.5, max(dets_diff)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto'  )
        plt.setp(ax9.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')

#        ax10 = fig.add_subplot(3,4,12)
#        plt.scatter(raw_run, bin_diff_summed)
                
        ax11 = fig.add_subplot(236, sharey=ax5, title = 'Diff norm')
        norm_diff = np.expand_dims(bin_diff_summed,1)
        plt.imshow(t_resolved_diff / norm_diff, cmap = 'magma', extent=[min(dets_diff)-0.5, max(dets_diff)+0.5,
        max(raw_run)+0.5, min(raw_run)-0.5], aspect = 'auto' )
        plt.setp(ax11.get_yticklabels(), visible=False)
#        plt.xlabel('Detector')
#        plt.ylabel('Run')
        
        plt.tight_layout()
        plt.show()
        
        ws_junk = [ws, ws_mon, ws_pg, ws_mi, ws_diff]
        DeleteWorkspaces(ws_junk)
        
# Register algorithm with Mantid
AlgorithmFactory.subscribe(irisDiagnostics)