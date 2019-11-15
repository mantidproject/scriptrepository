from __future__ import print_function
from mantid import *
from mantid.simpleapi import *
from scipy.optimize import curve_fit
import numpy as np

class Reduce():

    def __init__ (self, sample, energies, PF,
                  he_pressure = 0.9,
                  he_path_length = 0.06,
                  polmon_distance = 25.38,
                  polmon_delay = 100,
                  name_format = 'LET{0}_{1:<3.2f}meV_OneToOne.nxspe',
                  label = '',
                  mask = 'PLET_184_msk.xml',
                  NSF_first = True,
                  separate = False
                  ):

        self.energies = energies
        self.sample_runs = sample
        self.PF = PF
        self.he_pressure = he_pressure          # in bar
        self.he_path_length = he_path_length    # in metres
        self.polmon_distance = polmon_distance  # in metres
        self.polmon_delay = polmon_delay        # in microsecs
        self.NSF_first = NSF_first
        self.separate = separate
        self.label = label
        self.mask = mask
        self.name_format = name_format

#---------------------------------------------------------------------------

    def generate_dummy(self, run, ei):
        """Clears out a workspace with the correct dimensions to be used
           as an empty workspace that can be cloned and populated"""
        print("      ... Generating dummy workspace...")
        Dummy = LoadNXSPE(self.name_format.format(run, ei))
        DummyWorkspace = mtd['Dummy']*0.0
        self.DummyWorkspace = mtd['DummyWorkspace']
        DeleteWorkspace(Dummy)

#---------------------------------------------------------------------------

    def get_PF_from_quartz(self): 
        """Calculates PF (product of polariser and flipper efficiencies) using the out of plane
           angular dependence of the flipping ratio of quartz through the 3He analyser."""

        for ei in self.energies:
            wl = 9.045 / np.sqrt(ei)
            self.generate_dummy(self.sample_runs[0],ei)

            NSF_quartz_total = CloneWorkspace(self.DummyWorkspace)
            SF_quartz_total  = CloneWorkspace(self.DummyWorkspace)

# array of out-of-plane angles
            gamma = np.deg2rad(np.linspace(-30,30,num=256))
            
            print("****************************************")        

# read in NSF and SF quartz and totalize
            for run in self.sample_runs[::2]:
                if self.NSF_first:
                    NSF_run = run
                    SF_run = run + 1
                else:
                    NSF_run = run + 1
                    SF_run = run

                print(("Loading quartz runs {0} (NSF) and {1} (SF) at {2:<3.2f}meV".format(NSF_run,SF_run, ei)))
                NSF_quartz = LoadNXSPE(self.name_format.format(NSF_run, ei))
                SF_quartz  = LoadNXSPE(self.name_format.format(SF_run, ei))
                
                NSF_quartz_total += NSF_quartz
                SF_quartz_total  += SF_quartz 

            if not "Masking" in mtd:
                LoadMask('let',InputFile=self.mask,RefWorkspace=NSF_quartz_total)

# integrate over the elastic line and mask
            NSF_quartz_total = Integration(NSF_quartz_total,RangeLower=-ei*0.02,RangeUpper=ei*0.02)
            SF_quartz_total  = Integration(SF_quartz_total, RangeLower=-ei*0.02,RangeUpper=ei*0.02)
            MaskDetectors(NSF_quartz_total, MaskedWorkspace='Masking')
            MaskDetectors(SF_quartz_total,  MaskedWorkspace='Masking')

# LET_gamma_grouping sums over 2-theta to leave workspaces as function of gamma
            NSF_quartz_total = GroupDetectors('NSF_quartz_total', MapFile='LET_gamma_grouping.xml')
            SF_quartz_total  = GroupDetectors('SF_quartz_total',  MapFile='LET_gamma_grouping.xml')

            FAP = (NSF_quartz_total - SF_quartz_total) / (NSF_quartz_total + SF_quartz_total)
            FAP = Transpose(FAP)
            FAP_tofit  = FAP.extractY()[0]
            FAP_tofite = FAP.extractE()[0]

# neutron polarization vs. gamma            
            def FAP_gamma(gamma, PF, PHe):
                return PF*np.tanh(7.33*wl*self.he_path_length*self.he_pressure*PHe*(1.0/np.cos(gamma)))

# Only fit to values where our secant approximations holds well and with physically meaningful numbers
# of counts
            indices = np.intersect1d(list(np.where(np.abs(gamma) > 0.1)[0]), list(np.where(np.abs(gamma < 0.3))[0]))
            indices = np.intersect1d(list(np.where(FAP_tofit > 0.)[0]), indices)

            popt, pcov = curve_fit(FAP_gamma, gamma[indices], FAP_tofit[indices],p0=[0.9, 0.5], sigma=FAP_tofite[indices])

            PF = popt[0]
            PFe = np.sqrt(np.diag(pcov))[0]
            PHe = popt[1]
            PHee = np.sqrt(np.diag(pcov))[1]

# output fit to workspace as check
            fit = FAP_gamma(gamma,PF,PHe)
            fit_check = CreateWorkspace(np.linspace(1,256,num=256),np.array(fit))

            print("****************************************") 
            print(("PF   = {0:1.3f} +/- {1:1.3f}".format(PF, PFe)))
            print(("P_He = {0:1.3f} +/- {1:1.3f}".format(PHe, PHee)))
            print("****************************************") 

#-------------------------------------------------------------------------------------
    def get_helium_from_ROI(self, min_t, max_t, roi="MaskWorkspace"):
        """Using PF and the flipping ratio from a region of interest (ROI) in the data
           (i.e. a Bragg peak from the Al sample can) determine PHe0 and T1.
           The minimum and maximum TOF need to be provided as well as the Mask workspace
           containing the ROI."""

        if not roi in mtd:
            print("ERROR: ROI workspace not found")
            return

        ei = self.energies[0]
        PF = self.PF[0]
        wl = 9.045 / np.sqrt(ei)

        first = True
        t0    = 0
        PHes  = []
        PHese = []
        times = []
        
        print("\nStarting get_helium_from_ROI...")
        print(("Using {0} meV rep with PF={1}".format(ei,PF)))
        print("****************************************************")        
        for run in self.sample_runs[::2]:
            if self.NSF_first:
                NSF_run = run
                SF_run = run + 1
            else:
                NSF_run = run + 1
                SF_run = run

            print(("Calculating P_He from ROI for runs NSF:{0} and SF:{1} at {2:<3.2f}meV ".format(NSF_run,SF_run,ei)))

# load data and extract dummy monitor spectra from ROI
            NSF_roi = Load("LET000{0}.nxs".format(NSF_run))
            SF_roi  = Load("LET000{0}.nxs".format(SF_run))
            NSF_roi = NormaliseByCurrent(NSF_roi)
            SF_roi  = NormaliseByCurrent(SF_roi)
            NSF_roi = Rebin(NSF_roi,10)
            SF_roi  = Rebin(SF_roi,10)
            MaskDetectors(NSF_roi,MaskedWorkspace=roi)
            MaskDetectors(SF_roi,MaskedWorkspace=roi)
            NSF_monitors = SumSpectra(NSF_roi)
            SF_monitors  = SumSpectra(SF_roi)
            
            start_time = NSF_monitors.getSampleDetails().startTime().to_datetime64()
            end_time   =  SF_monitors.getSampleDetails().endTime().to_datetime64()
            
            time = start_time + (end_time - start_time)/2.0

            if first:
                t0 = time
                first = False
                
# integrate between user defined time limits passed                
            NSF_int = Integration(NSF_monitors,RangeLower=min_t,RangeUpper=max_t)                                      
            SF_int  = Integration(SF_monitors, RangeLower=min_t,RangeUpper=max_t)
            
# calculate the helium polarization, and append times and PHe to arrays                                                                    
            FAP = (NSF_int - SF_int) / (NSF_int + SF_int)
            A   = FAP.readY(0)[0] / PF
            Ae_fractional = FAP.readE(0)[0] / FAP.readY(0)[0]
            PHe = np.abs(np.arctanh(A) / (7.33*wl*self.he_path_length*self.he_pressure))
            
            if not np.isfinite(PHe):
                print("ERROR: Bad P_He value...")
                return
            
            PHes.append(PHe)
            PHese.append(Ae_fractional*PHe)
            times.append(float((time - t0))*1e-9/60./60.)
            
# fit the helium polarization vs time
        def exp_decay(t,T1,P0):
            return P0 * np.exp(t/T1)
        
        popt, pcov = curve_fit(exp_decay, times, PHes, p0=[-20, 0.5], sigma=PHese)
        self.T1    = popt[0]
        self.T1e   = np.sqrt(np.diag(pcov))[0]
        self.PHe0  = popt[1]
        self.PHe0e = np.sqrt(np.diag(pcov))[1]
        self.t0    = t0
        self.times = times
        
# output debugging workspaces
        he_T1_check = CreateWorkspace(np.array(times),np.array(PHes),np.array(PHese))
        he_fit      = exp_decay(times, self.T1, self.PHe0)
        he_fit_ws   = CreateWorkspace(np.array(times),np.array(he_fit))

        print(("\nInitial Cell polarization P0={0:1.3f}+/-{1:1.3f} with lifetime T1={2:3.2f}+/-{3:1.2f}".format(self.PHe0,self.PHe0e,-self.T1,self.T1e)))
        print("************************************************************************") 

#-------------------------------------------------------------------------------------

    def get_helium_parameters(self):
        """Using PF and the flipping ratio in the straight-thru monitor, determine PHe0 and T1"""

        ei = self.energies[0]
        PF = self.PF[0]
        wl = 9.045 / np.sqrt(ei)
        TOF = 251.9 * wl * self.polmon_distance + self.polmon_delay

# Now we have PF we can use the cell transmission monitor in a sample run to fit the lifetime
# and polarisation of the 3He cell, calculate the flipping ratio for each run, and apply the pol
# corrections
        first = True
        t0    = 0
        PHes  = []
        PHese = []
        times = []

        print("\nStarting get_helium_parameters...")
        print(("Using {0} meV rep with PF={1}".format(ei,PF)))
        print("****************************************************")        
        for run in self.sample_runs[::2]:
            if self.NSF_first:
                NSF_run = run
                SF_run = run + 1
            else:
                NSF_run = run + 1
                SF_run = run

            print(("Calculating P_He for runs NSF:{0} and SF:{1} at {2:<3.2f}meV with TOF {3:>6.0f} us".format(NSF_run,SF_run,ei,TOF)))

#load monitors and calculate flipping ratios
            NSF_monitors = LoadNexusMonitors("LET000{0}.nxs".format(NSF_run))
            SF_monitors  = LoadNexusMonitors("LET000{0}.nxs".format(SF_run))
            NSF_monitors = NormaliseByCurrent(NSF_monitors)
            SF_monitors  = NormaliseByCurrent(SF_monitors)
            
            start_time = NSF_monitors.getSampleDetails().startTime().to_datetime64()
            end_time   =  SF_monitors.getSampleDetails().endTime().to_datetime64()
            
            time = start_time + (end_time - start_time)/2.0
            
            if first:
                t0 = time
                first = False

# Integrate the monitors over the appropriate TOF                
            NSF_int = Integration(NSF_monitors,RangeLower=TOF-100,RangeUpper=TOF+100,StartWorkspaceIndex=7)                                      
            SF_int  = Integration(SF_monitors, RangeLower=TOF-100,RangeUpper=TOF+100,StartWorkspaceIndex=7)
                                                                    
            FAP = (NSF_int - SF_int) / (NSF_int + SF_int)
            A   = FAP.readY(0)[0] / PF
            Ae_fractional = FAP.readE(0)[0] / FAP.readY(0)[0]
            PHe = np.abs(np.arctanh(A) / (7.33*wl*self.he_path_length*self.he_pressure))
            
            if not np.isfinite(PHe):
                print("ERROR: Bad P_He value...")
                return
            
            PHes.append(PHe)
            PHese.append(Ae_fractional*PHe)
            times.append(float((time - t0))*1e-9/60./60.)
            
# fit the helium polarization vs time
        def exp_decay(t,T1,P0):
            return P0 * np.exp(t/T1)
            
        popt, pcov = curve_fit(exp_decay, times, PHes, p0=[-20, 0.5], sigma=PHese)
        self.T1    = popt[0]
        self.T1e   = np.sqrt(np.diag(pcov))[0]
        self.PHe0  = popt[1]
        self.PHe0e = np.sqrt(np.diag(pcov))[1]
        self.t0    = t0
        self.times = times
        
# output debugging workspaces
        he_T1_check = CreateWorkspace(np.array(times),np.array(PHes),np.array(PHese))
        he_fit      = exp_decay(times, self.T1, self.PHe0)
        he_fit_ws   = CreateWorkspace(np.array(times),np.array(he_fit))

        print(("\nInitial Cell polarization P0={0:1.3f}+/-{1:1.3f} with lifetime T1={2:3.2f}+/-{3:1.2f}".format(self.PHe0,self.PHe0e,-self.T1,self.T1e)))
        print("************************************************************************") 
        
#--------------------------------------------------------------------------

    def correct_data(self):
        """Using PF, PHe0, and T1, calculate the flipping ratio for all (gamma, t) and correct runs,
           summing and averaging at the end."""

        PF_iter = iter(self.PF)

        for ei in self.energies:

            PF = next(PF_iter)
            print("****************************************")
            print(("\ncorrect_data: Correcting {0:<3.2f}meV rep with PF={1} \n".format(ei,PF)))
            NSF_out = "PLET_{0}_{1:<3.2f}meV_NSF".format(self.label,ei)
            SF_out  = "PLET_{0}_{1:<3.2f}meV_SF".format(self.label,ei)   

            self.generate_dummy(self.sample_runs[0],ei)
            NSF_total       = CloneWorkspace(self.DummyWorkspace)
            SF_total        = CloneWorkspace(self.DummyWorkspace)
            Scharpf_ws      = CloneWorkspace(self.DummyWorkspace)
            transmission_ws = CloneWorkspace(self.DummyWorkspace)

# Make wavelength and gamma arrays, same shape as the DummyWorkspace
            y_shape = NSF_total.extractY().shape   
            delta_e_binning = NSF_total.extractX()
            final_energy_binning = np.abs(delta_e_binning[:y_shape[0],:y_shape[1]] - ei)
            wl = 9.045 / np.sqrt(final_energy_binning)
            dummy_gammaspace      = np.zeros(y_shape).T
            dummy_gammaspace[:,:] = np.tile(np.deg2rad(np.linspace(-30,30,num=256)), 384)
            gammaspace            = dummy_gammaspace.T

# Calculate cell opacity
            opacity = 7.33 * wl * self.he_path_length * self.he_pressure
            itime = 0

# Cycle through runs - uses same "times" array as get_helium_parrameters
            for run in self.sample_runs[::2]:
                if self.NSF_first:
                    NSF_run = run
                    SF_run = run + 1
                else:
                    NSF_run = run + 1
                    SF_run = run

                print(("Correcting runs NSF:{0} and SF:{1} at {2:<3.2f}meV".format(NSF_run,SF_run,ei)))
                NSF_One2One = LoadNXSPE(self.name_format.format(NSF_run, ei))
                SF_One2One  = LoadNXSPE(self.name_format.format(SF_run, ei))   

# Calculate FAP, cell transmission and Scharpf correction factor for these runs            
                time = self.times[itime]
                itime += 1
                PHe = self.PHe0 * np.exp(time / self.T1)

                FAP = PF * np.tanh(opacity*(1.0/np.cos(gammaspace)*PHe))
                print(("After {0:1.2f} hours, cell polarization is {1:1.3f}".format(time,PHe)))
                flipping_ratio = (1.0 + FAP) / (1.0 - FAP)
                transmission = np.exp(-opacity) * np.cosh(opacity * PHe)
                Scharpf = (1.0 / (flipping_ratio - 1.0))

# Populate Scharpf and transmission workspaces
                for i in range(y_shape[0]):
                    Scharpf_ws.setY(i, Scharpf[i,:y_shape[1]])
                    Scharpf_ws.setX(i, delta_e_binning[i,:])
                for i in range(y_shape[0]):
                    transmission_ws.setY(i, transmission[i,:y_shape[1]])
                    transmission_ws.setX(i, delta_e_binning[i,:])
                 
# Apply corrections, totalise and rename workspace to save it
                Diff     = (NSF_One2One - SF_One2One) * Scharpf_ws
                NSF_corr = (NSF_One2One + Diff) / transmission_ws
                SF_corr  = (SF_One2One  - Diff) / transmission_ws
                NSF_total.setYUnit('')
                SF_total.setYUnit('')
                NSF_total = NSF_total + NSF_corr
                SF_total  = SF_total  + SF_corr
    
            NSF_total /= len(self.sample_runs[::2])             
            SF_total  /= len(self.sample_runs[::2])
            RenameWorkspace(NSF_total,OutputWorkspace=NSF_out)
            RenameWorkspace(SF_total,OutputWorkspace=SF_out)
            
            print("****************************************")

#----------------------------------------------------------------------------

    def components(self):

        self.separate = True
        print("Combining components...")
        for ei in self.energies:

            NSF_in  = "PLET_{0}_{1:<3.2f}meV_NSF".format(self.label,ei)
            SF_in   = "PLET_{0}_{1:<3.2f}meV_SF".format(self.label,ei)            
            coh_out = "PLET_{0}_{1:<3.2f}meV_coh".format(self.label,ei)
            inc_out = "PLET_{0}_{1:<3.2f}meV_inc".format(self.label,ei)

            NSF = mtd[NSF_in]
            SF  = mtd[SF_in]

            coh = NSF - 0.5 * SF 
            inc = 1.5 * SF

            RenameWorkspace(coh,OutputWorkspace=coh_out)
            RenameWorkspace(inc,OutputWorkspace=inc_out)

#---------------------------------------------------------------------------

    def rings_output(self, rings_map='LET_rings_153.map'):

        print("Outputting rings nxspe files...")

        if self.separate:
            ext1 = "coh"
            ext2 = "inc"
        else:
            ext1 = "NSF"
            ext2 = "SF"

        for ei in self.energies:
            NSF_in  = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext1)
            SF_in   = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext2)
            NSF_out  = "PLET_{0}_{1:<3.2f}meV_{2}_rings.nxspe".format(self.label,ei,ext1)
            SF_out   = "PLET_{0}_{1:<3.2f}meV_{2}_rings.nxspe".format(self.label,ei,ext2)
            NSF = mtd[NSF_in]
            SF  = mtd[SF_in]
            if not "Masking" in mtd:
                LoadMask('let',InputFile=self.mask,RefWorkspace=NSF)
            MaskDetectors(NSF,MaskedWorkspace="Masking")
            grouped_NSF = GroupDetectors(NSF,MapFile=rings_map,PreserveEvents=False,Behaviour='Average')
            SaveNXSPE(grouped_NSF,Filename=NSF_out,Efixed=ei, KiOverKfScaling=False)
            MaskDetectors(SF,MaskedWorkspace="Masking")
            grouped_SF  = GroupDetectors(SF, MapFile=rings_map,PreserveEvents=False,Behaviour='Average')
            SaveNXSPE(grouped_SF, Filename=SF_out, Efixed=ei, KiOverKfScaling=False)

#---------------------------------------------------------------------------

    def one2one_output(self):

        print("Outputting one2one nxspe files...")

        if self.separate:
            ext1 = "coh"
            ext2 = "inc"
        else:
            ext1 = "NSF"
            ext2 = "SF"

        for ei in self.energies:
            NSF_in  = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext1)
            SF_in   = "PLET_{0}_{1:<3.2f}meV_{2}".format(self.label,ei,ext2)
            NSF_out  = "PLET_{0}_{1:<3.2f}meV_{2}_One2One.nxspe".format(self.label,ei,ext1)
            SF_out   = "PLET_{0}_{1:<3.2f}meV_{2}_One2One.nxspe".format(self.label,ei,ext2)
            NSF = mtd[NSF_in]
            SF  = mtd[SF_in]
            if not "Masking" in mtd:
                LoadMask('let',InputFile=self.mask,RefWorkspace=NSF)
            MaskDetectors(NSF,MaskedWorkspace="Masking")
            SaveNXSPE(NSF,Filename=NSF_out,Efixed=ei, KiOverKfScaling=False)
            MaskDetectors(SF,MaskedWorkspace="Masking")
            SaveNXSPE(SF, Filename=SF_out, Efixed=ei, KiOverKfScaling=False)






