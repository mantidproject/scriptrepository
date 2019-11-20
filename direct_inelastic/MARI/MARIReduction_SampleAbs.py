#pylint: disable=invalid-name
""" Sample MARI reduction script """
from __future__ import print_function
import os,sys
from numpy import *
from mantid import *
from Direct.AbsorptionShapes import *
from Direct.ReductionWrapper import *
from mantid.simpleapi import *
from mantid.kernel import funcinspect
from mantid.dataobjects import EventWorkspace
import six
import types
from Direct.PropertyManager import PropertyManager
import Direct
import numpy as np

class MARIReduction(ReductionWrapper):
    @MainProperties
    def def_main_properties(self):
        """Define main properties used in reduction. These are the property
           a user usually wants to change

        MARI Instrument scientist beware!!!!
        -- the properties set up here may be overridden in iliad_mari (below ) if you use it, or
            in section __name__=='__main__' below if you do not use iliad_mari
        """
        prop = {}
        # if energy is specified as a list (even with single value e.g. ei=[81])
        # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is
        # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
        #
        #prop['incident_energy'] = 50
        #prop['energy_bins'] = [-20,0.1,49]
        prop['incident_energy'] = [20, 5.5]
        prop['energy_bins'] = [-1, 0.005, 0.97]
        #
        # the range of files to reduce. This range ignored when deployed from autoreduction,
        # unless you going to sum these files.
        # The range of numbers or run number is used when you run reduction from PC.

        # If you "save" a run without ending it, you have to give the file name
        #prop['sample_run'] = ['MAR25360.n001']
        # Otherwise just give the run number
        #prop['sample_run'] = [25362, 25363, 25364, 25365]
        prop['sample_run'] = 25781
        prop['wb_run'] = 25779

        prop['sum_runs'] = False # set to true to sum everything provided to sample_run
        #                        # list
        # Absolute units reduction properties. Set prop['monovan_run']=None to do relative units
        prop['monovan_run'] = None
        prop['sample_mass'] = 0
        prop['sample_rmm'] = 0
        return prop

    @AdvancedProperties
    def def_advanced_properties(self):
        """Set up advanced properties, describing reduction.
           These are the properties, usually provided by an instrument
           scientist

           separation between simple and advanced properties depends
           on scientist, experiment and user.   All are necessary for reduction
           to work properly

        MARI Instrument scientist beware!!!!
        -- the properties set up here may be overridden in iliad_mari (below ) if you use it, or
            in section __name__=='__main__' below if you do not use iliad_mari
        """
        prop = {}
        prop['normalise_method'] = 'current'
        prop['map_file'] = "mari_res2013.map"
        prop['monovan_mapfile'] = "mari_res2013.map"

        # Next lines are for removing detector artifacts which should not be needed
        #prop['remove_streaks'] = True
        #prop['fakewb'] = True
        #
        #prop['hardmaskOnly']=maskfile # disable diag, use only hard mask
        #prop['hard_mask_file'] = "mari_mask2019.msk"
        prop['det_cal_file'] = ''
        # Comment out the next line if you want to use the data run for background masking
        #prop['mask_run'] = 25035
        #prop['use_hard_mask_only'] = True
        prop['save_format'] = 'nxspe'
        #
        #prop['wb_integr_range'] = [2,10]
        prop['data_file_ext'] = '.nxs' # if two input files with the same name and
                                       # different extension found, what to prefer.
        prop['load_monitors_with_workspace'] = False
        # change this to correct value and verify that motor_log_names refers correct and existing
        # log name for crystal rotation to write correct psi value into nxspe files
        prop['motor_offset']=None
        prop['check_background']=False
        prop['bkgd-range-min']=18000
        prop['bkgd-range-max']=19000
        # Uncomment two following properties to correct for absorption
        # on sample or sample container during the experiment.
        # 1) Define the sample material and sample shape:
        #
        #prop['correct_absorption_on'] = Cylinder(['Fe'],{'Height':10,'Radius':2})
        #
        #  The shapes currently available are:
        #  Cylinder, FlatPlate, HollowCylinder and Sphere
        # Each class accepts two parameters, namely dictionary or list describing
        # sample material, as accepted by SetSampleMaterial algorithm
        # and Shape parameters as accepted by SetSample algorithm.

        # Go to iPython window, type from AbsorptionShapes import * and type
        # help(Cylinder),
        # help(Sphere) or help(anAbsorptionShape) to learn more about different
        # ways of setting the material and shape parameters:
        ##----------------------------
        # 2) Optionally provide additional parameters defining speed and accuracy of the
        #    absorption corrections algorithm in abs_corr_info dictionary.
        #
        #prop['abs_corr_info'] = {'is_mc':True}
        #
        #Two generic algorithms are currently available to
        #    correct for absorption:
        # a) MonteCarloAbsorption (invoked by 'is_mc':True key of abs_corr_info
        #    property and by default)
        # and
        # b) AbsorptionCorrections (invoked by 'is_fast':True key of
        #    abs_corr_info property). This algorithm has number of optimizations
        #    for the special shapes cases and need further scientific validations

        # All other properties, to be provided in the input dictionary of this
        # property are the non-sample related properties of the
        # correspondent correction algorithm.

        return prop

    @iliad
    def reduce(self,input_file=None,output_directory=None):
        """ Method executes reduction over single file
          Modify only if custom pre or post-processing is needed, namely:
        """

        """
            Define custom preprocessing procedure, to be applied to the whole
            run, obtained from the instrument before diagnostics is applied
        1) Get access to the pointer to the input workspace
            (the workspace will be loaded  and calibrated at this point if it has not been done yet)
            using sample_run property class:
            run_ws = PropertyManager.sample_run.get_worksapce()

        2) Get access to any property defined in Instrument_Properties.xml file
           or redefined in the reduction script:
           properties = self.reducer.prop_man

        e.g:
        RunNumber = properties.sample_run
        ei = properties.incident_energy

            Perform custom preprocessing procedure (with the values you specified)
        preprocessed_ws = custom_preprocessing_procedure(run_ws,RunNumber,ei,...)
        3) Store preprocessed workspace in the sample_run property for further analysis from preprocessing
            for the case where the workspace name have changed in the Analysis Data Service
            (this situation is difficult to predict so better always do this operation)
        PropertyManager.sample_run.synchronize_ws(preprocessed_ws)
        """
        output = ReductionWrapper.reduce(self,input_file,output_directory)
        # Auto-reduction returns workspace list, so for compatibility with auto-reduction
        # we better process any output as reduction list
        """ Defined custom post-processing procedure, in the way, similar to
            the preprocessing procedure, using the WS pointer, returned by the reduction procedure
            above. If the run is reduced in multirep mode, the WS is the list of the reduced
            workspace pointers, so the procedure should be applied to each workspace.

            At this stage standard saving had already been done, so save your results
            if you need them in a future:
            for i in range(0:len(WS)):
                SaveNexus(WS[i],Filename = 'CustomFilenameN%d.nxs'.format(i))
        """

        if not isinstance(output,list):
            output = [output]
        for ws in output:
            ei = ws.getEFixed(1)
            q_min = 0.04*sqrt(ei)
            q_max = 1.3*sqrt(ei)
            q_bins = str(q_min)+','+str(q_max/285.)+','+str(q_max)
            wsn = ws.name()
            SofQW3(InputWorkspace=ws, OutputWorkspace=wsn+'_SQW', QAxisBinning=q_bins, Emode='Direct')
            Transpose(InputWorkspace=wsn+'_SQW', OutputWorkspace=wsn+'_SQW')
        return output

    def run_reduction(self, out_ws_name=None):
        """" Reduces runs one by one or sum all them together and reduce after this

            if wait_for_file time is > 0, it will until  missing files appear on the
            data search path
        """
        try:
            _,r = funcinspect.lhs_info('both')
            out_ws_name = r[0]
        except:
            pass

        if not hasattr(PropertyManager.wb_run, '_old_get_workspace'):
            PropertyManager.wb_run._old_get_workspace = PropertyManager.wb_run.get_workspace

        old_wb_get_workspace = PropertyManager.wb_run._old_get_workspace
        if self.reducer.prop_man.fakewb is True:
            def new_wb_get_workspace(self):
                ws = old_wb_get_workspace()
                if ((self._run_number is not None and self._run_number != 25035)
                    or ('25035' not in self._ws_name)) or ws.run().hasProperty('faked'):
                    return ws
                print("*** Faking Whitebeam run")
                x = ws.extractX()
                y = ws.extractY()
                e = ws.extractE()
                for ifake0, ifake1, ireal0, ireal1 in [[404, 440, 143, 179], [663, 699, 143, 179],
                                [441, 470, 182, 211], [700, 729, 182, 211], [276, 371, 535, 630],
                                [373, 378, 632, 637], [381, 386, 640, 645], [389, 394, 648, 653]]:
                    for ifake, ireal in np.array([list(range(ifake0-1, ifake1)), list(range(ireal0-1, ireal1))]).T:
                        y[ifake, :] = y[ireal, :]
                        e[ifake, :] = e[ireal, :]
                # Masking is messed up if we use this fake white beam so put masks here directly...
                for ifake in [351, 617, 846]:
                    y[ifake,:] = y[ifake,:] * 0
                    e[ifake,:] = e[ifake,:] * 0
                for isp in range(ws.getNumberHistograms()):
                    ws.setY(isp, y[isp, :])
                    ws.setE(isp, e[isp, :])
                AddSampleLog(ws, 'faked', 'already_faked')
                return ws
        else:
            def new_wb_get_workspace(self):
                ws = old_wb_get_workspace()
                return ws

        PropertyManager.wb_run.get_workspace = types.MethodType(new_wb_get_workspace, PropertyManager.wb_run)

        if not hasattr(PropertyManager.sample_run, '_old_get_workspace'):
            PropertyManager.sample_run._old_get_workspace = PropertyManager.sample_run.get_workspace

        old_get_workspace = PropertyManager.sample_run._old_get_workspace
        if self.reducer.prop_man.remove_streaks is True:
            def new_get_workspace(self):
                ws = old_get_workspace()
                if isinstance(ws, EventWorkspace) and not ws.run().hasProperty('unstreaked'):
                    print('*** Removing Streaks')
                    wsn = ws.name()

                    t0 = 221.75                          # ToF of first streak in us
                    stp = 85.33333333333333              # ToF between streaks in us (==256/3)
                    w1 = 150                             # Number of ToF bins too look for streak around expected position
                    w2 = 10                              # Number of ToF bins around streaks to calculate background level
                    spikes_tof = np.arange(t0, 20000, stp)
                    spikes_tof = np.round(spikes_tof * 4) / 4

                    SumSpectra(wsn, IncludeMonitors=False, OutputWorkspace=wsn+'s')
                    wsr = Rebin(wsn+'s', '1,0.25,20000',PreserveEvents=False, OutputWorkspace=wsn+'s')
                    xx = (np.array(wsr.extractX()).T)[:,0]
                    yy = (np.array(wsr.extractY()).T)[:,0]
                    ee = (np.array(wsr.extractE()).T)[:,0]

                    bad = []
                    ymax = np.max(yy)
                    for spk in spikes_tof:
                        ix = np.where(xx == spk)[0][0]
                        yy2 = yy[(ix-w1):(ix+w1)]
                        iy = np.where(yy2 == np.max(yy2))[0][0] + ix - w1
                        yv = yy[iy]
                        yy3 = yy[(iy-w2):(iy+w2)]
                        mv = np.mean(yy3[np.where(yy3 != yv)])
                        if yv > mv * 1.5 and yv > ymax/500:
                            bad.append(iy)

                    badtof = xx[bad]
                    for id in range(ws.getNumberHistograms()):
                         ev = ws.getEventList(id)
                         for tof in badtof:
                            ev.maskTof(tof-0.075, tof+0.225)
                    AddSampleLog(ws, 'unstreaked', 'unstreaked')
                    DeleteWorkspace(wsr)

                return ws

        else:
            def new_get_workspace(self):
                ws = old_get_workspace()
                return ws

        PropertyManager.sample_run.get_workspace = types.MethodType(new_get_workspace, PropertyManager.sample_run)

        # if this is not None, we want to run validation not reduction
        if self.validate_run_number:
            self.reducer.prop_man.log\
                ("**************************************************************************************",'warning')
            self.reducer.prop_man.log\
                ("**************************************************************************************",'warning')
            rez,mess=self.build_or_validate_result()
            if rez:
                self.reducer.prop_man.log("*** SUCCESS! {0}".format(mess))
                self.reducer.prop_man.log\
                    ("**************************************************************************************",'warning')

            else:
                self.reducer.prop_man.log("*** VALIDATION FAILED! {0}".format(mess))
                self.reducer.prop_man.log\
                    ("**************************************************************************************",'warning')
                raise RuntimeError("Validation against old data file failed")
            self.validate_run_number=None
            return rez,mess

        sam_run = self.reducer.prop_man.sample_run

        setattr(Direct.diagnostics, 'normalise_background', mari_normalise_background)
        if self.reducer.sum_runs:
        ### sum runs provided
            if out_ws_name is None:
                return self.sum_and_reduce()
            else:
                red_ws = self.sum_and_reduce()
                if len(red_ws) > 1:
                    ws_list = []
                    for id, ws_out in enumerate(red_ws):
                        ws_list.append('{0}_{1}_sum_SQW'.format(out_ws_name, id))
                        RenameWorkspace(InputWorkspace=ws_out.name()+'_SQW', OutputWorkspace=ws_list[-1])
                        ws_list.append('{0}_{1}_sum'.format(out_ws_name, id))
                        RenameWorkspace(InputWorkspace=ws_out, OutputWorkspace=ws_list[-1])
                    GroupWorkspaces(InputWorkspaces=ws_list, OutputWorkspace=out_ws_name)
                else:
                    RenameWorkspace(InputWorkspace=red_ws[0].name()+'_SQW', OutputWorkspace=out_ws_name+'_sum_SQW')
                    RenameWorkspace(InputWorkspace=red_ws[0], OutputWorkspace=out_ws_name+'_sum')
                return red_ws
        else:
        ### reduce list of runs one by one
            runfiles = PropertyManager.sample_run.get_run_file_list()
            #if hasattr(runfiles, '__len__') and len(runfiles) > 1:
            #    runfiles = [runfiles[-1]]
            if out_ws_name is None:
                ws_refs = []
                for file_name in runfiles:
                    ws_refs.append(self.reduce(file_name))
                return ws_refs if len(runfiles) > 1 else ws_refs[0]
            else:
                results = []
                nruns = len(runfiles)
                for num, file_name in enumerate(runfiles):
                    red_ws = self.reduce(file_name)
                    if isinstance(red_ws, list):
                        for ws in red_ws:
                            results.append(ws)
                        if len(red_ws) > 1:
                            ws_list = []
                            for id, ws_out in enumerate(red_ws):
                                print('--------------------')
                                print(ws_out.name())
                                print('--------------------')
                                ws_list.append('{0}_{1}_SQW'.format(out_ws_name, id))
                                RenameWorkspace(InputWorkspace=ws_out.name()+'_SQW', OutputWorkspace=ws_list[-1])
                                ws_list.append('{0}_{1}'.format(out_ws_name, id))
                                RenameWorkspace(InputWorkspace=ws_out, OutputWorkspace=ws_list[-1])
                            GroupWorkspaces(InputWorkspaces=ws_list, OutputWorkspace=out_ws_name)
                        else:
                            RenameWorkspace(InputWorkspace=red_ws[0].name()+'_SQW', OutputWorkspace=out_ws_name+'_SQW')
                            RenameWorkspace(InputWorkspace=red_ws[0], OutputWorkspace=out_ws_name)
                    else:
                        if nruns == 1:
                            if red_ws.name() != out_ws_name:
                                RenameWorkspace(InputWorkspace=red_ws, OutputWorkspace=out_ws_name)
                                RenameWorkspace(InputWorkspace=red_ws.name()+'_SQW', OutputWorkspace=out_ws_name+'_SQW')
                            results.append(mtd[out_ws_name])
                        else:
                            OutWSName = '{0}#{1}of{2}'.format(out_ws_name,num+1,nruns)
                            if red_ws.name() != out_ws_name:
                                RenameWorkspace(InputWorkspace=red_ws, OutputWorkspace=OutWSName)
                                RenameWorkspace(InputWorkspace=red_ws.name()+'_SQW', OutputWorkspace=OutWSName+'_SQW')
                            results.append(mtd[OutWSName])
                if len(results) == 1:
                    return results[0]
                else:
                    return results
    def do_preprocessing(self,reducer,ws):
        """ Custom function, applied to each run or every workspace, the run is divided to
            in multirep mode
            Applied after diagnostics but before any further reduction is invoked.
            Inputs:
            self    -- initialized instance of the instrument reduction class
            reducer -- initialized instance of the reducer
                       (DirectEnergyConversion class initialized for specific reduction)
            ws         the workspace, describing the run or partial run in multirep mode
                       to preprocess

            By default, does nothing.
            Add code to do custom preprocessing.
            Must return pointer to the preprocessed workspace
        """
        return ws
      #
    def do_postprocessing(self,reducer,ws):
        """ Custom function, applied to each reduced run or every reduced workspace,
            the run is divided into, in multirep mode.
            Applied after reduction is completed but before saving the result.

            Inputs:
            self    -- initialized instance of the instrument reduction class
            reducer -- initialized instance of the reducer
                       (DirectEnergyConversion class initialized for specific reduction)
            ws         the workspace, describing the run or partial run in multirep mode
                       after reduction to postprocess
            By default, does nothing.
            Add code to do custom postprocessing.
            Must return pointer to the postprocessed workspace.
            The postprocessed workspace should be consistent with selected save method.
            (E.g. if you decide to convert workspace units to wavelength, you can not save result as nxspe)
        """
        return ws
    #
    def set_custom_output_filename(self):
        """define custom name of output files if standard one is not satisfactory

          In addition to that, example of accessing complex reduction properties
          Simple reduction properties can be accessed as e.g.: value= prop_man.sum_runs
        """
        def custom_name(prop_man):
            """Sample function which builds filename from
              incident energy and run number and adds some auxiliary information
              to it.
            """
            # Note -- properties have the same names as the list of advanced and
            # main properties
            # Note: the properties are stored in prop_man class and accessed as
            # below.
            ei = PropertyManager.incident_energy.get_current()
            # sample run is more then just list of runs, so we use
            # the formalization below to access its methods
            if self.reducer.prop_man.filename_prefix:
                return reduced_filename(0, ei, False, self.reducer.prop_man.filename_prefix)
            else:
                runs_list = PropertyManager.sample_run.get_run_list()
                return reduced_filename(runs_list, ei, self.reducer.prop_man.sum_runs)
        # Uncomment this to use custom filename function
        return lambda : custom_name(self.reducer.prop_man)
        # Uncomment this to use standard file name generating function
        #return None

    def __init__(self,web_var=None):
        """ sets properties defaults for the instrument with Name"""
        ReductionWrapper.__init__(self,'MAR',web_var)
        object.__setattr__(self.reducer.prop_man, 'remove_streaks', False)
        object.__setattr__(self.reducer.prop_man, 'fakewb', False)
        object.__setattr__(self.reducer.prop_man, 'filename_prefix', '')
        Mt = MethodType(self.do_preprocessing, self.reducer)
        DirectEnergyConversion.__setattr__(self.reducer,'do_preprocessing',Mt)
        Mt = MethodType(self.do_postprocessing, self.reducer)
        DirectEnergyConversion.__setattr__(self.reducer,'do_postprocessing',Mt)

#------------------------------------------------------------------------------

# Defines a function to return the data file names
def reduced_filename(runs, ei, is_sum, prefix=None):
    runs = [runs] if not isinstance(runs, list) else runs
    is_sum = is_sum if len(runs) > 1 else False
    if not prefix:
        prefix = 'MAR{}to{}sum'.format(runs[0], runs[-1]) if is_sum else 'MAR{}'.format(runs[0])
    return '{}_Ei{:<3.2f}meV'.format(prefix, ei)

def iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs=False,**kwargs):
    """Helper function, which allow to run MARIReduction in old iliad way
     inputs:
        runno       -- one or list of run numbers to process
        ei            -- incident energy or list of incident energies
        wbvan      --  white beam vanadium run number or file name of the vanadium
        monovan  -- monochromatic vanadium run number or file name
        sam_mass-- mass of the sample under investigation
        sam_rmm -- rmm of sample under investigation
        sum_runs -- if true, all runs provided in runno list should be added together
        **kwargs -- list of any reduction properties, found in MARI_Parameters.xml file
                         written in the form property=value
        NOTE: to avoid duplication, all default parameters are set up within def_advanced properites
                  and def_main properties functions. They of course may be overwritten here.
    """
    rd = MARIReduction()
    # set up advanced and main properties, specified in code above
    rd.def_advanced_properties()
    rd.def_main_properties()
    prop_man = rd.reducer.prop_man

#    if not hasattr(runno, '__len__') or isinstance(runno, six.string_types):
#        runno = [runno]
    if sum_runs and len(runno)==1:
        sum_runs = False

    #assign input arguments:
    prop_man.incident_energy = ei
    prop_man.sum_runs = sum_runs
    prop_man.sample_run = runno
    prop_man.wb_run = wbvan

    multirun = False
    if hasattr(ei, '__len__') and len(ei) > 1:
        prop_man.energy_bins=[-1, 1./400., 0.97]
        multirun = True if sum_runs else False
    elif ei != 'auto' and not hasattr(ei, '__len__'):
        prop_man.energy_bins=[-1*ei, ei/400., 0.97*ei]

    if ( sam_rmm!=0 and sam_mass!=0 ) :
        prop_man.sample_mass=sam_mass
        prop_man.sample_rmm=sam_rmm
        prop_man.monovan_run=monovan
    else:
        prop_man.monovan_run=None

    outws = None
    for key,val in list(kwargs.items()):
        if key == 'save_file_name':
            if isinstance(runno, (list, tuple)) or isinstance(ei,(list, tuple)) :
                print("**************************************************************************************")
                print("*** WARNING: you can not set up single file name for list of files or list of energies")
                print("*** change ''set_custom_output_filename'' function, which returns lamda function used ")
                print("*** to calculate file name as function of each incident energy and run number.")
                print("**************************************************************************************")
                continue
        if key == 'wait_for_file':
             rd.wait_for_file = kwargs['wait_for_file']
             continue
        if key == 'OutputWorkspace':
             outws = kwargs['OutputWorkspace']
             continue
        if key == 'dos_background':
             continue

        setattr(prop_man,key,val);

    rd.reducer.prop_man = prop_man
    #rd.reducer.prop_man.save_file_name='mar'+str(runno)+'_ei'+str(int(round(ei)))
    return rd.run_reduction(outws)


class Runs(object):
    """Helper class for iliad_dos - a list of runs and associated metadata"""
    def __init__(self, run_nums, wbvan, ei, monovan=0, sam_mass=0, sam_rmm=0, sum_runs=True, **kwargs):
        self.runs = run_nums if hasattr(run_nums, '__iter__') and not isinstance(run_nums, six.string_types) else [run_nums]
        self.ei, self.wbvan, self.monovan, self.sam_mass, self.sam_rmm, self.sum_runs = (ei, wbvan, monovan, sam_mass, sam_rmm, sum_runs)
        self.kwargs = kwargs
        self.prefix = self.kwargs['filename_prefix'] if 'filename_prefix' in self.kwargs else None
        self.outputworkspace = self.kwargs.pop('OutputWorkspace', None)
        self.recalc = self.kwargs.pop('recalc', False)

    def run_iliad(self):
        if self.sam_mass == 0 or self.sam_rmm == 0:
            self.monovan = None
        if self.outputworkspace:
            if self.recalc or self.outputworkspace not in mtd.getObjectNames():
                return iliad_mari(self.runs, self.ei, self.wbvan, self.monovan, self.sam_mass, self.sam_rmm, self.sum_runs,
                                  OutputWorkspace=self.outputworkspace, **self.kwargs)
        else:
            return iliad_mari(self.runs, self.ei, self.wbvan, self.monovan, self.sam_mass, self.sam_rmm, self.sum_runs, **self.kwargs)

    def load_reduce(self, wd):
        ws = []
        for ei in self.ei if hasattr(self.ei, '__iter__') else [self.ei]:
            filename = reduced_filename(self.runs, ei, self.sum_runs, self.prefix)
            if filename in mtd.getObjectNames():
                ws.append(mtd[filename])
                continue
            for ext in PropertyManager.save_format.save_formats:
                try:
                    ws1 = Load('{}/{}.{}'.format(wd, filename, ext), OutputWorkspace=filename)
                except ValueError:
                    pass
                else:
                    ws.append(ws1)
                    continue
        return ws if len(ws) > 0 else self.run_iliad()


def reduce_runs(runs_dict, wbvan, ei, monovan, **kwargs):
    """Parses a dictionary of runs / samples / temperatures"""
    load_reduce = kwargs.pop('load_reduce', False)
    use_subdirs = kwargs.pop('use_sub_directories', False)
    wd0 = config['defaultsave.directory']
    new_dict = {}
    for sam in list(runs_dict.keys()):
        if use_subdirs:
            wd = '{}/{}'.format(wd0, sam)
            if not os.path.isdir(wd):
                os.mkdir(wd)
            config['defaultsave.directory'] = wd
        new_dict[sam] = {}
        ws_list = []
        (ei0, monovan0, sam_mass, sam_rmm) = (runs_dict[sam].pop(ky, df)
            for ky, df in list(zip(['ei', 'monovan', 'sam_mass', 'sam_rmm'], [ei, monovan, 0, 0])))
        for tt in list(runs_dict[sam].keys()):
            new_dict[sam][tt] = {}
            runobj = Runs(runs_dict[sam][tt]['data'], wbvan, ei0, monovan0, sam_mass, sam_rmm, **kwargs)
            new_dict[sam][tt]['data'] = runobj.load_reduce(wd) if load_reduce else runobj.run_iliad()
            #ws_list = ws_list + [ws.name() for ws in new_dict[sam][tt]['data']]
            #ws_list = ws_list + [ws.name()+'_SQW' for ws in new_dict[sam][tt]['data']]
            if 'background' in runs_dict[sam][tt]:
                runobj = Runs(runs_dict[sam][tt]['background'], wbvan, ei0, monovan0, sam_mass, sam_rmm, **kwargs)
                new_dict[sam][tt]['background'] = runobj.load_reduce(wd) if load_reduce else runobj.run_iliad()
                #ws_list = ws_list + [ws.name() for ws in new_dict[sam][tt]['background']]
                #ws_list = ws_list + [ws.name()+'_SQW' for ws in new_dict[sam][tt]['background']]
            #GroupWorkspaces(InputWorkspaces=ws_list, OutputWorkspace='{}_{}K_reduced'.format(sam, tt))
    return new_dict

def _parseqe(qe, ei):
    if isinstance(qe, list):
        if isinstance(qe[0], six.string_types) and len(qe)==len(ei):
            return qe
        elif isinstance(qe[0], list) and len(qe)==len(ei):
            return [','.join(v) for v in qe]
        else:
            return [','.join(qe)] * len(ei)
    else:
        return [qe] * len(ei)

def iliad_dos(runno, wbvan, ei=None, monovan=None, sam_mass=0, sam_rmm=0, sum_runs=False, **kwargs):
    """Reduces a set of data (and optionally background) runs and calculates the phonon density of states
         in the incoherent approximation from the data (or background subtracted data).
       inputs:
         runno - either a list of run numbers (in which case the next 5 parameters must be:
                    ei, wbvan, monovan, sam_mass, sam_rmm just like in the iliad_mari function)
                    in this case you must also specify the temperature keyword with the sample temperature
                    in this case you can also specify the sum_runs parameter like in iliad_mari (default: False)

                 or runno can be a python dictionary of dictionaries with the following structure:
                    runno = {'sample_name': { temperature: sample_dict, 'ei':ei, 'monovan': n, 'sam_mass': n, 'sam_rmm': y }, ... }
                    (e.g. a dictionary with keys which are sample names containing another dictionary with the
                    sample temperature as keys whose values is another dictionary with the following keys:
                      'data' - a list of data run numbers. These runs will be summed and reduced
                      'background' - an optional list of background run numbers. These will also be summed and reduced
                                     and subtracted from the data
                      'recalc' - by default this routine does not recalculate the reduction if it sees that the output
                                     workspaces are present in the *Analysis Data Service*. If this key is present and
                                     set to True, then it will force a recalculation of the reduction.
                      'ssf' - a per sample and per temperature self-shielding factor for background subtraction
                      'msd' - a per sample and per temperature mean-square displacement factor for DOS calculation
                    In addition to the sample temperature in the samples dict, you can also provide the following optional keys:
                      'ei' - the incident energ(y)(ies) of the runs. If you don't provide the ei(s) in the keyword arguments to
                                     iliad_dos, it must be provided on a per dataset basis.
                      'monovan' - the run number of a vanadium calibration run with the same spectrometer
                                     setting as the sample and data for absolute units normalisation.
                      'sam_mass' - if 'monovan' is set you must provide this key, which is the sample mass in g (ignored if monovan not set)
                      'sam_rmm' - if 'monovan' is set you must provide this key, which is the sample molar mass (ignored if monovan not set)
         wbvan - the white beam vanadium run number (mandatory)
         ei - either a number or a list of the incident energies in the measurement. If you provide the ei here it will be assumed
                    that all runs have this ei. If you have measured different samples / temperatures with different ei's you have
                    to use the runno dictionary input and give the ei in each samples' dictionary.
         monovan - the monochromatic vanadium run number for absolute units calibration (assuming all runs have the same
                    ei, otherwise this should also be in the samples' dictionaries). **Note that if you must also define the sample
                    mass and molar mass in the sample's dictionary otherwise this option will be ignored.**

       In addition this function understands the following keyword arguments (and will pass on other keyword args to iliad_mari):

         ssf - the global self shielding factor for background subtraction (default: 1.) This is overriden by any SSF defined in the runno dict
         msd - the global mean square displacement for DOS calculation (default: 0.) This is overriden by any MSD defined in the runno dict
         qrange - a string or list of string or list of lists of two numbers denoting the |Q| range to sum over for the DOS calculation.
            if it is a list of strings or list of list it must be the same size as the number of incident energies and corresponds to those.
            (default: Qmax/3 to Qmax at the elastic line)
         ebins - a string or list of string or list of lists of three numbers denoting the energy transfer bins for the DOS calculation.
            if it is a list of strings or list of list it must be the same size as the number of incident energies and corresponds to those.
            (default: Emax/10 to Emax*0.95 in steps of Emax/100)
         temperature - if runno is not a dictionary, you *must* specify the sample temperature using this keyword argument
         background - if runno is not a dictionary, you can specify the list of background runs here
         load_reduce - if this is set to True, the function will try to load in the reduced data files rather than recalculate them
                    Note that this option will override any 'recalc' keys in the runno dict (if you want to force recalculation
                    set this to False or omit this keyword altogether).
         save_text - if True this will save the calculated DOS as 3-column x,y,error text files.
         nsmooth - if set, this will apply an n-point moving average filter to the calculated DOS creating another file/workspace
                    nsmooth should be an odd number greater than 2. (Default: None - do not apply smoothing).
         save_folder - if set the function will save the reduce data to this folder instead of the Mantid default folder
         use_sub_directories - if set and is True then for each sample create a new subdirectory and save its file there

         E.g.:
            iliad_dos([25000, 25001], ei=[120, 10], wbvan=25035, background=[25004, 25005], temperature=5)
         will run the reduction for one set of data files with background subtraction and calculate the DOS at 5K.

            iliad_dos({'sam1':
                              {5: {'data'=[25000,25001], 'background'=[25004,25005]},
                               300: {'data'=[25002,25003], 'background=[25006,25007]},
                               'sam_mass':10, 'sam_rmm':177.77},
                       'sam2':
                              {10: {'data'=[25010,25011], 'background'=[25014,25015]},
                               600: {'data'=[25012,25013], 'background=[25016,25017]},
                               'sam_mass':8, 'sam_rmm':187.77},
                      }, ei=[120,10], wbvan=25035, monovan=25008)
         will run the reduction for two sets of samples (one at 5K and 300K, one at 10K and 600K), and calculate the
         density of states for the four sets of measurements, normalising to absolute units. All runs are with Ei=120 and 10meV.
    """
    # Parses the input
    save_text = kwargs.pop('save_text', False)
    nsmooth = kwargs.pop('nsmooth', None)
    save_folder = kwargs.pop('save_folder', None)
    use_subdirs = kwargs['use_sub_directories'] if 'use_sub_directories' in kwargs else False
    global_ssf = kwargs.pop('ssf', 1.0)
    global_msd = kwargs.pop('msd', 0.0)
    global_qrange = kwargs.pop('qrange', 'Qmax/3, Qmax')
    global_ebins = kwargs.pop('ebins', 'Emax/10, Emax/100, Emax*0.95')
    global_ei = ei
    oldwd = config['defaultsave.directory']
    wd0 = save_folder if save_folder is not None else oldwd
    config['defaultsave.directory'] = wd0

    # Runs the reduction
    if isinstance(runno, dict):
        runs_dict = runno
    else:
        if not hasattr(runno, '__len__') or isinstance(runno, six.string_types):
            runno = [runno]
        if sum_runs and len(runno)==1:
            sum_runs = False

        if 'temperature' not in kwargs:
            raise ValueError('No sample temperature given')
        temperature = kwargs.pop('temperature')
        if ei is None:
            raise ValueError('Incident energy not defined')
        if sum_runs:
            runs_dict = {None: {temperature: {'data':runno}}}
            if 'background' in kwargs:
                runs_dict[None][temperature]['background'] = kwargs.pop('background')
        else:
            runs_dict = {'MAR{}'.format(run): {temperature: {'data':run}} for run in runno}
            if monovan and sam_mass:
                for ky in list(runs_dict.keys()):
                    runs_dict[ky]['sam_mass'] = sam_mass
            if monovan and sam_rmm:
                for ky in list(runs_dict.keys()):
                    runs_dict[ky]['sam_rmm'] = sam_rmm
            if 'background' in kwargs:
                background = kwargs.pop('background')
                for idx, run in enumerate(runno):
                    runs_dict['MAR{}'.format(run)][temperature]['background'] = background[idx]
    ws_dict = reduce_runs(runs_dict, wbvan, ei, monovan, **kwargs)

    # Calculates the DOS (with optional background subtraction)
    for sam in list(ws_dict.keys()):
        if use_subdirs:
            wd = '{}/{}'.format(wd0, sam)
            if not os.path.isdir(wd):
                os.mkdir(wd)
            config['defaultsave.directory'] = wd
        for tt in list(ws_dict[sam].keys()):
            def_ei = runs_dict[sam][tt]['ei'] if 'ei' in runs_dict[sam][tt] else global_ei
            if not hasattr(def_ei, '__len__'):
                def_ei = [def_ei]
            ws_ei = [ws.getEFixed(1) for ws in ws_dict[sam][tt]['data']]
            id_ei = [np.argsort([np.abs(ei1-ei0) for ei1 in ws_ei])[0] for ei0 in def_ei]
            data_ws = [ws_dict[sam][tt]['data'][id_ei[ii]] for ii in range(len(ws_ei))]
            msd = runs_dict[sam][tt]['msd'] if 'msd' in runs_dict[sam][tt] else global_msd
            # Calculates the sample DOS (without background subtraction)
            qstr = _parseqe(runs_dict[sam][tt]['qrange'] if 'qrange' in runs_dict[sam][tt] else global_qrange, def_ei)
            estr = _parseqe(runs_dict[sam][tt]['ebins'] if 'ebins' in runs_dict[sam][tt] else global_ebins, def_ei)
            for ws, ei, qq, ee in list(zip(data_ws, def_ei, qstr, estr)):
                if ws.name()+'_SQW' not in mtd.getObjectNames():
                    q_min, q_max = tuple([v*sqrt(ei) for v in [0.04, 1.3]])
                    ws_sqw = SofQW3(ws, '{},{},{}'.format(q_min, q_max/285., q_max), EMode='Direct', OutputWorkspace=ws.name()+'_SQW')
                else:
                    ws_sqw = mtd[ws.name()+'_SQW']
                ws_dos = ComputeIncoherentDOS(ws_sqw, tt, msd, qq, ee, OutputWorkspace='{}_{}K_Ei{}_data_DOS'.format(sam, tt, ei))
                if save_text:
                    SaveAscii(ws_dos, ws_dos.name()+'.txt', Separator='Space')
                if nsmooth > 2:
                    SmoothData(ws_dos, nsmooth, OutputWorkspace=ws_dos.name()+'_smooth')
                    if save_text:
                        SaveAscii(ws_dos.name()+'_smooth', ws_dos.name()+'_smooth.txt', Separator='Space')
            if 'background' in list(ws_dict[sam][tt].keys()):
                bkg_ei = [ws.getEFixed(1) for ws in ws_dict[sam][tt]['background']]
                id_bkg_ei = [np.argsort([np.abs(ei1-ei0) for ei1 in bkg_ei])[0] for ei0 in def_ei]
                bkg_ws = [ws_dict[sam][tt]['background'][id_bkg_ei[ii]] for ii in range(len(bkg_ei))]
                ssf = runs_dict[sam][tt]['ssf'] if 'ssf' in runs_dict[sam][tt] else global_ssf
                sub_ws = [data_ws[ii] - ssf*bkg_ws[ii] for ii in range(len(bkg_ws))]
                for ws, ei, qq, ee in list(zip(sub_ws, def_ei, qstr, estr)):
                    SaveNXSPE(ws, '{}_{}K_Ei{:.2f}meV_subtracted.nxspe'.format(sam, tt, ei))
                    q_min, q_max = tuple([v*sqrt(ei) for v in [0.04, 1.3]])
                    ws_sqw = SofQW3(ws, '{},{},{}'.format(q_min, q_max/285., q_max), EMode='Direct', OutputWorkspace=ws.name()+'_SQW')
                    ws_dos = ComputeIncoherentDOS(ws_sqw, tt, msd, qq, ee, OutputWorkspace='{}_{}K_Ei{}_subtracted_DOS'.format(sam, tt, ei))
                    if save_text:
                        SaveAscii(ws_dos, ws_dos.name()+'.txt', Separator='Space')
                    if nsmooth > 2:
                        ws_dos_smooth = SmoothData(ws_dos, nsmooth, OutputWorkspace=ws_dos.name()+'_smooth')
                        if save_text:
                            SaveAscii(ws_dos_smooth, ws_dos_smooth.name()+'.txt', Separator='Space')
    config['defaultsave.directory'] = oldwd

def mari_normalise_background(background_int, white_int, second_white_int=None):
    """Normalize the background integrals"""
    if second_white_int is None:
        if background_int.getNumberHistograms() == 919 and white_int.getNumberHistograms() == 918:
            background_int = CropWorkspace(background_int, StartWorkspaceIndex=1)
        background_int =  Divide(LHSWorkspace=background_int,RHSWorkspace=white_int,WarnOnZeroDivide='0')
    else:
        hmean = 2.0*white_int*second_white_int/(white_int+second_white_int)
        background_int =  Divide(LHSWorkspace=background_int,RHSWorkspace=hmean,WarnOnZeroDivide='0')
        DeleteWorkspace(hmean)

if __name__ == "__main__" or __name__ == "__builtin__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    #map_mask_dir = '/usr/local/mprogs/InstrumentFiles/maps'
    # folder where input data can be found
    #data_dir = r'\\isis\inst$\NDXMARI\Instrument\data\cycle_14_2'
    #config.appendDataSearchDir(map_mask_dir)
    #config.appendDataSearchDir(data_dir)

    root=os.path.dirname(os.path.realpath(__file__))
    #data_dir = os.path.join(root,r'data')

    #config.appendDataSearchDir(root)
    #config.appendDataSearchDir(data_dir)
    #config['defaultsave.directory']=root
###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = MARIReduction()
    rd.def_advanced_properties()
    rd.def_main_properties()

#### uncomment rows below to generate web variables and save then to transfer to ###
    ## web services.
    run_dir = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(run_dir,'reduce_vars.py')
    rd.save_web_variables(file)

#### Set up time interval (sec) for reducer to check for input data file.  ####
    #  If this file is not present and this value is 0,reduction fails
    #  if this value >0 the reduction waits until file appears on the data
    #  search path checking after time specified below.
    rd.wait_for_file = 0  # waiting time interval in seconds


####get reduction parameters from properties above, override what you want locally ###
   # and run reduction.  Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g.
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or
   ## rd.reducer.prop_man.sum_runs = False
   #
###### Run reduction over all run numbers or files assigned to ######
    # sample_run variable

    # return output workspace only if you are going to do
    # something with it here.  Running range of runs will return the array
    # of workspace pointers.
    #red_ws = rd.run_reduction()
    # usual way to go is to reduce workspace and save it internally

    rd.run_reduction()

###### Test absorption corrections to find optimal settings for corrections algorithm         ######
#     corr = rd.eval_absorption_corrections()
