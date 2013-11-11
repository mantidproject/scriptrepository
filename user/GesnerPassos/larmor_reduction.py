from mantid.simpleapi import * 
import ISISCommandInterface as ici
import SANSUtility as su


def hasCan(reducer):
    return reducer.background_subtracter is not None


def bankReduction(reducer):
    """ The reduction comprises: 
         - Reducing sample
         - Reducing can
         - sample - can
    """
    
    reducer.pre_process()
    outputname = reducer.output_wksp
    
    reduced_sample = _reduceSANSData(outputname, reducer)
    trans_log = _getTransmissionLogValue(reduced_sample)

    if hasCan(reducer):
        
        reduced_sample = RenameWorkspace(reduced_sample, OutputWorkspace=str(reduced_sample) + '_sam_tmp')
        
        reducer._process_can = True
        
        reduced_can = _reduceSANSData(outputname + '_can_tmp' , reducer)
        
        trans_can_log = _getTransmissionLogValue(reduced_can)
        
        reducer._process_can = False
        
        reduced_sample = Minus(reduced_sample, reduced_can, OutputWorkspace=outputname)
        
        _setLogValue(reduced_sample, 'TransmissionCan', trans_can_log)

    _setLogValue(reduced_sample, 'Transmission', trans_log)
    _setLogValue(reduced_sample, 'UserFile', reducer.user_settings.filename)

    reduced_sample = step_tideUpResult(reduced_sample, reducer)
    
    reducer.output_wksp = outputname

    return reduced_sample



def _reduceSANSData(outname, reducer):
    """ The core of the reduction"""
    reducer.output_wksp = outname

    # crop, mask and adjust units
    cropped = step_crop2detector(reducer.output_wksp, reducer)
    masked = step_maskDetectors(cropped, reducer)
    sample_wav = step_fromTOF2Wavelenght(masked, reducer)
    
    # calculate transmissions
    monitor_norm = step_normalize2Monitor(sample_wav, reducer)    
    trans = step_calculateTransmission(reducer)

    # geometry corrections
    sample_wav = step_rescale2AbsoluteUnit(sample_wav, reducer)
    sample_ws = step_applyGeometryCorrection(sample_wav, reducer)

    # convert to Q
    reduced = step_convert2QISIS(sample_ws, monitor_norm, trans, reducer)
    
    if trans:
        AddSampleLog(reduced, LogName='Transmission',LogText=str(trans)+'_unfitted')

    return mtd[reduced]

####
# Auxiliary methods to get and set logs
####


def _getTransmissionLogValue(ws, key='Transmission'):
    if ws.getRun().hasProperty(key):
        return ws.getRun().getLogData(key).value
    return ""

def _setLogValue(ws, key, value):
    if value:
        AddSampleLog(ws, LogName=key, LogText=value)


###############################
#  REDUCTION STEPS
#  Wrappers to the real steps called. You can see these 
#  reduction steps inside the file: isis_reducer.py 
#  And the implementation of these steps is inside isis_reduction_steps
###############################

def step_crop2detector(ws, reducer):
    """ execute isis_reduction_steps.CropDetBank """
    logger.notice(str(inspect.stack()[0][3]))
    cropper = reducer.crop_detector
    cropper.execute(reducer, ws)
    return str(ws)


def step_maskDetectors(ws, reducer):
    """ execute isis_reduction_steps.Mask_ISIS"""
    logger.notice(str(inspect.stack()[0][3]))
    masker = reducer.mask
    masker.execute(reducer, ws)
    return str(ws)

def step_fromTOF2Wavelenght(ws, reducer):
    """ execute isis_reduction_steps.UnitsConvert('Wavelength')"""
    logger.notice(str(inspect.stack()[0][3]))
    converter = reducer.to_wavelen
    converter.execute(reducer, ws)
    return str(ws)

def step_rescale2AbsoluteUnit(ws, reducer):
    """ execute isis_reduction_steps.AbsoluteUnitsISIS"""
    logger.notice(str(inspect.stack()[0][3]))
    rescale = reducer._corr_and_scale
    rescale.execute(reducer, ws)
    return str(ws)

def step_normalize2Monitor(ws, reducer):
    """ execute isis_reduction_steps.NormalizeToMonitor"""
    logger.notice(str(inspect.stack()[0][3]))
    norm = reducer.norm_mon
    norm.execute(reducer, ws)
    return reducer.norm_mon.output_wksp

def step_calculateTransmission(reducer):
    """ execute isis_reduction_steps.TransmissionCalc(loader=None)"""
    logger.notice(str(inspect.stack()[0][3]))
    transcalc = reducer.transmission_calculator
    transcalc.execute(reducer, None)
    return transcalc.output_wksp


def step_applyGeometryCorrection(ws, reducer):
    """ execute sans_reduction_steps.SampleGeomCor(self._sample_run.geometry)
    
    It is defined in scripts/reduction/instruments/sans/sans_reduction_steps.py
    """
    logger.notice(str(inspect.stack()[0][3]))
    correcter = reducer.geometry_correcter
    correcter.execute(reducer, ws)
    return str(ws)


def step_convert2QISIS(sample_ws, monitor_norm, trans, reducer):
    """ execute isis_reduction_steps.ConvertToQISIS(isis_reduction_steps.CalculateNormISIS(
                            [self.norm_mon, self.transmission_calculator])
    
    It is the main command that will eventually execute the Q1D algorithm
    """			    
    logger.notice(str(inspect.stack()[0][3]))
    converter = reducer.to_Q
    converter.execute(reducer, sample_ws)
    return str(sample_ws)


def step_tideUpResult(ws, reducer):
    """ execute sans_reduction_steps.StripEndNans"""    
    logger.notice(str(inspect.stack()[0][3]))
    tide = reducer._rem_nans
    name = str(ws)
    tide.execute(reducer, name)
    return mtd[name]


###############################################
##  REDUCTION TESTS
#################################################


###########
# required
############
ici.LARMOR()
ici.MaskFile('~/MaskLarmor.txt')

ici.AssignSample('73') 

##########
# optional
##########
# TransmissionSample (transmission, direct)
ici.TransmissionSample('73', '53')
# AssignCan
ici.AssignCan('70')
# TransmissionCan (transmission, direct)
ici.TransmissionCan('70','72')

############
# required
############
reduced = bankReduction(ici.ReductionSingleton())
