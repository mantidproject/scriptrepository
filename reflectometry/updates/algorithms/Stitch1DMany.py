"""*WIKI* 

Stitches single histogram [[MatrixWorkspace|Matrix Workspaces]] together outputting a stitched Matrix Workspace. This algorithm is a wrapper over [[Stitch1D]].

The algorithm expects  pairs of StartOverlaps and EndOverlaps values. The order in which these are provided determines the pairing. 
There should be N entries in each of these StartOverlaps and EndOverlaps lists, where N = 1 -(No of workspaces to stitch). 
StartOverlaps and EndOverlaps are in the same units as the X-axis for the workspace and are optional.

The workspaces must be histogrammed. Use [[ConvertToHistogram]] on workspaces prior to passing them to this algorithm.
*WIKI*"""
#from mantid.simpleapi import *

from mantid.simpleapi import *
from mantid.api import *
from mantid.kernel import *
import numpy as np

class Stitch1DMany(PythonAlgorithm):

    def category(self):
        return "Reflectometry\\Developmental"
    
    def version(self):
        return 2

    def name(self):
        return "Stitch1D"

    def PyInit(self):
    
        input_validator = StringMandatoryValidator()
    
        self.declareProperty(name="InputWorkspaces", defaultValue="", direction=Direction.Input, validator=input_validator, doc="Input workspaces")
        self.declareProperty(WorkspaceProperty("OutputWorkspace", "", Direction.Output), "Output stitched workspace")
        
        self.declareProperty(FloatArrayProperty(name="StartOverlaps", values=[]), doc="Overlap in Q.")
        self.declareProperty(FloatArrayProperty(name="EndOverlaps", values=[]), doc="End overlap in Q.")
        self.declareProperty(FloatArrayProperty(name="Params", validator=FloatArrayMandatoryValidator()), doc="Rebinning Parameters. See Rebin for format.")
        self.declareProperty(name="ScaleRHSWorkspace", defaultValue=True, doc="Scaling either with respect to workspace 1 or workspace 2.")
        self.declareProperty(name="UseManualScaleFactor", defaultValue=False, doc="True to use a provided value for the scale factor.")
        self.declareProperty(name="ManualScaleFactor", defaultValue=1.0, doc="Provided value for the scale factor.")
        self.declareProperty(name="OutScaleFactor", defaultValue=-2.0, direction = Direction.Output, doc="The actual used value for the scaling factor.")
        
    def __workspace_from_split_name(self, list_of_names, index):
        return mtd[list_of_names[index].strip()]
        
    def __workspaces_from_split_name(self, list_of_names):
        workspaces = list()
        for name in list_of_names:
            workspaces.append(mtd[name.strip()])
        return workspaces
            
    '''
    If the property value has been provided, use that. If default, then create an array of Nones so that Stitch1D can still be correctly looped over.
    '''
    def __input_or_safe_default(self, property_name, n_entries):
        property = self.getProperty(property_name)
        property_value = property.value
        if property.isDefault:
            property_value = list()
            for i in range(0, n_entries):
                property_value.append(None)
        return property_value
                
       
    def __do_stitch_workspace(self, lhs_ws, rhs_ws, start_overlap, end_overlap, params, scale_rhs_ws, use_manual_scale_factor, manual_scale_factor):     
        out_name = lhs_ws.name() + rhs_ws.name()
        
        alg = self.createChildAlgorithm("Stitch1D")
        alg.initialize()
        alg.setProperty("LHSWorkspace", lhs_ws)
        alg.setProperty("RHSWorkspace", rhs_ws)
        if start_overlap:
            alg.setProperty("StartOverlap", start_overlap)
        if end_overlap:
            alg.setProperty("EndOverlap", end_overlap)
        alg.setProperty("Params", params)
        alg.setProperty("ScaleRHSWorkspace", scale_rhs_ws)
        alg.setProperty("UseManualScaleFactor", use_manual_scale_factor)
        if manual_scale_factor:
            alg.setProperty("ManualScaleFactor", manual_scale_factor)
        alg.setProperty("OutputWorkspace", "from_sub_alg" + out_name)
        alg.execute()
        out_ws = alg.getProperty("OutputWorkspace").value
        scale_factor = alg.getProperty("OutScaleFactor").value
        
        return (out_ws, scale_factor)
    
    def __check_workspaces_are_common(self, input_workspace_names):
        workspaces = self.__workspaces_from_split_name(input_workspace_names)
        exemplar = workspaces[0]
        for i in range(1, len(workspaces)):
            test_ws = workspaces[i]
            if type(exemplar) != type(test_ws):
                raise RuntimeError("Input Workspaces must all be of the same type.")
            if isinstance(test_ws, WorkspaceGroup):
                if test_ws.size() != exemplar.size():
                    raise RuntimeError("Group Workspaces as InputWorkspaces must have the same number of sub-workspaces.")
    
    def __are_processing_groups(self, input_workspace_names):
        test_ws = self.__workspace_from_split_name(input_workspace_names, 0)   
        return isinstance(test_ws, WorkspaceGroup)
            
    def PyExec(self):
    
        inputWorkspaces = self.getProperty("InputWorkspaces").value
        inputWorkspaces = inputWorkspaces.split(',')
        numberOfWorkspaces = len(inputWorkspaces)
        
        startOverlaps = self.__input_or_safe_default('StartOverlaps', numberOfWorkspaces-1)
        endOverlaps = self.__input_or_safe_default('EndOverlaps', numberOfWorkspaces-1)
        
        # Just forward the other properties on.
        scaleRHSWorkspace = self.getProperty('ScaleRHSWorkspace').value
        useManualScaleFactor = self.getProperty('UseManualScaleFactor').value
        manualScaleFactor = self.getProperty('ManualScaleFactor').value
        params = self.getProperty("Params").value
        
        
        if not numberOfWorkspaces > 1:
            raise ValueError("Too few workspaces to stitch")
        if not (len(startOverlaps) == len(endOverlaps)):
            raise ValueError("StartOverlaps and EndOverlaps are different lengths")
        if not (len(startOverlaps) == (numberOfWorkspaces- 1)):
            raise ValueError("Wrong number of StartOverlaps, should be %i not %i" % (numberOfWorkspaces - 1, startOverlaps))
        self.__check_workspaces_are_common(inputWorkspaces)
        
        scaleFactor = None
        comma_separator = ","
        no_separator = str()
        
        # Identify and process as group workspaces
        if self.__are_processing_groups(inputWorkspaces):
            workspace_groups = self.__workspaces_from_split_name(inputWorkspaces)
                
            out_group_separator = no_separator
            out_group_workspaces = str()
                
            n_sub_workspaces = workspace_groups[0].size() 
            for i in range(n_sub_workspaces):
                
                to_process = str()
                out_name = str()
                separator = no_separator
                
                for j in range(0, numberOfWorkspaces, 1):
                        
                    to_process += separator + workspace_groups[j][i].name()
                    out_name += workspace_groups[j][i].name()
                    separator=comma_separator
                
                startOverlaps = self.getProperty("StartOverlaps").value
                endOverlaps = self.getProperty("EndOverlaps").value
                
                stitched, scaleFactor = Stitch1DMany(InputWorkspaces=to_process, OutputWorkspace=out_name, StartOverlaps=startOverlaps, EndOverlaps=endOverlaps, 
                                                         Params=params, ScaleRHSWorkspace=scaleRHSWorkspace, UseManualScaleFactor=useManualScaleFactor,  
                                                         ManualScaleFactor=manualScaleFactor)
                       
                out_group_workspaces += out_group_separator + out_name
                out_group_separator = comma_separator
             
            out_workspace_name = self.getPropertyValue("OutputWorkspace")        
            out_group = GroupWorkspaces(InputWorkspaces=out_group_workspaces, OutputWorkspace=out_workspace_name)
            self.setProperty("OutputWorkspace", out_group)   
    
        else:
                
            # Iterate forward through the workspaces
            if scaleRHSWorkspace:
                lhsWS = self.__workspace_from_split_name(inputWorkspaces, 0)   
    
                for i in range(1, numberOfWorkspaces, 1):
                    rhsWS = self.__workspace_from_split_name(inputWorkspaces, i)
                    lhsWS, scaleFactor = self.__do_stitch_workspace(lhsWS, rhsWS, startOverlaps[i-1], endOverlaps[i-1], params, scaleRHSWorkspace,  useManualScaleFactor, manualScaleFactor)
                self.setProperty('OutputWorkspace', lhsWS)
                
            # Iterate backwards through the workspaces.
            else:
                rhsWS = self.__workspace_from_split_name(inputWorkspaces, -1) 
                for i in range(0, numberOfWorkspaces-1, 1):
                    lhsWS = self.__workspace_from_split_name(inputWorkspaces, i)
                    rhsWS, scaleFactor = Stitch1D(LHSWorkspace=lhsWS, RHSWorkspace=rhsWS, StartOverlap=startOverlaps[i-1], EndOverlap=endOverlaps[i-1], Params=params, ScaleRHSWorkspace=scaleRHSWorkspace, UseManualScaleFactor=useManualScaleFactor,  ManualScaleFactor=manualScaleFactor)            
                self.setProperty('OutputWorkspace', rhsWS)
        
        self.setProperty('OutScaleFactor', scaleFactor)
        return None
        

#############################################################################################

import mantid
ver_str = mantid.__version__.split('.')
major = int(ver_str[0])
minor = int(ver_str[1])
if major <= 3 and minor <= 1:
    AlgorithmFactory.subscribe(Stitch1DMany)
else:
    logger.warning("Custom Stitch1DMany algorithm from script repository will not be used.")
    
