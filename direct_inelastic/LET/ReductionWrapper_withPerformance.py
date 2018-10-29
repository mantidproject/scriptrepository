# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#pylint: disable=invalid-name
from __future__ import (absolute_import, division, print_function)
from mantid.simpleapi import *
from mantid import config,api
from mantid.kernel import funcinspect

from Direct.PropertyManager import PropertyManager
# this import is used by children
from Direct.ReductionWrapper import *
import os
import re
import time
import platform
import subprocess
from six import iteritems
from abc import abstractmethod
import resource
import math

# R0921 abstract class not referenced -- wrong, client references it.
# pylint: disable=too-many-instance-attributes, R0921


class ReductionWrapper_withPerformance(ReductionWrapper):
    """ Abstract class provides interface to direct inelastic reduction
        allowing it to be run  from Mantid, web services, or system tests
        using the same interface and the same run file placed in different
        locations.
    """

    def __init__(self,instrumentName,web_var=None):
        """ sets properties defaults for the instrument with Name
          and define if wrapper runs from web services or not
        """
        super(ReductionWrapper_withPerformance, self).__init__(instrumentName,web_var)        
    
      # The variables responsible for measuring the global performance of a reduction script
        self._measure_perfrormance = False
        self._performance_results_file = 'dummy_performance.txt';

 
    def run_reduction(self):
        """" Reduces runs one by one or sum all them together and reduce after this

            if wait_for_file time is > 0, it will until  missing files appear on the
            data search path
        """
        class time_logger:
            """" helper class to run using Python 'with' statement
                 and log reduction execution time and optionally
                 time spent on reducing every input file.
            """
            def __init__(self,filename):
                "filename -- short name of the log file"
                self._start_time = time.time()
                self._tick_time = self._start_time
                filepath = os.path.dirname(os.path.realpath(__file__))
                self._log_file = os.path.join(filepath,filename)
                self.fh = -1;
            def __enter__(self):
                self.fh = open(self._log_file,"w")
                lt= time.localtime(self._start_time )        
                res_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024);
                ch_memory = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/(1024)
                self.fh.write("*** Test started:  {0}/{1}/{2} at {3}:{4}\n".format(lt.tm_mday,lt.tm_mon,lt.tm_year,lt.tm_hour,lt.tm_min))
                self.fh.write("*** Self Memory: {0}Kb; Kids memory {1}Kb\n".format(res_memory,ch_memory))
                pv = subprocess.check_output(['free','-m'])
                pvs = pv.split('\n')                
                self.fh.write("***      {0}\n".format(pvs[0]))                
                self.fh.write("***      {0}\n".format(pvs[1]))
                self.fh.write("***      {0}\n".format(pvs[2]))
                
                self.fh.flush()
                return self
            def __exit__(self, type, value, traceback):
                fin_time = time.time()
                lt = time.localtime(fin_time)
                self.fh.write("*** Test finished: {0}/{1}/{2} at {3}:{4}\n".format(lt.tm_mday,lt.tm_mon,lt.tm_year,lt.tm_hour,lt.tm_min))                
                self.fh.write("*** Total execution time: {0:.2f}(sec)\n".format(fin_time -self._start_time))
                self.fh.close()
            def tick(self,fileID):
                start_time = self._tick_time;
                end_time  = time.time()
                res_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024)           
                ch_memory = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/(1024)
                self.fh.write("*** ---> File {0} processed in {1:.2f}sec\n".format(fileID,end_time-start_time))                
                self.fh.write("*** Self Memory: {0}Kb; Kids memory {1}Kb\n".format(math.ceil(res_memory),math.ceil(ch_memory)))
                try:
                    pv = subprocess.check_output(['free','-m'])
                    pvs = pv.split('\n')                
                    self.fh.write("***      {0}\n".format(pvs[1]))
                    self.fh.write("***      {0}\n".format(pvs[2]))
                except:
                    #ClearCache(True,True,True,True,True,True,True)                    
                    #self.fh.write("***      Can not launch subprocess to evaluate free memory. Clearing all Mantid Caches\n")
                    self.fh.write("***      Can not launch subprocess to evaluate free memory.\n")
                    
                self.fh.flush()
                self._tick_time = end_time
                
        try:
            _,r = funcinspect.lhs_info('both')
            out_ws_name = r[0]
# no-exception-type(s) specified. Who knows what exception this internal procedure rises...
#pylint: disable=W0702
        except:
            out_ws_name = None
        host = platform.node()
        host = host.replace('.','_')
        host = host.replace('-','_')
        inst  = config.getInstrument()
		count = 1
		inst_name = inst.name()
        log_file_name = "{0}_performance_{1}_test{2}.txt".format(inst_name,host,count)
        filepath = os.path.dirname(os.path.realpath(__file__))		
		ff = os.path.join(filepath,log_file_name)
		while os.path.isfile(ff):
			count = count+1;
			log_file_name = "{0}_performance_{1}_test{2}.txt".format(inst_name,host,count)
			ff = os.path.join(filepath,log_file_name)
			
        print (' ******************* storing performance data to file: ',log_file_name)
        if self.reducer.sum_runs:
# --------### sum runs provided ------------------------------------###
            with time_logger(log_file_name) as log:
                if out_ws_name is None:
                    self.sum_and_reduce()
                    return None
                else:
                    red_ws = self.sum_and_reduce()
                    RenameWorkspace(InputWorkspace=red_ws,OutputWorkspace=out_ws_name)
                    return mtd[out_ws_name]
        else:
# --------### reduce list of runs one by one ----------------------------###
            runfiles = PropertyManager.sample_run.get_run_file_list()
            with time_logger(log_file_name) as log:            
                if out_ws_name is None:
                    for file_name in runfiles:
                        self.reduce(file_name)
                        log.tick(file_name)
                    return None
                else:
                    results = []
                    nruns = len(runfiles)
                    for num,file_name in enumerate(runfiles):
                        red_ws = self.reduce(file_name)
                        log.tick(file_name)
                        if isinstance(red_ws,list):
                            for ws in red_ws:
                                results.append(ws)
                        else:
                            if nruns == 1:
                                if red_ws.name() != out_ws_name:
                                    RenameWorkspace(InputWorkspace=red_ws,OutputWorkspace=out_ws_name)
                                results.append(mtd[out_ws_name])
                            else:
                                OutWSName = '{0}#{1}of{2}'.format(out_ws_name,num+1,nruns)
                                if red_ws.name() != out_ws_name:
                                    RenameWorkspace(InputWorkspace=red_ws,OutputWorkspace=OutWSName)
                                results.append(mtd[OutWSName])
                #end
                if len(results) == 1:
                    return results[0]
                else:
                    return results
                #end if
            #end if
        #end


if __name__ == "__main__":
    pass
