''' Temporary solutions to the messy problem of importing F2Py libraries into
the Indirect scripts depending on platform and numpy package version.

We also deal with importing the mantidplot module outside of MantidPlot here.
'''

import numpy
import platform
import sys

def import_mantidplot():
    ''' Currently, all scripts in the PythonAlgorithms directory are imported 
    during system tests.  Unfortunately, these tests are run outside of 
    MantidPlot and so are incompatible with scripts that import the 
    "mantidplot" module.  As a result, an error message is dumped to the
    results log for each PythonAlgorithm in the directory that imports
    mantidplot, for each and every test.
    
    Here, we silently catch all ImportErrors so that this does not occur.
    
    @returns the mantidplot module.
    '''
    try:
        import mantidplot
        return mantidplot
    except ImportError:
        # Not a problem since we are only in a system test anyway, and these
        # scripts are not needed there.
        return None

def _os_env():
    return platform.system() + platform.architecture()[0]

def _lib_suffix():
    if platform.system() == "Windows":
        suffix = "win"
    elif platform.system() == "Linux":
        suffix = "lnx"
    else:
        return ""
    return "_" + suffix + platform.architecture()[0][0:2]

def _numpy_ver():
    return numpy.version.short_version

def _linux_distro_name():
    return platform.linux_distribution()[0]

def _linux_distro_version():
    return platform.linux_distribution()[1]

def unsupported_message():
    sys.exit('F2Py functionality not currently available on your platform.')

def is_supported_f2py_platform():
    ''' We check for numpy version, as if Linux we check its distro and version
    as well.
    
    @returns True if we are currently on a platform that supports the F2Py
    libraries, else False.
    '''
    if _os_env().startswith("Windows") and _numpy_ver() == "1.6.2":
        return True
    if _os_env() == "Linux64bit" and \
       _linux_distro_name()[0:24] == "Ubuntu" and \
       _linux_distro_version() == "16.04" and \
       _numpy_ver() == "1.11.0":
        return True
    return False

def import_f2py(lib_base_name):
    ''' Until we can include the compilation process of Indirect F2Py modules 
    into the automated build of Mantid, we are forced to compile the libraries
    separately on every platform we wish to support.
    
    Here, we provide a centralised method through which we can import these
    modules, which hopefully makes the other Indirect scripts a lot less messy.
    
    @param lib_base_name :: is the prefix of the library name.  For example, 
    the CEfit_lnx64.so and CEfit_win32.pyd libraries share the same base name 
    of "CEfit".
    
    @returns the imported module.
    '''
    # Only proceed if we are indeed on one of the supported platforms.
    assert is_supported_f2py_platform()
    
    lib_name = lib_base_name + _lib_suffix()
    
    return __import__(lib_name)

def run_f2py_compatibility_test():
    ''' Convenience method that raises an exception should a user try to run
    the F2Py libraries on an incompatible platform.
    '''
    if not is_supported_f2py_platform():
        raise RuntimeError("F2Py programs NOT available on this platform.")
