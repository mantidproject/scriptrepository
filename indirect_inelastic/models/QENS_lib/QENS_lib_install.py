#run this in Mantid scripting window
from __future__ import print_function
import subprocess

"""
You will need to download/checkout the QENS library: https://github.com/QENSlibrary/QENSmodels
The path below needs to be the full path to where you have saved the above.
Run this script and the QENSmodels should be available as part of Mantid
"""
path_to_QENS_lib = "C:\\Users\BTR75544\work\QENSlib\QENSmodels"

# install phase
print(subprocess.Popen("python -m pip install -U --no-deps "+path_to_QENS_lib ,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate())

