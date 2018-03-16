from __future__ import (absolute_import, division, print_function)
import os
from muesr.settings import config
"""
This is to set up the python paths for Vesta to be called from scripts
VESTA: http://jp-minerals.org/vesta/en/download.html

Save/copy this file to the same directory as your vesta executable
and then run this script
"""

head,tail = os.path.split(os.path.realpath(__file__))

config.VESTAExec =os.path.join(head,"VESTA.exe")
config.store()