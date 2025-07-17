# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from Engineering.texture.correction.correction_model import TextureCorrectionModel

# Create an example Reference Workspace

exp_name = "Example"
root_dir = fr"C:\Users\Name\Engineering_Mantid\User\{exp_name}"
instr = "ENGINX"


model = TextureCorrectionModel()
LoadEmptyInstrument(InstrumentName=instr, OutputWorkspace="")

model.create_reference_ws(exp_name)

# either set or load sample shape
#set:
shape_xml = ""
SetSampleShape(model.reference_ws, shape_xml)

#load:
shape_file = ""
LoadSampleShape(model.reference_ws, shape_file)

# set material
SetSampleMaterial(model.reference_ws, "Fe")

# save reference file
model.save_reference_file(exp_name, None, root_dir) # just set group as None here