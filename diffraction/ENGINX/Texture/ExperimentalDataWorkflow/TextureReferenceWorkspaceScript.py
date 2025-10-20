# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from Engineering.texture.correction.correction_model import TextureCorrectionModel

# Create an example Reference Workspace

# set experiment name
exp_name = "Example"

# set the directory where your workflow files should be saved
save_root = r"C:\Users\fedid12345\Engineering_Mantid"
root_dir = fr"{save_root}\User\{exp_name}"
instr = "ENGINX"

# set shape info to either be a shape xml string or a file to an stl
example_shape_info = """
<hollow-cylinder id="A">
  <centre-of-bottom-base x="-0.01315" y="-0.01315" z="-0.00756" /> 
  <axis x="0.0" y="0.0" z="1.0" />
  <inner-radius val="0.0145" />
  <outer-radius val="0.0223" />
  <height val="0.01512" />
</hollow-cylinder>
"""

sample_material = "Zr"


model = TextureCorrectionModel()
model.create_reference_ws(exp_name, instr)

# if it ends with .stl assume we have been given the file path
model.set_sample_info(model.reference_ws, example_shape_info, sample_material)

# save reference file
model.save_reference_file(exp_name, None, root_dir) # just set group as None here