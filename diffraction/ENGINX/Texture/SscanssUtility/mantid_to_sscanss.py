# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from scipy.spatial.transform import Rotation

def convert_to_sscanss_frame(angs, euler_scheme):
    # Define your intrinsic YZY rotation in ZXY system
    r_zxy = Rotation.from_euler(euler_scheme, angs, degrees=True)

    # Define M: matrix to convert vectors from XYZ to ZXY
    M = np.array([
        [0, 0, 1],  # X in ZXY = Z in XYZ
        [1, 0, 0],  # Y in ZXY = X in XYZ
        [0, 1, 0]   # Z in ZXY = Y in XYZ
    ])
    M_inv = M.T  # since M is orthonormal

    # Apply the similarity transform to express R in XYZ frame
    r_xyz = Rotation.from_matrix(M @ r_zxy.as_matrix() @ M.T)

    # Now extract Euler angles in XYZ axes (extrinsic or intrinsic as needed)
    return -r_xyz.as_euler("xyz", degrees=True)

angle_file =  r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\SteelSpatialResolve\Sscanss\Euler_angles.txt"
save_file = angle_file.replace(".txt",".angles").replace(".angles", "_sscanss.angles")
N_SCAN_POINTS = 2
euler_scheme = "YZY"

with open(angle_file, "r") as f:
    lines = f.readlines()
    
header = ["xyz\n",]
angle_lines = lines
print(angle_lines)
angles = [[float(y) for y in x.rstrip('\n').split('\t')] for x in angle_lines]
sscanss_angles = [convert_to_sscanss_frame(angs) for angs in angles]
sscanss_lines = []
for angs in sscanss_angles:
    for _ in range(N_SCAN_POINTS):
        sscanss_lines.append(f"{np.round(angs[0],2)}\t{np.round(angs[1],2)}\t{np.round(angs[2],2)}\n")



with open(save_file, "w") as f:
    f.writelines(header+sscanss_lines)