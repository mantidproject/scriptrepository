# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

# script to covert a file with flattened matrices that have been generated in sscanss (and
# thus us in the sscanss reference frame where beam = X, detector = Y, roof = Z) into a
# matrix that is in the mantid reference frame

# Just set the txt file path and the tell it the number of scan points there were and you
# will get a _mantid_point_n.txt file created for each point


#~~~~~~~~~~~~~~~~~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

txt_file = r"C:\Users\kcd17618\Documents\dev\TextureCommisioning\Day3\ZrRing\Sscanss\Split\Zirc_ring_pose_matrices.txt"
NUM_POINTS = 3


#~~~~~~~~~~~~~~~~~~ Script Execution ~~~~~~~~~~~~~~~~~~~~~~~~

with open(txt_file, "r") as f:
    goniometer_strings = [line.replace("\t", ",") for line in f.readlines()]

transformed_strings = []


def convert_from_sscanss_frame(r_zxy):
    # Define M: matrix to convert vectors from XYZ to ZXY
    M = np.array([
        [0, 0, 1],  # X in ZXY = Z in XYZ
        [1, 0, 0],  # Y in ZXY = X in XYZ
        [0, 1, 0]   # Z in ZXY = Y in XYZ
    ])
    M_inv = M.T  # since M is orthonormal

    # Apply the similarity transform in reverse express R in XYZ frame
    return M.T @ r_zxy @ M


for gs in goniometer_strings:
    or_vals = gs.split(",")
    trans_vals = or_vals[9:]
    run_mat = np.asarray(or_vals[:9], dtype=float).reshape((3, 3)).T
    
    mantid_mat = convert_from_sscanss_frame(run_mat)
    new_string = ",".join([str(x) for x in mantid_mat.reshape(-1)]+trans_vals)
    transformed_strings.append(new_string)
    
num_scans = len(goniometer_strings)//NUM_POINTS

for scan_ind in range(NUM_POINTS):
    save_file = txt_file.replace(".txt", f"_mantid_point_{scan_ind}.txt")

    with open(save_file, "w") as f:
        f.writelines(transformed_strings[scan_ind*num_scans:(scan_ind+1)*num_scans])