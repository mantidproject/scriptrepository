import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
from mantid.simpleapi import *
from mantid.api import AnalysisDataService as ADS
from Engineering.common.calibration_info import CalibrationInfo
from Engineering.EnggUtils import GROUP
from Engineering.common.texture_sample_viewer import plot_sample_directions
from mantidqt.plotting.sample_shape import plot_sample_only
import os

# ========== CONFIGURATION ==========

instr = "ENGINX"
group = GROUP("Texture30")
groupfile_root = r"C:\Users\kcd17618\Documents\dev\mantid\mantid\scripts\Engineering\calib"
shape = None
euler_scheme = "YZY"
angle_res = 5
angles = [0, 0, 0]  # Make it mutable

ax_transform = np.array([(1,0,0), 
                         (0,0,1), 
                         (0,1,0)])
dir_names = ["RD", "ND", "TD"]
dir_cols = ["red", "green", "blue"]
projection = "azim"

save_scheme = "euler" #euler, matrix, sscanss
save_path = r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data\i_euler_angles.txt"

# ========== INITIAL SETUP ==========

def get_default_shape():
    return """
     <cuboid id='some-cuboid'> \
    <height val='0.01'  /> \
    <width val='0.01' />  \
    <depth  val='0.01' />  \
    <centre x='0.0' y='0.0' z='0.0'  />  \
    </cuboid>  \
    <algebra val='some-cuboid' /> \
    """
    
def ring(r = 1, res = 100, offset=(0,0,0), axis = "z"):    
    u = np.linspace(0, 2 * np.pi, res)
    a = r*np.cos(u) + offset[0]
    b = r*np.sin(u) + offset[1]
    c = np.zeros_like(a) + offset[2]
    match axis:
        case ax if ax.lower() == "z":
            return np.concatenate((a[None,:],b[None,:],c[None,:]), axis = 0)
        case ax if ax.lower() == "y":
            return np.concatenate((b[None,:],c[None,:],a[None,:]), axis = 0)
        case ax if ax.lower() == "x":
            return np.concatenate((c[None,:],a[None,:],b[None,:]), axis = 0)
            
            
def write_to_euler_angle_txt(arr, fp):
    with open(fp, "w") as f:
        f.writelines(["\t".join([str(np.round(a,2)) for a in r])+"\n" for r in arr])
        
def write_matrix_txt_file(arr, fp, rot_scheme):
    with open(fp, "w") as f:
        f.writelines(["\t".join([str(a) for a in Rotation.from_euler(rot_scheme, r, degrees=True).as_matrix().reshape(-1)])+"\n" for r in arr])
        
def get_alpha_beta_from_cart(q_sample_cart: np.ndarray) -> np.ndarray:
    """
    get spherical angles from cartesian coordinates
    alpha is angle from positive x towards positive z
    beta is angle from positive y
    """
    q_sample_cart = np.clip(q_sample_cart.copy(), -1.0, 1.0)  # numerical inaccuracies outside this range will give nan in the trig funcs
    q_sample_cart = np.where(q_sample_cart[1] < 0, -q_sample_cart, q_sample_cart)  # invert the southern points
    alphas = np.arctan2(q_sample_cart[2], q_sample_cart[0])
    betas = np.arccos(q_sample_cart[1])
    return np.concatenate([alphas[:, None], betas[:, None]], axis=1)
    
def spherical_to_cartesian(phi, theta):
    """
    Convert spherical angles to 3D Cartesian unit vector.
    phi: azimuthal angle [0, 2pi]
    theta: polar angle [0, pi/2]
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.cos(theta)
    z = np.sin(theta) * np.sin(phi)
    return np.stack((x, y, z), axis=-1)
    

def ster_proj(alphas: np.ndarray, betas: np.ndarray) -> np.ndarray:
    betas = np.pi - betas  # this formula projects onto the north pole, and beta is taken from the south
    r = np.sin(betas) / (1 - np.cos(betas))
    out = np.zeros((len(alphas), 2))
    out[:, 0] = r * np.cos(alphas)
    out[:, 1] = r * np.sin(alphas)
    return out
    
def azim_proj(alphas: np.ndarray, betas: np.ndarray) -> np.ndarray:
    betas = betas / (np.pi / 2)
    xs = (betas * np.cos(alphas))[:, None]
    zs = (betas * np.sin(alphas))[:, None]
    out = np.concatenate([xs, zs], axis=1)
    return out

        

ws = LoadEmptyInstrument(InstrumentName=instr)
calib_info = CalibrationInfo(group=group)
group_ws = GroupDetectors(InputWorkspace="ws", MapFile=os.path.join(groupfile_root, calib_info.get_group_file()),
                          OutputWorkspace="group_ws")

spec_info = group_ws.spectrumInfo()
det_pos = np.asarray([
    spec_info.position(i) / np.linalg.norm(spec_info.position(i))
    for i in range(1, group_ws.getNumberHistograms())
])
detQs_lab = det_pos - np.array((0,0,1))
detQs_lab = detQs_lab / np.linalg.norm(detQs_lab, axis=1)[:, None]

shape = get_default_shape() if not shape else shape
SetSampleShape(ws, shape)
shape_mesh = ws.sample().getShape().getMesh()
extent = (np.linalg.norm(shape_mesh, axis=(1, 2)).max() / 2) * 1.2

current_rot_index = [0]

# ========== PLOTTING FUNCTION ==========

fig = plt.figure()
lab_ax = fig.add_subplot(1, 2, 1, projection='3d')
lab_ax.view_init(elev=-140, azim=-55, roll = 135)
proj_ax = fig.add_subplot(1, 2, 2)
gon_ax = fig.add_subplot(3, 3, 1, projection='3d')
gon_ax.view_init(elev=-140, azim=-55, roll = 135)

nGon = len(euler_scheme)
gon_colors = ["hotpink", "orange", "purple"]
ge = ((nGon/2)+1)*0.866
gon_res = 100
gon_steps = np.linspace(0, gon_res, 360, endpoint=False).astype(int)
gon_ax_dict = {"x": np.array([1,0,0]), 
                "y": np.array([0,1,0]), 
                "z": np.array([0,0,1])}
                
                
                
temp_arr = []
angle_hashes = []
temp_arr_points = []
global temp_ind
temp_ind = [0]

def update_plot():
    lab_ax.clear()
    proj_ax.clear()
    gon_ax.clear()

    R = Rotation.from_euler(euler_scheme, angles, degrees=True)
    
    gRs = [Rotation.identity()]
    for i in range(nGon):
        r_step = Rotation.from_euler(euler_scheme[i], angles[i], degrees=True)
        gRs.append(gRs[-1] * r_step)
        
    gVecs = []
        
        
    ws.run().getGoniometer().setR(R.as_matrix())

    rot_mesh = R.apply(shape_mesh.reshape((-1, 3))).reshape(shape_mesh.shape)
    rot_pos = R.inv().apply(detQs_lab) @ ax_transform
    
    cart_pos = get_alpha_beta_from_cart(rot_pos.T)
    pf_xy = ster_proj(*cart_pos.T) if projection == "ster" else azim_proj(*cart_pos.T)
    
    for i, a in enumerate(euler_scheme):
        gon_ring = (gRs[i].apply(ring(1 + ((nGon-i)/2), res = gon_res, axis = a).T).T)
        gon_ax.plot(*gon_ring, color = "grey")
        pos_ind = gon_steps[angles[i]-1]
        gon_ax.plot(*gon_ring[:,:pos_ind], color = gon_colors[i])
        gon_ax.plot(*gon_ring[:,pos_ind:], color = "grey")
        
        vec = gRs[i].apply(gon_ax_dict[a.lower()])
        gVecs.append(vec)
        gon_ax.quiver(*np.zeros(3), *vec*2, color = gon_colors[i],ls = ("-", "--")[int(i != current_rot_index[0])])
        
    gPole = R.inv().apply(np.array(gVecs)) @ ax_transform
    cart_gPole = get_alpha_beta_from_cart(gPole.T)
    gPole_xy = ster_proj(*cart_gPole.T) if projection == "ster" else azim_proj(*cart_gPole.T)
    
    gon_ax.set_xlim([-ge, ge])
    gon_ax.set_ylim([-ge, ge])
    gon_ax.set_zlim([-ge, ge])
    gon_ax.set_axis_off()

    # 3D plot
    fig.sca(lab_ax)
    plot_sample_only(fig, rot_mesh*0.5, 0.5, "grey")
    plot_sample_directions(fig, "ws", ax_transform, dir_names)
    lab_ax.set_xlim([-extent, extent])
    lab_ax.set_ylim([-extent, extent])
    lab_ax.set_zlim([-extent, extent])
    [lab_ax.quiver(*np.zeros(3), *dQ*1.25*extent, arrow_length_ratio=0.05, color="grey", alpha=0.25) for dQ in detQs_lab]
    [lab_ax.scatter(*dQ*1.25*extent, color="dodgerblue", s=2,) for i, dQ in enumerate(detQs_lab)]
    lab_ax.set_axis_off()
    
    
    

    # 2D plot
    labels = ["1/2 : -/+", "4/5 : -/+", "7/8 : -/+"]
    for i, gP in enumerate(gPole_xy):
        pc = gon_colors[i]
        fc = "None" if i != current_rot_index[0] else pc
        if np.isclose(np.linalg.norm(gP),1):
            proj_ax.plot((gP[1], -gP[1]), (gP[0],-gP[0]),color = pc, 
            ls = ("-", "--")[int(i != current_rot_index[0])], label = labels[i])
        else:
            proj_ax.scatter(gP[1], gP[0], s=30, edgecolor=pc, facecolor = fc, label = labels[i])
    proj_ax.scatter(pf_xy[:, 1], pf_xy[:, 0], s=20, c='dodgerblue')
    # Get current handles and labels
    handles, labels = proj_ax.get_legend_handles_labels()
    
    for i, tpf in enumerate(temp_arr_points):
        fc = "None" if i != temp_ind[0] else "dodgerblue"
        proj_ax.scatter(tpf[:, 1], tpf[:, 0], s=20, facecolor=fc, edgecolor = 'dodgerblue', alpha = 0.2)
        
    
    # Reverse both
    proj_ax.legend(handles[::-1], labels[::-1], loc=1)

    proj_ax.set_aspect('equal')
    proj_ax.set_xlim(-1.1, 1.1)
    proj_ax.set_ylim(-1.1, 1.1)
    proj_ax.set_title("Experimental Coverage")
    [proj_ax.quiver(*np.array((-1, -1)), *bv, color=dir_cols[-1+i], scale=5) for i, bv in enumerate(np.eye(2))]
    circle = plt.Circle((0, 0), 1, color="grey", fill=False, linestyle="-")
    proj_ax.add_patch(circle)
    proj_ax.annotate(dir_names[0], (-0.95, -0.8))
    proj_ax.annotate(dir_names[2], (-0.8, -0.95))
    proj_ax.set_axis_off()

    fig.canvas.draw_idle()
    return pf_xy


# ========== INTERACTIVITY ==========

def on_key(event):
    delta = angle_res
    
    if event.key == '1':
        angles[0] = (angles[0] - delta)%360
        current_rot_index[0] = 0
    elif event.key == '2':
        angles[0] = (angles[0] + delta)%360
        current_rot_index[0] = 0
    elif event.key == '4':
        angles[1] = (angles[1] - delta)%360
        current_rot_index[0] = 1
    elif event.key == '5':
        angles[1] = (angles[1] + delta)%360
        current_rot_index[0] = 1
    elif event.key == '7':
        angles[2] = (angles[2] - delta)%360
        current_rot_index[0] = 2
    elif event.key == '8':
        angles[2] = (angles[2] + delta)%360
        current_rot_index[0] = 2
    if event.key == 'enter':
        ang_hash = hash(tuple(angles))
        if ang_hash not in angle_hashes:
            angs = tuple(angles)
            temp_arr.append(angs)
            angle_hashes.append(ang_hash)
            pf_xy = update_plot()
            temp_arr_points.append(pf_xy)
            temp_ind[0] = temp_ind[0] + 1 if len(temp_arr) != 1 else 0
    if event.key == 'backspace':
        if len(temp_arr) > 0:
            ang_hash = hash(tuple(angles))
            temp_arr.pop(temp_ind[0])
            angle_hashes.pop(temp_ind[0])
            temp_arr_points.pop(temp_ind[0])
            temp_ind[0] = temp_ind[0] - 1
    if event.key == 'right':
        if temp_ind[0] < len(temp_arr)-1:
            temp_ind[0]  = temp_ind[0] + 1
    if event.key == 'left':
        if temp_ind[0] > 0:
            temp_ind[0]  = temp_ind[0] - 1
    print(f"Updated angles: {angles}")
    update_plot()

def on_close(event):
    if save_scheme == "euler":
        write_to_euler_angle_txt(temp_arr, save_path)
        logger.notice("Figure closed, orientations saved as euler angles")
    if save_scheme == "matrix":
        write_to_euler_angle_txt(temp_arr, save_path, euler_scheme)
        logger.notice("Figure closed, orientations saved as euler angles")

fig.canvas.mpl_connect('key_press_event', on_key)
fig.canvas.mpl_connect('close_event', on_close)

# Initial plot
update_plot()
plt.show()



