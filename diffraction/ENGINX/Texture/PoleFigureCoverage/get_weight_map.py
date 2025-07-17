# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from coverage_funcs import *
from os import path
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path
from importlib import reload
import sys
reload(sys.modules["coverage_funcs"])


# This script might help you highlight the regions of interest in a pole figure to
# produce a weighting map that can be used within the coverage generator script

# questions: andy.bridger@stfc.ac.uk

#------------------------------------------------------
# ~~~~~~~~~~  Parameters of interest  ~~~~~~~~~~~~~~~~~
#------------------------------------------------------

# define your intrinsic direction names and colours (the second will be the projection direction)
dir_names = ["RD", "ND", "TD"]
dir_cols = ["red", "green", "blue"]

# define the resolution of the alpha-beta search space
grid_res = 256

#define the projection of the pole figure
projection = "azim" # ster or azim

# directory to save the files in once the figure has been closed
save_dir = r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data"

#------------------------------------------------------
# ~~~~~~~~~~  Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------------------------


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


#------------------------------------------------------
# ~~~~~~~~~~  Execution Code  ~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------------------------

# create the alpha-beta meshgrid
X,Y = np.meshgrid(np.linspace(0, 2*np.pi, grid_res, endpoint = True), 
                    np.linspace(0, np.pi/2, grid_res//4, endpoint = True))
ang_grid = np.stack((X, Y), axis=-1).reshape((grid_res//4, grid_res, 2))

# flatten mesh grid
grid_flat = ang_grid.reshape(-1, 2)
grid_phi, grid_theta = grid_flat[:, 0], grid_flat[:, 1]

# get the alpha-beta grid in pole figure plotting cartesian coordinates
grid_vecs = spherical_to_cartesian(grid_phi, grid_theta)  # shape (X*Y, 3)

# project the coordinates     
ster_xy = ster_proj(*grid_flat.T) if projection == "ster" else azim_proj(*grid_flat.T)

# mask array for saving the deired weights
paint_mask = np.zeros(len(ster_xy), dtype=int)

# define the radius for drawing
brush_radius = 0.1
brush_r2 = brush_radius ** 2

# create figure
fig, ax = plt.subplots()
pts = ax.scatter(ster_xy[:, 1], ster_xy[:, 0], s=2, c='blue') # we want Dir1 along ax-Y and Dir3 along ax-X
ax.set_aspect('equal')
ax.set_title("Left-click to paint (1), Ctrl+click to paint (2)")
[ax.quiver(*np.array((-1,-1)), *bv, color = dir_cols[-1+i], scale = 5) for i, bv in enumerate(np.eye(2))]
ax.annotate(dir_names[0], (-0.95, -0.8))
ax.annotate(dir_names[2], (-0.8, -0.95))
ax.set_axis_off()

# Function to update colors based on mask
def update_colors():
    colors = np.where(paint_mask == 0, 'blue',
             np.where(paint_mask == 1, 'yellow', 'red'))
    pts.set_color(colors)

def on_move(event):
    if not event.button or not event.inaxes:
        return

    x, y = event.xdata, event.ydata  # we have flipped x and y in the plot so Dir3 (ster_xy[:,1]) corresponds to plot x, etc
    d2 = (ster_xy[:, 1] - x) ** 2 + (ster_xy[:, 0] - y) ** 2

    # Determine paint value
    if event.button == 3 or event.key == 'control':
        paint_val = 2
    elif event.button == 1:
        paint_val = 1
    else:
        return

    # Apply paint
    paint_mask[d2 < brush_r2] = paint_val
    update_colors()
    fig.canvas.draw_idle()

def save_weight_map(dir):
    np.save(path.join(dir, "ang_grid.npy"), ang_grid)
    np.save(path.join(dir, "grid_vecs.npy"), grid_vecs)
    np.save(path.join(dir, "ster_xy.npy"), ster_xy)
    np.save(path.join(dir, "weight_mask.npy"), paint_mask.reshape(ang_grid.shape[:2]))        

def on_close(event):
    save_weight_map(save_dir)
    logger.notice("Figure closed, files saved")


fig.canvas.mpl_connect('motion_notify_event', on_move)
fig.canvas.mpl_connect('close_event', on_close)
update_colors()
plt.show()
