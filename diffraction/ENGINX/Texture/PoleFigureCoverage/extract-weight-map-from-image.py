# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import minimize
from importlib import reload
import sys
import os
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull
from scipy.interpolate import griddata

# This script might help you load a pole figure that you have acquired by other means
# and produce a weighting map to use within the coverage generator script

# To get it to work might require some user effort, depending how different your input pf 
# is from those used to develop the script. 

# Broadly, we are looking for images of spherical pole figures with white backgrounds outside
# the boarder. If there are white regions within the figure, use mask_type = "hull".

# If you know what colour map has been used you can pass that as best_guess_cmap, otherwise 
# if a colorbar is available, you can provide an image of just the cropped colorbar by setting
# best_guess_cmap = "custom" and providing the image path

# questions: andy.bridger@stfc.ac.uk

# Define image path
im_path = r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data\TextureImages\example2.png"

# define type of masking either "hull" or "thresh"(see above)
mask_type = "hull"
# thresh_val defines how white is considered background where white-white is (255,255,255)
thresh_val = 230

# this will scale the image up, after the circle has been fit - this is to account for any
# inscribed border around the figure
border_scale = 1.05

# define cmap (matplotlib name or "custom")
best_guess_cmap = "jet"
custom_cbar_path = None #r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data\TextureImages\mtex_cbar.png"
cbar_ax = 1 # 1 for vertical, 0 for horizontal

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

def points_in_hull(points, hull, tolerance=1e-12):
    A = hull.equations[:, :-1]  # shape (F, D)
    b = hull.equations[:, -1]   # shape (F,)
    
    # Compute dot products and apply inequalities
    inside = np.all(np.dot(points, A.T) + b <= tolerance, axis=1)
    return inside

def get_convex_hull_mask(im, tmask):
    tc = np.asarray(np.where(tmask ==1)).T
    hull = ConvexHull(tc)
    all_coords = np.asarray(np.where(np.any(im, axis = -1))).T
    return points_in_hull(all_coords, hull).reshape(im.shape[:2])
    
    
def get_threshold_mask(im, val = 230):
    return 1- np.all(im[:,:,:3]>val, axis = -1).astype(int)
    
    
def get_rgb_tree_from_cmap(cmap_name, num_samples):
    cmap = plt.get_cmap(cmap_name, num_samples)
    scalar_samples = np.linspace(0, 1, num_samples)
    samples = cmap(scalar_samples)
    return cKDTree(samples), scalar_samples
    
def get_rgb_tree_from_custom_cbar(custom_cbar, cbar_ax):
    print(custom_cbar.shape)
    cbar = custom_cbar.mean(axis = cbar_ax)
    cbar_vals = np.clip(cbar / 255.0, 0, 1)
    scalar_samples = np.linspace(1, 0, len(cbar_vals))
    return cKDTree(cbar_vals), scalar_samples
    

def rgb_to_scalar(rgb_array, cmap_name='viridis', num_samples=1000, custom_cbar = None, cbar_ax = 0):
    # Normalize input RGB values to [0, 1]
    rgb_array = np.clip(rgb_array / 255.0, 0, 1)
    
    if cmap_name == "custom":
        tree, scalar_samples = get_rgb_tree_from_custom_cbar(custom_cbar, cbar_ax)
    else:
        tree, scalar_samples = get_rgb_tree_from_cmap(cmap_name, num_samples)

    # Find the closest match for each input RGB
    _, indices = tree.query(rgb_array)

    # Map indices back to scalar values
    scalar_values = scalar_samples[indices]

    return scalar_values
    
    

def fit_pf_image(im, 
                 mask_type = "hull", 
                 thersh_val = 230, 
                 border_scale = 1,
                 plot_circle_fit=False, 
                 best_guess_cmap = "jet",
                 plot_output =False, 
                 custom_cbar = None, 
                 cbar_ax = 0):
                     
    # get a mask of the 'non-white-background' around the pole figure
    tmask = get_threshold_mask(im, thresh_val)
    # if mask type is hull it will make a convex hull mask - this is for when there are white
    # regions within the pole figure itself 
    approx_mask = get_convex_hull_mask(im, tmask) if mask_type == "hull" else tmask
    
    # get the pixel coordinates 
    coords = np.asarray(np.where(approx_mask ==1)).T
    
    # we are going to get the distance, r, of each point from the proposed centre.
    # each pixel can be considered as an area element, so if we sort the rs we expect
    # a circle to have a distribution where a given r = sqrt(A/pi)
    
    As = np.linspace(0,np.pi,len(coords))
    t_rs = np.sqrt(As/np.pi)

    def cost(params):
        x, y, a, b = params
        norm_coords = (coords - np.array((x,y)))/np.array((a,b))
        rs = np.sort(np.linalg.norm(norm_coords, axis = 1))
        return np.abs(rs-t_rs).sum()
        
    com = coords.mean(axis = 0)
    scale = np.abs(coords-com).max(axis =0)

    initial_guess = [*com, *scale]

    result = minimize(cost, initial_guess, method='L-BFGS-B', bounds=[
        (0, im.shape[0]),  # x
        (0, im.shape[1]),  # y
        (1, im.shape[0]),  # a
        (1, im.shape[1])   # b
    ])

    
    if plot_circle_fit:
        norm_coords = (coords - np.array((x_opt,y_opt)))/np.array((a_opt,b_opt))

        rs = np.linalg.norm(norm_coords, axis = 1)

        plt.figure()
        plt.plot(ps, np.sort(rs))
        plt.plot(ps, np.sqrt(ps/np.pi))
        plt.show()
    
    x_opt, y_opt, a_opt, b_opt = result.x
    all_coords = np.asarray(np.where(np.any(im, axis = -1))).T

    norm_all_coords = border_scale*(all_coords - np.array((x_opt,y_opt)))/np.array((a_opt,b_opt))
    norm_all_coords[:,0] = -norm_all_coords[:,0] # image has origin at top left
    
    val_i = np.zeros((im.shape[0],im.shape[1]))
    vals = im[tmask==1,:]
    
    val_i[tmask==1] = rgb_to_scalar(vals, best_guess_cmap, 1000, custom_cbar, cbar_ax)
    val_i = val_i.reshape(-1)



    if plot_output:
        plt.figure()
        plt.imshow(val_i.reshape([vals.shape[0], vals.shape[1]]))
        plt.show()
        
    return norm_all_coords, val_i
    
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

#load image
image = Image.open(im_path)
im = np.asarray(image)

#load colorbar
custom_cbar = np.asarray(Image.open(custom_cbar_path)) if custom_cbar_path else None

# resize the provided pole figure so it is between 1 and -1 and centered on the origin
norm_all_coords, val_i = fit_pf_image(im, mask_type, thresh_val, border_scale,
                                      best_guess_cmap= best_guess_cmap,
                                      custom_cbar = custom_cbar, 
                                      cbar_ax = cbar_ax
                                      )

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

# interpolate the value of these points from the provided image
interp_mask = griddata(norm_all_coords, val_i, ster_xy, method='linear')

fig, ax = plt.subplots()
pts = ax.scatter(ster_xy[:, 1], ster_xy[:, 0], s=2, c=interp_mask) # we want Dir1 along ax-Y and Dir3 along ax-X
ax.set_aspect('equal')
[ax.quiver(*np.array((-1,-1)), *bv, color = dir_cols[-1+i], scale = 5) for i, bv in enumerate(np.eye(2))]
ax.annotate(dir_names[0], (-0.95, -0.8))
ax.annotate(dir_names[2], (-0.8, -0.95))
ax.set_axis_off()
fig.show()

def save_weight_map(dir):
    np.save(os.path.join(dir, "ang_grid.npy"), ang_grid)
    np.save(os.path.join(dir, "grid_vecs.npy"), grid_vecs)
    np.save(os.path.join(dir, "ster_xy.npy"), ster_xy)
    np.save(os.path.join(dir, "weight_mask.npy"), interp_mask.reshape(ang_grid.shape[:2]))       
   
save_weight_map(save_dir) 

