# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import os
from Engineering.common.calibration_info import CalibrationInfo
from Engineering.EnggUtils import GROUP
from scipy.spatial.transform import Rotation
from coverage_funcs import *
from Engineering.common.texture_sample_viewer import plot_sample_directions
from importlib import reload
import sys
reload(sys.modules["coverage_funcs"])

# Script that might help you identify a set of sample orientations to use to get desired
# pole figure coverage. You essentially generate a map of regions of interest within the pole figure
# either by drawing the regions or by loading an existing pole figure (other scripts).

# You then state the starting configuration of your sample intrinsic axes, intended detector grouping,
# and positioner constraints (i.e. goniometer axis and accessible angles). The script then uses
# a genetic algorithm to search for a set of orientations with the specified number of runs
# which covers the desired regions of the pole figure. 

# questions: andy.bridger@stfc.ac.uk

#------------------------------------------------------
# ~~~~~~~~~~  Parameters of interest  ~~~~~~~~~~~~~~~~~
#------------------------------------------------------

instr = "ENGINX"  # currently just ENGINX here
group = GROUP("Texture30")   # set this to view the coverage with the correct calibration grouping

euler_scheme = "YZY"   # scheme for both the goniometer if output file is euler and for generating orientation search space
                       # search space for 

point_confidence = 0.5  # this controls the size of the region around each experimental pole figure point
                      # this is a factor applied to the steradian solid angle of the individual grouping spectra

ax_transform = np.array([(1,0,0),   # defines the intrinsic directions in the 0,0,0 goniometer state
                         (0,0,1),
                         (0,1,0)])
dir_names = ["RD", "ND", "TD"]  # names of these three directions
dir_cols = ["red", "green", "blue"]   # colours for these three directions

mask_weight_boost = (0,2,2)       # controls the weighting values, taking a,b,c:
                                   # w_x' = ((bw_x)^c)+a where w_x = 0,1,2, from the weighting map provided by the user

N_ANGLES = 8        # (n x n x n) search space for angles passed into euler scheme
angle_range_1 = np.linspace(-np.pi/2, np.pi/2, N_ANGLES)
angle_range_2 = np.linspace(-np.pi/2, np.pi/2, N_ANGLES)
angle_range_3 = np.linspace(-np.pi, np.pi, N_ANGLES)

# genetic algorithm parameters
num_runs = 8 # number of experimental runs you want to optimise pole figure coverage for
pop_size = 256  # population of sets of experimental runs
n_generations = 128  # number of iterations to evolve the population over
pool_size = pop_size   # the number of parents to pair up for the next generation
mutation_rate_lookup = np.linspace(0.5, 0.1, n_generations)**2  # probability of a given orientation being randomly replaced

projection = "azim" # ster or azim
output_file_type = "matrix" # euler, matrix, sscanss

#------------------------------------------------------
# ~~~~~~~~~~  Configuration Data  ~~~~~~~~~~~~~~~~~~~~~
#------------------------------------------------------

groupfile_root = r"C:\Users\kcd17618\Documents\dev\mantid\mantid\scripts\Engineering\calib"
save_path = r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data"
map_files = r"C:\Users\kcd17618\Documents\MantidScripts\PoleFigureCoverageEstimate\Data"

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

def plot_fitness_evolution(fitnesses, n_generations):
    plt.figure()
    ci = np.linspace(0.2, 1, n_generations)
    for i, f in enumerate(fitnesses):
        plt.plot(f, c = (ci[i], 0, 0))
    plt.show()
    
def plot_gene_evolution(pops, full_masks, n_generations):
    gene_abundance = np.zeros((full_masks.shape[0], n_generations))

    plt.figure()
    for igen, pop in enumerate(pops):
        ugenes = np.unique(pop)
        abundance = []
        for ug in ugenes:
            gene_abundance[ug, igen] = np.sum(pop == ug)
    plt.imshow(gene_abundance)
    plt.show()
    
def write_to_euler_angle_txt(arr, fp):
    with open(fp, "w") as f:
        f.writelines(["\t".join([str(np.round(a,2)) for a in np.rad2deg(r)])+"\n" for r in arr])
        
def write_matrix_txt_file(arr, fp, rot_scheme):
    with open(fp, "w") as f:
        f.writelines(["\t".join([str(a) for a in Rotation.from_euler(rot_scheme, r).as_matrix().reshape(-1)])+"\n" for r in arr])
        
        

#------------------------------------------------------
# ~~~~~~~~~~  Actual Code  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------------------------

# load the instrument geometry
ws = LoadEmptyInstrument(InstrumentName = instr)

#setup the calibration info
calib_info = CalibrationInfo(group = group)

#get the calibration detector grouping
group_ws = GroupDetectors(InputWorkspace = "ws", MapFile = os.path.join(groupfile_root, calib_info.get_group_file()),
               OutputWorkspace = "group_ws")

# get the solid angle coverage of each of these spectra 
solid_ang = SolidAngle(InputWorkspace = "group_ws", OutputWorkspace = "solid_ang")

# scale based on point confidence
patch_r = solid_ang.extractY()[1:]*point_confidence
               
# calculate the detector Q vectors in the lab frame
spec_info = group_ws.spectrumInfo()

det_pos = np.asarray([spec_info.position(spec_ind)/np.linalg.norm(spec_info.position(spec_ind)) 
                    for spec_ind in range(1, group_ws.getNumberHistograms())])

detQs_lab = det_pos - np.array((0,0,1))
detQs_lab = detQs_lab/np.linalg.norm(detQs_lab, axis = 1)[:,None]   

# get the detQs in sample axes frame
#detQs = np.asarray([dQ @ ax_transform for dQ in detQs_lab])

# Load Weighting Map etc
ang_grid = np.load(os.path.join(map_files, "ang_grid.npy"))
grid_vecs = np.load(os.path.join(map_files, "grid_vecs.npy"))
ster_xy = np.load(os.path.join(map_files, "ster_xy.npy"))
weight_mask = np.load(os.path.join(map_files, "weight_mask.npy"))[None, :,:]
weight_mask = ((weight_mask*mask_weight_boost[1])**mask_weight_boost[2])+mask_weight_boost[0]
grid_res = ang_grid.shape[1]


# Define orientation search space
X, Y, Z = np.meshgrid(angle_range_1, angle_range_2, angle_range_3, indexing='ij')
euler_angles = np.stack((X, Y, Z), axis=-1).reshape(-1, 3)

# get the detQs for each orientation 
rot_pos = []
for angles in euler_angles:
    R = Rotation.from_euler(euler_scheme, angles)
    rot_pos.append(R.inv().apply(detQs_lab) @ ax_transform) # rotate in lab frame then transform ax to sample
rot_pos = np.moveaxis(np.asarray(rot_pos), 0, 1)



# Get detector masks for each set of angles
cos_angles = rot_pos @ grid_vecs.T  # rot_pos is (group, n_orientations, 3) (target_pf_vals, 3)


theta_r = np.arccos(1 - patch_r / (2 * np.pi))  # convert solid angle to angular radius
cos_patch_r = np.cos(theta_r)[:, :, None]       # shape (ndet, 1, 1)

full_masks = np.abs(cos_angles) >= cos_patch_r
full_masks = np.any(full_masks.reshape(full_masks.shape[0], full_masks.shape[1], grid_res//4, grid_res), axis=0)

# Find a solution using genetic algorithm

# define the GA functions

def calc_fitness(mask_set, pop, weights):
    return np.tensordot(combine_masks(mask_set, pop), weights[0], axes=([1, 2], [0, 1]))
    
def combine_masks(mask_set, pop):    
    return (mask_set[pop].sum(axis = 1) > 0).astype(int)
    
def breed(pool, mutation_rate=0.05):
    parents1 = pool[::2]
    parents2 = pool[1::2]

    mask = np.random.rand(*parents1.shape) > 0.5
    offspring = np.where(mask, parents1, parents2)

    mutation_mask = np.random.rand(*offspring.shape) < mutation_rate
    offspring[mutation_mask] = np.random.randint(0, full_masks.shape[0], size=mutation_mask.sum())

    return offspring, calc_fitness(full_masks, offspring, weight_mask)

# generate an initial population
population = np.random.choice(full_masks.shape[0], (pop_size, num_runs), replace = True) # allow repeated genes for speed
# calc initial fitness
pop_fitness = calc_fitness(full_masks, population, weight_mask)

fitnesses = []
pops = []

best_weight = pop_fitness.max()
best = population[pop_fitness.argmax()]
for gen in range(n_generations):
    print(gen)

    weights = pop_fitness/np.sum(pop_fitness)

    ranks = np.argsort(pop_fitness)
    
    # evaluate here to save reranking 
    if pop_fitness[ranks[-1]] > best_weight:
        best = population[ranks[-1]]
    
    #track population 
    pops.append(population.copy())
    fitnesses.append(pop_fitness[ranks])
    
    # get breeding pool
    pool = population[np.random.choice(population.shape[0], pool_size, p = weights[ranks], replace = True)]
    
    # update the population and the fitnesses with new offspring
    population[ranks[:pool_size//2]], pop_fitness[ranks[:pool_size//2]] = breed(pool, mutation_rate_lookup[gen])

# check the final generation
if pop_fitness[ranks[-1]] > best_weight:
    best = population[ranks[-1]]


# uncomment these to view how the solution evolves

#plot_fitness_evolution(fitnesses, n_generations)

#plot_gene_evolution(pops, full_masks, n_generations)

# extract the best solution detQs
cover_vecs = np.concatenate(rot_pos[:,best,:], axis = 0)

# get the binary pole figure coverage of these detQs 
best_coverage_m = combine_masks(full_masks, best[None,:])
# and the weighted coverage
best_coverage_w = best_coverage_m*(weight_mask[0])
# covert these into cartesian for plotting a pole figure
best_ab = get_alpha_beta_from_cart(cover_vecs.T)
# project the pole figure
pf_xy = ster_proj(*best_ab.T) if projection == "ster" else azim_proj(*best_ab.T)

# plot the experimental pole figure sample points
fig, ax = plt.subplots()
pts = ax.scatter(pf_xy[:, 1], pf_xy[:, 0], s=2, c='blue')
ax.set_aspect('equal')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_title("Experimental Coverage")
[ax.quiver(*np.array((-1,-1)), *bv, color = dir_cols[-1+i], scale = 5) for i, bv in enumerate(np.eye(2))]
ax.annotate(dir_names[0], (-0.95, -0.8))
ax.annotate(dir_names[2], (-0.8, -0.95))
ax.set_axis_off()
fig.show()

# see the difference in the coverage
fig, axs = plt.subplots(1,3)
axs[0].scatter(ster_xy[:, 1], ster_xy[:, 0], s=2, c=weight_mask[0])
axs[0].set_title("Weighted Target")
axs[1].scatter(ster_xy[:, 1], ster_xy[:, 0], s=2, c=best_coverage_m)
axs[1].set_title("Effective Coverage")
axs[2].scatter(ster_xy[:, 1], ster_xy[:, 0], s=2, c=best_coverage_w)
axs[2].set_title("Weighted Coverage")
for ax in axs:
    ax.set_aspect('equal')
    ax.set_axis_off()
    [ax.quiver(*np.array((-1,-1)), *bv, color = dir_cols[-1+i], scale = 5) for i, bv in enumerate(np.eye(2))]
    ax.annotate(dir_names[0], (-0.95, -0.8), fontsize = 5)
    ax.annotate(dir_names[2], (-0.8, -0.95), fontsize = 5)

plt.show()

# output the orientation data
if output_file_type == "euler":
    write_to_euler_angle_txt(euler_angles[np.sort(best)], os.path.join(save_path, "euler_angles.txt"))
if output_file_type == "matrix":
    write_matrix_txt_file(euler_angles[np.sort(best)], os.path.join(save_path, "rotation_matrices.txt"), euler_scheme)




