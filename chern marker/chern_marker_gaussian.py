import numpy as np
from spinv import local_chern_marker, make_finite, onsite_disorder
from spinv.example_models import haldane_pythtb
import time
import sys
import os
import pythtb
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# GLOBALS

# size of the system
NX = 50
NY = 50
# model masses (use linspace bethween 0 and 3 later)
MASS = 2.0
# disorder strength list (bethween 0 and 10)
W = [2.5, 3.0, 3.5]

def onsite_disorder_gauss(source_model, w : float, spinstates : int = 2, seed : int = None):
    """
    Add onsite (Anderson) disorder to the specified model. The disorder amplitude per site is taken randomly in [-w/2, w/2].

        Args:
        - source_model : the model to add disorder to
        - w : disorder amplitude
        - spinstates : spin of the model
        - seed : seed for random number generator

        Returns:
        - model : the disordered model
    """
    # Quick return for no disorder
    if w == 0:
        return source_model

    if seed is not None:
        np.random.seed(seed)

    # Number of orbitals in the supercell model = norbs (original) x num
    norbs = source_model.get_num_orbitals()
    
    # Onsite energies per unit cell (2 is by convention with TBModels)
    # disorder = 0.5 * w * (2 * np.random.rand(norbs // spinstates) - 1.0)
    disorder = np.random.normal(loc=0.0, scale=0.5 * w, size=(norbs // spinstates))
    disorder = np.repeat(disorder, spinstates)
    onsite = source_model._site_energies + disorder

    newmodel = pythtb.tb_model(dim_k = 0, dim_r = source_model._dim_r, lat = source_model._lat, orb = source_model._orb, nspin = source_model._nspin)
    newmodel.set_onsite(onsite, mode = 'set')

    # Cycle over the rows of the hopping matrix
    hoppings = source_model._hoppings
    for k in range(len(hoppings)):
        if np.absolute(hoppings[k][0]) < 1e-10: continue
        newmodel.set_hop(hoppings[k][0], hoppings[k][1], hoppings[k][2], mode = "add")
    
    return newmodel


def lcm_calc(mass, w, Nx, Ny):
    """
    Constructs a Haldane model for given [mass: float] parameter
    Performs heavy local_chern_marker() calcluation for each of the disorder strength 
    parameters in the [w_list: float list]
    Returns a dictionary of {w: lcm (Nx*Ny array)}
    Times the procedure
    """

    # Timing of the task
    start_time = time.time()

    print(f'started calculations for mass: {mass}')

    # Create Haldane models through PythTB and TBmodels packages
    hmodel_pbc_pythtb = haldane_pythtb(delta = mass, t = -1, t2 = -1/3, phi = np.pi / 2, L=1)
    # Cut the model to make a sample of finite size (defines))
    hmodel_obc_pythtb = make_finite(model=hmodel_pbc_pythtb, nx_sites=Nx, ny_sites=Ny)
    # Initiallise chern marker matrix dictionary

    print(f"computing w = {w/2} for mass {mass}")
    # Add Anderson disorder within [-w/2, w/2]. The argument spinstates specifies the spin of the model
    hmodel_pythtb_disorder = onsite_disorder_gauss(hmodel_obc_pythtb, w=w, spinstates=1, seed=181)
    # hmodel_pythtb_disorder = onsite_disorder(model=hmodel_obc_pythtb, w=w, spinstates=1, seed=181)
    print("disorder added")

    # Compute the local Chern markers for TBmodels and PythTB
    chern_matrix, eval, evec = local_chern_marker(model=hmodel_pythtb_disorder, nx_sites=Nx, ny_sites=Ny, return_eigen=True)
    # chern_matrix = local_chern_marker(model=hmodel_pythtb_disorder, nx_sites=Nx, ny_sites=Ny)
    print("chern marker computed")
    
    # Timing of the task
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Total Elapsed time:", elapsed_time, "seconds", f"for mass {mass}")
    
    # return chern_matrices
    # return chern_matrix
    return chern_matrix, eval, evec

def plot_lcm(matrix, w):

    small_size = 16
    medium_size = 20
    bigger_size = 22

    plt.rc('font', size=medium_size +5)          # controls default text sizes
    plt.rc('axes', titlesize=bigger_size-2+5)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size-2+5)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small_size +5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small_size +5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small_size +5)    # legend fontsize
    plt.rc('figure', titlesize=bigger_size +5)  # fontsize of the figure title

    # Define colors
    #colors = [(0, 'green'), (0.33, 'blue'), (0.5, 'white'), (0.67, 'red'), (1, 'yellow')]
    colors = [(0, '#A63919'), (0.33, '#8C4324'), (0.5, 'white'), (0.67, '#292984'), (1, '#010590')]
    # Create colormap
    custom_cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)

    plt.figure()
    plt.gcf().set_facecolor('none')
    plt.title(f"w = {w/2}", color='#51230F')
    #plt.imshow(matrix, cmap='seismic', origin='lower', extent=(0, matrix.shape[1], 0, matrix.shape[0]), vmin=-np.max(np.abs(matrix)), vmax=np.max(np.abs(matrix)))
    plt.imshow(matrix, cmap=custom_cmap, origin='lower', extent=(0, matrix.shape[1], 0, matrix.shape[0]), vmin=-2, vmax=2)
    plt.xlabel('X', color='#51230F')
    plt.ylabel('Y', color='#51230F')
    plt.tick_params(axis='x', colors='#51230F')  # Change x-axis tick color to red
    plt.tick_params(axis='y', colors='#51230F') # Change y-axis tick color to blue

    cbar = plt.colorbar(label='Chern Marker')
    numticks = 4
    cbar.locator = plt.MaxNLocator(numticks)
    cbar.update_ticks()
    plt.show()
    # plt.savefig(path + f'w{(W)}'+ f'M{MASS}' + '.png', bbox_inches='tight', format = 'png', dpi=1200, transparent=True)

    return 0

if __name__ == "__main__":
    for w in W:
        lcm_matrix, eval, evec = lcm_calc(MASS, w, NX, NY)
        # lcm_matrix = lcm_calc(MASS, W, NX, NY)
        plot_lcm(lcm_matrix, w)
