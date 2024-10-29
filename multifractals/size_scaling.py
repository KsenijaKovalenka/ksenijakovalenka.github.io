import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# GLOBAL PARAMS
# Q = [0]
Q = np.arange(-10,11,1)
# Q = np.arange(0,11,1)
# Q = np.arange(-5,6,1)
NUMV = 4
PATH = "./data/M2W4/"
FILE_PREFIX = "_2_4.0_"
RUN_N = [2000, 1527, 2000]
SAVE_FIGS = False
OUTFILE = "output_spectra2"
SIZES = [100]
SIZES = [50, 80, 100]


def linear(x, a, b):
    return a + b*x 

def coord_to_index(x, y):

    # rewrite for Haldane and 3dTB separately

    index_1 = 2*x + 2*NX*y
    index_2 = index_1 + 1

    return int(index_1), int(index_2)

def read_data(filename):
    # Define the parameters
    n = NX**2 * 2

    # Initialize the array to store vectors
    vector = np.zeros((NUMV, n), dtype=np.complex128)

    # Open the file for reading
    with open(filename, 'r') as file:
        for k in range(NUMV):
            for i in range(n):
                # Read the real and imaginary parts from the file
                line = file.readline().strip()
                re, imag = map(float, line.split())
                vector[k, i] = complex(re, imag)
            
            # Read and skip the empty line separator
            file.readline()

    return vector

def find_moment(wavefunction, q_val, l, nx):

    q = q_val
    box_prob = 0
    size = int(nx / l)
    box_prob_list = np.empty((size**2))
    moment = 0

    # for each box
    # get the box coorinates
    b = 0
    for m in range(size):
        x_list = np.arange(m*l, (m+1)*l)
        for n in range(size):
            y_list = np.arange(n*l, (n+1)*l)
            # converst to list of indicies and run over them to compute the box probability
            indicies = np.array([], dtype=int)
            for x in x_list:
                for y in y_list:
                    index1, index2 = coord_to_index(x, y)
                    indicies = np.append(indicies, index1)
                    indicies = np.append(indicies, index2)
            # calculate the box probability for a set of indicies
            for i in indicies:
                box_prob += np.real(wavefunction[i]* np.conjugate(wavefunction[i]))
            box_prob_list[b] = box_prob
            b += 1
            box_prob = 0
    # sum over all boxes to get the moment
    for prob in box_prob_list:
        moment += prob**q
    return moment

def do_stats(moments):

    moments_mean = np.real(np.mean(moments))
    moments_error = np.std(moments)

    return moments_mean, moments_error

def fit_scaling(x, y, sigma_y):

    # Initial guess
    x0    = np.array([0.0, 0.0])
    if (qval == 0):
        sigma_y = np.ones(len(SIZES)) / 5

    params, cov_mat = curve_fit(linear, x, y, x0, sigma_y)

    # plot the fit
    plt.errorbar(x, y, sigma_y)
    plt.plot(x, linear(x, params[0], params[1]))
    plt.show()

    # extract parameters
    grad = params[1]
    grad_err = np.sqrt(cov_mat[1][1])

    if (qval == 0):
        grad_err = 0

    return grad, grad_err

def plot_tau(x, y, yerr):

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr)

    plt.xlabel('q', fontsize=18)
    plt.ylabel(r'$\tau$(q)', fontsize=18)
    plt.xlim([-11, 11])  # Set x-axis limits
    plt.xticks(np.arange(-10, 11, 1))  # Set x-axis ticks spaced by 1
    # Add vertical and horizontal lines
    plt.axvline(x=0, color='black', linestyle='--', alpha= 0.5)
    plt.axhline(y=-2, color='black', linestyle='--', alpha= 0.5)
    # Add support of a measure point at q=0
    if np.any(x == 0):
        plt.scatter(0, y[np.where(x == 0)], color='purple', zorder=5)  
    plt.show()
    
    return 0

def plot_dimension(x, y, yerr):

    # Find the index of the value 1
    index_to_remove = np.where(x == 1)[0]

    # Remove elements at that index
    new_Q = np.delete(x, index_to_remove)
    new_grad = np.delete(y, index_to_remove)
    new_error = np.delete(yerr, index_to_remove)
    dimension = new_grad / (new_Q - 1)
    dimension_error = new_error / (new_Q - 1)

    plt.scatter(new_Q, dimension)
    plt.errorbar(new_Q, dimension, np.abs(dimension_error))

    plt.xlabel('q', fontsize=18)  # Add x-axis label
    plt.ylabel(r'$D_q$', fontsize=18)  # Add y-axis label
    plt.xlim([-11, 11])  # Set x-axis limits from -10 to 10
    plt.xticks(np.arange(-10, 11, 1))  # Set x-axis ticks spaced by 1

    ax = plt.gca()  # Get current axes
    ticks = ax.get_xticklabels()
    for tick in ticks:
        if tick.get_text() == '1':  # Find the tick corresponding to 1
            tick.set_color('red')  # Set its color to red
        else:
            tick.set_color('black')  # Ensure other ticks are black
    
    # Add vertical and horizontal lines
    plt.axvline(x=0, color='black', linestyle='--', alpha= 0.5)
    plt.axhline(y=2, color='black', linestyle='--', alpha= 0.5)

    # Color the point at x = 0 in purple
    if np.any(new_Q == 0):
        plt.scatter(0, dimension[np.where(new_Q == 0)], color='purple', zorder=5)

    plt.show()
    
    return dimension, dimension_error

# DO THE BOX COUNTING
taus = np.array([])
tau_errors = np.array([])

for qval in Q:

    print(f"qval: {qval}")

    if (qval < 0):
        L = 5
    else: 
        L = 1
    
    means_arr =  np.array([])
    stds_arr =  np.array([])
    lambda_arr =  np.array([])
    for SIZE_I, SIZE in enumerate(SIZES):
        # READ OUT DATA
        N_RUNS = RUN_N[SIZE_I]
        NX = SIZE

        evecs = np.zeros((NUMV*N_RUNS, NX**2 *2), dtype=np.complex128)

        for run in range(N_RUNS):
            file = PATH + str(NX) + "/" + str(NX) + FILE_PREFIX + str(run + 1) + ".dat"
            print('reading' + file)
            evec = read_data(file)
            evecs[run*NUMV:(run+1)*NUMV, :] = evec
        print(f"read evecs for nx = {NX}")
        
        moment_arr =  np.array([])

        for evec in evecs:
            # sort out the readout here
            moment_val = find_moment(evec, qval, L, NX)
            moment_arr = np.append(moment_arr,moment_val)

        lambda_arr = np.append(lambda_arr, L/NX)
        mean_moment, std_moment = do_stats(moment_arr)
        # append to means
        means_arr = np.append(means_arr, mean_moment)
        stds_arr = np.append(stds_arr, std_moment)

    tau, tau_err = fit_scaling(np.log(lambda_arr), 
                            np.log(means_arr), 
                            (stds_arr/np.sqrt(N_RUNS))/means_arr)

    taus = np.append(taus, tau)
    tau_errors = np.append(tau_errors, tau_err)
# plt.show()


# PLOTTING THE TWO MAIN FIRGURES
plot_tau(Q, taus, tau_errors)
dim, dim_err = plot_dimension(Q, taus, tau_errors)

if SAVE_FIGS == True:
    np.savez(OUTFILE, arr1=Q, arr2=taus, arr3=tau_errors, arr4=dim, arr5=dim_err)

# # PLOTTING DoS
# from DOS_PLOTS import plot_dos
    
# plot_dos(evecs[0], NX)