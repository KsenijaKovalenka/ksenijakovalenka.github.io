import numpy as np
import matplotlib.pyplot as plt
from pythtb import * # import TB model class
import scipy.optimize as optimization

nx = 50

def coord_to_index(x, y):

    # rewrite for Haldane and 3dTB separately

    index_1 = 2*x + 2*nx*y
    index_2 = index_1 + 1

    return int(index_1), int(index_2)

def read_data(filename):
    # Define the parameters
    n = nx**2 * 2
    numv = 4

    # Initialize the array to store vectors
    vector = np.zeros((numv, n), dtype=np.complex128)

    # Open the file for reading
    with open(filename, 'r') as file:
        for k in range(numv):
            for i in range(n):
                # Read the real and imaginary parts from the file
                line = file.readline().strip()
                re, imag = map(float, line.split())
                vector[k, i] = complex(re, imag)
            
            # Read and skip the empty line separator
            file.readline()

    return vector

def find_moments(wavefunction, q_val, l_list):

    q = q_val

    moment = 0
    box_prob = 0
    moment_list = np.array([])
    lambda_list = np.array([])
    for l in l_list:
        # print(f'l = {l}')
        size = int(nx / l)
        box_prob_list = np.empty((size**2))
        lambda_list = np.append(lambda_list, l/nx)
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
        
        for prob in box_prob_list:
            # print("box_prob: ")
            # print(prob)
            moment += prob**q
        # print('moment')
        # print(moment)
        moment_list = np.append(moment_list, moment)
        moment = 0
    return lambda_list, moment_list

def fit_func(x, a, b):
    return a + b*x 

def find_spectra(Q):
    gradients = np.array([])

    for qval in Q:
        # print(f"qval: {qval}")
        lambda_matrix =  np.array([])
        moment_matrix =  np.array([])

        for enum, evec in enumerate(evecs):
            # print(f"evec {enum}")
            # sort outh the readout here
            lambdas, moments = find_moments(evec, qval, l_list)
            lambda_matrix = np.append(lambda_matrix,lambdas)
            moment_matrix = np.append(moment_matrix,moments)

        # moment_matrix_try = np.reshape(moment_matrix, (3,40))
        # lambda_matrix_try = np.reshape(lambda_matrix, (3,40))
        moment_matrix_try = np.empty((np.size(l_list),n_runs*numv))
        lambda_matrix_try = np.empty((np.size(l_list),n_runs*numv))
        
        for indl, l in enumerate(l_list):
            ind = np.where(lambda_matrix == l/nx)[0]
            moment_matrix_try[indl] = moment_matrix[ind]
            lambda_matrix_try[indl] = lambda_matrix[ind]

        mean = np.mean(moment_matrix_try, axis = 1)
        mean = np.real(mean)
        std = np.std(moment_matrix_try, axis = 1)

        std = std/mean
        mean = np.log(mean)

        xdata = np.log(lambda_matrix_try[:,0])
        ydata = mean
        # Initial guess.
        x0    = np.array([0.0, 0.0])
        sigma = std
        if (qval == 0):
            sigma = np.array([0.5, 0.5, 0.5])

        params = optimization.curve_fit(fit_func, xdata, ydata, x0, sigma)
        # print("qval")
        # print(qval)
        # plt.errorbar(xdata, mean, std)
        # plt.plot(xdata, fit_func(xdata, params[0][0], params[0][1]))
        # plt.show()
        gradients = np.append(gradients, params[0][1])
    return gradients

def plot_spectra(Q, gradients):
    # PLOTTING THE TWO MAIN FIRGURES--------------------------------------------------

    plt.scatter(Q, gradients)
    plt.xlabel('q', fontsize=18)  # Add x-axis label
    plt.ylabel(r'$\tau$(q)', fontsize=18)  # Add y-axis label
    plt.xlim([-11, 11])  # Set x-axis limits from -10 to 10
    plt.xticks(np.arange(-10, 11, 1))  # Set x-axis ticks spaced by 1
    # Add vertical line at x = 0
    plt.axvline(x=0, color='black', linestyle='--', alpha= 0.5)

    # Add horizontal line at y = 3
    plt.axhline(y=-2, color='black', linestyle='--', alpha= 0.5)
    plt.scatter(0, gradients[np.where(Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it
    plt.show()

    # Find the index of the value 4
    index_to_remove = np.where(Q == 1)[0]

    # Remove the element at that index
    new_Q = np.delete(Q, index_to_remove)
    new_grad = np.delete(gradients, index_to_remove)
    dimension = new_grad / (new_Q - 1)
    plt.scatter(new_Q, dimension)
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
    # Add vertical line at x = 0
    plt.axvline(x=0, color='black', linestyle='--', alpha= 0.5)

    # Add horizontal line at y = 3
    plt.axhline(y=2, color='black', linestyle='--', alpha= 0.5)

    # Color the point at x = 0 in purple
    plt.scatter(0, dimension[np.where(new_Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it
    plt.show()

    return new_Q, dimension

# GLOBALS
l_list = np.array([5, 10, 25])   #for 50by50
numv = 4
Q = np.arange(-10,11,1)
mass_values = np.array([1.34, 2.49, 1.76, 1.34, 2.49])
w_values = np.array([0.49, 0.73, 3.90, 9.01, 9.01])

for i in range(len(mass_values)):
    m = mass_values[i]
    w = w_values[i]
    print(f"mass:{m} w:{w}")
    path = "./disorder_data/differentpoints/"
    n_runs = 10
    evecs = np.zeros((4*n_runs, nx**2 *2), dtype=np.complex128)
    for run in range(n_runs):

        file = path + f"m{m:.2f}w{w:.2f}_hald_evecs_" + str(run + 1) + ".dat"
        evec = read_data(file)
        evecs[run*4:(run+1)*4, :] = evec

    gradients = find_spectra(Q)
    new_Q, dimension = plot_spectra(Q, gradients)

    # save the spectra
    outpath = "./outputs/differentpoints/"
    outfile = "MASS" + str(m) + "W" + str(w)
    np.savez(outpath + outfile, arr1=Q, arr2=gradients, arr3=new_Q, arr4=dimension)