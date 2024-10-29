import numpy as np
import matplotlib.pyplot as plt
from pythtb import * # import TB model class
import scipy.optimize as optimization

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

nx = 150
# path = "./disorder_data/M2W4/"
path = "./disorder_data/150by150/"
n_runs = 10
evecs = np.zeros((4*n_runs, nx**2 *2), dtype=np.complex128)
for run in range(n_runs):

    file = path + "150by150m2w4.0_hald_evecs_" + str(run + 1) + ".dat"
    evec = read_data(file)
    evecs[run*4:(run+1)*4, :] = evec

# l_list = np.array([5, 10, 25])
l_list = np.array([10, 15, 25, 30, 50])
numv = 4
# Q = [-5, 0, 5]
Q = np.arange(-10,11,1)
# Q = np.arange(0,11,1)
# Q = np.arange(-5,6,1)

def coord_to_index(x, y):

    # rewrite for Haldane and 3dTB separately

    index_1 = 2*x + 2*nx*y
    index_2 = index_1 + 1

    return int(index_1), int(index_2)

def edge_check(m, n, size):
    flag = False
    if (m == 0 or m == size - 1):
        flag = True
    elif (n == 0 or n == size - 1):
        flag = True

    return flag
    
def find_moments_edge(wavefunction, q_val, l_list):

    q = q_val

    moment = 0
    box_prob = 0
    moment_list_edge = np.array([])
    moment_list_bulk = np.array([])
    lambda_list = np.array([])
    for l in l_list:
        # print(f'l = {l}')
        size = int(nx / l)
        box_prob_list_edge = np.empty((size*4 - 4))
        box_prob_list_bulk = np.empty((size**2 - (size*4 - 4)))
        lambda_list = np.append(lambda_list, l/nx)
        # for each box
        # get the box coorinates
        b = 0
        v = 0
        for m in range(size):
            x_list = np.arange(m*l, (m+1)*l)
            for n in range(size):
                y_list = np.arange(n*l, (n+1)*l)
                # converst to list of indicies and run over them to compute the box probability
                edge = edge_check(m, n, size)
                indicies = np.array([], dtype=int)
                for x in x_list:
                    for y in y_list:
                        index1, index2 = coord_to_index(x, y)
                        indicies = np.append(indicies, index1)
                        indicies = np.append(indicies, index2)
                # calculate the box probability for a set of indicies
                for i in indicies:
                    box_prob += np.real(wavefunction[i]* np.conjugate(wavefunction[i]))
                if (edge):
                    box_prob_list_edge[b] = box_prob
                    b += 1
                    box_prob = 0
                else:
                    box_prob_list_bulk[v] = box_prob
                    v += 1
                    box_prob = 0

        for prob in box_prob_list_edge:
            moment += prob**q
        moment_list_edge = np.append(moment_list_edge, moment)
        moment = 0

        for prob in box_prob_list_bulk:
            moment += prob**q
        moment_list_bulk = np.append(moment_list_bulk, moment)
        moment = 0
        # print('edge:')
        # print(moment_list_edge)
        # print('bulk:')
        # print(moment_list_bulk)

    return lambda_list, moment_list_edge, moment_list_bulk

def plot_spectra2(Q, gradients1, gradients2):
# PLOTTING THE TWO MAIN FIRGURES--------------------------------------------------

    plt.scatter(Q, gradients1)
    plt.scatter(Q, gradients2)
    plt.xlabel('q', fontsize=18)  # Add x-axis label
    plt.ylabel(r'$\tau$(q)', fontsize=18)  # Add y-axis label
    plt.xlim([-11, 11])  # Set x-axis limits from -10 to 10
    plt.xticks(np.arange(-10, 11, 1))  # Set x-axis ticks spaced by 1
    # Add vertical line at x = 0
    plt.axvline(x=0, color='black', linestyle='--', alpha= 0.5)

    # Add horizontal line at y = 3
    plt.axhline(y=-2, color='black', linestyle='--', alpha= 0.5)
    plt.scatter(0, gradients1[np.where(Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it
    plt.scatter(0, gradients2[np.where(Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it
    plt.show()

    # Find the index of the value 4
    index_to_remove = np.where(Q == 1)[0]

    # Remove the element at that index
    new_Q = np.delete(Q, index_to_remove)
    new_grad1 = np.delete(gradients1, index_to_remove)
    new_grad2 = np.delete(gradients2, index_to_remove)
    dimension1 = new_grad1 / (new_Q - 1)
    dimension2 = new_grad2 / (new_Q - 1)
    plt.scatter(new_Q, dimension1)
    plt.scatter(new_Q, dimension2)
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
    plt.scatter(0, dimension1[np.where(new_Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it

    plt.show()

    # PLOTTING THE TWO MAIN FIRGURES END--------------------------------------------------
    return 0
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
    # plt.scatter(0, gradients[np.where(Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it
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
    # plt.scatter(0, dimension[np.where(new_Q == 0)], color='purple', zorder=5)  # Ensure Q has 0 to color it

    plt.show()

    # PLOTTING THE TWO MAIN FIRGURES END--------------------------------------------------
    return 0

def linear(x, a, b):
    return a + b*x 

def fit_moments(lambda_matrix, moment_matrix, removezero = False):
    moment_matrix_try = np.empty((np.size(l_list),n_runs*numv))
    lambda_matrix_try = np.empty((np.size(l_list),n_runs*numv))
    
    for indl, l in enumerate(l_list):
        ind = np.where(lambda_matrix == l/nx)[0]
        moment_matrix_try[indl] = moment_matrix[ind]
        lambda_matrix_try[indl] = lambda_matrix[ind]

    mean = np.mean(moment_matrix_try, axis = 1)
    mean = np.real(mean)
    std = np.std(moment_matrix_try, axis = 1)

    if (removezero):
        mean = mean[:-1]
        std = std[:-1]

    std = std/mean
    mean = np.log(mean)

    xdata = np.log(lambda_matrix_try[:,0])
    if (removezero):
        xdata = np.log(lambda_matrix_try[:2,0])
    ydata = mean
    # Initial guess.
    x0    = np.array([0.0, 0.0])
    sigma = std
    if (qval == 0):
        # sigma = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
        sigma = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
        if (removezero):
            sigma = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    params = optimization.curve_fit(linear, xdata, ydata, x0, sigma)
    print("qval")
    print(qval)
    plt.errorbar(xdata, mean, std)
    plt.plot(xdata, linear(xdata, params[0][0], params[0][1]))
    plt.show()
    return params[0][1]


# MAIN CODE ------------------------------------------------------
gradients_edge = np.array([])
gradients_bulk = np.array([])

for qval in Q:
    print(f"qval: {qval}")
    lambda_matrix =  np.array([])
    moment_matrix_edge =  np.array([])
    moment_matrix_bulk =  np.array([])

    for enum, evec in enumerate(evecs):
        # print(f"evec {enum}")
        # sort outh the readout here
        lambdas, moments_edge, moments_bulk = find_moments_edge(evec, qval, l_list)
        lambda_matrix = np.append(lambda_matrix,lambdas)
        moment_matrix_edge = np.append(moment_matrix_edge, moments_edge)
        moment_matrix_bulk = np.append(moment_matrix_bulk, moments_bulk)
    # print(lambda_matrix)
    # print("edge")
    # print(moment_matrix_edge)
    gradient_edge = fit_moments(lambda_matrix, moment_matrix_edge)
    # print("bulk")
    # print(moment_matrix_bulk)
    gradient_bulk = fit_moments(lambda_matrix, moment_matrix_bulk, removezero = False)
    gradients_edge = np.append(gradients_edge, gradient_edge)
    gradients_bulk = np.append(gradients_bulk, gradient_bulk)


plot_spectra(Q, gradients_edge)
plot_spectra(Q, gradients_bulk)

wavefunction_vis = evecs[0]
field = np.empty((nx,nx))

# visualise the wavefunction itself
totalprob = 0
for y in range(nx):
    for x in range(nx):
        index1, index2 = coord_to_index(x, y)
        prob = wavefunction_vis[index1]*np.conjugate(wavefunction_vis[index1]) + \
            wavefunction_vis[index2]*np.conjugate(wavefunction_vis[index2])
        field[x, y] = prob
        totalprob += prob


# Create a heatmap using matplotlib
plt.imshow(field, cmap='viridis', interpolation='nearest')
plt.colorbar()

# Display the heatmap
plt.show()

# lambda_list = np.array([5, 10, 25])
lambda_list = np.array([10, 15, 25, 30, 50])
box_field = np.empty((nx,nx))
moment = 0
box_prob = 0
box_prob_list = np.array([])
moment_list = np.array([])
lambda_list = np.array([])
for l in l_list:
    # print(f'l = {l}')
    size = int(nx / l)
    lambda_list = np.append(lambda_list, l/nx)
    # for each box
    # get the box coorinates
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
                box_prob += wavefunction_vis[i]* np.conjugate(wavefunction_vis[i])
            box_prob_list = np.append(box_prob_list, box_prob)
            edge = edge_check(m, n, size)
            if (edge):
                for x in x_list:
                    for y in y_list:
                        box_field[x, y] = 0.0
                box_prob = 0
            else:
                for x in x_list:
                    for y in y_list:
                        box_field[x, y] = box_prob
                box_prob = 0

    plt.imshow(box_field, cmap='viridis', interpolation='nearest')
    plt.colorbar()

    # Display the heatmap
    plt.show()