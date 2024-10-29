import numpy as np
import matplotlib.pyplot as plt
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

path = "./SPECTRALFUNC/"

n_runs = 6
evecs = np.zeros((4*n_runs, nx**2 *2), dtype=np.complex128)
for run in range(n_runs):

    file = path + "m3_hald_evecs_" + str(run + 1) + ".dat"
    evec = read_data(file)
    evecs[run*4:(run+1)*4, :] = evec

for evec in evecs:

    wavefunction_vis = evec
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


    # Create a heatmap using matplotlib. viridis, hot, flag
    plt.imshow(field, cmap='hot', interpolation='nearest')
    plt.colorbar()

    # Display the heatmap
    plt.show()