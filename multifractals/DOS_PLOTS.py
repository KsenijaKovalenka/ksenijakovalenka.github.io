import matplotlib.pyplot as plt
import numpy as np

def coord_to_index(x, y, nx):

    # rewrite for Haldane and 3dTB separately

    index_1 = 2*x + 2*nx*y
    index_2 = index_1 + 1

    return int(index_1), int(index_2)

def plot_dos(evec, nx):
    wavefunction_vis = evec
    field = np.empty((nx,nx))

    # visualise the wavefunction itself
    totalprob = 0
    for y in range(nx):
        for x in range(nx):
            index1, index2 = coord_to_index(x, y, nx)
            prob = np.real(wavefunction_vis[index1]*np.conjugate(wavefunction_vis[index1]) + \
                wavefunction_vis[index2]*np.conjugate(wavefunction_vis[index2]))
            field[x, y] = prob
            totalprob += prob

    plt.imshow(field, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.show()

def plot_boxes(evec, nx, l_list):
    box_field = np.empty((nx,nx))
    box_prob = 0
    box_prob_list = np.array([])
    lambda_list = np.array([])

    for l in l_list:
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
                        index1, index2 = coord_to_index(x, y, nx)
                        indicies = np.append(indicies, index1)
                        indicies = np.append(indicies, index2)
                # calculate the box probability for a set of indicies
                for i in indicies:
                    box_prob += evec[i]* np.conjugate(evec[i])
                box_prob_list = np.append(box_prob_list, box_prob)
                for x in x_list:
                    for y in y_list:
                        box_field[x, y] = np.real(box_prob)
                box_prob = 0

        plt.imshow(box_field, cmap='viridis', interpolation='nearest')
        plt.colorbar()
        plt.show()