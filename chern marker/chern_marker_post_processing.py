import numpy as np
from spinv import local_chern_marker, make_finite, onsite_disorder
from spinv.example_models import haldane_pythtb
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


#CALCULATING INTEGRATED LCM AND (OPTIONALLY) PLOTTING THE RESULTS

# GLOBALS
INTEGRATION_WIDTH = 3
PLOT = False
# size of the system
NX = 50
NY = 50
# model masses (use linspace bethween 0 and 3 later)

# MASS_LIST = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]
# MASS_LIST = np.array([1.275, 1.35 , 1.425, 1.5  , 1.575, 1.65 , 1.725, 1.8  ,
#        1.875, 1.95 , 2.025, 2.1  , 2.175, 2.25 , 2.325, 2.4  , 2.475,
#        2.55 , 2.625])
# MASS_LIST = np.array([1.275, 2.625])
# W_LIST = np.array([0.0, 9.5])
# disorder strength list (bethween 0 and 10)
# W_LIST = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
# W_LIST = np.linspace(0.0, 9.5, 20)
# path for saving figures
MASS_LIST = np.array([1.275     , 1.30961538, 1.34423077, 1.37884615, 1.41346154,
       1.44807692, 1.48269231, 1.51730769, 1.55192308, 1.58653846,
       1.62115385, 1.65576923, 1.69038462, 1.725     , 1.75961538,
       1.79423077, 1.82884615, 1.86346154, 1.89807692, 1.93269231,
       1.96730769, 2.00192308, 2.03653846, 2.07115385, 2.10576923,
       2.14038462, 2.175     , 2.20961538, 2.24423077, 2.27884615,
       2.31346154, 2.34807692, 2.38269231, 2.41730769, 2.45192308,
       2.48653846, 2.52115385, 2.55576923, 2.59038462, 2.625     ])
W_LIST = np.array([0.        , 0.24358974, 0.48717949, 0.73076923, 0.97435897,
       1.21794872, 1.46153846, 1.70512821, 1.94871795, 2.19230769,
       2.43589744, 2.67948718, 2.92307692, 3.16666667, 3.41025641,
       3.65384615, 3.8974359 , 4.14102564, 4.38461538, 4.62820513,
       4.87179487, 5.11538462, 5.35897436, 5.6025641 , 5.84615385,
       6.08974359, 6.33333333, 6.57692308, 6.82051282, 7.06410256,
       7.30769231, 7.55128205, 7.79487179, 8.03846154, 8.28205128,
       8.52564103, 8.76923077, 9.01282051, 9.25641026, 9.5       ])
INPUT_PATH = "./inputs/repetitions/"
OUTPUT_PATH = "./outputs/"
FIG_PATH = "./figures/"
RUN_N = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

def plot_lcm(matrix, path, w, mass):

    small_size = 16
    medium_size = 20
    bigger_size = 22

    plt.rc('font', size=small_size)          # controls default text sizes
    plt.rc('axes', titlesize=bigger_size-2)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size-2)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small_size)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small_size)    # legend fontsize
    plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title

    # Define colors
    colors = [(0, 'green'), (0.33, 'blue'), (0.5, 'white'), (0.67, 'red'), (1, 'yellow')]
    #colors = [(0, '#A63919'), (0.33, '#8C4324'), (0.5, 'white'), (0.67, '#292984'), (1, '#010590')]
    # Create colormap
    custom_cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)

    plt.figure()
    plt.gcf().set_facecolor('none')
    plt.title(f"w = {w/2}")
    #plt.imshow(matrix, cmap='seismic', origin='lower', extent=(0, matrix.shape[1], 0, matrix.shape[0]), vmin=-np.max(np.abs(matrix)), vmax=np.max(np.abs(matrix)))
    plt.imshow(matrix, cmap=custom_cmap, origin='lower', extent=(0, matrix.shape[1], 0, matrix.shape[0]), vmin=-2, vmax=2)
    plt.xlabel('X')
    plt.ylabel('Y')
    cbar = plt.colorbar(label='Chern Marker')
    numticks = 4
    cbar.locator = plt.MaxNLocator(numticks)
    cbar.update_ticks()
    plt.show()
    # plt.savefig(path + f'w{(w)}'+ f'M{mass}' + '.png', bbox_inches='tight', format = 'png', dpi=1200, transparent=True)

    return 0

if __name__ == "__main__":
    big_boy = np.empty((np.size(RUN_N), len(MASS_LIST), len(W_LIST)))
    for ir, r in enumerate(RUN_N):

        topology_class = np.zeros((len(MASS_LIST), len(W_LIST)))
        for mass_index, mass in enumerate(MASS_LIST):
            for disorder_index, w in enumerate(W_LIST):

                filename = f'w{w}' + f'M{mass}' + '.npy'
                #lcm_data = np.load(INPUT_PATH + filename)
                lcm_data = np.load(INPUT_PATH + f'run{r}/' + filename)

                l = INTEGRATION_WIDTH
                m = len(lcm_data[0])
                topology_class[mass_index, disorder_index] = sum(sum(lcm_data[0:l])) + sum(sum(lcm_data[m-l:m])) + sum(sum(lcm_data[l:m-l, 0:l])) + sum(sum(lcm_data[l:m-l, m-l:m]))

                if PLOT:
                    if os.path.exists(FIG_PATH):
                        plot_lcm(lcm_data, FIG_PATH, w, mass)
                    else: 
                        print("Could not find the path for figures :(")
        big_boy[ir, :, :] = topology_class
        # plotting the topology classification
        # heatmap
        small_size = 16
        medium_size = 20
        bigger_size = 22

        plt.rc('font', size=small_size)          # controls default text sizes
        plt.rc('axes', titlesize=bigger_size-2)     # fontsize of the axes title
        plt.rc('axes', labelsize=medium_size-2)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=small_size)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=small_size)    # fontsize of the tick labels
        plt.rc('legend', fontsize=small_size)    # legend fontsize
        plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title

        # Define colors
        #colors = [(0, 'green'), (0.33, 'blue'), (0.5, 'white'), (0.67, 'red'), (1, 'yellow')]
        colors = [(0, '#A63919'), (0.33, '#8C4324'), (0.5, 'white'), (0.67, '#292984'), (1, '#010590')]
        # Create colormap
        custom_cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)

        plt.figure()
        plt.gcf().set_facecolor('none')
        #plt.title("Integrated LCM plot")
        #plt.imshow(matrix, cmap='seismic', origin='lower', extent=(0, matrix.shape[1], 0, matrix.shape[0]), vmin=-np.max(np.abs(matrix)), vmax=np.max(np.abs(matrix)))
        #plt.imshow(topology_class, cmap="cool", origin='lower', extent=(0, topology_class.shape[1], 0, topology_class.shape[0]), vmin=0, vmax=2000)
        # print(topology_class)
        plt.imshow(topology_class, cmap="cool", origin='lower', vmin=0, vmax=2000)
        # plt.imshow(topology_class, cmap=custom_cmap, origin='lower', vmin=0, vmax=2000)
        # Remove ticks from both x and y axes
        # Set ticks and labels
        # Set ticks and labels with 5 ticks each
        x_indices = np.linspace(0, len(W_LIST) - 1, 5).astype(int)
        y_indices = np.linspace(0, len(MASS_LIST) - 1, 5).astype(int)

        plt.xticks(ticks=x_indices, labels=np.round(W_LIST[x_indices], 2), rotation=45)
        plt.yticks(ticks=y_indices, labels=np.round(MASS_LIST[y_indices], 3))
        plt.xlabel('disorder strength x2')
        plt.ylabel('mass')
        cbar = plt.colorbar(label='Integrated LCM')
        numticks = 4
        cbar.locator = plt.MaxNLocator(numticks)
        cbar.update_ticks()
        plt.show()
        #plt.savefig(path + f'w{(w)}'+ f'M{mass}' + '.png', bbox_inches='tight', format = 'png', dpi=1200, transparent=True)
    
averaged_phase_diagram = np.mean(big_boy, axis = 0)
# heatmap
small_size = 16
medium_size = 20
bigger_size = 22

plt.rc('font', size=small_size)          # controls default text sizes
plt.rc('axes', titlesize=bigger_size-2)     # fontsize of the axes title
plt.rc('axes', labelsize=medium_size-2)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small_size)    # fontsize of the tick labels
plt.rc('legend', fontsize=small_size)    # legend fontsize
plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title

plt.figure()
plt.gcf().set_facecolor('none')
# print(averaged_phase_diagram)
plt.imshow(averaged_phase_diagram, cmap="cool", origin='lower', vmin=0, vmax=2000)
x_indices = np.linspace(0, len(W_LIST) - 1, 5).astype(int)
y_indices = np.linspace(0, len(MASS_LIST) - 1, 5).astype(int)

plt.xticks(ticks=x_indices, labels=np.round(W_LIST[x_indices], 2), rotation=45)
plt.yticks(ticks=y_indices, labels=np.round(MASS_LIST[y_indices], 3))
plt.xlabel('disorder strength x2')
plt.ylabel('mass')
plt.title('Averaged Phase Diagram')
cbar = plt.colorbar(label='Integrated LCM')
numticks = 4
cbar.locator = plt.MaxNLocator(numticks)
cbar.update_ticks()
plt.show()
#plt.savefig(path + f'w{(w)}'+ f'M{mass}' + '.png', bbox_inches='tight', format = 'png', dpi=1200, transparent=True)


    # # separate line for each of the masses
    # for i, mass in enumerate(MASS_LIST):
    #     plt.plot(W_LIST, topology_class[i, :], label=f"Mass: {mass}")
    #     plt.xlabel("disorder strength")
    #     plt.ylabel("integrated LCM")
    # plt.legend(fontsize=9) 
    # plt.show()

    # from matplotlib.colors import LinearSegmentedColormap
    # from mpl_toolkits.mplot3d import Axes3D


    # # Create a meshgrid
    # x, y = np.meshgrid(W_LIST, MASS_LIST)

    # # Create a figure and a 3D Axes
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # plt.tight_layout()
    # plt.gcf().set_facecolor('none')
    # # Plot the wireframe
    # wire = ax.plot_wireframe(x, y, topology_class, color='#882E15', linewidth=0.5)
    # # Remove background and grid
    # ax.xaxis.pane.fill = False  # Remove background pane for x-axis
    # ax.yaxis.pane.fill = False  # Remove background pane for y-axis
    # ax.zaxis.pane.fill = False  # Remove background pane for z-axis
    # ax.grid(False)              # Turn off the grid
    # # Adjust z-axis ticks and labels
    # z_ticks = np.linspace(0, 2000, num=4)
    # ax.set_zticks(z_ticks)
    # ax.set_zticklabels([f'{tick:.0f}' for tick in z_ticks], rotation=25, ha='right')

    # ax.view_init(36, 56)
    # #ax.set_axis_off()
    # ax.tick_params(axis='x', colors='#51230F')  # Change x-axis tick color to red
    # ax.tick_params(axis='y', colors='#51230F') # Change y-axis tick color to blue
    # ax.tick_params(axis='z', colors='#51230F') # Change z-axis tick color to green
    # ax.xaxis.line.set_color('#51230F')   # Change x-axis line color to red
    # ax.yaxis.line.set_color('#51230F')  # Change y-axis line color to blue
    # ax.zaxis.line.set_color('#51230F') # Change z-axis line color to green

    # plt.savefig('wire.png', bbox_inches='tight', format = 'png', dpi=1200, transparent=True)



    # # # Show plot
    # # plt.show()

    