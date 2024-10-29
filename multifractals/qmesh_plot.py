import numpy as np
import matplotlib.pyplot as plt

path = "./qmesh_data/"
q = -4

data = np.loadtxt( path  + f'out_q{q}.txt')  # Use exact delimiter

# Split the data into four separate arrays (columns)
m, w, grad, dim = data.T  # Transpose to unpack the columns


# add 0s to the missing entries (so far)
for i in range(400 - len(m)):
    m = np.append(m, 0)
    w = np.append(w, 0)
    grad = np.append(grad, 0)
    dim = np.append(dim, 0)

m = np.reshape(m, (20, 20))
m = m.T 

plt.imshow(m, cmap='viridis', interpolation='nearest', origin='lower')
plt.colorbar()
plt.show()

w = np.reshape(w, (20, 20))
w = w.T 

plt.imshow(w, cmap='viridis', interpolation='nearest', origin='lower')
plt.colorbar()
plt.show()

grad = np.reshape(grad, (20, 20))
grad = grad.T

plt.imshow(grad, cmap='viridis', interpolation='nearest', origin='lower')
plt.colorbar()
plt.show()

dim = np.reshape(dim, (20, 20))
dim = dim.T

plt.imshow(dim, cmap='viridis', interpolation='nearest', origin='lower')
plt.colorbar()
plt.show()

mass = m[:,0]
wdis = w[0]

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
plt.imshow(dim, cmap='viridis', origin='lower')
x_indices = np.linspace(0, len(wdis) - 1, 5).astype(int)
y_indices = np.linspace(0, len(mass) - 1, 5).astype(int)

plt.xticks(ticks=x_indices, labels=np.round(wdis[x_indices], 2), rotation=45)
plt.yticks(ticks=y_indices, labels=np.round(mass[y_indices], 3))
plt.xlabel('disorder strength x2')
plt.ylabel('mass')
plt.title(f'q = {q}')
cbar = plt.colorbar(label='Dimension')
numticks = 4
cbar.locator = plt.MaxNLocator(numticks)
cbar.update_ticks()
plt.show()