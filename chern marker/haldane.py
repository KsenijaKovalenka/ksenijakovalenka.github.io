#!/usr/bin/env python

# Haldane my_model from Phys. Rev. Lett. 61, 2015 (1988)

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function
from pythtb import * # import TB my_model class
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

# define lattice vectors
lat=[[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
# define coordinates of orbitals
orb=[[1./3.,1./3.],[2./3.,2./3.]]

# make two dimensional tight-binding Haldane my_model
my_model=tb_model(2,2,lat,orb)

# set my_model parameters
delta= 2.0
t= -1.0
t2 = t / 3. * np.exp((1.j)*np.pi/2.)
t2c= t2.conjugate()
# # Add random disorder
W = 5.0
np.random.seed()
eps1 = np.random.rand() * W - W / 2.0 #THIS PART IS BULLSHIT, WON'T WORK (see notes)
eps2 = np.random.rand() * W - W / 2.0
my_model.set_onsite([delta,-delta])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t, 0, 1, [ 0, 0])
my_model.set_hop(t, 1, 0, [ 1, 0])
my_model.set_hop(t, 1, 0, [ 0, 1])
# add second neighbour complex hoppings
my_model.set_hop(t2 , 0, 0, [ 1, 0])
my_model.set_hop(t2 , 1, 1, [ 1,-1])
my_model.set_hop(t2 , 1, 1, [ 0, 1])
my_model.set_hop(t2c, 1, 1, [ 1, 0])
my_model.set_hop(t2c, 0, 0, [ 1,-1])
my_model.set_hop(t2c, 0, 0, [ 0, 1])

# print tight-binding my_model
my_model.display()

# generate list of k-points following a segmented path in the BZ
# list of nodes (high-symmetry points) that will be connected
path=[[0.,0.],[2./3.,1./3.],[.5,.5],[1./3.,2./3.], [0.,0.]]
# labels of the nodes
label=(r'$\Gamma $',r'$K$', r'$M$', r'$K^\prime$', r'$\Gamma $')

knum = 101
# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=my_model.k_path(path,knum)
# inputs:
#   path: see above
#   101: number of interpolated k-points to be plotted
# outputs:
#   k_vec: list of interpolated k-points
#   k_dist: horizontal axis position of each k-point in the list
#   k_node: horizontal axis position of each original node

# obtain eigenvalues to be plotted
evals=my_model.solve_all(k_vec)

# figure for bandstructure

fig, ax = plt.subplots()
# specify horizontal axis details
# set range of horizontal axis
ax.set_xlim(k_node[0],k_node[-1])
# put tickmarks and labels at node positions
ax.set_xticks(k_node)
ax.set_xticklabels(label)
# add vertical lines at node positions
for n in range(len(k_node)):
  ax.axvline(x=k_node[n],linewidth=0.5, color='k')
# put title
ax.set_title("Haldane my_model band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

# plot first band
ax.plot(k_dist,evals[0])
# plot second band
ax.plot(k_dist,evals[1])

# make an PDF figure of a plot
fig.tight_layout()
#fig.savefig("haldane_band.pdf")

n = 50

# cutout finite my_model first along direction x with no PBC
tmp1_my_model=my_model.cut_piece(n,0,glue_edgs=False)
# cutout also along y direction with no PBC
fin_my_model=tmp1_my_model.cut_piece(n,1,glue_edgs=False)


# print tight-binding my_model
fin_my_model.display()

# solve finite my_models
(evals,evecs)=fin_my_model.solve_all(eig_vectors=True)

# pick index of state in the middle of the gap
ed=fin_my_model.get_num_orbitals()//2

# draw one of the edge states in both cases
(fig,ax)=fin_my_model.visualize(0,1,eig_dr=evecs[ed,:],draw_hoppings=True)
ax.set_title("Edge state for finite my_model without periodic direction")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.tight_layout()
fig.savefig("edge_state.pdf")
plt.show()

# visualise the Hamiltonian
# Print orbital positions and Hamiltonian
orbital_positions = fin_my_model.get_orb()
hamiltonian = fin_my_model._gen_ham()

# print("Orbital positions (real space):")
# for i, pos in enumerate(orbital_positions):
#     print(f"Orbital {i}: Position {pos}")

# print("\nHamiltonian matrix:")
# print(hamiltonian)
count = 0
for thing in hamiltonian:
  for thing2 in thing:
    if (abs(thing2) > 0): count = count + 1
print(count)


# Plotting the Hamiltonian for visualization (optional)
import matplotlib.pyplot as plt
plt.imshow(hamiltonian.real, cmap='viridis')
plt.colorbar()
plt.title('Hamiltonian Matrix')
plt.xlabel('Orbital index')
plt.ylabel('Orbital index')
plt.show()

comp = hamiltonian.imag
plt.imshow(comp, cmap='viridis')
plt.colorbar()
plt.title('Hamiltonian Matrix')
plt.xlabel('Orbital index')
plt.ylabel('Orbital index')
plt.show()

# Diagonalize the Hamiltonian
eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)

# Print the results
print("Eigenvalues:", eigenvalues)

# # Create COO matrix
# csr = sp.csr_matrix(hamiltonian)


# # Extract CSR components
# ao = csr.data
# jao = csr.indices
# iao = csr.indptr

# # Print the CSR components
# print("Values (ao):")
# print(ao[0:10])
# print("Column indices (jao):")
# print(jao[0:10])
# print("Row pointers (iao):")
# print(iao[0:10])

