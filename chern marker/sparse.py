import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla


# Parameters
n = 5
N_sites = n**3
t = 1.0 + 0.0j  # Hopping term
W = 16.5  # Disorder strength

# COO format data arrays
data = []
row = []
col = []
print("matrix construction")
# Construct the Hamiltonian matrix with nearest-neighbor hoppings
for ix in range(1, n+1):
    for iy in range(1, n+1):
        for iz in range(1, n+1):
            i = (iz - 1) * n * n + (iy - 1) * n + ix - 1
            
            # Add nearest neighbors in 3D, avoiding going over the edge
            if ix < n:
                j = i + 1
                data.append(t)
                row.append(i)
                col.append(j)
                data.append(np.conj(t))
                row.append(j)
                col.append(i)
            if iy < n:
                j = i + n
                data.append(t)
                row.append(i)
                col.append(j)
                data.append(np.conj(t))
                row.append(j)
                col.append(i)
            if iz < n:
                j = i + n * n
                data.append(t)
                row.append(i)
                col.append(j)
                data.append(np.conj(t))
                row.append(j)
                col.append(i)

# Add random disorder
# np.random.seed(12345)
for ix in range(1, n+1):
    for iy in range(1, n+1):
        for iz in range(1, n+1):
            i = (iz - 1) * n * n + (iy - 1) * n + ix - 1
            eps = np.random.rand() * W - W / 2.0
            # someting_else = 5
            # data.append(someting_else)
            data.append(eps)
            row.append(i)
            col.append(i)

# Create COO matrix
coo = sp.coo_matrix((data, (row, col)), shape=(N_sites, N_sites))

# Convert COO to CSR format
csr = coo.tocsr()


# Extract CSR components
ao = csr.data
jao = csr.indices
iao = csr.indptr

# Print the CSR components
print("Values (ao):")
print(ao[0:11])
print("Column indices (jao):")
print(jao[0:11])
print("Row pointers (iao):")
print(iao[0:11])


# # # Convert CSR to dense format
# hamiltonian = csr.toarray()

# # Plot the Hamiltonian for visualization
# plt.imshow(np.abs(hamiltonian), cmap='viridis')
# plt.colorbar()
# plt.title('Hamiltonian Matrix')
# plt.xlabel('Orbital index')
# plt.ylabel('Orbital index')
# plt.show()
# print("eigenvalue calculation")



from timeit import default_timer as timer
time1 = timer()
# Find the lowest eigenvalues with ARPACK
num_eigenvalues = 4
eigenvalues, eigenvectors = spla.eigsh(csr, k=num_eigenvalues, which='SM', tol=0.01, maxiter = 500000)
# tol=1E-2
time2 = timer()

# Print the lowest eigenvalues
print("Lowest eigenvalues:")
print(eigenvalues)
print("this took", time2 - time1)
# # Create the vector (1, 0, 0, ...)
# vector = np.zeros(N_sites, dtype=np.complex128)
# vector[0] = 1.0
# vector[3] = 1.0
# vector[9] = 2.0
# vector[10] = 1.0
# vector[24] = 7.0

# # Multiply CSR matrix by vector
# result_vector = csr.dot(vector)

# # Print the result
# print("Result of matrix-vector multiplication:")
# print(result_vector)

# lineval, linevec = np.linalg.eigh(hamiltonian)
# print("from dense matrix:")
# print(lineval)
