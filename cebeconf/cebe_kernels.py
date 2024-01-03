import numpy as np

# Gaussian and Laplacian kernels
def kernel(option,sigma,dT, dQ):
    if option == 'L':
        dij=np.sum(np.abs(dT-dQ))
        val = np.exp(-dij / sigma)
    elif option == 'G':
        dij=np.sqrt(np.sum(np.abs(dT-dQ)**2))
        val = np.exp( -dij**2   / (2*sigma**2) )
    return val
