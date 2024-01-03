import numpy as np

def kernel(option,sigma,dT, dQ):
    '''
    Calculates kernel matrix element

            Input:
                    option (str): Kernel definition, 'L' or 'G'
                    sigma (float): Kernel width
                    dT(np.array, float): descriptor for atom T
                    dQ(np.array, float): descriptor for atom Q

            Returns:
                    val (float): Kernel matrix element
    '''
    if option == 'L':
        dij=np.sum(np.abs(dT-dQ))
        val = np.exp(-dij / sigma)
    elif option == 'G':
        dij=np.sqrt(np.sum(np.abs(dT-dQ)**2))
        val = np.exp( -dij**2   / (2*sigma**2) )
    return val
