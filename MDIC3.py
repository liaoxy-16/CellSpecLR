import numpy as np


def MDIC3(G,E):
    T = np.dot(G, E)   # Compute the dot product or matrix multiplication between two arrays.
    #print('T:\n',T)
    Tp = np.linalg.pinv(T)  # Compute the pseudoinverse of a matrix.
    #print('Tp:\n', Tp)
    M = np.dot(Tp, E)
    #print('M:\n', M)

    return M







