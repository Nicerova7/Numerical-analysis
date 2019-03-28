#Newton divided diferences
#input : X array of nodes
#      : Y array of value of those nodes
#ouput : C vector that contains coefficients of Newton interpolation polynomial

import numpy as np

def newpoly(X,Y):

    n = len(X)
    D = np.zeros((n,n))
    D[:,1] = np.tranpose(Y)

    for j in range(2, n):
        for k in range(j, n):
            D[k,j] = (D[k,j-1]-D[k-1,j-1])/(X[k]-X[k-j+1])

    C = D[n,n]

def main():


main()
