#Newton divided diferences
#input : X array of nodes
#      : Y array of value of those nodes
#ouput : C vector that contains coefficients of Newton interpolation polynomial

import numpy as np

def newpoly(X,Y):

    n = len(X)
    D = np.zeros((n,n))
    D[:,0] = np.transpose(Y)

    for j in range(1, n):
        for k in range(j, n):
            D[k,j] = (D[k,j-1]-D[k-1,j-1])/(X[k]-X[k-j])

    C = D[n-1,n-1]

    for k in range(n-1-1,-1,-1):
        C = np.convolve(C,np.poly(X[k]))
        m = len(C) - 1 #cuz in python array starts at 0
        C[m] +=  D[k,k]

    return C,D

def main():

    #X'
    X = [0, 2, 3, 5, 6, 1]
    #Y' (x^3)
    Y = [0, 8, 27, 125, 216, 1]
    R = newpoly(X,Y)
    print("C: \n",R[0])
    print("D: \n",R[1])

main()
