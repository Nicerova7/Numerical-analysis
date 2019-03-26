#Lagrange interpolation
#input : X array of nodes
#      : Y array of value of those nodes
#ouput : C matrix that contains the coefficents of the Lagrange interpolator polynomial
#      : L matrix that contais the coefficients of the Lagrange interpolant

import numpy as np

def lagran(X,Y):

    w = len(X)
    n = w-1
    L = np.zeros((w,w))

    for k in range(0, n+1):
        V = 1
        for j in range(0, n+1):
            if (k != j):
                V = np.convolve(V,np.poly(X[j]))/(X[k]-X[j])
        L[k,:] = V

    C = np.dot(Y,L)
    return C,L

def main():

    #X
    X = [-1.5, -0.75, 0, 0.75, 1.5]
    #Y
    Y = [-14.1014, -0.9315, 0, 0.9315, 14.1014]
    R = lagran(X,Y)
    print("C: \n",R[0])
    print("L: \n",R[1])
main()
