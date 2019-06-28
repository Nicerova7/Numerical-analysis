## Regresive difference method to solve heat equation

"""
 du(x,t)/dt - alpha^2*d^2u(x,y)/dx^2  = 0
 0 < x < l
 0 < t < T
  
 Boundary Conditions
  
 u(0,t) = u(l,t) = 0 ,  0 < t < T

 Initial conditions
 
 u(x,0) = f(x)  ,   0 <= x <= l

INPUT: l,T,alpha
       m >= 3
       N >= 1

OUTPUT Wij aproximation to u(xi,tj) to all i = 1, ... , m-1 and
                                           j = 1, ... , N

"""
import numpy as np

def diffRegressive(l,T,alpha,m,N):

    h = l/m
    k = T/N
    lambd = alpha**2*k/(h**2)

    w = np.zeros([m])
    l = np.zeros([m])
    u = np.zeros([m])
    W = np.zeros([m+1,N])
    z = np.zeros([m])

    for i in range(1,m):  w[i-1] = f(i*h) # Initial values
    print("W initial values (1:m)")
    print(w)
   
 
    # Solve a tridiagonal linear system by Crout factorization
    l[0] = 1 + 2*lambd
    u[0] = -lambd/l[0]

    for i in range(2,m-1):
        l[i-1] = 1 + 2*lambd + lambd*u[i-2]
        u[i-1] = -lambd/l[i-1]

    l[m-2] = 1 + 2*lambd + lambd*u[m-3]

    for j in range(1,N+1):
        t = j*k
        z[0] = w[0]/l[0]

        for i in range(2,m):
            z[i-1] = (w[i-1] + lambd*z[i-2])/l[i-1]

        w[m-2] = z[m-2]

        for i in range(m-2,0,-1):
            w[i-1] = z[i-1] - u[i-1]*w[i]

        W[0,j-1] = f(0)  # restriction
        W[1:,j-1] = w
            
    return W

def f(x):
    return np.sin(np.pi*x)
    
    
def main():

    l = 1.0  # limit of x
    T = 1.0  # limit of time
    alpha = 1.0 # alpha coeficient
    m = 10   # number of nodes
    N = 100  # number of nodes

    W = diffRegressive(l,T,alpha,m,N)
    print("W in another t (0:m)")
    print(W[:,49])
    
main()
