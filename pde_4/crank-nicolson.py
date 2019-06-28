## Crank-Nicolson method to parabolic eq

"""
   du(x,t)/dt  - alpha^2*d^2u(x,t)/dt^2 = 0
   0 < x < l
   0 < t < T

   Boundary Conditions
   
   u(0,t) = u(l,t) = 0 ,  0 < t < T

   Initial Conditions
   u(x,0) = f(x) ,   0 <= x <= l

INPUT: l,T, max values
       alpha
       m >= 3
       N >= 1

OUTPUT: aproximation wij to u(xi,tj) to all  i = 1, ... , m-1
                                             j = 1, ... , N
"""
import numpy as np

def crankNicolson(l,T,alpha,m,N):

    h = l/m
    k = T/N
    lambd = alpha**2*k/(h**2)

    w = np.zeros([m])
    l = np.zeros([m])
    u = np.zeros([m])
    z = np.zeros([m])
    W = np.zeros([m+1,N])
    
    for i in range(1,m): w[i-1] = f(i*h) # initial values

    l[0] = 1 + lambd
    u[0] = -lambd/(2*l[0])

    for i in range(2,m-1):
        l[i-1] = 1 + lambd + lambd*u[i-2]/2
        u[i-1] = -lambd/(2*l[i-1])

    l[m-2] = 1 + lambd + lambd*u[m-3]/3

    for j in range(1,N+1):

        z[0] = ((1-lambd)*w[0] + lambd*w[1]/2)/l[0]

        for i in range(2,m):
            z[i-1] = ((1-lambd)*w[i-1] + lambd*(w[i] + w[i-2] + z[i-2])/2 )/l[i-1]

        w[m-2] = z[m-2]

        for i in range(m-2,0,-1):
            w[i-1] = z[i-1] - u[i-1]*w[i]

        W[0,j-1] = f(0) # restriction
        W[1:,j-1] = w

    return W
        
def f(x):
    return np.sin(np.pi*x)
    
def main():
    l = 1.0
    T = 0.5
    alpha = 1.0
    m = 10
    N = 50

    w = crankNicolson(l,T,alpha,m,N)
    print("W in another t (0:m)")
    print(w[:,49])

main()
