# Finite differences method to solve Poisson equation.

"""
 d^2u(x,y)/dx^2  + d^2u(x,y)/dy^2 = f(x,y)  

 a <= x <= b
 c <= y <= d

 Boundary Conditions:
 
 u(x,y) = g(x,y)   If x = a   or   x = b  and  c <= y <= d
  
 and

 u(x,y) = g(x,y)   If y = c   or   y = d and   a <= x <= b

INPUT: a, b, c, d
        m,n >= 3
        tol -> tolerance
        N   -> Number of iterations. 

OUTPU: w_i,j approximation to u(xi,yi) to all i = 1, ... , n-1 and
                                       to all j = 1, ... , m-1 

       Or a message that the number of iterations was exceeded

"""
import numpy as np


def finDifPoisson(a,b,c,d,m,n,tol,N):

    h = (b-a)/n
    k = (d-c)/m

    x = np.zeros([n-1])
    y = np.zeros([n-1])
    w = np.zeros([n-1,m-1])
    
    for i in range(1,n): x[i-1] = a + i*h
    for j in range(1,m): y[j-1] = c + j*k

    lambd = h**2/k**2
    mu = 2*(1+lambd)
    l = 0

    # Gauss-Siedel
    while(l <= N):
        z = (-h**2*f(x[0],y[m-2]) + g(a,y[m-2]) + lambd*g(x[0],d) +
             lambd*w[0,m-3] + w[1,m-2])/mu

        norm = abs(z - w[0,m-2])
        w[0,m-2] = z

        for i in range(2,n-1): # 2 .... n-2
            z = (-h**2*f(x[i-1],y[m-2]) + lambd*g(x[i-1],d) + w[i-2,m-2] +
                 w[i,m-2] + lambd*w[i-1,m-3])/mu
            if (abs(w[i-1,m-2] - z) > norm):
                norm = abs(w[i-1,m-2] - z)
            w[i-1,m-2] = z

        z = (-h**2*f(x[n-2],y[m-2]) + g(b,y[m-2]) + lambd*g(x[n-2],d) +
             lambd*w[n-3,m-2] + lambd*w[n-2,m-3])/mu

        if(abs(w[n-2,m-2] - z) > norm):
            norm = abs(w[n-2,m-2] - z)
        w[n-2,m-2] = z

        for j in range(m-2,1,-1):
            z = (-h**2*f(x[0],y[j-1]) + g(a,y[j-1]) + lambd*w[0,j] +
                 lambd*w[0,j-2] + w[1,j-1])/mu
            if( abs(w[0,j-1] - z) > norm):
                norm = abs(w[0,j-1] - z)
            w[0,j-1] = z

            for i in range(2,n-1):
                z = (-h**2*f(x[i-1],y[j-1]) + w[i-2,j-1] + lambd*w[i-1,j] +
                    w[i,j-1] + lambd*w[i-1,j-2])/mu
                if(abs(w[i-1,j-1] - z) > norm):
                    norm = abs(w[i-1,j-1] - z)
                w[i-1,j-1] = z

            z = (-h**2*f(x[n-2],y[j-1]) + g(b,y[j-1]) + w[n-3,j-1] +
                     lambd*w[n-2,j] + lambd*w[n-2,j-2])/mu

            if(abs(w[n-2,j-1] - z) > norm):
                norm = abs(w[n-2,j-1] - z)
            w[n-2,j-1] = z

        z = (-h**2*f(x[0],y[0]) + g(a,y[0]) + lambd*g(x[0],c) +
             lambd*w[0,1] + w[1,0])/mu

        if(abs(w[0,0] - z) > norm):
            norm = abs(w[0,0] - z)
        w[0,0] = z

        for i in range(2,n-1):
            z = (-h**2*f(x[i-1],y[0]) + lambd*g(x[i-1],c) + w[i-2,0] +
                 lambd*w[i-1,1] + w[i,1])/mu
            if(abs(w[i-1,0] - z) > norm):
                norm = abs(w[i-1,0] - z)
            w[i-1,0] = z

        z = (-h**2*f(x[n-2],y[0]) + g(b,y[0]) + lambd*g(x[n-2],c) +
             w[n-3,0] + lambd*w[n-2,1])/mu

        if(abs(w[n-2,0] - z) > norm):
            norm = abs(w[n-2,0] - z)
        w[n-2,0] = z
                
        if(norm <= tol):
            print(l)
            return x,j,w

        l += 1
    
    print("Max number of iterations was exceeded")
    
        
def f(x,y):
    return x*np.exp(y)

def g(x,y):
    if(x == 0):
        return 0
    if(x == 2):
        return 2*np.exp(y)
    if(y == 0):
        return x
    if(y == 1):
        return np.exp(1)*x

def main():
    a = 0.0
    b = 2.0
    c = 0.0
    d = 1.0
    m = 5
    n = 6
    tol = 1E-10
    N = 100

    x,y,w = finDifPoisson(a,b,c,d,m,n,tol,N)
    print(w)

main()
