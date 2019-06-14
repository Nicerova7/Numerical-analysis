# Finite difference method lineal case

"""
y'' = p(x)y' + q(x)y + r(x)   ,  a <= x <= b
                                 y(a) = alpha , y(b) = beta

INPUT : p(x), q(x), r(x) , a , b , alpha , beta , n

OUTPUT : wi as approximation of y
"""

import numpy as np
import matplotlib.pyplot as plt

def dif(aa,bb,alpha,beta,n):
    
    a = np.zeros([n+2])  # 2 cuz we need w_0 and w_n+1
    b = np.zeros([n+2])  
    c = np.zeros([n+2])  
    d = np.zeros([n+2])
    
    # n number of x points
    h = (bb-aa)/(n+2)
    print(h)
    x = aa + h
    a[1] = 2.0 + (h**2)*q(x)
    b[1] = -1.0 + (h/2)*p(x)
    d[1] = -(h**2)*r(x) + (1 + (h/2)*p(x))*alpha
    
    for i in range(2,n):
        x = aa + i*h
        a[i] = 2.0 + (h**2)*q(x)
        b[i] = -1.0 + (h/2)*p(x)
        c[i] = -1.0 - (h/2)*p(x)
        d[i] = -(h**2)*r(x)
	
    x = bb-h
    a[n] = 2.0 + (h**2)*q(x)
    c[n] = -1.0 - (h/2)*p(x)
    d[n] = -(h**2)*r(x) + (1.0 - (h/2)*p(x))*beta
    
    l = np.zeros([n+2])
    u = np.zeros([n+2])
    z = np.zeros([n+2])
    
    # Crout algorithm
    l[1] = a[1]
    u[1] = b[1]/a[1]
    z[1] = d[1]/l[1]
    
    for i in range(2,n):
        l[i] = a[i]-c[i]*u[i-1]
        u[i] = b[i]/l[i]
        z[i] = (d[i] - c[i]*z[i-1])/l[i]
	
	
    l[n] = a[n] - c[n]*u[n-1]
    z[n] = (d[n] -c[n]*z[n-1])/l[n]
    
    w = np.zeros([n+2])
    
    w[0] = alpha
    w[n+1] = beta
    w[n] = z[n]
	
    for i in range(n-1,0,-1):
        w[i] = z[i] - u[i]*w[i+1]
	    
    print(w)
    return w
    
	
def p(x):
    return -2/x
	
def q(x):
    return 2/(x**2)

def r(x):
    return np.sin(np.log(x))/(x**2)
	
def main():
    a = 1.0
    b = 2.0
    alpha = 1.0
    beta = 2.0
    n = 8
    
    w = dif(a,b,alpha,beta,n)
    
    x = np.linspace(a,b,n+2) # add x_0 and x_n+1
    plt.plot(x,w)
    plt.show()

main()
