# Shooting Method
"""
# INPUT : y'' = p(x)y' + q(x)y + r(x) ,  y(a) = alfa, y(b) = beta

         p(x)
         q(x)
         r(x)
         a,b
         alfa,beta 

# OUTPUT : approximations w1,i to y(xi)
                          w2,i to y'(xi)  to all i = 0, ... ,N
"""

import numpy as np
import rungekutta4-system as rk4

def shootingMethod(p,q,r,a,b,alfa,beta,n):

    # n points -> n-1 sub intervals
    
    h = (b-a)/(n-1)
    u = np.zeros([2,n+1])
    v = np.zeros([2,n+1])
    v[1,0] = 1
    x = np.linspace(a,b,n)

    for i in range(0,n-1):
        #Apply Runge -  Kutta to u1,i+1  u2,i+1  v1,i+1  v2,i+1

    w = np.zeros([2,n+1])

    for i in range(1,n-1):
        w[0,i] = u[0,i] + w[1,0]*v[0,i]
        w[1,i] = u[1,i] + w[1,0]*v[1,i]
        x = a+i*h

    return (x,w)

def fun1(x, y, z):

    
def fun2(x, y, z):

    

def main():

    shootingMethod()
    
main()
