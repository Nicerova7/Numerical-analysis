# Shooting Method
"""
#INPUT : y'' = p(x)y' + q(x)y + r(x) ,  y(a) = alfa, y(b) = beta

         p(x)
         q(x)
         r(x)
         a,b
         alfa,beta         
"""

import numpy as np

def shootingMethod(p,q,r,a,b,alfa,beta,n):

    # n points
    
    h = (b-a)/(n-1)
    u = np.zeros([2,n+1])
    v = np.zeros([2,n+1])
    x = np.linspace(a,b,n)

    for i in range(0,n-1):
        #Apply Runge -  Kutta to u1,i+1  u2,i+1  v1,i+1  v2,i+1

    w = np.zeros([2,n+1])

    for i in range(1,n-1):
        w[1,i] = u[1,i] + w[2,0]*v[1,i]
        w[2,i] = u[2,i] + w[2,0]*v[2,i]
        x = a+i*h

    return (x,w)

def main():

    shootingMethod()
    
main()
