# Euler to equation system of two odes.

"""
INPUT:  dy/dx = f1(x,y,z)  : differential equation 1
        dz/dx = f2(x,y,z)  : differential equation 2
        y(x(0)) = y0 : initial condition 1
        z(x(0)) = z0 : initial condition 2
        xf           : variable final which define range
        n            : number of points
OUTPUT: Aproximation solution of differential equations system whit its graph
""" 

import numpy as np
import matplotlib.pyplot as plt

def euler(x0,y0,z0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])
    z = np.zeros([n,1])

    y[0] = y0
    z[0] = z0
    
    for i in range(1,n):
        y[i] = h*(f1(x[i-1],y[i-1],z[i-1])) + y[i-1]
        z[i] = h*(f2(x[i-1],y[i-1],z[i-1])) + z[i-1]
        
    plt.plot(x,y,label='Y aprox')
    plt.plot(x,z,label='Z aprox')
    plt.xlabel("X values")
    plt.ylabel("Y,Z values")
    plt.legend()
    plt.title("Sytem of First Order Differential Equations")
    plt.show()

def f1(x,y,z):
    return 32*y + 66*z + 2*x/3 + 2/3

def f2(x,y,z):
    return -66*y - 133*z - 1*x/3 - 1/3

def main():
    
    x0 = 0.0
    y0 = 1/3
    z0 = 1/3
    xf = 0.5
    n  = 21  # 21 points, 20 spaces so h = 0.025
    euler(x0,y0,z0,xf,n,)

main()

"""
Notice that the error in the end of the range of the X values is so big.
"""
    
