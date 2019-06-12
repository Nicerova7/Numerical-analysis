# Runge Kutta 4th order to edo system of two equations.

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

def solveSystem(f1,f2,x0,y0,z0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])
    z = np.zeros([n,1])

    y[0] = y0
    z[0] = z0
    
    for i in range(1,n):
        k1 = h*f1(x[i-1],y[i-1],z[i-1])
        l1 = h*f2(x[i-1],y[i-1],z[i-1])
        k2 = h*f1(x[i-1]+h/2,y[i-1]+k1/2,z[i-1]+l1/2)
        l2 = h*f2(x[i-1]+h/2,y[i-1]+k1/2,z[i-1]+l1/2)
        k3 = h*f1(x[i-1]+h/2,y[i-1]+k2/2,z[i-1]+l2/2)
        l3 = h*f2(x[i-1]+h/2,y[i-1]+k2/2,z[i-1]+l2/2)
        k4 = h*f1(x[i-1]+h, y[i-1]+k3,z[i-1]+l3)
        l4 = h*f2(x[i-1]+h, y[i-1]+k3,z[i-1]+l3)
        
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6        
        z[i] = z[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6
        
    return np.transpose(y),np.transpose(z)
