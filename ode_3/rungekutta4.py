# Runge-Kutta 4th order

# INPUT: y' = f(x,y)  : diferential equation
#        y(x(0)) = y0 : initial condition
#        xf           : variable final which define range
#        n            : number of points
# OUTPUT: Aproximation solution of ordinary diferential equation whit its
#         respective graph


import numpy as np
import matplotlib.pyplot as plt

def rg4(x0,y0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])

    y[0] = y0
    for i in range(1,n):
        k1 = h*f(x[i-1],y[i-1])
        k2 = h*f(x[i-1]+h/2,y[i-1]+k1/2)
        k3 = h*f(x[i-1]+h/2,y[i-1]+k2/2)
        k4 = h*f(x[i-1]+h, y[i-1]+k3)
        y[i] = y[i-1]+(k1+2*k2+2*k3+k4)/6
        
    plt.plot(x,y,'*')
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.title("Solution Runge-Kutta 4th Order")
    plt.show()

def f(x,y):
    return -y+np.sin(x) # modify function to use this file.

def main():

    x0 = 0.0
    y0 = 1
    xf = 10 
    n  = 101 
    rg4(x0,y0,xf,n)

main()
                                      
    
