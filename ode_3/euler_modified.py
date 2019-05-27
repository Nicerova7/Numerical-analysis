# Euler modified method

# Description : In this method we use the method euler with a midpoint
#               to calculate the f(x,y), with a midpoint the error decreases.
#               This lead to the family of Runge-Kutta methods.

# INPUT: y' = f(x,y)  : diferential equation
#        y(x(0)) = y0 : initial condition
#        xf           : variable final which define range
#        n            : number of points
# OUTPUT: Aproximation solution of ordinary diferential equation whit its
#         respective graph


import numpy as np
import matplotlib.pyplot as plt

def euler(x0,y0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])

    y[0] = y0
    for i in range(1,n):
        y[i] = h*(f(x[i-1]+h/2,y[i-1]+f(x[i-1],y[i-1])*h/2))+y[i-1]

    plt.plot(x,y,'o')
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.title("Solution aprox")
    plt.show()

def f(x,y):
    return -y+np.sin(x) # modify function to use
    
def main():

    x0 = 0
    y0 = 1
    xf = 10
    n  = 50
    euler(x0,y0,xf,n,)

main()
