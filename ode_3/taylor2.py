# Taylor second order method

# INPUT: y' = f(x,y)  : diferential equation
#        y(x(0)) = y0 : initial condition
#        xf           : variable final which define range
#        n            : number of points
# OUTPUT: Aproximation solution of ordinary diferential equation whit its
#         respective graph


import numpy as np
import matplotlib.pyplot as plt

def taylor2(f,x0,y0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])

    y[0] = y0
    for i in range(1,n):
        y[i] = y[i-1] + h*f(x[i-1],y[i-1]) + (h**2)/2 * (dfx(f,x[i-1],y[i-1])+f(x[i-1],y[i-1])*dfy(f,x[i-1],y[i-1]))

    plt.plot(x,y,'o')
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.title("Solution aprox")
    plt.show()

def function(x,y):
    return -y+np.sin(x) # modified function to use

def dfx(f,x,y):
    h = 0.00001
    return (f(x+h,y)-f(x,y))/h

def dfy(f,x,y):
    h = 0.00001
    return (f(x,y+h)-f(x,y))/h
    
def main():

    x0 = 0
    y0 = 1
    xf = 10
    n = 101
    taylor2(function,x0,y0,xf,n,)

main()
    
