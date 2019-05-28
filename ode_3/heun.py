# Heun method

# INPUT: y' = f(x,y)  : diferential equation
#        y(x(0)) = y0 : initial condition
#        xf           : variable final which define range
#        n            : number of points
# OUTPUT: Aproximation solution of ordinary diferential equation whit its
#         respective graph


import numpy as np
import matplotlib.pyplot as plt

def heun(x0,y0,xf,n):

    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])

    y[0] = y0
    for i in range(1,n):
        y[i] = (h/2)*(f(x[i-1],y[i-1])+f(x[i],y[i-1]+h*f(x[i-1],y[i-1])))+y[i-1]
        # y_i+1 = y_i + h/2 * (f(x_i,y_i) + f(x_i+1,y_i+h*f(x_i,y_i)))

    plt.plot(x,y,'o')
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.title("Solution aprox")
    plt.show()

def f(x,y):
    return -y+np.sin(x) # modified function to use
    
def main():

    x0 = 0
    y0 = 1
    xf = 10
    n = 101
    heun(x0,y0,xf,n,)

main()
    
