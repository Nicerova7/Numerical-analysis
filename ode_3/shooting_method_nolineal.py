# Shooting Method no lineal
"""
# INPUT : y'' = f(x,y,y')  ,  a <= x <= b  
                           ,  y(a) = alpha
                           ,  y(b) = beta

          tol = tolerance
          M   = Number max of iterations

# OUTPUT : approximations w1,i to y(xi)
                          w2,i to y'(x)   to all i = 0, ..., N
           OR message indicating that the maximum number of iterations
              was exceeded.
"""

import numpy as np
import rungekutta4mod as rk4
import matplotlib.pyplot as plt

def shootingNoLineal(a,b,alpha,beta,n,tol,M):

    # n intervals.
    
    h = (b-a)/n
    k = 1
    TK = (beta-alpha)/(b-a)
    u1 = 0
    u2 = 1
    
    w = np.zeros([2,n+1])
    u = np.zeros([2,n+1])
    x = np.linspace(a,b,n+1)
    
    while(k <= M):
        w = rk4.solveSystem(fun1,f,a,alpha,TK,b,n+1)
        u = rk4.solveSystem(fun1,f2(u1,u2),a,0,1,b,n+1)
        u1 = u[0,n]
        u2 = u[1,n]
        if(np.abs(w[0,n]-beta)<tol):
            return x,w
        k = k+1
        TK = TK-(w[0,n]-beta)/u1

    print("Maximum number of iterations exceeded")
    return [0],0
        
def fun1(x,y,z):
    return z

def f(x,y,z):
    # f(x,y,y')
    return (1/8)*(32 + 2*x**3 - y*z)  # modify here to change function

# df(x,y,y')/dy
def dfy(x,y,z): 
    dy = 0.0001
    return (f(x,y+dy,z)-f(x,y,z))/dy

# df(x,y,y')/dy'
def dfyp(x,y,z):
    dzz = 0.0001
    return (f(x,y,z+dzz)-f(x,y,z))/dzz

def f2(u1,u2):
    def f22(x,y,z):
        return dfy(x,y,z)*u1 + dfyp(x,y,z)*u2
    return f22
    
def main():

    a = 1.0
    b = 3.0
    alpha = 17
    beta = 43/3
    n = 20
    tol = 0.00001
    M = 33

    x,w = shootingNoLineal(a,b,alpha,beta,n,tol,M)

    if(len(x)!=1):
        plt.plot(x,w[0],label='Y')
        plt.plot(x,w[1],label='Y\'')
        plt.xlabel("X values")
        plt.ylabel("Y values")
        plt.legend()
        plt.title("Shooting method to no lineal case")
        plt.show()
    
main()
