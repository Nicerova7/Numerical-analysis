# Shooting Method
"""
# INPUT : y'' = p(x)y' + q(x)y + r(x) ,  y(a) = alfa, y(b) = beta

         p(x)
         q(x)
         r(x)
         a,b
         alpha,beta 

# OUTPUT : approximations w1,i to y(xi)
                          w2,i to y'(xi)  to all i = 0, ... ,N
"""

import numpy as np
import rungekutta4mod as rk4
import matplotlib.pyplot as plt

def shootingMethod(p,q,r,a,b,alpha,beta,n):

    # n subintervals
    
    h = (b-a)/n
    u = np.zeros([2,n+1])
    v = np.zeros([2,n+1])
    x = np.linspace(a,b,n+1)

    #Apply Runge -  Kutta to u1,i+1  u2,i+1  v1,i+1  v2,i+1
    u = rk4.solveSystem(fun1,fun2(p,q,r),a,alpha,0,b,n+1)
    v = rk4.solveSystem(fun1,fun2(p,q,r),a,0,1,b,n+1)
                        
    w = np.zeros([2,n+1])
    w[0,0] = alpha
    w[1,0] = (beta-u[0,n])/v[0,n]
    
    for i in range(1,n):
        w[0,i] = u[0,i] + w[1,0]*v[0,i]
        w[1,i] = u[1,i] + w[1,0]*v[1,i]

    return(x,w)


"""
u'' = p(x)u' + q(x)u + r(x)

u'  = n  -> n(a)  = alpha
u'' = n' -> n'(a) = 0

"""
def fun1(x, u, n):
    return n

"""
v'' = p(x)u' + q(x)u + r(x)

v'  = n  -> n(a)  = 0
v'' = n' -> n'(a) = 1
"""
def fun2(p, q, r):
    def fun22(x,v,n):
        return p(x)*n + q(x)*v + r(x)
    return fun22


# Modify p(x) q(x) and r(x) to use.

def p(x):
    return 0

def q(x):
    return -4

def r(x):
    return 0

def main():

    a = 0.0
    b = np.pi
    alpha = -2.0
    beta = 10.0
    n = 40
    x,w = shootingMethod(p,q,r,a,b,alpha,beta,n) # To use in general

    plt.plot(x,w[0],label='Y')
    plt.plot(x,w[1],label='Y\'')
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.legend()
    plt.title("Solution with shooting method")
    plt.show()
    
main()
