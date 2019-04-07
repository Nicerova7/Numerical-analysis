#Simpson rule
import numpy as np

def simpson(f,a,b,n):

    h = (b-a)*1.0/n
    x = np.linspace(a,b,n+1)
    C = np.ones([n+1,1])
    C[1:n:2] = 4
    C[2:n-1:2] = 2
    z = (h/3)*np.dot(f(x),C)
    return z

def func(x):

    z = x**2 #
    return z

print(simpson(func,1,2,100))
