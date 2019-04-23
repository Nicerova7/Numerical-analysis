#Simpson 3/8 rule
import numpy as np

def simpson(f,a,b,n):

    h = (b-a)*1.0/n
    x = np.linspace(a,b,n+1)
    C = np.ones([n+1,1])     
    C[1:n-1:3] = 3
    C[2:n:3] = 3
    C[3:n-2:3] = 2
    z = sum((3*h/8)*np.dot(f(x),C)) # sum to make sure if we have funct const
    return z

def func(x):

    z = x**2 #
    return z

print(simpson(func,1,2,100))
