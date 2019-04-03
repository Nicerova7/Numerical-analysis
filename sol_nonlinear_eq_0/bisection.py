# Bisection
import numpy as np

tolerance = 0.000000001

def bisection(f,a,b):

    def close_enough(v1,v2):
        if(np.abs(v1-v2)<tolerance):
            return True
        else:
            return False

    c = (a+b)/2.0
    nextt = f(c)
    if (close_enough(nextt,0.0)):
        return c
    else:
        fa = f(a)
        fb = f(b)
        fc = f(c)
        if (fa*fc<0):
            return (bisection(f,a,c))
        if (fc*fb<0):
            return (bisection(f,c,b))

def func(x):
    z = x**4-2*x**3-4*x**2+4*x+4
    return z

print(bisection(func,-2,-1))
