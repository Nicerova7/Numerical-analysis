# False position method
import numpy as np

tolerance = 0.00001

def false_position(f,a,b):

    def close_enough(v1,v2):
        if(np.abs(v1-v2)<tolerance):
            return True
        else:
            return False

    c = (a*f(b)-b*f(a))/(f(b)-f(a))*1.0
    nextt = f(c)
    if (close_enough(nextt,0.0)):
        return c
    else:
        fa = f(a)
        fb = f(b)
        fc = f(c)
        if (fa*fc<0):
            return (false_position(f,a,c))
        if (fc*fb<0):
            return (false_position(f,c,b))

def func(x):
    z = x**4-2*x**3-4*x**2+4*x+4
    return z

print(false_position(func,-2,-1))
