# Fixed points method
import numpy as np

tolerance = 0.00001

def fixed_points(f,first_guess):

    def close_enough(v1,v2):
        if (np.abs(v1-v2)/(v1+np.finfo(np.float32).eps) < tolerance): #relative error
            return True
        else:
            return False

    def try_guess(guess):
        nextt = f(guess)
        if (close_enough(guess,nextt)):
            return nextt
        else:
            return(try_guess(nextt))

    return try_guess(first_guess)

def func(x):
    return (2-np.exp(1)**x+x**2)/3

print(fixed_points(func,1/6))
