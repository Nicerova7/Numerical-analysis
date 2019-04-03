# Recursive method of Newton
import numpy as np

i = 0

def close_enough(v1,v2):
    tolerance = 1E-5
    if (np.abs(v1 - v2) < tolerance):
        return True
    else:
        return False

def fixed_point(f,first_guess,g):

    def try_guess(guess):
        global i
        fx = g(guess) #only print
        nextt = f(guess)
        i += 1
        print("{0:2d} {1:8.5f} {2:8.5f}".format(i,guess,fx))
        if(close_enough(nextt,guess)):
            return nextt
        else:
            return(try_guess(nextt))

    return try_guess(first_guess)

def deriv(g):

    dx = 0.00001
    return lambda x: (g(x+dx)-g(x))/dx

def trans_newton(g):

    return lambda x: x-(g(x)/deriv(g)(x))

def newton_method(g,guess):

    print('Solution aprox x = {0:5.5f}'.format(fixed_point(trans_newton(g),guess,g)))

def f1 (x):

    z = np.pi*0.1*20*(x-300)+np.pi*0.1*0.8*5.67*10E-8*(x**4-300*4)-100
    return z

def main():
    print('-'*25)
    print(" {0:3s} {1:6s} {2:8s}".format('i','x','f(x)'))
    newton_method(f1,2)
    print('-'*25)

main()
