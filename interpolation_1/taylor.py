#Taylor interpolation around XO
#input : XO value for interpolation
#      : X value to evaluate
#      : n degree of taylor aprox
#ouput : Y value of interpolation function in X

import numpy as np
import math

def taylor(n,XO,Y,X):

    t = [None]*(n+1)
    for i in range(len(t)):
        t[i] = ((X-XO)**i)/math.factorial(i)
    
    
    return np.dot(t,Y)

def main():

    n = 4        #degree
    XO = np.pi/2
    X  = np.pi/3 #point to evaluate

    #f sin(x)
    #y derivate values of f in x0: f, f' ,f'',...
    Y = [1,0,-1,0,1]
    
    R = taylor(n,XO,Y,X)
    
    print("Aproximate is: \n",R)

main()
