# Integration of Gauss-Legendre

# INPUT: n degree of a polinomial
# OUPUT: integrate value of -1 to 1

import numpy as np

def gaussian(n):
    
    if n == 1:
        x = [0]
        w = [2]
        return w,x
    elif n == 2:
        x = [-1/np.sqrt(3),1/np.sqrt(3)]
        w = [1,1]
        return w,x
    elif n == 3:
        x = [-1/np.sqrt(3/5), 0, 1/np.sqrt(3/5)]
        w = [5/9, 8/9, 5/9]
        return w,x
    else:
        print("Not possible")
 

def function(x):
    return x**2
    
def main():

    [w,x] = gaussian(2)
    vfunc = np.vectorize(function) # to apply function to each value in array
    int_gauss = np.dot(w,vfunc(x))
    print(int_gauss)

main()
