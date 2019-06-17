# Power method also known as the Von Mises Iteration

"""
 A*v = lambda*v
 
 INPUT: Matrix A 
        Vector x random...
        error

 OUTPUT: Vector v and scalar lambda

"""

import numpy as np

def power(A,x,err):

	y = A.dot(x)
	c = np.max(np.abs(y))
	y = (1/c)*y
	if(np.sum(np.abs(y-x)) >= err):
		return power(A,y,err)
	else:
		return (c,y)


A = np.array([[4,2,2,1],[2,-3,1,1],[2,1,3,1],[1,1,1,2]],dtype='d')
X0 = np.array([2,3,-1,12],dtype='f').T

lambd,v = power(A,X0,1E-7)
print("Lambda is: "+str(lambd))
print("v vector is: "+str(v))
